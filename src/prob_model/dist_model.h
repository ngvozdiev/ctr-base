#ifndef CTR_DIST_MODEL_H
#define CTR_DIST_MODEL_H

#include <fftw3.h>
#include <stddef.h>
#include <chrono>
#include <complex>
#include <cstdint>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "ncode/common.h"
#include "ncode/strutil.h"
#include "ncode/net/net_common.h"
#include "ncode/thread_runner.h"
#include "../common.h"

namespace ctr {
class AggregateHistory;
} /* namespace ctr */

namespace ctr {

struct ConvolutionTimingData {
  // Time to bin data before doing an FFT.
  std::chrono::microseconds binning;

  // Time to do all FFTs, one per aggregate.
  std::chrono::microseconds fft;

  // Time to do a single inverse FFT.
  std::chrono::microseconds ifft;

  std::string ToString() const;
};

using Complex = std::complex<double>;

class FFTRunner {
 public:
  FFTRunner(size_t fft_size);

  ~FFTRunner();

  const std::vector<Complex>& ComputeForward(const std::vector<double>& input);

  std::vector<double> ComputeBackward(const std::vector<Complex>& input);

 private:
  void InitPlans();

  size_t fft_size_;

  // This is where the computation is performed.
  std::vector<Complex> out_;

  // The plans.
  fftw_plan plan_fw_;
  fftw_plan plan_bw_;

  // Initializing the plans is relatively expensive, so they are not initalized
  // until a call to Compute* is made.
  bool init_;

  DISALLOW_COPY_AND_ASSIGN(FFTRunner);
};

// A discrete distribution.
class Distribution {
 public:
  // Convolves a number of distributions and returns their sum. The resulting
  // distribution will have 'resolution' levels. If the last argument is not
  // null will populate it with timing data.
  static std::unique_ptr<Distribution> Sum(
      const std::vector<Distribution*>& distributions, size_t resolution,
      FFTRunner* fft_runner, ConvolutionTimingData* timing_data);

  Distribution(const std::vector<double>& probabilities, uint64_t base);

  Distribution(const std::vector<uint64_t>& values);

  uint64_t max_element() const { return max_element_; }

  uint64_t min_element() const { return base_; }

  uint64_t bin_size() const { return bin_size_; }

  // Bins this distribution. Only valid if the current bin size is 1. Should
  // have cumulative probabilities cached.
  std::unique_ptr<Distribution> Bin(size_t new_bin_size);

  // Returns the probabilities as a map.
  std::map<uint64_t, double> GetProbabilities() const;

  // Returns a series or x,y points that can be used in a histogram or a CDF
  // plot (if cumulative is true).
  std::vector<std::pair<double, double>> GetDataSeries(bool cumulative) const;

  // Returns the p-th percentile of this distribution. The percentile value will
  // be linearly interpolated between the values of the closest cumulative
  // probabilities.
  double Percentile(double p) const;

  // Samples 'count' values from this distribution.
  std::vector<uint64_t> Sample(size_t count, std::mt19937* rnd) const;

  // The probabilities array --- the value of element i is the probability of
  // min_element + i.
  const std::vector<double>& probabilities_raw() const {
    return probabilities_;
  }

 private:
  // Caches cumulative probabilities. Should be called before attempting to bin.
  void CacheCumulativeProbabilities();

  Distribution(uint64_t base, uint64_t bin_size, uint64_t max_element);

  std::vector<double> probabilities_;
  uint64_t base_;
  uint64_t bin_size_;
  uint64_t max_element_;

  // A cumulative sum of probabilities.
  std::vector<double> probabilities_cumulative_;

  // Protects the cumulative probabilities.
  std::mutex cumulative_probs_mu_;

  DISALLOW_COPY_AND_ASSIGN(Distribution);
};

// A single query --- when summing up a set of aggregates, can they fit within a
// given rate?
struct ProbModelQuery {
  // How to process the query.
  enum QueryType {
    // Combines all bins from aggregates and runs them through a simulated
    // queue.
    QUEUE_SIMULATION = 1,

    // Treats each aggrgate's bins as independent random variables, convolves
    // them to get a distribution of their sum and checks if the rate from the
    // query is statistically likely to cause queues.
    CONVOLUTION = 2,

    // Does both.
    BOTH = 3,
  };

  QueryType type;

  // A list aggregates and fractions. All the aggregate's bins will be
  // multiplied by the fraction.
  std::vector<std::pair<AggregateId, double>> aggregates;

  // Rate to test.
  nc::net::Bandwidth rate;

  std::string ToString() const {
    std::vector<std::string> pieces;
    for (const auto& aggregate : aggregates) {
      pieces.emplace_back(nc::StrCat(std::to_string(aggregate.first.src()), " ",
                                     std::to_string(aggregate.first.dst()), " ",
                                     std::to_string(aggregate.second)));
    }
    return nc::Join(pieces, ",");
  }
};

struct ProbModelReply {
  // True if the aggregates from the query fit the rate from the query.
  bool fit;

  // Max simulated queue size. If this is above the threshold in ProbModelConfig
  // fit will be false.
  std::chrono::milliseconds max_simulated_queue_size;

  // The min rate that will fit the traffic with low probability for queues.
  nc::net::Bandwidth optimal_rate;
};

struct ProbModelConfig {
  // Initial quantization to apply to values. Each value will be equal to
  // floor(bin_size / initial_quantization).
  size_t initial_quantization = 1000;

  // How many levels a distribution should have.
  size_t distribution_levels = 2048;

  // Probability the convolved distribution to exceed the queried rate.
  double exceed_probability = 0.005;

  // In order to detect mutually-dependent spikes in multiple aggregates will
  // also run a simple queue simulation and check if the queue is ever above
  // this threshold.
  std::chrono::milliseconds simulated_queue_threshold =
      std::chrono::milliseconds(10);
};

// Per-query timing data.
struct ProbModelTimingData {
  // Time to perform queue simulation.
  std::chrono::microseconds queue_simulation;

  // Total processing time for the query.
  std::chrono::microseconds total;

  // Fractional aggregates will need to be split and rebinned. This is the total
  // time it took to do that.
  std::chrono::microseconds split_aggregates;

  // Timing to do the convolution.
  ConvolutionTimingData convolution_timing;

  std::string ToString() const;
};

// The model associates with each aggregate an AggregateHistory instance. It can
// be queried about whether or not a set of aggregates, when combined, will fit
// in a given rate.
class ProbModel {
 public:
  ProbModel(const ProbModelConfig& config)
      : config_(config), batch_processor_(kParallelThreads) {
    for (size_t i = 0; i < kParallelThreads; ++i) {
      runners_.emplace_back(
          nc::make_unique<FFTRunner>(config.distribution_levels));
    }
  }

  void AddAggregate(const AggregateId& aggregate_id,
                    const AggregateHistory* history);

  void ClearAggregates();

  // Answers to multiple queries. Item i in the return vector corresponds to
  // query i in 'queries'. If the last argument is not null, will populate it
  // with timing data for each query.
  std::vector<ProbModelReply> Query(
      const std::vector<ProbModelQuery>& queries,
      std::vector<ProbModelTimingData>* timing_data = nullptr);

 private:
  // Number of threads to process queries on.
  static constexpr size_t kParallelThreads = 4;

  // Returns the max and the min bins.
  std::pair<double, double> MaxMinBins(const ProbModelQuery& query) const;

  // Sums up the bins from the query. Also populates the last argument with the
  // total of all max values of all histories of aggregates in the query.
  std::vector<uint64_t> SumUpBins(const ProbModelQuery& query) const;

  // Returns the common bin size for all aggregates in the query.
  std::chrono::milliseconds GetBinSize(const ProbModelQuery& query) const;

  // Runs the sum of the aggregates in the query through a simulated queue of
  // the same service rate as the one in the query. Each bin is a single
  // simulated 'packet'. Returns the max queue size in milliseconds.
  std::chrono::milliseconds SimulateQueue(
      const ProbModelQuery& query, const std::vector<uint64_t>& bins_summed,
      double bytes_per_bin) const;

  // The rate that should (most of the time) fit the traffic.
  nc::net::Bandwidth OptimalRate(const ProbModelQuery& query,
                                 FFTRunner* fft_runner,
                                 ProbModelTimingData* timing_data);

  // Performs a single query.
  ProbModelReply SingleQuery(const ProbModelQuery& query, FFTRunner* fft_runner,
                             ProbModelTimingData* timing_data);

  class AggregateState {
   public:
    AggregateState(const AggregateHistory* history, size_t quantization);

    Distribution* distribution();

    const AggregateHistory* history() const { return history_; }

   private:
    // The history.
    const AggregateHistory* history_;

    // The distribution, created on demand.
    std::unique_ptr<Distribution> distribution_;

    // Quantization for the distribution.
    size_t quantization_;

    // Protects the distribution.
    std::mutex distribution_mu_;
  };

  // Configuration
  const ProbModelConfig config_;

  std::map<AggregateId, AggregateState> aggregate_states_;

  // Processes batches of queries.
  nc::ThreadBatchProcessor<ProbModelQuery> batch_processor_;

  // Per-thread FFT runner.
  std::vector<std::unique_ptr<FFTRunner>> runners_;
};

}  // namespace fubar
#endif
