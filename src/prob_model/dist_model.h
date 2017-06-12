#ifndef FUBAR_DIST_MODEL_H
#define FUBAR_DIST_MODEL_H

#include <fftw3.h>
#include <stddef.h>
#include <chrono>
#include <complex>
#include <cstdint>
#include <map>
#include <memory>
#include <condition_variable>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "ncode_common/src/common.h"
#include "ncode_common/src/net/net_common.h"
#include "ncode_common/src/thread_runner.h"

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
  size_t fft_size_;

  // This is where the computation is performed.
  std::vector<Complex> out_;

  // The plans.
  fftw_plan plan_fw_;
  fftw_plan plan_bw_;

  DISALLOW_COPY_AND_ASSIGN(FFTRunner);
};

// A discrete distribution.
class Distribution {
 public:
  // Convolves a number of distributions and returns their sum. The resulting
  // distribution will have 'resolution' levels. If the last argument is not
  // null will populate it with timing data.
  static Distribution Sum(const std::vector<const Distribution*>& distributions,
                          size_t resolution, FFTRunner* fft_runner,
                          ConvolutionTimingData* timing_data);

  Distribution(const std::vector<double>& probabilities, uint64_t base);

  Distribution(const std::vector<uint64_t>& values);

  uint64_t max_element() const { return max_element_; }

  uint64_t min_element() const { return base_; }

  uint64_t bin_size() const { return bin_size_; }

  // Bins this distribution. Only valid if the current bin size is 1. Should
  // have cumulative probabilities cached.
  Distribution Bin(size_t new_bin_size) const;

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

  // Caches cumulative probabilities. Should be called before attempting to bin.
  void CacheCumulativeProbabilities();

 private:
  Distribution(uint64_t base, uint64_t bin_size, uint64_t max_element);

  std::vector<double> probabilities_;
  uint64_t base_;
  uint64_t bin_size_;
  uint64_t max_element_;

  // A cumulative sum of probabilities.
  std::vector<double> probabilities_cumulative_;
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

  void AddAggregate(AggregateId aggregate_id, const AggregateHistory* history);

  // Answers to multiple queries. Item i in the return vector corresponds to
  // query i in 'queries'. If the last argument is not null, will populate it
  // with timing data for each query.
  std::vector<ProbModelReply> Query(
      const std::vector<ProbModelQuery>& queries,
      std::vector<ProbModelTimingData>* timing_data = nullptr);

 private:
  // Number of threads to process queries on.
  static constexpr size_t kParallelThreads = 4;

  // Sums up the bins from the query. Also populates the last argument with the
  // total of all max values of all histories of aggregates in the query.
  std::vector<uint64_t> SumUpBins(const ProbModelQuery& query,
                                  double* max_sum) const;

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
                                 ProbModelTimingData* timing_data) const;

  // Performs a single query.
  ProbModelReply SingleQuery(const ProbModelQuery& query, FFTRunner* fft_runner,
                             ProbModelTimingData* timing_data) const;

  struct AggregateState {
    AggregateState(const AggregateHistory* history, size_t quantization);

    const AggregateHistory* history;
    Distribution distribution;
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
