#include <complex>
#include <ctime>
#include <fftw3.h>
#include "dist_model.h"

#include "ncode_common/src/substitute.h"
#include "../common.h"

namespace ctr {

FFTRunner::FFTRunner(size_t fft_size) : fft_size_(fft_size) {
  out_.resize(fft_size);
  fftw_complex* out_ptr = reinterpret_cast<fftw_complex*>(out_.data());
  plan_fw_ =
      fftw_plan_dft_1d(fft_size, out_ptr, out_ptr, FFTW_FORWARD, FFTW_MEASURE);
  plan_bw_ =
      fftw_plan_dft_1d(fft_size, out_ptr, out_ptr, FFTW_BACKWARD, FFTW_MEASURE);
}

FFTRunner::~FFTRunner() {
  fftw_destroy_plan(plan_fw_);
  fftw_destroy_plan(plan_bw_);
}

const std::vector<Complex>& FFTRunner::ComputeForward(
    const std::vector<double>& input) {
  CHECK(fft_size_ >= input.size()) << fft_size_ << " vs " << input.size();
  for (size_t i = 0; i < input.size(); ++i) {
    out_[i] = {input[i], 0.0};
  }

  for (size_t i = input.size(); i < fft_size_; ++i) {
    out_[i] = {0, 0};
  }

  fftw_execute(plan_fw_);
  return out_;
}

std::vector<double> FFTRunner::ComputeBackward(
    const std::vector<Complex>& input) {
  CHECK(input.size() == out_.size());
  std::copy(input.begin(), input.end(), out_.begin());

  fftw_execute(plan_bw_);
  std::vector<double> out(fft_size_);
  for (size_t i = 0; i < fft_size_; ++i) {
    out[i] = out_[i].real() / fft_size_;
  }
  return out;
}

std::string ConvolutionTimingData::ToString() const {
  return nc::Substitute("binning: $0µs, FFT: $1µs, IFFT: $2µs", binning.count(),
                        fft.count(), ifft.count());
}

std::string ProbModelTimingData::ToString() const {
  return nc::Substitute(
      "queue sim: $0µs, aggregate split: $1µs, $2, total: $3µs",
      queue_simulation.count(), split_aggregates.count(),
      convolution_timing.ToString(), total.count());
}

Distribution::Distribution(const std::vector<double>& probabilities,
                           uint64_t base)
    : probabilities_(probabilities), base_(base), bin_size_(1) {
  CHECK(!probabilities.empty());
  CHECK(probabilities.back() != 0);
  max_element_ = base_ + probabilities.size() - 1;
}

Distribution::Distribution(const std::vector<uint64_t>& values) : bin_size_(1) {
  CHECK(!values.empty());
  uint64_t min_element = std::numeric_limits<uint64_t>::max();
  uint64_t max_element = 0;

  std::map<uint64_t, uint64_t> counts;
  for (uint64_t value : values) {
    if (value > max_element) {
      max_element = value;
    }

    if (value < min_element) {
      min_element = value;
    }

    ++counts[value];
  }

  base_ = min_element;
  max_element_ = max_element;
  probabilities_.resize(max_element - min_element + 1, 0.0);
  for (const auto& value_and_count : counts) {
    uint64_t value = value_and_count.first;
    double count = value_and_count.second;
    double probability = count / values.size();
    probabilities_[value - min_element] = probability;
  }
}

Distribution::Distribution(uint64_t base, uint64_t bin_size,
                           uint64_t max_element)
    : base_(base), bin_size_(bin_size), max_element_(max_element) {}

Distribution Distribution::Bin(size_t new_bin_size) const {
  CHECK(bin_size_ == 1);
  CHECK(!probabilities_cumulative_.empty());
  Distribution new_dist(base_, new_bin_size, max_element_);

  std::vector<double>& new_probabilities = new_dist.probabilities_;
  size_t est_count = probabilities_.size() / new_bin_size + 1;
  new_probabilities.resize(est_count, 0);

  size_t offset = new_bin_size - 1;
  double sum_so_far = 0;
  for (size_t i = 0; i < est_count - 1; ++i) {
    size_t index = new_bin_size * i + offset;
    double total = probabilities_cumulative_[index];

    new_probabilities[i] = total - sum_so_far;
    sum_so_far = total;
  }

  new_probabilities[est_count - 1] =
      probabilities_cumulative_.back() - sum_so_far;
  return new_dist;
}

std::map<uint64_t, double> Distribution::GetProbabilities() const {
  std::map<uint64_t, double> out;
  for (size_t i = 0; i < probabilities_.size(); ++i) {
    double p = probabilities_[i];
    if (p > 0) {
      uint64_t value = base_ + i * bin_size_;
      out[value] = p;
    }
  }

  return out;
}

std::vector<std::pair<double, double>> Distribution::GetDataSeries(
    bool cumulative) const {
  std::vector<std::pair<double, double>> data;
  double total = 0;
  for (size_t i = 0; i < probabilities_.size(); ++i) {
    double p = probabilities_[i];
    if (p <= 0) {
      continue;
    }

    uint64_t value = base_ + i * bin_size_;
    if (cumulative) {
      total += p;
      data.emplace_back(value, total);
    } else {
      data.emplace_back(value, p);
    }
  }

  return data;
}

double Distribution::Percentile(double p) const {
  CHECK(p >= 0 && p <= 1);

  double lower_p = 0;
  double upper_p = 0;
  ssize_t upper_bound_index = -1;
  for (size_t i = 0; i < probabilities_.size(); ++i) {
    upper_p += probabilities_[i];
    if (upper_p >= p) {
      upper_bound_index = i;
      break;
    }

    lower_p += probabilities_[i];
  }

  if (upper_bound_index == -1) {
    return base_ + (probabilities_.size() - 1) * bin_size_;
  }

  if (upper_bound_index == 0) {
    return base_;
  }

  size_t lower_bound_index = upper_bound_index - 1;
  uint64_t upper_value = base_ + upper_bound_index * bin_size_;
  uint64_t lower_value = base_ + lower_bound_index * bin_size_;

  CHECK(upper_p > lower_p);
  CHECK(upper_value > lower_value);

  double f = (upper_p - p) / (upper_p - lower_p);
  return lower_value + (1 - f) * (upper_value - lower_value);
}

std::vector<uint64_t> Distribution::Sample(size_t count,
                                           std::mt19937* rnd) const {
  std::uniform_real_distribution<double> distribusion(0.0, 1.0);

  std::vector<uint64_t> out;
  for (size_t v_index = 0; v_index < count; ++v_index) {
    double r = distribusion(*rnd);
    for (size_t i = 0; i < probabilities_.size(); ++i) {
      double p = probabilities_[i];
      if (p <= 0) {
        continue;
      }

      uint64_t value = base_ + i * bin_size_;
      if (r <= p) {
        out.emplace_back(value);
        break;
      }

      r -= p;
    }
  }

  return out;
}

void Distribution::CacheCumulativeProbabilities() {
  probabilities_cumulative_.resize(probabilities_.size());
  double total = 0;
  for (size_t i = 0; i < probabilities_.size(); ++i) {
    total += probabilities_[i];
    probabilities_cumulative_[i] = total;
  }
}

Distribution Distribution::Sum(
    const std::vector<const Distribution*>& distributions, size_t resolution,
    FFTRunner* fft_runner, ConvolutionTimingData* timing_data) {
  using namespace std::chrono;

  CHECK(!distributions.empty());
  CHECK(resolution > 1);

  uint64_t min_value = 0;
  uint64_t max_value = 0;
  for (const Distribution* distribution : distributions) {
    min_value += distribution->min_element();
    max_value += distribution->max_element();
  }

  high_resolution_clock::time_point start;
  high_resolution_clock::time_point end;
  microseconds total_fft = microseconds::zero();
  microseconds total_bin = microseconds::zero();

  std::vector<Complex> current_product;
  size_t bin_size = (max_value - min_value) / (resolution - 1) + 1;
  for (const Distribution* distribution : distributions) {
    if (timing_data != nullptr) {
      start = high_resolution_clock::now();
    }
    Distribution binned = distribution->Bin(bin_size);
    if (timing_data != nullptr) {
      end = high_resolution_clock::now();
      total_bin += duration_cast<microseconds>(end - start);
    }

    if (timing_data != nullptr) {
      start = high_resolution_clock::now();
    }
    const std::vector<Complex>& fft =
        fft_runner->ComputeForward(binned.probabilities_raw());
    if (timing_data != nullptr) {
      end = high_resolution_clock::now();
      total_fft += duration_cast<microseconds>(end - start);
    }

    if (current_product.empty()) {
      current_product = fft;
      continue;
    }

    CHECK(current_product.size() == fft.size());
    for (size_t i = 0; i < fft.size(); ++i) {
      current_product[i] *= fft[i];
    }
  }

  Distribution new_dist(min_value, bin_size, max_value);
  if (timing_data != nullptr) {
    timing_data->binning = total_bin;
    timing_data->fft = total_fft;

    start = high_resolution_clock::now();
  }
  new_dist.probabilities_ = fft_runner->ComputeBackward(current_product);
  if (timing_data != nullptr) {
    end = high_resolution_clock::now();
    timing_data->ifft = duration_cast<microseconds>(end - start);
  }

  return new_dist;
}

static std::vector<uint64_t> GetQuantizedValues(const AggregateHistory& history,
                                                size_t quantization,
                                                double fraction) {
  std::vector<uint64_t> values;
  for (size_t bin : history.bins()) {
    values.emplace_back(bin * fraction / quantization);
  }

  return values;
}

ProbModel::AggregateState::AggregateState(const AggregateHistory* history,
                                          size_t quantization)
    : history(history),
      distribution(GetQuantizedValues(*history, quantization, 1.0)) {
  distribution.CacheCumulativeProbabilities();
}

void ProbModel::AddAggregate(uint64_t aggregate_id,
                             const AggregateHistory* history) {
  aggregate_states_.emplace(
      std::piecewise_construct, std::forward_as_tuple(aggregate_id),
      std::forward_as_tuple(history, config_.initial_quantization));
}

std::chrono::milliseconds ProbModel::GetBinSize(
    const ProbModelQuery& query) const {
  using namespace std::chrono;

  milliseconds bin_size = milliseconds::max();
  for (const auto& aggregate_and_fraction : query.aggregates) {
    uint64_t aggregate_id = aggregate_and_fraction.first;

    const AggregateState& aggregate_state =
        nc::FindOrDie(aggregate_states_, aggregate_id);
    const AggregateHistory* history = aggregate_state.history;

    if (bin_size == milliseconds::max()) {
      bin_size = history->bin_size();
    } else {
      CHECK(bin_size == history->bin_size()) << "Bin size mismatch "
                                             << bin_size.count() << " vs "
                                             << history->bin_size().count();
    }
  }

  CHECK(bin_size != milliseconds::max());
  return bin_size;
}

std::vector<uint64_t> ProbModel::SumUpBins(const ProbModelQuery& query,
                                           double* max_sum) const {
  std::vector<uint64_t> bins_total;

  *max_sum = 0;
  for (const auto& aggregate_and_fraction : query.aggregates) {
    uint64_t aggregate_id = aggregate_and_fraction.first;
    double fraction = aggregate_and_fraction.second;
    const AggregateState& aggregate_state =
        nc::FindOrDie(aggregate_states_, aggregate_id);
    const AggregateHistory* history = aggregate_state.history;

    double max_bin = 0;
    const std::vector<uint64_t>& bins = history->bins();
    for (size_t i = 0; i < bins.size(); ++i) {
      double bin_scaled = bins[i] * fraction;
      max_bin = std::max(max_bin, bin_scaled);

      bins_total.resize(std::max(i + 1, bins_total.size()), 0);
      bins_total[i] += bin_scaled;
    }

    *max_sum += max_bin;
  }

  return bins_total;
}

std::chrono::milliseconds ProbModel::SimulateQueue(
    const ProbModelQuery& query, const std::vector<uint64_t>& bins_total,
    double bytes_per_bin) const {
  using namespace std::chrono;

  double time_to_drain_byte_picos =
      (pow(10, 12.0) / static_cast<double>(query.rate.bps())) * 8.0;

  milliseconds max_time_to_drain = milliseconds::zero();
  double leftover = 0;
  for (double bin : bins_total) {
    leftover += (bin - bytes_per_bin);
    leftover = std::max(0.0, leftover);

    uint64_t time_to_drain_picos = time_to_drain_byte_picos * leftover;
    nanoseconds time_to_drain(time_to_drain_picos / 1000);
    if (max_time_to_drain < time_to_drain) {
      max_time_to_drain = duration_cast<milliseconds>(time_to_drain);
    }
  }

  return max_time_to_drain;
}

nc::net::Bandwidth ProbModel::OptimalRate(
    const ProbModelQuery& query, FFTRunner* fft_runner,
    ProbModelTimingData* timing_data) const {
  using namespace std::chrono;

  std::vector<std::unique_ptr<Distribution>> partial_splits;
  std::vector<const Distribution*> to_sum;

  high_resolution_clock::time_point start;
  high_resolution_clock::time_point end;
  microseconds split_total = microseconds::zero();

  for (const auto& aggregate_and_fraction : query.aggregates) {
    uint64_t aggregate_id = aggregate_and_fraction.first;
    double fraction = aggregate_and_fraction.second;

    const AggregateState& aggregate_state =
        nc::FindOrDie(aggregate_states_, aggregate_id);
    if (fraction != 1) {
      // This will be slow -- the distribution in aggregate_state is the one
      // that has the entire aggregate, not just a fraction of it. If we need
      // only a fraction of the aggregate we have to regenerate the distribution
      // here. If most of the aggregates are not being split this is fine, but
      // if they are this will trigger very often.
      if (timing_data != nullptr) {
        start = high_resolution_clock::now();
      }

      const AggregateHistory& history = *aggregate_state.history;
      auto new_dist = nc::make_unique<Distribution>(
          GetQuantizedValues(history, config_.initial_quantization, fraction));
      new_dist->CacheCumulativeProbabilities();

      to_sum.emplace_back(new_dist.get());
      partial_splits.emplace_back(std::move(new_dist));

      if (timing_data != nullptr) {
        end = high_resolution_clock::now();
        split_total += duration_cast<microseconds>(end - start);
      }

    } else {
      to_sum.emplace_back(&aggregate_state.distribution);
    }
  }

  if (timing_data != nullptr) {
    timing_data->split_aggregates = split_total;
  }

  Distribution sum = Distribution::Sum(
      to_sum, config_.distribution_levels, fft_runner,
      timing_data == nullptr ? nullptr : &timing_data->convolution_timing);
  double bytes_per_bin = sum.Percentile(1 - config_.exceed_probability);

  milliseconds bin_size = GetBinSize(query);
  double bins_per_second = 1000.0 / bin_size.count();
  return nc::net::Bandwidth::FromBitsPerSecond(
      config_.initial_quantization * bins_per_second * bytes_per_bin * 8);
}

ProbModelReply ProbModel::SingleQuery(const ProbModelQuery& query,
                                      FFTRunner* fft_runner,
                                      ProbModelTimingData* timing_data) const {
  using namespace std::chrono;
  high_resolution_clock::time_point start_outer;
  high_resolution_clock::time_point end_outer;

  if (timing_data != nullptr) {
    start_outer = high_resolution_clock::now();
  }

  // The rate is in bps. Need to figure out how many bytes per bin this is.
  milliseconds bin_size = GetBinSize(query);
  double rate_bins_in_second = 1000.0 / bin_size.count();
  CHECK(rate_bins_in_second > 0);
  double rate_bytes_per_bin = (query.rate.bps() / 8.0) / rate_bins_in_second;
  CHECK(rate_bytes_per_bin > 0);

  double max_sum_bin;
  std::vector<uint64_t> bins_summed = SumUpBins(query, &max_sum_bin);

  // max_sum_bin now contains the worst possible value a bin can have -- if all
  // peaks of the distributions of all aggregates are in sync. If this fits we
  // can shortcut the evaluation.
  if (max_sum_bin <= rate_bytes_per_bin) {
    return {true, std::chrono::milliseconds::max(), nc::net::Bandwidth::Zero()};
  }

  std::chrono::milliseconds max_queue_size = std::chrono::milliseconds::max();
  if (query.type == ProbModelQuery::QUEUE_SIMULATION ||
      query.type == ProbModelQuery::BOTH) {
    high_resolution_clock::time_point start;
    high_resolution_clock::time_point end;

    if (timing_data != nullptr) {
      start = high_resolution_clock::now();
    }

    max_queue_size = SimulateQueue(query, bins_summed, rate_bytes_per_bin);
    if (timing_data != nullptr) {
      end = high_resolution_clock::now();
      timing_data->queue_simulation = duration_cast<microseconds>(end - start);
    }

    if (max_queue_size > config_.simulated_queue_threshold) {
      if (timing_data != nullptr) {
        end_outer = high_resolution_clock::now();
        timing_data->total =
            duration_cast<microseconds>(end_outer - start_outer);
      }
      return {false, max_queue_size, nc::net::Bandwidth::Zero()};
    }
  }

  if (query.type == ProbModelQuery::CONVOLUTION ||
      query.type == ProbModelQuery::BOTH) {
    nc::net::Bandwidth bw = OptimalRate(query, fft_runner, timing_data);
    if (timing_data != nullptr) {
      end_outer = high_resolution_clock::now();
      timing_data->total = duration_cast<microseconds>(end_outer - start_outer);
    }
    return {query.rate >= bw, max_queue_size, bw};
  }

  LOG(FATAL) << "Bad query type";
  return {false, max_queue_size, nc::net::Bandwidth::Zero()};
}

std::vector<ProbModelReply> ProbModel::Query(
    const std::vector<ProbModelQuery>& queries,
    std::vector<ProbModelTimingData>* timing_data) {
  std::vector<ProbModelReply> out(queries.size());
  if (timing_data != nullptr) {
    timing_data->clear();
    timing_data->resize(queries.size());
  }

  batch_processor_.RunInParallel(
      queries, [this, &out, timing_data](const ProbModelQuery& query, size_t i,
                                         size_t thread_index) {
        FFTRunner* fft_runner = runners_[thread_index].get();
        ProbModelTimingData* timing =
            timing_data == nullptr ? nullptr : &(*timing_data)[i];
        out[i] = SingleQuery(query, fft_runner, timing);
      });
  return out;
}

}  // namespace ctr
