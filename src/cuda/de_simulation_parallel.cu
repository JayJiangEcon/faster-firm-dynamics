// FP64 micro simulation with source-ordered macro and moment reductions.
// Uses GPU-solved policies; optional independent CPU reference comparisons.

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct SolvedInputs {
 std::vector<std::int32_t> policy,defaults;
 std::vector<double> price,profit;
};
thread_local const SolvedInputs* solved_inputs=nullptr;
thread_local std::ostream* report_stream=nullptr;
thread_local std::ostream* error_stream=nullptr;
thread_local std::string dump_directory;
std::ostream& report(){return report_stream?*report_stream:std::cout;}
std::ostream& errors(){return error_stream?*error_stream:std::cerr;}
thread_local std::vector<double> last_moments;
thread_local std::vector<double> last_macro_output;
thread_local bool reuse_macro = false;
thread_local std::vector<double> cached_output, cached_default, cached_spread;
namespace {

constexpr int kExogenous = 625;
constexpr int kCapital = 30;
constexpr int kDebt = 30;
constexpr int kEndogenous = kCapital * kDebt;
constexpr int kA = 5;
constexpr int kZ = 5;
constexpr int kActiveExogenous = kZ * kZ;
constexpr int kFirms = 1000;
constexpr int kTimes = 302;
// Micro moments use released t=51..56 and at most t+2, so firm paths need only
// through zero-based t=57.  Macro moments still require the full 302-period
// aggregate path.
constexpr int kMicroLast = 57;
constexpr int kMacroLast = 300;
constexpr int kStats = 13;
constexpr int kMoments = 17;
// Match the author's promoted default-real literals exactly.
constexpr double kDelta = 0x1.99999a0000000p-4;
constexpr double kRate = 0x1.47ae140000000p-5;

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

void cuda_check(cudaError_t status, const char* expression) {
  if (status != cudaSuccess) {
    fail(std::string(expression) + ": " + cudaGetErrorString(status));
  }
}

#define CUDA_CHECK(expr) cuda_check((expr), #expr)

template <typename T>
std::vector<T> read_exact(const std::string& path, std::size_t count) {
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input) fail("cannot open " + path);
  const auto expected = static_cast<std::streamoff>(count * sizeof(T));
  if (input.tellg() != expected) fail("unexpected byte count for " + path);
  input.seekg(0);
  std::vector<T> result(count);
  input.read(reinterpret_cast<char*>(result.data()), expected);
  if (!input) fail("short read from " + path);
  return result;
}

template <typename T>
class DeviceBuffer {
 public:
  explicit DeviceBuffer(std::size_t count) : count_(count) {
    pointer_=static_cast<T*>(::device_storage.acquire(bytes()));
  }
  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;
  ~DeviceBuffer() {
    if (pointer_) ::device_storage.release(pointer_,bytes());
  }
  T* get() { return pointer_; }
  const T* get() const { return pointer_; }
  std::size_t bytes() const { return count_ * sizeof(T); }

 private:
  T* pointer_ = nullptr;
  std::size_t count_ = 0;
};

struct Options {
  std::string fixture;
  bool validate = true;
  double sigma_pi = std::numeric_limits<double>::quiet_NaN();
  int repeats = 3;
  int threads = 256;
};

double parse_double(const char* value, const char* name) {
  const double result = std::stod(value);
  if (!std::isfinite(result)) fail(std::string(name) + " must be finite");
  return result;
}

int parse_positive(const char* value, const char* name) {
  const long result = std::stol(value);
  if (result <= 0 || result > std::numeric_limits<int>::max()) {
    fail(std::string(name) + " must be a positive integer");
  }
  return static_cast<int>(result);
}

Options parse_options(int argc, char** argv) {
  Options result;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--production") {
      result.validate = false; result.repeats = 1;
    } else if (argument == "--fixture" && index + 1 < argc) {
      result.fixture = argv[++index];
    } else if (argument == "--sigma-pi" && index + 1 < argc) {
      result.sigma_pi = parse_double(argv[++index], "--sigma-pi");
    } else if (argument == "--repeats" && index + 1 < argc) {
      result.repeats = parse_positive(argv[++index], "--repeats");
    } else if (argument == "--threads" && index + 1 < argc) {
      result.threads = parse_positive(argv[++index], "--threads");
    } else {
      fail("unknown or incomplete argument: " + argument);
    }
  }
  if (!std::isfinite(result.sigma_pi)) fail("--sigma-pi is required");
  if (result.fixture.empty()) fail("--fixture is required");
  if (result.threads > 1024) fail("--threads cannot exceed 1024");
  return result;
}

int blocks_for(std::size_t count, int threads) {
  return static_cast<int>((count + threads - 1) / threads);
}

__device__ __forceinline__ int full_exogenous(
    int lag_a, int current_a, int lag_z, int current_z) {
  return 125 * lag_a + 25 * current_a + 5 * lag_z + current_z;
}

__global__ void micro_path_kernel(
    const std::int32_t* __restrict__ aggregate_path,
    const std::int32_t* __restrict__ idiosyncratic_path,
    const double* __restrict__ profit_shocks,
    const std::int32_t* __restrict__ policy,
    const std::int32_t* __restrict__ default_policy,
    const double* __restrict__ transition,
    const double* __restrict__ profit,
    const double* __restrict__ price,
    const double* __restrict__ capital_grid,
    const double* __restrict__ debt_grid,
    double sigma_pi,
    std::int32_t* __restrict__ endogenous_path,
    double* __restrict__ capital_path,
    double* __restrict__ debt_path,
    double* __restrict__ investment_path,
    double* __restrict__ measured_profit_path,
    double* __restrict__ forecast_profit_path,
    double* __restrict__ spread_path) {
  const int firm = blockIdx.x * blockDim.x + threadIdx.x;
  if (firm >= kFirms) return;
  const std::size_t base = static_cast<std::size_t>(firm) * kTimes;
  endogenous_path[base] = kEndogenous / 2 - 1;
  endogenous_path[base + 1] = kEndogenous / 2 - 1;

  for (int time = 1; time <= kMicroLast; ++time) {
    const int endogenous = endogenous_path[base + time];
    const int lag_a = aggregate_path[time - 1] - 1;
    const int current_a = aggregate_path[time] - 1;
    const int lag_z = idiosyncratic_path[base + time - 1] - 1;
    const int current_z = idiosyncratic_path[base + time] - 1;
    const int exogenous =
        full_exogenous(lag_a, current_a, lag_z, current_z);
    const std::size_t state =
        static_cast<std::size_t>(endogenous) * kExogenous + exogenous;
    const int action = default_policy[state] == 0 ? policy[state] - 1 : 0;
    endogenous_path[base + time + 1] = action;

    const double capital = capital_grid[endogenous / kDebt];
    const double debt = debt_grid[endogenous % kDebt];
    const double next_capital = capital_grid[action / kDebt];
    capital_path[base + time] = capital;
    debt_path[base + time] = debt;
    investment_path[base + time] = fma(-(1.0-kDelta),capital,next_capital);
    measured_profit_path[base + time] =
        fma(sigma_pi,profit_shocks[base+time],profit[state]);

    double forecast = 0.0;
    for (int next_a = 0; next_a < kA; ++next_a) {
      for (int next_z = 0; next_z < kZ; ++next_z) {
        const int next =
            full_exogenous(current_a, next_a, current_z, next_z);
        forecast = fma(profit[static_cast<std::size_t>(action) * kExogenous + next],
            transition[static_cast<std::size_t>(next) * kExogenous + exogenous],forecast);
      }
    }
    forecast_profit_path[base + time + 1] = forecast;
    const double q = price[static_cast<std::size_t>(action) * kExogenous +
                           exogenous];
    double spread = q > 0.0 ? 1.0 / q - (1.0 + kRate) : 0.0;
    if (spread > 0.5) spread = 0.5;
    spread_path[base + time] = spread;
  }
}

__global__ void initialize_distribution(double* distribution) {
  const int state = blockIdx.x * blockDim.x + threadIdx.x;
  if (state < kActiveExogenous * kEndogenous) {
    distribution[state] = 1.0 /
        static_cast<double>(kActiveExogenous * kEndogenous);
  }
}

__global__ void macro_aggregate_and_push_kernel(
    int lag_a, int current_a,
    const double* __restrict__ distribution,
    const std::int32_t* __restrict__ policy,
    const std::int32_t* __restrict__ default_policy,
    const double* __restrict__ transition_z,
    const double* __restrict__ output_net_cache,
    const double* __restrict__ price,
    double* __restrict__ next_distribution,
    double* __restrict__ aggregate,
    int time) {
  const int active_state = blockIdx.x * blockDim.x + threadIdx.x;
  if (active_state >= kActiveExogenous * kEndogenous) return;
  const int row = active_state % kActiveExogenous;
  const int endogenous = active_state / kActiveExogenous;
  const int lag_z = row / kZ;
  const int current_z = row % kZ;
  const int exogenous =
      full_exogenous(lag_a, current_a, lag_z, current_z);
  const std::size_t policy_state =
      static_cast<std::size_t>(endogenous) * kExogenous + exogenous;
  const double mass = distribution[active_state];
  if (!(mass > 0.0)) return;
  const int is_default = default_policy[policy_state];
  const int action = is_default == 0 ? policy[policy_state] - 1 : 0;

  const int capital_index = endogenous / kDebt;
  const double output = is_default == 0 ?
      2.0 * output_net_cache[current_a + 5 * current_z +
                             25 * capital_index] : 0.0;
  const double q =
      price[static_cast<std::size_t>(action) * kExogenous + exogenous];
  atomicAdd(&aggregate[3 * time], output * mass);
  atomicAdd(&aggregate[3 * time + 1],
            static_cast<double>(is_default) * mass);
  atomicAdd(&aggregate[3 * time + 2],
            (1.0 / q - (1.0 + kRate)) * mass);
  for (int next_z = 0; next_z < kZ; ++next_z) {
    const int next_row = current_z * kZ + next_z;
    const int destination = action * kActiveExogenous + next_row;
    const double probability = transition_z[current_z + kZ * next_z];
    atomicAdd(&next_distribution[destination], mass * probability);
  }
}

__global__ void normalize_distribution(
    double* distribution, const double* total) {
  const int state = blockIdx.x * blockDim.x + threadIdx.x;
  if (state < kActiveExogenous * kEndogenous) {
    distribution[state] = distribution[state] / total[0];
  }
}

__global__ void sum_distribution_kernel(
    const double* __restrict__ distribution, double* total) {
  const int state = blockIdx.x * blockDim.x + threadIdx.x;
  if (state < kActiveExogenous * kEndogenous) {
    atomicAdd(total, distribution[state]);
  }
}

void fill_statistics(
    int time, int firm,
    const std::vector<double>& capital,
    const std::vector<double>& debt,
    const std::vector<double>& investment,
    const std::vector<double>& measured_profit,
    const std::vector<double>& forecast_profit,
    const std::vector<double>& spread,
    double* result) {
  const std::size_t base = static_cast<std::size_t>(firm) * kTimes;
  auto at = [base](const std::vector<double>& values, int index) {
    return values[base + index];
  };
  result[0] =
      (at(measured_profit, time + 1) - at(forecast_profit, time + 1)) /
          at(capital, time) -
      (at(measured_profit, time) - at(forecast_profit, time)) /
          at(capital, time - 1);
  result[1] =
      (at(measured_profit, time + 2) - at(forecast_profit, time + 2)) /
          at(capital, time + 1) -
      (at(measured_profit, time) - at(forecast_profit, time)) /
          at(capital, time - 1);
  result[2] = at(measured_profit, time - 1) / at(capital, time - 1) -
      at(measured_profit, time - 2) / at(capital, time - 2);
  result[3] = at(measured_profit, time - 1) / at(capital, time - 1) -
      at(measured_profit, time - 3) / at(capital, time - 3);
  result[4] = at(investment, time - 1) / at(capital, time - 1) -
      at(investment, time - 2) / at(capital, time - 2);
  result[5] = at(investment, time - 1) / at(capital, time - 1) -
      at(investment, time - 3) / at(capital, time - 3);
  result[6] = at(debt, time) / at(capital, time - 1) -
      at(debt, time - 1) / at(capital, time - 2);
  result[7] = at(debt, time) / at(capital, time - 1) -
      at(debt, time - 2) / at(capital, time - 3);
  result[8] =
      (at(measured_profit, time + 1) - at(forecast_profit, time + 1)) /
      at(capital, time);
  result[9] = at(measured_profit, time) / at(capital, time);
  result[10] = at(investment, time) / at(capital, time);
  result[11] = at(debt, time + 1) / at(capital, time);
  result[12] = at(spread, time);
}

bool valid_micro_observation(
    int time, int firm, const std::vector<double>& capital) {
  const std::size_t base = static_cast<std::size_t>(firm) * kTimes;
  for (int offset = -3; offset <= 2; ++offset) {
    if (capital[base + time + offset] <= 1.0e-4) return false;
  }
  return true;
}

std::vector<double> compute_moments(
    const std::vector<double>& capital,
    const std::vector<double>& debt,
    const std::vector<double>& investment,
    const std::vector<double>& measured_profit,
    const std::vector<double>& forecast_profit,
    const std::vector<double>& spread,
    const std::vector<double>& aggregate_output,
    const std::vector<double>& aggregate_default,
    const std::vector<double>& aggregate_spread) {
  std::vector<double> firm_mean(kFirms * kStats, 0.0);
  std::vector<int> firm_count(kFirms, 0);
  double values[kStats];
  for (int firm = 0; firm < kFirms; ++firm) {
    for (int time = 50; time <= 55; ++time) {
      if (!valid_micro_observation(time, firm, capital)) continue;
      ++firm_count[firm];
      for (int stat = 0; stat < 8; ++stat) values[stat] = 0.0;
      fill_statistics(time, firm, capital, debt, investment, measured_profit,
                      forecast_profit, spread, values);
      for (int stat = 0; stat < 8; ++stat) values[stat] = 0.0;
      for (int stat = 0; stat < kStats; ++stat) {
        firm_mean[firm * kStats + stat] += values[stat];
      }
    }
    // The released code forms 0/0 for firms with no retained observation, but
    // those firms never enter either later pass.  Leave their unused row zero.
    if (firm_count[firm] > 0) {
      for (int stat = 0; stat < kStats; ++stat) {
        firm_mean[firm * kStats + stat] /=
            static_cast<double>(firm_count[firm]);
      }
    }
  }

  double time_mean[6][kStats]{};
  int time_count[6]{};
  for (int firm = 0; firm < kFirms; ++firm) {
    for (int time = 50; time <= 55; ++time) {
      if (!valid_micro_observation(time, firm, capital)) continue;
      ++time_count[time - 50];
      fill_statistics(time, firm, capital, debt, investment, measured_profit,
                      forecast_profit, spread, values);
      for (int stat = 0; stat < kStats; ++stat) {
        time_mean[time - 50][stat] +=
            values[stat] - firm_mean[firm * kStats + stat];
      }
    }
  }
  for (int time = 0; time < 6; ++time) {
    if (time_count[time] == 0) fail("time has no valid micro observations");
    for (int stat = 0; stat < kStats; ++stat) {
      time_mean[time][stat] /= static_cast<double>(time_count[time]);
    }
  }

  double mean[kStats]{};
  double covariance[kStats][kStats]{};
  std::uint64_t count = 0;
  for (int firm = 0; firm < kFirms; ++firm) {
    for (int time = 50; time <= 55; ++time) {
      if (!valid_micro_observation(time, firm, capital)) continue;
      ++count;
      fill_statistics(time, firm, capital, debt, investment, measured_profit,
                      forecast_profit, spread, values);
      for (int stat = 0; stat < kStats; ++stat) {
        values[stat] = (values[stat]-time_mean[time-50][stat])-firm_mean[firm*kStats+stat];
        mean[stat] += values[stat];
      }
      for (int first = 0; first < kStats; ++first) {
        for (int second = 0; second < kStats; ++second) {
          covariance[first][second] = std::fma(values[first],values[second],covariance[first][second]);
        }
      }
    }
  }
  for(int first=0;first<kStats;++first)mean[first]/=static_cast<double>(count);
  for (int first = 0; first < kStats; ++first) {
    for (int second = 0; second < kStats; ++second) {
      covariance[first][second] =
          std::fma(-mean[first],mean[second],covariance[first][second]/static_cast<double>(count));
    }
  }

  std::vector<double> moments(kMoments, 0.0);
  moments[0] = 0.25 * (covariance[0][2] + covariance[0][3] +
                       covariance[1][2] + covariance[1][3]);
  moments[1] = 0.25 * (covariance[0][4] + covariance[0][5] +
                       covariance[1][4] + covariance[1][5]);
  moments[2] = 0.25 * (covariance[0][6] + covariance[0][7] +
                       covariance[1][6] + covariance[1][7]);
  moments[3] = covariance[8][8];
  moments[4] = covariance[9][9];
  moments[5] = covariance[9][10];
  moments[6] = covariance[9][11];
  moments[7] = covariance[9][12];
  moments[8] = covariance[10][10];
  moments[9] = covariance[10][11];
  moments[10] = covariance[10][12];
  moments[11] = covariance[11][11];
  moments[12] = covariance[11][12];
  moments[13] = covariance[12][12];

  double macro_mean[3]{};
  double macro_covariance[3][3]{};
  constexpr int kMacroCount = 251;
  for (int time = 49; time <= 299; ++time) {
    const double macro[3] = {
        aggregate_spread[time], aggregate_default[time],
        std::log(aggregate_output[time]) -
            std::log(aggregate_output[time - 1])};
    for (int first = 0; first < 3; ++first) {
      macro_mean[first] += macro[first];
      for (int second = 0; second < 3; ++second) {
        macro_covariance[first][second]=std::fma(macro[first],macro[second],macro_covariance[first][second]);
      }
    }
  }
  for(int first=0;first<3;++first)macro_mean[first]/=static_cast<double>(kMacroCount);
  for (int first = 0; first < 3; ++first) {
    for (int second = 0; second < 3; ++second) {
      macro_covariance[first][second] =
          std::fma(-macro_mean[first],macro_mean[second],macro_covariance[first][second]/static_cast<double>(kMacroCount));
    }
  }
  moments[14] = macro_mean[0];
  moments[15] = macro_mean[1];
  moments[16] = macro_covariance[2][2];
  return moments;
}

template <typename T>
std::uint64_t bit_mismatches(
    const std::vector<T>& left, const std::vector<T>& right) {
  if (left.size() != right.size()) fail("comparison size mismatch");
  std::uint64_t count = 0;
  for (std::size_t index = 0; index < left.size(); ++index) {
    if (std::memcmp(&left[index], &right[index], sizeof(T)) != 0) ++count;
  }
  return count;
}

double maximum_error(
    const std::vector<double>& left, const std::vector<double>& right) {
  if (left.size() != right.size()) fail("comparison size mismatch");
  double result = 0.0;
  for (std::size_t index = 0; index < left.size(); ++index) {
    result = std::max(result, std::abs(left[index] - right[index]));
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    const std::size_t states =
        static_cast<std::size_t>(kExogenous) * kEndogenous;
    const std::size_t paths = static_cast<std::size_t>(kTimes) * kFirms;
    const std::size_t active_states =
        static_cast<std::size_t>(kActiveExogenous) * kEndogenous;
    auto aggregate_path = read_exact<std::int32_t>(
        options.fixture + "/RCC_Asimpos.bin", kTimes);
    auto idiosyncratic_path = read_exact<std::int32_t>(
        options.fixture + "/RCC_zfirmpos.bin", paths);
    auto profit_shocks = read_exact<double>(
        options.fixture + "/RCC_pifirmshocks.bin", paths);
    auto policy = solved_inputs ? solved_inputs->policy : read_exact<std::int32_t>(options.fixture+"/RCC_kbpolind.bin",states);
    auto default_policy = solved_inputs ? solved_inputs->defaults : read_exact<std::int32_t>(options.fixture+"/RCC_dfpolind.bin",states);
    auto transition = read_exact<double>(
        options.fixture + "/RCC_pr_mat_exog_theta.bin",
        static_cast<std::size_t>(kExogenous) * kExogenous);
    auto transition_z = read_exact<double>(
        options.fixture + "/RCC_pr_mat_z.bin", kZ * kZ);
    auto profit = solved_inputs ? solved_inputs->profit : read_exact<double>(options.fixture+"/RCC_pimat.bin",states);
    auto price = solved_inputs ? solved_inputs->price : read_exact<double>(options.fixture+"/RCC_qmat.bin",states);
    auto capital_grid = read_exact<double>(options.fixture + "/RCC_k0.bin", kCapital);
    auto debt_grid = read_exact<double>(options.fixture + "/RCC_b0.bin", kDebt);
    auto output_net_cache = read_exact<double>(
        options.fixture + "/RCC_ynetcache.bin", kA * kZ * kCapital);
    std::vector<std::int32_t> target_endogenous;
    std::vector<double> target_capital,target_debt,target_investment,target_measured_profit,target_forecast_profit,target_spread,target_output,target_default,target_aggregate_spread,target_moments;
    if (options.validate) {
    target_endogenous = read_exact<std::int32_t>(
        options.fixture + "/RCC_sim_endogfirmpos.bin", paths);
    for (auto& state : target_endogenous) {
      if (state <= 0) fail("CPU endogenous simulation target is not populated");
      --state;
    }
    target_capital = read_exact<double>(
        options.fixture + "/RCC_sim_kfirm.bin", paths);
    target_debt = read_exact<double>(
        options.fixture + "/RCC_sim_bfirm.bin", paths);
    target_investment = read_exact<double>(
        options.fixture + "/RCC_sim_ifirm.bin", paths);
    target_measured_profit = read_exact<double>(
        options.fixture + "/RCC_sim_pifirm.bin", paths);
    target_forecast_profit = read_exact<double>(
        options.fixture + "/RCC_sim_pifcstfirm.bin", paths);
    target_spread = read_exact<double>(
        options.fixture + "/RCC_sim_sprfirm.bin", paths);
    target_output = read_exact<double>(
        options.fixture + "/RCC_sim_Ysim.bin", kTimes);
    target_default = read_exact<double>(
        options.fixture + "/RCC_sim_dfsim.bin", kTimes);
    target_aggregate_spread = read_exact<double>(
        options.fixture + "/RCC_sim_spreadsim.bin", kTimes);
    target_moments = read_exact<double>(
        options.fixture + "/RCC_GMMmodelmom.bin", kMoments);
    for (int firm = 0; firm < kFirms; ++firm) {
      const std::size_t base = static_cast<std::size_t>(firm) * kTimes;
      for (int time = kMicroLast + 1; time < kTimes; ++time) {
        target_capital[base + time] = 0.0;
        target_debt[base + time] = 0.0;
        target_investment[base + time] = 0.0;
        target_measured_profit[base + time] = 0.0;
        target_spread[base + time] = 0.0;
      }
      for (int time = kMicroLast + 2; time < kTimes; ++time) {
        target_endogenous[base + time] = 0;
        target_forecast_profit[base + time] = 0.0;
      }
    }

    }

    DeviceBuffer<std::int32_t> d_aggregate_path(kTimes);
    DeviceBuffer<std::int32_t> d_idiosyncratic_path(paths);
    DeviceBuffer<double> d_profit_shocks(paths);
    DeviceBuffer<std::int32_t> d_policy(states);
    DeviceBuffer<std::int32_t> d_default_policy(states);
    DeviceBuffer<double> d_transition(
        static_cast<std::size_t>(kExogenous) * kExogenous);
    DeviceBuffer<double> d_transition_z(kZ * kZ);
    DeviceBuffer<double> d_profit(states);
    DeviceBuffer<double> d_price(states);
    DeviceBuffer<double> d_capital_grid(kCapital);
    DeviceBuffer<double> d_debt_grid(kDebt);
    DeviceBuffer<double> d_output_net_cache(kA * kZ * kCapital);
    DeviceBuffer<std::int32_t> d_endogenous_path(paths);
    DeviceBuffer<double> d_capital_path(paths);
    DeviceBuffer<double> d_debt_path(paths);
    DeviceBuffer<double> d_investment_path(paths);
    DeviceBuffer<double> d_measured_profit_path(paths);
    DeviceBuffer<double> d_forecast_profit_path(paths);
    DeviceBuffer<double> d_spread_path(paths);
    DeviceBuffer<double> d_distribution(active_states);
    DeviceBuffer<double> d_next_distribution(active_states);
    DeviceBuffer<double> d_total(1);
    DeviceBuffer<double> d_aggregate(3 * kTimes);

#define COPY_TO_DEVICE(target, source) \
    CUDA_CHECK(cudaMemcpy((target).get(), (source).data(), (target).bytes(), \
                          cudaMemcpyHostToDevice))
    COPY_TO_DEVICE(d_aggregate_path, aggregate_path);
    COPY_TO_DEVICE(d_idiosyncratic_path, idiosyncratic_path);
    COPY_TO_DEVICE(d_profit_shocks, profit_shocks);
    COPY_TO_DEVICE(d_policy, policy);
    COPY_TO_DEVICE(d_default_policy, default_policy);
    COPY_TO_DEVICE(d_transition, transition);
    COPY_TO_DEVICE(d_transition_z, transition_z);
    COPY_TO_DEVICE(d_profit, profit);
    COPY_TO_DEVICE(d_price, price);
    COPY_TO_DEVICE(d_capital_grid, capital_grid);
    COPY_TO_DEVICE(d_debt_grid, debt_grid);
    COPY_TO_DEVICE(d_output_net_cache, output_net_cache);
#undef COPY_TO_DEVICE

    auto clear_paths = [&]() {
      CUDA_CHECK(cudaMemset(d_endogenous_path.get(), 0, d_endogenous_path.bytes()));
      CUDA_CHECK(cudaMemset(d_capital_path.get(), 0, d_capital_path.bytes()));
      CUDA_CHECK(cudaMemset(d_debt_path.get(), 0, d_debt_path.bytes()));
      CUDA_CHECK(cudaMemset(d_investment_path.get(), 0, d_investment_path.bytes()));
      CUDA_CHECK(cudaMemset(d_measured_profit_path.get(), 0,
                            d_measured_profit_path.bytes()));
      CUDA_CHECK(cudaMemset(d_forecast_profit_path.get(), 0,
                            d_forecast_profit_path.bytes()));
      CUDA_CHECK(cudaMemset(d_spread_path.get(), 0, d_spread_path.bytes()));
    };
    auto launch_micro = [&]() {
      micro_path_kernel<<<blocks_for(kFirms, options.threads), options.threads>>>(
          d_aggregate_path.get(), d_idiosyncratic_path.get(),
          d_profit_shocks.get(), d_policy.get(), d_default_policy.get(),
          d_transition.get(), d_profit.get(), d_price.get(),
          d_capital_grid.get(), d_debt_grid.get(), options.sigma_pi,
          d_endogenous_path.get(), d_capital_path.get(), d_debt_path.get(),
          d_investment_path.get(), d_measured_profit_path.get(),
          d_forecast_profit_path.get(), d_spread_path.get());
    };

    if (options.validate) {
    clear_paths();
    launch_micro();
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(cudaStreamPerThread));
    }
    cudaEvent_t start = nullptr;
    cudaEvent_t stop = nullptr;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    CUDA_CHECK(cudaEventRecord(start));
    for (int repeat = 0; repeat < options.repeats; ++repeat) {
      clear_paths();
      launch_micro();
    }
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float micro_ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&micro_ms, start, stop));

    std::vector<double> aggregate_output(kTimes, 0.0);
    std::vector<double> aggregate_default(kTimes, 0.0);
    std::vector<double> aggregate_spread(kTimes, 0.0);
    const auto macro_start=std::chrono::steady_clock::now();
    if (reuse_macro) {
      if (cached_output.size()!=kTimes) fail("Missing conditional-profile macro cache");
      aggregate_output=cached_output; aggregate_default=cached_default; aggregate_spread=cached_spread;
    } else {
    std::vector<double> distribution(22500,1.0/22500.0),next_distribution(22500);
    for(int time=1;time<=kMacroLast;++time) {
      int lag_a=aggregate_path[time-1]-1,current_a=aggregate_path[time]-1;
      std::fill(next_distribution.begin(),next_distribution.end(),0.0);
      // Preserve the source accumulation order independently for each destination.
      for(int lz=0;lz<5;++lz)for(int z=0;z<5;++z)for(int j=0;j<900;++j) {
        int e=125*lag_a+25*current_a+5*lz+z,st=j*625+e;
        double mass=distribution[j*25+5*lz+z];if(!(mass>0))continue;
        int a=default_policy[st]?0:policy[st]-1;
        for(int nz=0;nz<5;++nz){int dest=a*25+5*z+nz;next_distribution[dest]=std::fma(mass,transition_z[z+5*nz],next_distribution[dest]);}
      }
      double y=0,d=0,sp=0;
      for(int z=0;z<5;++z)for(int lz=0;lz<5;++lz)for(int j=0;j<900;++j) {
        int e=125*lag_a+25*current_a+5*lz+z,st=j*625+e;
        double mass=distribution[j*25+5*lz+z];if(!(mass>0))continue;
        int df=default_policy[st],a=df?0:policy[st]-1;
        double output=df?0:2.0*output_net_cache[current_a+5*z+25*(j/30)];
        y=std::fma(output,mass,y);d=std::fma(double(df),mass,d);
        sp=std::fma(1.0/price[a*625+e]-(1.0+kRate),mass,sp);
      }
      aggregate_output[time]=y;aggregate_default[time]=d;aggregate_spread[time]=sp;
      double total=0;for(double v:next_distribution)total+=v;
      for(double&v:next_distribution)v/=total;distribution.swap(next_distribution);
    }
      cached_output=aggregate_output; cached_default=aggregate_default; cached_spread=aggregate_spread;
    }
    const auto macro_stop=std::chrono::steady_clock::now();
    const double macro_ms = std::chrono::duration<double, std::milli>(
        macro_stop - macro_start).count();

    std::vector<std::int32_t> endogenous_path(paths);
    std::vector<double> capital_path(paths), debt_path(paths),
        investment_path(paths), measured_profit_path(paths),
        forecast_profit_path(paths), spread_path(paths);
#define COPY_TO_HOST(target, source) \
    CUDA_CHECK(cudaMemcpy((target).data(), (source).get(), (source).bytes(), \
                          cudaMemcpyDeviceToHost))
    const auto copy_start = std::chrono::steady_clock::now();
    COPY_TO_HOST(endogenous_path, d_endogenous_path);
    COPY_TO_HOST(capital_path, d_capital_path);
    COPY_TO_HOST(debt_path, d_debt_path);
    COPY_TO_HOST(investment_path, d_investment_path);
    COPY_TO_HOST(measured_profit_path, d_measured_profit_path);
    COPY_TO_HOST(forecast_profit_path, d_forecast_profit_path);
    COPY_TO_HOST(spread_path, d_spread_path);
#undef COPY_TO_HOST
    const auto copy_stop = std::chrono::steady_clock::now();
    const auto moment_start = std::chrono::steady_clock::now();
    auto moments = compute_moments(
        capital_path, debt_path, investment_path, measured_profit_path,
        forecast_profit_path, spread_path, aggregate_output,
        aggregate_default, aggregate_spread);
    const auto moment_stop = std::chrono::steady_clock::now();

    if(const char* dir=dump_directory.empty()?std::getenv("RCC_GPU_DUMP"):dump_directory.c_str()) {
      auto save=[&](const char* name,const std::vector<double>&v){std::ofstream f(std::string(dir)+"/"+name,std::ios::binary);f.write((char*)v.data(),v.size()*8);};
      save("moments.bin",moments);save("output.bin",aggregate_output);save("default.bin",aggregate_default);save("spread.bin",aggregate_spread);
      save("investment.bin",investment_path);save("profit.bin",measured_profit_path);save("forecast.bin",forecast_profit_path);
    }
    last_moments=moments;
    last_macro_output=aggregate_output;
    if (!options.validate) {
      report()<<std::setprecision(17)<<"{\"family\":\"simulation\",\"reference_checked\":false,\"sigma_pi\":"<<options.sigma_pi
              <<",\"macro_wall_ms\":"<<macro_ms<<",\"moments\":[";
      for (size_t i=0;i<moments.size();++i) { if(i)report()<<','; if(std::isfinite(moments[i]))report()<<moments[i];else report()<<"null"; }
      report()<<"]}\n";
      CUDA_CHECK(cudaEventDestroy(start)); CUDA_CHECK(cudaEventDestroy(stop));
      return 0;
    }
    double maximum_moment_error = maximum_error(moments, target_moments);
    int device = 0;
    cudaDeviceProp properties{};
    CUDA_CHECK(cudaGetDevice(&device));
    CUDA_CHECK(cudaGetDeviceProperties(&properties, device));
    const double copy_ms = std::chrono::duration<double, std::milli>(
        copy_stop - copy_start).count();
    const double moment_ms = std::chrono::duration<double, std::milli>(
        moment_stop - moment_start).count();
    report() << std::setprecision(17)
              << "{\"experiment\":\"de_hybrid_simulation_exact_gate_v2\""
              << ",\"scientific_role\":\"fixed_parameter_simulation_parity_gate\""
              << ",\"device\":\"" << properties.name << "\""
              << ",\"sigma_pi\":" << options.sigma_pi
              << ",\"mean_micro_kernel_ms\":"
              << static_cast<double>(micro_ms) / options.repeats
              << ",\"macro_wall_ms\":" << macro_ms
              << ",\"device_to_host_ms\":" << copy_ms
              << ",\"host_moment_ms\":" << moment_ms
              << ",\"endogenous_path_bit_mismatches\":"
              << bit_mismatches(endogenous_path, target_endogenous)
              << ",\"capital_path_bit_mismatches\":"
              << bit_mismatches(capital_path, target_capital)
              << ",\"debt_path_bit_mismatches\":"
              << bit_mismatches(debt_path, target_debt)
              << ",\"investment_path_bit_mismatches\":"
              << bit_mismatches(investment_path, target_investment)
              << ",\"measured_profit_bit_mismatches\":"
              << bit_mismatches(measured_profit_path, target_measured_profit)
              << ",\"forecast_profit_bit_mismatches\":"
              << bit_mismatches(forecast_profit_path, target_forecast_profit)
              << ",\"spread_path_bit_mismatches\":"
              << bit_mismatches(spread_path, target_spread)
              << ",\"maximum_capital_error\":"
              << maximum_error(capital_path, target_capital)
              << ",\"maximum_debt_error\":"
              << maximum_error(debt_path, target_debt)
              << ",\"maximum_investment_error\":"
              << maximum_error(investment_path, target_investment)
              << ",\"maximum_measured_profit_error\":"
              << maximum_error(measured_profit_path, target_measured_profit)
              << ",\"maximum_forecast_profit_error\":"
              << maximum_error(forecast_profit_path, target_forecast_profit)
              << ",\"maximum_spread_path_error\":"
              << maximum_error(spread_path, target_spread)
              << ",\"maximum_output_error\":"
              << maximum_error(aggregate_output, target_output)
              << ",\"maximum_default_error\":"
              << maximum_error(aggregate_default, target_default)
              << ",\"maximum_aggregate_spread_error\":"
              << maximum_error(aggregate_spread, target_aggregate_spread)
              << ",\"moment_bit_mismatches\":" << bit_mismatches(moments,target_moments)
              << ",\"output_bit_mismatches\":" << bit_mismatches(aggregate_output,target_output)
              << ",\"macro_default_bit_mismatches\":" << bit_mismatches(aggregate_default,target_default)
              << ",\"macro_spread_bit_mismatches\":" << bit_mismatches(aggregate_spread,target_aggregate_spread)
              << ",\"maximum_moment_error\":" << maximum_moment_error
              << ",\"moments\":[";
    for (int index = 0; index < kMoments; ++index) {
      if (index) report() << ',';
      report() << moments[index];
    }
    report() << "]}\n";
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
    return (bit_mismatches(moments,target_moments) ||
      bit_mismatches(endogenous_path,target_endogenous) || bit_mismatches(capital_path,target_capital) ||
      bit_mismatches(debt_path,target_debt) || bit_mismatches(investment_path,target_investment) ||
      bit_mismatches(measured_profit_path,target_measured_profit) || bit_mismatches(forecast_profit_path,target_forecast_profit) ||
      bit_mismatches(spread_path,target_spread) || bit_mismatches(aggregate_output,target_output) ||
      bit_mismatches(aggregate_default,target_default) || bit_mismatches(aggregate_spread,target_aggregate_spread)) ? 2 : 0;
  } catch (const std::exception& error) {
    errors() << "de_simulation_snapshot: " << error.what() << '\n';
    return 1;
  }
}
