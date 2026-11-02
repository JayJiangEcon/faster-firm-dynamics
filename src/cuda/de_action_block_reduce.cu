// CUDA gate for one complete DE debt-envelope action-improvement snapshot.
//
// The input is a separately executed current-CPU successor-local fixture at
// a declared fixed-parameter point. The native batch driver includes these action kernels.

#include <cuda_runtime.h>
#include <cstdlib>

#include <algorithm>
#include <chrono>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int kExogenous = 625;
constexpr int kCapital = 30;
constexpr int kDebt = 30;
constexpr int kEndogenous = kCapital * kDebt;
constexpr int kSuccessors = 25;
constexpr double kDelta = 0x1.99999a0000000p-4;
constexpr double kInvestmentPrice = 1.0;
constexpr double kRate = 0x1.47ae140000000p-5;
constexpr double kBeta = 0x1.ec4ec4f368855p-1;
constexpr double kTax = 0x1.99999a0000000p-3;
constexpr double kDebtAdjustment = 0.0;
constexpr double kTieGuard = 1.0e-12;

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
  std::vector<T> output(count);
  input.read(reinterpret_cast<char*>(output.data()), expected);
  if (!input) fail("short read from " + path);
  return output;
}

template <typename T>
class DeviceBuffer {
 public:
  explicit DeviceBuffer(std::size_t count) : count_(count) {
    pointer_=static_cast<T*>(::device_storage.acquire(count*sizeof(T)));
  }
  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;
  ~DeviceBuffer() {
    if (pointer_) ::device_storage.release(pointer_,count_*sizeof(T));
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
  bool production = false;
  bool profile = false;
  double eta_k = std::numeric_limits<double>::quiet_NaN();
  double gamma = std::numeric_limits<double>::quiet_NaN();
  double sigma_pi = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();
  double eta_f = std::numeric_limits<double>::quiet_NaN();
  double eta_d = std::numeric_limits<double>::quiet_NaN();
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
      result.production = true;
    } else if (argument == "--sigma-profile") {
      result.profile = true;
    } else if (argument == "--fixture" && index + 1 < argc) {
      result.fixture = argv[++index];
    } else if (argument == "--eta-k" && index + 1 < argc) {
      result.eta_k = parse_double(argv[++index], "--eta-k");
    } else if (argument == "--gamma" && index + 1 < argc) {
      result.gamma = parse_double(argv[++index], "--gamma");
    } else if (argument == "--sigma-pi" && index + 1 < argc) {
      result.sigma_pi = parse_double(argv[++index], "--sigma-pi");
    } else if (argument == "--phi" && index + 1 < argc) {
      result.phi = parse_double(argv[++index], "--phi");
    } else if (argument == "--eta-f" && index + 1 < argc) {
      result.eta_f = parse_double(argv[++index], "--eta-f");
    } else if (argument == "--eta-d" && index + 1 < argc) {
      result.eta_d = parse_double(argv[++index], "--eta-d");
    } else if (argument == "--repeats" && index + 1 < argc) {
      result.repeats = parse_positive(argv[++index], "--repeats");
    } else if (argument == "--threads" && index + 1 < argc) {
      result.threads = parse_positive(argv[++index], "--threads");
    } else {
      fail("unknown or incomplete argument: " + argument);
    }
  }
  for (double value : {result.eta_k, result.gamma, result.sigma_pi, result.phi, result.eta_f, result.eta_d})
    if (!std::isfinite(value)) fail("All six economic parameter options are required");
  if (result.profile && !result.production) fail("--sigma-profile requires --production");
  if (result.fixture.empty()) fail("--fixture is required");
  if (result.threads > 1024) fail("--threads cannot exceed 1024");
  return result;
}

int blocks_for(std::size_t count, int threads) {
  return static_cast<int>((count + threads - 1) / threads);
}

__device__ __forceinline__ int next_exogenous(int exogenous, int successor) {
  const int current_a = (exogenous % 125) / 25;
  const int current_z = exogenous % 5;
  const int next_a = successor / 5;
  const int next_z = successor % 5;
  return 125 * current_a + 25 * next_a + 5 * current_z + next_z;
}

__global__ void compact_transition_kernel(
    const double* __restrict__ transition,
    double* __restrict__ compact_transition) {
  const int linear = blockIdx.x * blockDim.x + threadIdx.x;
  if (linear >= kExogenous * kSuccessors) return;
  const int exogenous = linear / kSuccessors;
  const int successor = linear % kSuccessors;
  const int next = next_exogenous(exogenous, successor);
  // Fortran column-major transition(exogenous,next).
  compact_transition[linear] =
      transition[static_cast<std::size_t>(next) * kExogenous + exogenous];
}

__global__ void expected_value_kernel(
    const double* __restrict__ compact_transition,
    const double* __restrict__ value_old,
    double* __restrict__ expected_value,
    double* __restrict__ default_value) {
  const int state = blockIdx.x * blockDim.x + threadIdx.x;
  if (state >= kExogenous * kEndogenous) return;
  const int exogenous = state % kExogenous;
  const int action = state / kExogenous;
  double result = 0.0;
  for (int successor = 0; successor < kSuccessors; ++successor) {
    const int next = next_exogenous(exogenous, successor);
    result = fma(kBeta * compact_transition[exogenous * kSuccessors + successor],
        value_old[static_cast<std::size_t>(action) * kExogenous + next], result);
  }
  expected_value[state] = result;
  if (action == 0) default_value[exogenous] = result;
}

__device__ __forceinline__ double adjustment_cost(
    double capital, double capital_prime, double eta_k) {
  const double investment = fma(-(1.0-kDelta),capital,capital_prime);
  if (capital > 0.0) {
    double result = eta_k / 2.0;
    const double ratio = investment / capital;
    result = result * (ratio * ratio);
    result = result * capital;
    return result;
  }
  return 0.0;
}

__global__ void exact_adjustment_kernel(const double* k,double eta,double* ac) {
 int i=blockIdx.x*blockDim.x+threadIdx.x;
 if(i<900) ac[i]=adjustment_cost(k[i%30],k[i/30],eta);
}
struct Candidate { int action, split; double positive, negative; };

__device__ __forceinline__ double exact_dividend(
    double output_factor, double capital_root, double capital, double capital_prime,
    double adjustment, double price, double debt_prime, double debt,
    double phi, double eta_f, double eta_d) {
  const double investment = fma(-(1.0-kDelta),capital,capital_prime);
  double operating = fma(output_factor,capital_root,-adjustment)-phi;
  double result = fma(1.0-kTax,operating,-investment);
  double tax_base = fma(kRate,debt,kDelta*capital);
  result = fma(kTax,tax_base,result);
  result = fma(price,debt_prime,result)-debt;
  if (debt_prime>debt) {
    const double change=debt_prime-debt;
    result = fma(-(kDebtAdjustment/2.0),change*change,result);
  }
  if (result<0.0) result=fma(-eta_d,-result,result-eta_f);
  return result;
}

__global__ void shortlist_kernel(
    const double* __restrict__ compact_transition,
    const double* __restrict__ value_old,
    const double* __restrict__ expected_value,
    const double* __restrict__ default_value,
    const double* __restrict__ price,
    const double* __restrict__ capital_grid,
    const double* __restrict__ debt_grid,
    const double* __restrict__ output_net_cache,
    const double* __restrict__ adjustment_cache,
    const double* __restrict__ action_capital_cache,
    double eta_k, double phi, double eta_f, double eta_d,
    std::int32_t* __restrict__ policy,
    std::int32_t* __restrict__ default_policy,
    double* __restrict__ value, Candidate* candidates, int* counts, double* bins) {
  const int group = blockIdx.x;
  const int group_count = kExogenous * kCapital;
  if (group >= group_count) return;
  const int exogenous = group % kExogenous;
  const int current_capital_index = group / kExogenous;
  const int current_a = (exogenous % 125) / 25;
  const int current_z = exogenous % 5;
  const double capital = capital_grid[current_capital_index];
  // Fortran column-major ynetcache(A,z,k).
  const int output_index=current_a + 5 * current_z + 25 * current_capital_index;
  const double output_factor=output_net_cache[output_index];
  const double capital_root=output_net_cache[output_index+750];
  const double ynet=output_factor*capital_root;
  const double dividend_slope = kTax * kRate - 1.0;

  __shared__ double positive_values[870],negative_values[870];
  __shared__ int splits[870];
  __shared__ double positive_bin[31],negative_bin[31];
  __shared__ int warp_counts[8],total;
  for(int action=30+threadIdx.x;action<900;action+=blockDim.x){
    const int next_capital_index = action / kDebt;
    const int next_debt_index = action % kDebt;
    const double debt_prime = debt_grid[next_debt_index];
    const double action_cost = action_capital_cache[
        current_capital_index + kCapital * next_capital_index];
    const double action_intercept =
        (1.0 - kTax) * (ynet - phi) +
        kInvestmentPrice * (1.0 - kDelta) * capital +
        kTax * kDelta * capital + action_cost +
        price[static_cast<std::size_t>(action) * kExogenous + exogenous] *
            debt_prime;
    int split = 0;
    int low = 0;
    int high = kDebt - 1;
    while (low <= high) {
      const int middle = (low + high) / 2;
      const double raw = action_intercept + dividend_slope * debt_grid[middle];
      if (raw >= 0.0) {
        split = middle + 1;
        low = middle + 1;
      } else {
        high = middle - 1;
      }
    }
    const double continuation =
        expected_value[static_cast<std::size_t>(action) * kExogenous +
                       exogenous];
    const double positive = action_intercept + continuation;
    const double negative =
        (1.0 + eta_d) * action_intercept - eta_f + continuation;

    positive_values[action-30]=positive;negative_values[action-30]=negative;splits[action-30]=split;
  }
  __syncthreads();
  // Warp argmax preserves the earliest action's bits on an exact tie.
  int reduce_lane=threadIdx.x%32,reduce_warp=threadIdx.x/32;
  for(int bin=reduce_warp;bin<31;bin+=blockDim.x/32){double p=-DBL_MAX,n=-DBL_MAX;int pi=2147483647,ni=2147483647;
   for(int i=reduce_lane;i<870;i+=32)if(splits[i]==bin){if(positive_values[i]>p){p=positive_values[i];pi=i;}if(negative_values[i]>n){n=negative_values[i];ni=i;}}
   for(int step=16;step>0;step/=2){
    double op=__shfl_down_sync(0xffffffff,p,step),on=__shfl_down_sync(0xffffffff,n,step);int oi=__shfl_down_sync(0xffffffff,pi,step),oj=__shfl_down_sync(0xffffffff,ni,step);
    if(reduce_lane+step<32){if(op>p || (op==p && oi<pi)){p=op;pi=oi;}if(on>n || (on==n && oj<ni)){n=on;ni=oj;}}
   }
   if(reduce_lane==0){positive_bin[bin]=p;negative_bin[bin]=n;bins[bin*18750+group]=p;bins[(31+bin)*18750+group]=n;}
  }
  if(threadIdx.x==0)total=0;
  __syncthreads();
  int lane=threadIdx.x%32,warp=threadIdx.x/32;unsigned prior=lane==0?0u:((1u<<lane)-1u);
  // Stable cooperative compaction avoids writing the complete action table.
  for(int base=0;base<870;base+=blockDim.x){int i=base+threadIdx.x;bool keep=false;
   if(i<870)keep=positive_values[i]>=positive_bin[splits[i]]-kTieGuard || negative_values[i]>=negative_bin[splits[i]]-kTieGuard;
   unsigned ballot=__ballot_sync(0xffffffff,keep);if(lane==0)warp_counts[warp]=__popc(ballot);
   __syncthreads();int offset=total;for(int w=0;w<warp;++w)offset+=warp_counts[w];
   if(keep)candidates[(offset+__popc(ballot&prior))*18750+group]={i+30,splits[i],positive_values[i],negative_values[i]};
   __syncthreads();if(threadIdx.x==0)for(int w=0;w<int(blockDim.x)/32;++w)total+=warp_counts[w];
   __syncthreads();
  }
  if(threadIdx.x==0)counts[group]=total;
}

__global__ void compact_action_kernel(
    const double* __restrict__ compact_transition,
    const double* __restrict__ value_old,
    const double* __restrict__ expected_value,
    const double* __restrict__ default_value,
    const double* __restrict__ price,
    const double* __restrict__ capital_grid,
    const double* __restrict__ debt_grid,
    const double* __restrict__ output_net_cache,
    const double* __restrict__ adjustment_cache,
    const double* __restrict__ action_capital_cache,
    double eta_k, double phi, double eta_f, double eta_d,
    std::int32_t* __restrict__ policy,
    std::int32_t* __restrict__ default_policy,
    double* __restrict__ value, Candidate* candidates, int* counts, double* bins) {
  const int state=blockIdx.x*blockDim.x+threadIdx.x;
  if(state>=kExogenous*kEndogenous)return;
  const int exogenous=state%625;
  const int current_capital_index=(state/625)/30;
  const int debt_index=(state/625)%30;
  const int group=current_capital_index*625+exogenous;
  const double capital=capital_grid[current_capital_index];
  const int output_index=(exogenous%125)/25+5*(exogenous%5)+25*current_capital_index;
  const double output_factor=output_net_cache[output_index],capital_root=output_net_cache[output_index+750];
  const double dividend_slope=kTax*kRate-1.0;
  const double* positive_bin=bins+group;
  const double* negative_bin=positive_bin+31*18750;
    const double debt = debt_grid[debt_index];
    const double debt_shift = dividend_slope * debt;
    double approximate_best = -DBL_MAX;
    for (int split = debt_index + 1; split <= kDebt; ++split) {
      const double candidate = positive_bin[split*18750] + debt_shift;
      if (candidate > approximate_best) approximate_best = candidate;
    }
    for (int split = 0; split <= debt_index; ++split) {
      const double candidate =
          negative_bin[split*18750] + (1.0 + eta_d) * debt_shift;
      if (candidate > approximate_best) approximate_best = candidate;
    }

    double best = -1000000.0;
    int best_action = 0;
    for(int ci=0;ci<counts[group];++ci) {
      const Candidate c=candidates[ci*18750+group];
      const int action=c.action,split=c.split,next_capital_index=action/30,next_debt_index=action%30;
      const double positive=c.positive,negative=c.negative;
      const double capital_prime=capital_grid[next_capital_index],debt_prime=debt_grid[next_debt_index];
      bool eligible = false;
      if (split >= debt_index + 1) {
        const double regime = positive_bin[split*18750] + debt_shift;
        eligible = positive >= positive_bin[split*18750] - kTieGuard &&
            regime >= approximate_best - kTieGuard;
      } else {
        const double regime =
            negative_bin[split*18750] + (1.0 + eta_d) * debt_shift;
        eligible = negative >= negative_bin[split*18750] - kTieGuard &&
            regime >= approximate_best - kTieGuard;
      }
      if (!eligible) continue;

      const double dividend = exact_dividend(
          output_factor, capital_root, capital, capital_prime,
          adjustment_cache[current_capital_index +
                           kCapital * next_capital_index],
          price[static_cast<std::size_t>(action) * kExogenous + exogenous],
          debt_prime, debt, phi, eta_f, eta_d);
      double candidate_value = dividend;
      for (int successor = 0; successor < kSuccessors; ++successor) {
        const int next = next_exogenous(exogenous, successor);
        candidate_value = fma(kBeta * compact_transition[exogenous * kSuccessors + successor],
            value_old[static_cast<std::size_t>(action) * kExogenous + next], candidate_value);
      }
      // Ascending action order plus strict greater-than reproduces MAXLOC's
      // first-index rule for finite exact ties.
      if (candidate_value > best) {
        best = candidate_value;
        best_action = action;
      }
    }

    policy[state] = best_action + 1;
    int is_default = best >= default_value[exogenous] ? 0 : 1;
    if (debt_index == 0) is_default = 0;
    default_policy[state] = is_default;
    value[state] = is_default == 0 ? best : default_value[exogenous];
}

__global__ void reference_score_kernel(const double* tr,const double* old,const double* price,
 const double* capital_grid,const double* debt_grid,const double* operands,const double* ac,
 const std::int32_t* policy,double phi,double ef,double ed,double* dividends,double* scores) {
 int state=blockIdx.x*blockDim.x+threadIdx.x;if(state>=625*900)return;
 int exog=state%625,j=state/625,k=j/30,b=j%30,action=policy[state]-1,kp=action/30,bp=action%30;
 int oi=(exog%125)/25+5*(exog%5)+25*k;
 double v=exact_dividend(operands[oi],operands[oi+750],capital_grid[k],capital_grid[kp],
 ac[k+30*kp],price[action*625+exog],debt_grid[bp],debt_grid[b],phi,ef,ed);
 dividends[state]=v;
 for(int n=0;n<25;++n)v=fma(kBeta*tr[exog*25+n],old[action*625+next_exogenous(exog,n)],v);
 scores[state]=v;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    const std::size_t states =
        static_cast<std::size_t>(kExogenous) * kEndogenous;
    const std::size_t transition_states =
        static_cast<std::size_t>(kExogenous) * kExogenous;
    auto transition = read_exact<double>(
        options.fixture + "/RCC_pr_mat_exog_theta.bin", transition_states);
    auto value_old = read_exact<double>(
        options.fixture + "/RCC_Vthetaold.bin", states);
    auto price = read_exact<double>(options.fixture + "/RCC_qmat.bin", states);
    auto capital_grid = read_exact<double>(options.fixture + "/RCC_k0.bin", 30);
    auto debt_grid = read_exact<double>(options.fixture + "/RCC_b0.bin", 30);
    auto agrid=read_exact<double>(options.fixture+"/RCC_A0.bin",5);
    auto zgrid=read_exact<double>(options.fixture+"/RCC_z0.bin",5);
    std::vector<double> output_net_cache(1500);
    for(int k=0;k<30;++k)for(int z=0;z<5;++z)for(int a=0;a<5;++a) {
      int i=a+5*z+25*k;double az=agrid[a]*zgrid[z];
      output_net_cache[i]=(az*az)*0.5;output_net_cache[i+750]=std::pow(capital_grid[k],0.5);
    }
    std::vector<double> adjustment_cache(900);
    auto action_capital_cache = read_exact<double>(
        options.fixture + "/RCC_actionkbase.bin", 30 * 30);
    auto target_policy = read_exact<std::int32_t>(
        options.fixture + "/RCC_kbpolind.bin", states);
    auto target_default = read_exact<std::int32_t>(
        options.fixture + "/RCC_dfpolind.bin", states);
    auto target_value =
        read_exact<double>(options.fixture + "/RCC_Vtheta.bin", states);

    DeviceBuffer<double> d_transition(transition_states);
    DeviceBuffer<double> d_compact_transition(kExogenous * kSuccessors);
    DeviceBuffer<double> d_value_old(states);
    DeviceBuffer<double> d_expected_value(states);
    DeviceBuffer<double> d_default_value(kExogenous);
    DeviceBuffer<double> d_price(states);
    DeviceBuffer<double> d_capital_grid(30);
    DeviceBuffer<double> d_debt_grid(30);
    DeviceBuffer<double> d_output_net_cache(2 * 5 * 5 * 30);
    DeviceBuffer<double> d_adjustment_cache(30 * 30);
    DeviceBuffer<double> d_action_capital_cache(30 * 30);
    DeviceBuffer<std::int32_t> d_policy(states);
    DeviceBuffer<std::int32_t> d_default(states);
    DeviceBuffer<double> d_value(states);
    DeviceBuffer<Candidate> d_candidates(625*30*870);
    DeviceBuffer<int> d_counts(625*30);
    DeviceBuffer<double> d_bins(625*30*62);

    CUDA_CHECK(cudaMemcpy(d_transition.get(), transition.data(),
                          d_transition.bytes(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_value_old.get(), value_old.data(),
                          d_value_old.bytes(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_price.get(), price.data(), d_price.bytes(),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_capital_grid.get(), capital_grid.data(),
                          d_capital_grid.bytes(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_debt_grid.get(), debt_grid.data(),
                          d_debt_grid.bytes(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_output_net_cache.get(), output_net_cache.data(),
                          d_output_net_cache.bytes(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_adjustment_cache.get(), adjustment_cache.data(),
                          d_adjustment_cache.bytes(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_action_capital_cache.get(),
                          action_capital_cache.data(),
                          d_action_capital_cache.bytes(),
                          cudaMemcpyHostToDevice));
    exact_adjustment_kernel<<<4,256>>>(d_capital_grid.get(),options.eta_k,d_adjustment_cache.get());
    compact_transition_kernel<<<
        blocks_for(kExogenous * kSuccessors, options.threads),
        options.threads>>>(d_transition.get(), d_compact_transition.get());
    auto launch_expected = [&]() {
      expected_value_kernel<<<blocks_for(states, options.threads),
                              options.threads>>>(
          d_compact_transition.get(), d_value_old.get(),
          d_expected_value.get(), d_default_value.get());
    };
    launch_expected();
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    auto launch_action = [&]() {
      shortlist_kernel<<<
          blocks_for(kExogenous * kCapital, options.threads),
          options.threads>>>(
          d_compact_transition.get(), d_value_old.get(),
          d_expected_value.get(), d_default_value.get(), d_price.get(),
          d_capital_grid.get(), d_debt_grid.get(), d_output_net_cache.get(),
          d_adjustment_cache.get(), d_action_capital_cache.get(),
          options.eta_k, options.phi,
          options.eta_f, options.eta_d, d_policy.get(), d_default.get(),
          d_value.get(),d_candidates.get(),d_counts.get(),d_bins.get());
      compact_action_kernel<<<
          blocks_for(states, options.threads),
          options.threads>>>(
          d_compact_transition.get(), d_value_old.get(),
          d_expected_value.get(), d_default_value.get(), d_price.get(),
          d_capital_grid.get(), d_debt_grid.get(), d_output_net_cache.get(),
          d_adjustment_cache.get(), d_action_capital_cache.get(),
          options.eta_k, options.phi,
          options.eta_f, options.eta_d, d_policy.get(), d_default.get(),
          d_value.get(),d_candidates.get(),d_counts.get(),d_bins.get());
    };
    launch_action();
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    cudaEvent_t start = nullptr;
    cudaEvent_t stop = nullptr;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    CUDA_CHECK(cudaEventRecord(start));
    for (int repeat = 0; repeat < options.repeats; ++repeat) launch_action();
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float elapsed_ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&elapsed_ms, start, stop));

    CUDA_CHECK(cudaEventRecord(start));
    for (int repeat = 0; repeat < options.repeats; ++repeat) {
      launch_expected();
      launch_action();
    }
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float combined_ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&combined_ms, start, stop));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));

    std::vector<std::int32_t> policy(states);
    std::vector<std::int32_t> default_policy(states);
    std::vector<double> value(states);
    CUDA_CHECK(cudaMemcpy(policy.data(), d_policy.get(), d_policy.bytes(),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(default_policy.data(), d_default.get(),
                          d_default.bytes(), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(value.data(), d_value.get(), d_value.bytes(),
                          cudaMemcpyDeviceToHost));

    if (const char* directory = std::getenv("RCC_GPU_DUMP")) {
      auto save = [&](const char* name, const void* data, std::size_t bytes) {
        std::ofstream out(std::string(directory)+"/"+name,std::ios::binary);
        if (!out) fail("cannot write GPU dump");
        out.write(reinterpret_cast<const char*>(data),bytes);
      };
      save("policy.bin",policy.data(),policy.size()*sizeof(std::int32_t));
      save("default.bin",default_policy.data(),default_policy.size()*sizeof(std::int32_t));
      save("value.bin",value.data(),value.size()*sizeof(double));
      std::vector<double> ev(states);
      CUDA_CHECK(cudaMemcpy(ev.data(),d_expected_value.get(),d_expected_value.bytes(),cudaMemcpyDeviceToHost));
      save("expected.bin",ev.data(),ev.size()*sizeof(double));
      DeviceBuffer<std::int32_t> ref_policy(states);
      DeviceBuffer<double> ref_div(states),ref_score(states);
      CUDA_CHECK(cudaMemcpy(ref_policy.get(),target_policy.data(),ref_policy.bytes(),cudaMemcpyHostToDevice));
      reference_score_kernel<<<blocks_for(states,options.threads),options.threads>>>(
        d_compact_transition.get(),d_value_old.get(),d_price.get(),d_capital_grid.get(),
        d_debt_grid.get(),d_output_net_cache.get(),d_adjustment_cache.get(),ref_policy.get(),
        options.phi,options.eta_f,options.eta_d,ref_div.get(),ref_score.get());
      CUDA_CHECK(cudaDeviceSynchronize());
      CUDA_CHECK(cudaMemcpy(ev.data(),ref_div.get(),ref_div.bytes(),cudaMemcpyDeviceToHost));
      save("dividend.bin",ev.data(),ev.size()*sizeof(double));
      CUDA_CHECK(cudaMemcpy(ev.data(),ref_score.get(),ref_score.bytes(),cudaMemcpyDeviceToHost));
      save("score.bin",ev.data(),ev.size()*sizeof(double));

    }
    std::uint64_t policy_mismatches = 0;
    std::uint64_t capital_policy_mismatches = 0;
    std::uint64_t default_mismatches = 0;
    std::uint64_t value_bit_mismatches = 0;
    double maximum_value_error = 0.0;
    for (std::size_t state = 0; state < states; ++state) {
      if (policy[state] != target_policy[state]) {
        ++policy_mismatches;
        if ((policy[state] - 1) / kDebt !=
            (target_policy[state] - 1) / kDebt) {
          ++capital_policy_mismatches;
        }
      }
      if (default_policy[state] != target_default[state]) ++default_mismatches;
      if (std::memcmp(&value[state], &target_value[state], sizeof(double)) != 0) {
        ++value_bit_mismatches;
      }
      maximum_value_error = std::max(
          maximum_value_error, std::abs(value[state] - target_value[state]));
    }
    int device = 0;
    cudaDeviceProp properties{};
    CUDA_CHECK(cudaGetDevice(&device));
    CUDA_CHECK(cudaGetDeviceProperties(&properties, device));
    std::cout << std::setprecision(17)
              << "{\"experiment\":\"de_gpu_compact_exact_action_gate_v2\""
              << ",\"scientific_role\":\"action_snapshot_not_smm_evaluator\""
              << ",\"device\":\"" << properties.name << "\""
              << ",\"threads_per_block\":" << options.threads
              << ",\"repeats\":" << options.repeats
              << ",\"mean_action_kernel_ms\":"
              << static_cast<double>(elapsed_ms) / options.repeats
              << ",\"mean_expected_plus_action_ms\":"
              << static_cast<double>(combined_ms) / options.repeats
              << ",\"policy_mismatches\":" << policy_mismatches
              << ",\"capital_policy_mismatches\":"
              << capital_policy_mismatches
              << ",\"default_mismatches\":" << default_mismatches
              << ",\"value_bit_mismatches\":" << value_bit_mismatches
              << ",\"maximum_absolute_value_error\":"
              << maximum_value_error << "}\n";
    return (policy_mismatches || default_mismatches || value_bit_mismatches) ? 2 : 0;
  } catch (const std::exception& error) {
    std::cerr << "de_action_snapshot: " << error.what() << '\n';
    return 1;
  }
}
