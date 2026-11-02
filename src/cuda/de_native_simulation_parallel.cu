#include <thread>
#include <mutex>
#include <atomic>
#include <exception>
#include <cub/cub.cuh>
#include <cuda_runtime.h>
#include <vector>
#include <utility>
#include <stdexcept>

struct ReusableDeviceStorage {
 std::vector<std::pair<size_t,void*>> free;std::mutex guard;
 void* acquire(size_t n){std::lock_guard<std::mutex> lock(guard);for(size_t i=0;i<free.size();++i)if(free[i].first==n){void* p=free[i].second;free.erase(free.begin()+i);return p;}void* p=nullptr;cudaError_t e=cudaMalloc(&p,n);if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));return p;}
 void release(void* p,size_t n){std::lock_guard<std::mutex> lock(guard);free.push_back({n,p});}
 ~ReusableDeviceStorage(){for(auto& x:free)cudaFree(x.second);}
};
ReusableDeviceStorage device_storage;
#include <sstream>
// Experimental cold-start DE policy solver. All numerical inputs are fixed-parameter grids.
// Target arrays are opened only after solving, for independent comparison.
#define main snapshot_main
#include "de_action_block_reduce.cu"
#undef main
#undef CUDA_CHECK
namespace simulation {
#define main run_simulation
#include "de_simulation_parallel.cu"
#undef main
}
#undef CUDA_CHECK
#define CUDA_CHECK(expr) cuda_check((expr), #expr)

__global__ void initialize_policy(int* p,int* d,double* v,double* dist) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 p[s]=450;d[s]=0;v[s]=0;dist[s]=0;
}
__global__ void prices(const double* tr,const int* def,const double* k,const double* b,double gamma,double* q) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 int e=s%625,a=s/625;double bp=b[a%30],kp=k[a/30],v=1.0;
 if(bp>0) {
  double rec=((1.0-kTax)*gamma)*((1.0-kDelta)*kp)/bp;
  rec=fmax(0.0,fmin(1.0,rec));
  for(int n=0;n<25;++n)v=fma(tr[e*25+n]*double(def[a*625+next_exogenous(e,n)]),rec-1.0,v);
 }
 q[s]=v/(1.0+kRate);
}
__global__ void reach(const int* policy,const int* def,int* mask) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 if(s<625)atomicExch(mask+s,1);
 if(def[s])return;int e=s%625,a=policy[s]-1;
 for(int n=0;n<25;++n)atomicExch(mask+a*625+next_exogenous(e,n),1);
}
__global__ void current_dividend(const int* policy,const double* q,const double* k,const double* b,
 const double* op,const double* ac,double phi,double ef,double ed,double* out) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 int e=s%625,j=s/625,a=policy[s]-1,ki=j/30,oi=(e%125)/25+5*(e%5)+25*ki;
 out[s]=exact_dividend(op[oi],op[oi+750],k[ki],k[a/30],ac[ki+30*(a/30)],q[a*625+e],b[a%30],b[j%30],phi,ef,ed);
}
__global__ void weighted_transition(const double* tr,double* weighted){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<625*25)weighted[(i%25)*625+i/25]=kBeta*tr[i];
}
__global__ void default_values(const double* tr,const double* v,double* vd) {
 int e=blockIdx.x*blockDim.x+threadIdx.x;if(e>=625)return;
 double r=0;for(int n=0;n<25;++n)r=fma(tr[n*625+e],v[next_exogenous(e,n)],r);vd[e]=r;
}
__global__ void howard(const double* tr,const int* p,const int* d,const int* mask,const double* div,
 const double* vd,const double* old,double* out,bool full) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500 || (!full&&!mask[s]))return;
 int e=s%625,a=p[s]-1;if(d[s]) {out[s]=vd[e];return;}
 double v=div[s];for(int n=0;n<25;++n)v=fma(tr[n*625+e],old[a*625+next_exogenous(e,n)],v);out[s]=v;
}
__global__ void policy_error(const int* a,const int* b,int* err) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s<562500)atomicMax(err,abs(a[s]-b[s]));
}
__global__ void dist_initialize(double* old) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s<562500)old[s]=1.0/562500.0;
}
// Compact zero-mass sources once per group; preserve every positive source's order.
// A warp compacts each group in source order; ballots affect offsets, never sums.
__global__ void compact_positive_sources(const int* offsets,const int* sources,
 const double* old,int* compact,int* counts) {
 int linear=blockIdx.x*blockDim.x+threadIdx.x,g=linear/32,lane=linear%32;if(g>=22500)return;
 int begin=offsets[g],end=offsets[g+1],n=0;
 unsigned prior=lane==0?0u:((1u<<lane)-1u);
 for(int base=begin;base<end;base+=32){
  int i=base+lane,src=i<end?sources[i]:0;bool keep=i<end && old[src]>0.0;
  unsigned mask=__ballot_sync(0xffffffff,keep);
  if(keep)compact[begin+n+__popc(mask&prior)]=src;
  n+=__popc(mask);
 }
 if(lane==0)counts[g]=n;
}

__global__ void dist_step(const double* tr,const int* offsets,const int* sources,const int* active_counts,const double* old,double* out) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 int e=s%625,a=s/625,pair=(e/125)*5+(e%25)/5,g=pair*900+a;
 double v=out[s];
 for(int i=offsets[g];i<offsets[g]+active_counts[g];++i) {
  int src=sources[i];double w=old[src];if(w>0)v=fma(tr[e*625+src%625],w,v);
 }
 out[s]=v;
}
__global__ void dist_stats(const double* old,const double* out,double* stats) {
 if(blockIdx.x||threadIdx.x)return;double total=0,err=-1e300;
 for(int s=0;s<562500;++s) {total+=out[s];err=fmax(err,out[s]-old[s]);}
 stats[0]=total;stats[1]=err;
}
__global__ void dist_normalize(double* old,double* out,const double* stats) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s<562500){old[s]=out[s]/stats[0];out[s]=0;}
}
__global__ void default_error(const double* dist,const int* current,const int* old,double* stats) {
 if(blockIdx.x||threadIdx.x)return;double err=0;
 for(int s=0;s<562500;++s)err=fma(dist[s],double(abs(current[s]-old[s])),err);
 stats[0]=err;
}

__global__ void profit_table(const int* policy,const int* def,const double* k,const double* b,
 const double* op,const double* ac,double phi,double* profit) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 int e=s%625,j=s/625,ki=j/30,a=def[s]?0:policy[s]-1,oi=(e%125)/25+5*(e%5)+25*ki;
 double operating=fma(op[oi],op[oi+750],-ac[ki+30*(a/30)])-phi;
 double p=fma(1.0-kTax,operating,-kDelta*k[ki]);
 profit[s]=fma(kTax,fma(kRate,b[j%30],kDelta*k[ki]),p);
}

template<class T> void upload(DeviceBuffer<T>& d,const std::vector<T>& x) {CUDA_CHECK(cudaMemcpy(d.get(),x.data(),d.bytes(),cudaMemcpyHostToDevice));}
template<class T> std::vector<T> download(DeviceBuffer<T>& d) {std::vector<T>x(d.bytes()/sizeof(T));CUDA_CHECK(cudaMemcpy(x.data(),d.get(),d.bytes(),cudaMemcpyDeviceToHost));return x;}
template<class T> void save_array(const std::string& path,const std::vector<T>&x) {std::ofstream f(path,std::ios::binary);f.write((char*)x.data(),x.size()*sizeof(T));}
// Pack strictly positive mass within each contiguous tile, retaining linear order.
__global__ void pack_mass_tiles(const double* mass,const double* old,double* packed,int* counts,double* errors) {
 __shared__ int warp_counts[8];__shared__ double warp_errors[8];
 int lane=threadIdx.x%32,warp=threadIdx.x/32,index=blockIdx.x*256+threadIdx.x;
 double v=index<562500?mass[index]:0.0;
 double err=index<562500?fmax(v-old[index],0.0):0.0;
 unsigned mask=__ballot_sync(0xffffffff,v>0.0);
 if(lane==0)warp_counts[warp]=__popc(mask);
 for(int step=16;step>0;step/=2)err=fmax(err,__shfl_down_sync(0xffffffff,err,step));
 if(lane==0)warp_errors[warp]=err;
 __syncthreads();
 int offset=0;for(int w=0;w<warp;++w)offset+=warp_counts[w];
 unsigned prior=lane==0?0u:((1u<<lane)-1u);
 if(v>0.0)packed[blockIdx.x*256+offset+__popc(mask&prior)]=v;
 if(threadIdx.x==0){int count=0;double e=0;for(int w=0;w<8;++w){count+=warp_counts[w];e=fmax(e,warp_errors[w]);}counts[blockIdx.x]=count;errors[blockIdx.x]=e;}
}
__global__ void mass_tile_prefix(const int* counts,const double* errors,int* offsets,double* meta) {
 __shared__ int sums[2198];__shared__ double errs[2198];
 for(int i=threadIdx.x;i<2198;i+=blockDim.x){sums[i]=counts[i];errs[i]=errors[i];}
 __syncthreads();
 if(threadIdx.x==0){int total=0;double e=0;for(int i=0;i<2198;++i){int n=sums[i];sums[i]=total;total+=n;e=fmax(e,errs[i]);}meta[0]=double(total);meta[1]=e;}
 __syncthreads();
 for(int i=threadIdx.x;i<2198;i+=blockDim.x)offsets[i]=sums[i];
}
__global__ void scatter_mass_tiles(const double* tiles,const int* counts,const int* offsets,double* compact) {
 int i=threadIdx.x;if(i<counts[blockIdx.x])compact[offsets[blockIdx.x]+i]=tiles[blockIdx.x*256+i];
}

template<class T> class PinnedBuffer {
 T* p=nullptr;
 public:
 explicit PinnedBuffer(size_t n){CUDA_CHECK(cudaMallocHost(reinterpret_cast<void**>(&p),n*sizeof(T)));}
 ~PinnedBuffer(){if(p)cudaFreeHost(p);}
 T* get(){return p;}
};


__device__ void initialize_policy_body(int* p,int* d,double* v,double* dist) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 p[s]=450;d[s]=0;v[s]=0;dist[s]=0;
}
__device__ void prices_body(const double* tr,const int* def,const double* k,const double* b,double gamma,double* q) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 int e=s%625,a=s/625;double bp=b[a%30],kp=k[a/30],v=1.0;
 if(bp>0) {
  double rec=((1.0-kTax)*gamma)*((1.0-kDelta)*kp)/bp;
  rec=fmax(0.0,fmin(1.0,rec));
  for(int n=0;n<25;++n)v=fma(tr[e*25+n]*double(def[a*625+next_exogenous(e,n)]),rec-1.0,v);
 }
 q[s]=v/(1.0+kRate);
}
__device__ void reach_body(const int* policy,const int* def,int* mask) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 if(s<625)atomicExch(mask+s,1);
 if(def[s])return;int e=s%625,a=policy[s]-1;
 for(int n=0;n<25;++n)atomicExch(mask+a*625+next_exogenous(e,n),1);
}
__device__ void current_dividend_body(const int* policy,const double* q,const double* k,const double* b,
 const double* op,const double* ac,double phi,double ef,double ed,double* out) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 int e=s%625,j=s/625,a=policy[s]-1,ki=j/30,oi=(e%125)/25+5*(e%5)+25*ki;
 out[s]=exact_dividend(op[oi],op[oi+750],k[ki],k[a/30],ac[ki+30*(a/30)],q[a*625+e],b[a%30],b[j%30],phi,ef,ed);
}
__device__ void weighted_transition_body(const double* tr,double* weighted){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<625*25)weighted[(i%25)*625+i/25]=kBeta*tr[i];
}
__device__ void default_values_body(const double* tr,const double* v,double* vd) {
 int e=blockIdx.x*blockDim.x+threadIdx.x;if(e>=625)return;
 double r=0;for(int n=0;n<25;++n)r=fma(tr[n*625+e],v[next_exogenous(e,n)],r);vd[e]=r;
}
__device__ void howard_body(const double* tr,const int* p,const int* d,const int* mask,const double* div,
 const double* vd,const double* old,double* out,bool full) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500 || (!full&&!mask[s]))return;
 int e=s%625,a=p[s]-1;if(d[s]) {out[s]=vd[e];return;}
 double v=div[s];for(int n=0;n<25;++n)v=fma(tr[n*625+e],old[a*625+next_exogenous(e,n)],v);out[s]=v;
}
__device__ void policy_error_body(const int* a,const int* b,int* err) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s<562500)atomicMax(err,abs(a[s]-b[s]));
}
__device__ void dist_initialize_body(double* old) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s<562500)old[s]=1.0/562500.0;
}
__device__ void compact_positive_sources_body(const int* offsets,const int* sources,
 const double* old,int* compact,int* counts) {
 int linear=blockIdx.x*blockDim.x+threadIdx.x,g=linear/32,lane=linear%32;if(g>=22500)return;
 int begin=offsets[g],end=offsets[g+1],n=0;
 unsigned prior=lane==0?0u:((1u<<lane)-1u);
 for(int base=begin;base<end;base+=32){
  int i=base+lane,src=i<end?sources[i]:0;bool keep=i<end && old[src]>0.0;
  unsigned mask=__ballot_sync(0xffffffff,keep);
  if(keep)compact[begin+n+__popc(mask&prior)]=src;
  n+=__popc(mask);
 }
 if(lane==0)counts[g]=n;
}
__device__ void dist_step_body(const double* tr,const int* offsets,const int* sources,const int* active_counts,const double* old,double* out) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 int e=s%625,a=s/625,pair=(e/125)*5+(e%25)/5,g=pair*900+a;
 double v=out[s];
 for(int i=offsets[g];i<offsets[g]+active_counts[g];++i) {
  int src=sources[i];double w=old[src];if(w>0)v=fma(tr[e*625+src%625],w,v);
 }
 out[s]=v;
}
__device__ void dist_normalize_body(double* old,double* out,const double* stats) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s<562500){old[s]=out[s]/stats[0];out[s]=0;}
}
__device__ void pack_mass_tiles_body(const double* mass,const double* old,double* packed,int* counts,double* errors) {
 __shared__ int warp_counts[8];__shared__ double warp_errors[8];
 int lane=threadIdx.x%32,warp=threadIdx.x/32,index=blockIdx.x*256+threadIdx.x;
 double v=index<562500?mass[index]:0.0;
 double err=index<562500?fmax(v-old[index],0.0):0.0;
 unsigned mask=__ballot_sync(0xffffffff,v>0.0);
 if(lane==0)warp_counts[warp]=__popc(mask);
 for(int step=16;step>0;step/=2)err=fmax(err,__shfl_down_sync(0xffffffff,err,step));
 if(lane==0)warp_errors[warp]=err;
 __syncthreads();
 int offset=0;for(int w=0;w<warp;++w)offset+=warp_counts[w];
 unsigned prior=lane==0?0u:((1u<<lane)-1u);
 if(v>0.0)packed[blockIdx.x*256+offset+__popc(mask&prior)]=v;
 if(threadIdx.x==0){int count=0;double e=0;for(int w=0;w<8;++w){count+=warp_counts[w];e=fmax(e,warp_errors[w]);}counts[blockIdx.x]=count;errors[blockIdx.x]=e;}
}
__device__ void mass_tile_prefix_body(const int* counts,const double* errors,int* offsets,double* meta) {
 __shared__ int sums[2198];__shared__ double errs[2198];
 for(int i=threadIdx.x;i<2198;i+=blockDim.x){sums[i]=counts[i];errs[i]=errors[i];}
 __syncthreads();
 if(threadIdx.x==0){int total=0;double e=0;for(int i=0;i<2198;++i){int n=sums[i];sums[i]=total;total+=n;e=fmax(e,errs[i]);}meta[0]=double(total);meta[1]=e;}
 __syncthreads();
 for(int i=threadIdx.x;i<2198;i+=blockDim.x)offsets[i]=sums[i];
}
__device__ void scatter_mass_tiles_body(const double* tiles,const int* counts,const int* offsets,double* compact) {
 int i=threadIdx.x;if(i<counts[blockIdx.x])compact[offsets[blockIdx.x]+i]=tiles[blockIdx.x*256+i];
}
__device__ void profit_table_body(const int* policy,const int* def,const double* k,const double* b,
 const double* op,const double* ac,double phi,double* profit) {
 int s=blockIdx.x*blockDim.x+threadIdx.x;if(s>=562500)return;
 int e=s%625,j=s/625,ki=j/30,a=def[s]?0:policy[s]-1,oi=(e%125)/25+5*(e%5)+25*ki;
 double operating=fma(op[oi],op[oi+750],-ac[ki+30*(a/30)])-phi;
 double p=fma(1.0-kTax,operating,-kDelta*k[ki]);
 profit[s]=fma(kTax,fma(kRate,b[j%30],kDelta*k[ki]),p);
}
__device__ void expected_value_kernel_body(
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
__device__ void shortlist_kernel_body(
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
__device__ void compact_action_kernel_body(
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
__device__ void ordered_mass_sum_body(const double* mass,const double* meta,double* stats) {
 __shared__ double tile[1024];
 int n=int(meta[0]);double sum=0.0;
 for(int base=0;base<n;base+=1024){
  for(int i=threadIdx.x;i<1024;i+=blockDim.x)tile[i]=(base+i<n)?mass[base+i]:0.0;
  __syncthreads();
  if(threadIdx.x==0)for(int i=0;i<1024 && base+i<n;++i)sum+=tile[i];
  __syncthreads();
 }
 if(threadIdx.x==0){stats[0]=sum;stats[1]=meta[1];}
}
__device__ void pack_default_body(const double* mass,const int* current,const int* prior_def,double* packed,int* counts,double* errors) {
 __shared__ int warp_counts[8];__shared__ double warp_errors[8];
 int lane=threadIdx.x%32,warp=threadIdx.x/32,index=blockIdx.x*256+threadIdx.x;
 double v=index<562500?mass[index]*double(abs(current[index]-prior_def[index])):0.0;
 double err=0.0;
 unsigned mask=__ballot_sync(0xffffffff,v>0.0);
 if(lane==0)warp_counts[warp]=__popc(mask);
 for(int step=16;step>0;step/=2)err=fmax(err,__shfl_down_sync(0xffffffff,err,step));
 if(lane==0)warp_errors[warp]=err;
 __syncthreads();
 int offset=0;for(int w=0;w<warp;++w)offset+=warp_counts[w];
 unsigned prior=lane==0?0u:((1u<<lane)-1u);
 if(v>0.0)packed[blockIdx.x*256+offset+__popc(mask&prior)]=v;
 if(threadIdx.x==0){int count=0;double e=0;for(int w=0;w<8;++w){count+=warp_counts[w];e=fmax(e,warp_errors[w]);}counts[blockIdx.x]=count;errors[blockIdx.x]=e;}
}
__device__ __forceinline__ int hw_full_state(int group,int lane){
 int a=group/25,pair=group%25;return a*625+(pair/5)*125+(lane/5)*25+(pair%5)*5+lane%5;
}
__global__ void hw_reorder_weights(const double* tr,double* weights){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=15625)return;int n=i/625,e=i%625;
 int ef=((e/125)*5+(e%25)/5)*25+((e%125)/25)*5+e%5;weights[n*625+ef]=tr[n*625+e];
}
__device__ void hw_group_flags_body(const int* mask,int* flags){
 int g=blockIdx.x*blockDim.x+threadIdx.x;if(g<22500)flags[g]=mask[hw_full_state(g,0)];
}
__device__ void hw_group_inverse_body(const int* ids,int groups,int* inverse){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<groups)inverse[ids[i]]=i*25;
}
__device__ void hw_group_plan_body(const int* ids,int m,const int* inverse,const int* p,const int* d,
 const double* div,const double* old,int* next,int* exog,double* base,double* packed){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=m)return;
 int s=hw_full_state(ids[i/25],i%25),e=s%625,a=d[s]?0:p[s]-1;
 int pair=((e%125)/25)*5+e%5;
 next[i]=inverse[a*25+pair];exog[i]=((e/125)*5+(e%25)/5)*25+((e%125)/25)*5+e%5;base[i]=d[s]?0.0:div[s];packed[i]=old[s];
}
__device__ void hw_group_step_body(int m,const int* __restrict__ next,const int* __restrict__ exog,
 const double* __restrict__ weights,const double* __restrict__ base,const double* __restrict__ old,double* __restrict__ out){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=m)return;
 int a=next[i],e=exog[i];double v=base[i];
 #pragma unroll
 for(int n=0;n<25;++n)v=fma(weights[n*625+e],old[a+n],v);
 out[i]=v;
}
__device__ void hw_group_scatter_body(const int* ids,int m,const double* compact,double* v){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<m)v[hw_full_state(ids[i/25],i%25)]=compact[i];
}

struct BatchJob {
 double *hw_weighted,*hw_base,*hw_va,*hw_vb;
 int *hw_flags,*hw_ids,*hw_count,*hw_inverse,*hw_edges,*hw_exog;
 double* tr;
 double* ctr;
 double* weighted_tr;
 double* true_tr;
 double* kg;
 double* bg;
 double* op;
 double* ac;
 double* base;
 double* vin;
 double* vtemp;
 double* vout;
 double* q;
 double* div;
 double* ev;
 double* vd;
 double* dist;
 double* distold;
 double* stats;
 double* bins;
 double* packed_tiles;
 double* compact_mass;
 double* tile_errors;
 double* mass_meta;
 double* prof;
 int* policy;
 int* oldpolicy;
 int* def;
 int* olddef;
 int* mask;
 int* err;
 int* counts;
 int* offsets;
 int* sources;
 int* active_sources;
 int* active_counts;
 int* tile_counts;
 int* tile_offsets;
 Candidate* candidates;
 double eta_k,phi,gamma,eta_f,eta_d,sigma_pi;
 int active=1,policy_active=1,dist_active=0,copy_policy=0,outer=0,total_pol=0,pe=0,dist_iterations=0;
};
__global__ void batch_prices(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.active))return;prices_body(j.ctr,j.olddef,j.kg,j.bg,j.gamma,j.q);}
__global__ void batch_dividend(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.policy_active))return;current_dividend_body(j.oldpolicy,j.q,j.kg,j.bg,j.op,j.ac,j.phi,j.eta_f,j.eta_d,j.div);}
__global__ void batch_reach(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.policy_active))return;reach_body(j.oldpolicy,j.olddef,j.mask);}
__global__ void batch_default(BatchJob* jobs,int it){auto& j=jobs[blockIdx.y];if(!(j.policy_active))return;default_values_body(j.weighted_tr,(it%2?j.vtemp:j.vin),j.vd);}
__global__ void batch_howard(BatchJob* jobs,int it){auto& j=jobs[blockIdx.y];if(!(j.policy_active))return;howard_body(j.weighted_tr,j.oldpolicy,j.olddef,j.mask,j.div,j.vd,(it%2?j.vtemp:j.vin),(it%2?j.vin:j.vtemp),it==99);}
__global__ void batch_expected(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.policy_active))return;expected_value_kernel_body(j.ctr,j.vin,j.ev,j.vd);}
__global__ void batch_shortlist(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.policy_active))return;shortlist_kernel_body(j.ctr,j.vin,j.ev,j.vd,j.q,j.kg,j.bg,j.op,j.ac,j.base,j.eta_k,j.phi,j.eta_f,j.eta_d,j.policy,j.def,j.vout,j.candidates,j.counts,j.bins);}
__global__ void batch_action(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.policy_active))return;compact_action_kernel_body(j.ctr,j.vin,j.ev,j.vd,j.q,j.kg,j.bg,j.op,j.ac,j.base,j.eta_k,j.phi,j.eta_f,j.eta_d,j.policy,j.def,j.vout,j.candidates,j.counts,j.bins);}
__global__ void batch_error(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.policy_active))return;policy_error_body(j.policy,j.oldpolicy,j.err);}
__global__ void batch_compact(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.dist_active))return;compact_positive_sources_body(j.offsets,j.sources,j.distold,j.active_sources,j.active_counts);}
__global__ void batch_dist(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.dist_active))return;dist_step_body(j.true_tr,j.offsets,j.active_sources,j.active_counts,j.distold,j.dist);}
__global__ void batch_pack(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.dist_active))return;pack_mass_tiles_body(j.dist,j.distold,j.packed_tiles,j.tile_counts,j.tile_errors);}
__global__ void batch_prefix(BatchJob* jobs,bool defaults){auto& j=jobs[blockIdx.y];if(!(defaults?j.active:j.dist_active))return;mass_tile_prefix_body(j.tile_counts,j.tile_errors,j.tile_offsets,j.mass_meta);}
__global__ void batch_scatter(BatchJob* jobs,bool defaults){auto& j=jobs[blockIdx.y];if(!(defaults?j.active:j.dist_active))return;scatter_mass_tiles_body(j.packed_tiles,j.tile_counts,j.tile_offsets,j.compact_mass);}
__global__ void batch_sum(BatchJob* jobs,bool defaults){auto& j=jobs[blockIdx.y];if(!(defaults?j.active:j.dist_active))return;ordered_mass_sum_body(j.compact_mass,j.mass_meta,j.stats);}
__global__ void batch_normalize(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.dist_active))return;dist_normalize_body(j.distold,j.dist,j.stats);}
__global__ void batch_profit(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(true))return;profit_table_body(j.policy,j.def,j.kg,j.bg,j.op,j.ac,j.phi,j.prof);}
__global__ void batch_pack_default(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!(j.active))return;pack_default_body(j.dist,j.def,j.olddef,j.packed_tiles,j.tile_counts,j.tile_errors);}


__global__ void batch_hw_flags(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!j.policy_active)return;hw_group_flags_body(j.mask,j.hw_flags);}
__global__ void batch_hw_count(BatchJob* jobs,int* counts){int c=blockIdx.x;auto& j=jobs[c];if(threadIdx.x==0)counts[c]=j.policy_active?25*(*j.hw_count):0;}
__global__ void batch_hw_inverse(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!j.policy_active)return;hw_group_inverse_body(j.hw_ids,*j.hw_count,j.hw_inverse);}
__global__ void batch_hw_plan(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!j.policy_active)return;hw_group_plan_body(j.hw_ids,25*(*j.hw_count),j.hw_inverse,j.oldpolicy,j.olddef,j.div,j.vin,j.hw_edges,j.hw_exog,j.hw_base,j.hw_va);}
__global__ void batch_hw_step(BatchJob* jobs,int it){auto& j=jobs[blockIdx.y];if(!j.policy_active)return;hw_group_step_body(25*(*j.hw_count),j.hw_edges,j.hw_exog,j.hw_weighted,j.hw_base,it%2?j.hw_vb:j.hw_va,it%2?j.hw_va:j.hw_vb);}
__global__ void batch_hw_scatter(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!j.policy_active)return;hw_group_scatter_body(j.hw_ids,25*(*j.hw_count),j.hw_vb,j.vtemp);}
__global__ void begin_outer(BatchJob* jobs,int outer){auto& j=jobs[blockIdx.x];if(threadIdx.x||!j.active)return;j.outer=outer;j.policy_active=1;}
__global__ void clear_policy(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!j.policy_active)return;int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<562500)j.mask[i]=0;if(i==0)*j.err=0;}
__global__ void decide_policy(BatchJob* jobs){auto& j=jobs[blockIdx.x];if(threadIdx.x)return;j.copy_policy=0;if(!j.policy_active)return;j.pe=*j.err;++j.total_pol;j.copy_policy=j.pe!=0;if(j.pe==0)j.policy_active=0;}
__global__ void accept_policy(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!j.copy_policy)return;int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<562500){j.oldpolicy[i]=j.policy[i];j.vin[i]=j.vout[i];}}
__global__ void begin_dist(BatchJob* jobs){auto& j=jobs[blockIdx.x];if(threadIdx.x)return;j.dist_active=j.active;j.policy_active=0;j.dist_iterations=0;}
__global__ void decide_dist(BatchJob* jobs){auto& j=jobs[blockIdx.x];if(threadIdx.x||!j.dist_active)return;++j.dist_iterations;if(j.stats[1]<1e-6)j.dist_active=0;}
__global__ void decide_outer(BatchJob* jobs){auto& j=jobs[blockIdx.x];if(threadIdx.x||!j.active)return;if(j.stats[0]<1e-4 && j.outer>1)j.active=0;}
__global__ void accept_defaults(BatchJob* jobs){auto& j=jobs[blockIdx.y];if(!j.active)return;int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<562500)j.olddef[i]=j.def[i];}

#include <filesystem>
#include <memory>
template<class T> std::vector<T> take(T* p,size_t n){std::vector<T>x(n);CUDA_CHECK(cudaMemcpy(x.data(),p,n*sizeof(T),cudaMemcpyDeviceToHost));return x;}
template<class T> void put(T* p,const std::vector<T>&x){CUDA_CHECK(cudaMemcpy(p,x.data(),x.size()*sizeof(T),cudaMemcpyHostToDevice));}

__global__ void adjacency_keys(const int* policy,const int* def,int* keys,int* values){
 int u=blockIdx.x*blockDim.x+threadIdx.x;if(u>=562500)return;int e=u/900,k=u%900,state=k*625+e,action=def[state]?0:policy[state]-1,pair=((e%125)/25)*5+e%5;
 keys[u]=pair*900+action;values[u]=state;
}
__global__ void adjacency_offsets(const int* keys,int* offsets){
 int group=blockIdx.x*blockDim.x+threadIdx.x;if(group>22500)return;int lo=0,hi=562500;
 while(lo<hi){int mid=(lo+hi)/2;if(keys[mid]<group)lo=mid+1;else hi=mid;}offsets[group]=lo;
}

struct OwnedJob {
 void* hw_select_temp;size_t hw_select_bytes=0;
 int* adj_keys;int* adj_sorted;int* adj_values;void* sort_temp;size_t sort_bytes=0;
 BatchJob j;Options options;std::vector<std::pair<void*,size_t>> allocations;
 template<class T>T* alloc(size_t n){void* p=device_storage.acquire(n*sizeof(T));allocations.push_back({p,n*sizeof(T)});return static_cast<T*>(p);}
 ~OwnedJob(){for(auto& a:allocations)device_storage.release(a.first,a.second);}
 OwnedJob(const Options& o):options(o){
 j.hw_weighted=alloc<double>(15625);j.hw_base=alloc<double>(562500);j.hw_va=alloc<double>(562500);j.hw_vb=alloc<double>(562500);
 j.hw_flags=alloc<int>(22500);j.hw_ids=alloc<int>(22500);j.hw_count=alloc<int>(1);j.hw_inverse=alloc<int>(22500);j.hw_edges=alloc<int>(562500);j.hw_exog=alloc<int>(562500);
 cub::CountingInputIterator<int> hw_input(0);CUDA_CHECK(cub::DeviceSelect::Flagged(nullptr,hw_select_bytes,hw_input,j.hw_flags,j.hw_ids,j.hw_count,22500));hw_select_temp=alloc<unsigned char>(hw_select_bytes);
 j.tr=alloc<double>(390625);
 j.ctr=alloc<double>(15625);
 j.weighted_tr=alloc<double>(15625);
 j.true_tr=alloc<double>(390625);
 j.kg=alloc<double>(30);
 j.bg=alloc<double>(30);
 j.op=alloc<double>(1500);
 j.ac=alloc<double>(900);
 j.base=alloc<double>(900);
 j.vin=alloc<double>(562500);
 j.vtemp=alloc<double>(562500);
 j.vout=alloc<double>(562500);
 j.q=alloc<double>(562500);
 j.div=alloc<double>(562500);
 j.ev=alloc<double>(562500);
 j.vd=alloc<double>(625);
 j.dist=alloc<double>(562500);
 j.distold=alloc<double>(562500);
 j.stats=alloc<double>(2);
 j.bins=alloc<double>(1162500);
 j.packed_tiles=alloc<double>(562688);
 j.compact_mass=alloc<double>(562500);
 j.tile_errors=alloc<double>(2198);
 j.mass_meta=alloc<double>(2);
 j.prof=alloc<double>(562500);
 j.policy=alloc<int>(562500);
 j.oldpolicy=alloc<int>(562500);
 j.def=alloc<int>(562500);
 j.olddef=alloc<int>(562500);
 j.mask=alloc<int>(562500);
 j.err=alloc<int>(1);
 j.counts=alloc<int>(18750);
 j.offsets=alloc<int>(22501);
 j.sources=alloc<int>(562500);
 j.active_sources=alloc<int>(562500);
 j.active_counts=alloc<int>(22500);
 j.tile_counts=alloc<int>(2198);
 j.tile_offsets=alloc<int>(2198);

 adj_keys=alloc<int>(562500);adj_sorted=alloc<int>(562500);adj_values=alloc<int>(562500);
 CUDA_CHECK(cub::DeviceRadixSort::SortPairs(nullptr,sort_bytes,adj_keys,adj_sorted,adj_values,j.sources,562500,0,15));sort_temp=alloc<unsigned char>(sort_bytes);
 j.candidates=alloc<Candidate>(18750*870);j.eta_k=o.eta_k;j.phi=o.phi;j.gamma=o.gamma;j.eta_f=o.eta_f;j.eta_d=o.eta_d;j.sigma_pi=o.sigma_pi;
 put(j.tr,read_exact<double>(o.fixture+"/RCC_pr_mat_exog_theta.bin",390625));put(j.true_tr,read_exact<double>(o.fixture+"/RCC_pr_mat_exog.bin",390625));
 auto kh=read_exact<double>(o.fixture+"/RCC_k0.bin",30);put(j.kg,kh);put(j.bg,read_exact<double>(o.fixture+"/RCC_b0.bin",30));
 auto ah=read_exact<double>(o.fixture+"/RCC_A0.bin",5),zh=read_exact<double>(o.fixture+"/RCC_z0.bin",5);std::vector<double> oh(1500);
 for(int k=0;k<30;++k)for(int z=0;z<5;++z)for(int a=0;a<5;++a){int i=a+5*z+25*k;double az=ah[a]*zh[z];oh[i]=(az*az)*.5;oh[i+750]=std::pow(kh[k],.5);}put(j.op,oh);
 put(j.base,read_exact<double>(o.fixture+"/RCC_actionkbase.bin",900));
 compact_transition_kernel<<<62,256>>>(j.tr,j.ctr);weighted_transition<<<62,256>>>(j.ctr,j.weighted_tr);hw_reorder_weights<<<62,256>>>(j.weighted_tr,j.hw_weighted);exact_adjustment_kernel<<<4,256>>>(j.kg,o.eta_k,j.ac);
 initialize_policy<<<2198,256>>>(j.oldpolicy,j.olddef,j.vin,j.dist);
 }
 void howard_select(){cub::CountingInputIterator<int> input(0);CUDA_CHECK(cub::DeviceSelect::Flagged(hw_select_temp,hw_select_bytes,input,j.hw_flags,j.hw_ids,j.hw_count,22500));}
 void adjacency(){
 adjacency_keys<<<2198,256>>>(j.policy,j.def,adj_keys,adj_values);
 CUDA_CHECK(cub::DeviceRadixSort::SortPairs(sort_temp,sort_bytes,adj_keys,adj_sorted,adj_values,j.sources,562500,0,15));
 adjacency_offsets<<<88,256>>>(adj_sorted,j.offsets);dist_initialize<<<2198,256>>>(j.distold);
 }

};
int finish_job(OwnedJob& owner,const BatchJob& j,int index,std::ostringstream& report,std::ostringstream& errors){
 auto& o=owner.options;std::string dump;
 if(const char* p=getenv("RCC_BATCH_DUMP")){dump=std::string(p)+"/job"+std::to_string(index);std::filesystem::create_directories(dump);}
 auto ph=take(j.policy,562500),dh=take(j.def,562500);auto qh=take(j.q,562500),profh=take(j.prof,562500);long bad=0;
 auto compare=[&](auto& values,const char* file,const char* output){using T=typename std::decay_t<decltype(values)>::value_type;
   if(!o.production){auto expected=read_exact<T>(o.fixture+"/"+file,values.size());for(size_t i=0;i<values.size();++i)bad+=memcmp(&values[i],&expected[i],sizeof(T))!=0;}
   if(!dump.empty())save_array(dump+"/"+output,values);};
 compare(ph,"RCC_kbpolind.bin","policy.bin");compare(dh,"RCC_dfpolind.bin","policy_default.bin");compare(qh,"RCC_qmat.bin","q.bin");compare(profh,"RCC_pimat.bin","profit_table.bin");
 if(!o.production || !dump.empty()) {auto v=take(j.vout,562500),old=take(j.vin,562500),mass=take(j.dist,562500);
 compare(v,"RCC_Vtheta.bin","value.bin");compare(old,"RCC_Vthetaold.bin","old.bin");compare(mass,"RCC_dist.bin","dist.bin");}
 report<<"{\"family\":\"native_batch_candidate\",\"job\":"<<index<<",\"outer_iterations\":"<<j.outer<<",\"policy_iterations\":"<<j.total_pol
       <<",\"reference_checked\":"<<(o.production?"false":"true");if(!o.production)report<<",\"array_bit_mismatches\":"<<bad;report<<"}\n";
 simulation::SolvedInputs inputs{ph,dh,qh,profh};simulation::solved_inputs=&inputs;simulation::report_stream=&report;simulation::error_stream=&errors;
 std::vector<std::vector<double>> basis;int status=0;
 for(int lane=0;lane<(o.profile?3:1);++lane){
   double sigma=o.profile?(lane==0?0.0:(lane==1?1.0:-1.0)):o.sigma_pi;
   simulation::reuse_macro=(lane>0);simulation::dump_directory=dump;
   if(o.profile&&!dump.empty()){simulation::dump_directory=dump+"/lane"+std::to_string(lane);std::filesystem::create_directories(simulation::dump_directory);}
   std::ostringstream value;value<<std::setprecision(17)<<sigma;
   std::vector<std::string> text={"simulation","--fixture",o.fixture,"--sigma-pi",value.str()};if(o.production)text.push_back("--production");
   std::vector<char*> args;for(auto& x:text)args.push_back(x.data());status=simulation::run_simulation(args.size(),args.data());if(status)break;basis.push_back(simulation::last_moments);
 }
 simulation::solved_inputs=nullptr;simulation::report_stream=nullptr;simulation::error_stream=nullptr;simulation::dump_directory.clear();simulation::reuse_macro=false;
 if(status)return status;
 auto write_vector=[&](const std::vector<double>& values){report<<'[';for(size_t i=0;i<values.size();++i){if(i)report<<',';if(std::isfinite(values[i]))report<<values[i];else report<<"null";}report<<']';};
 report<<std::setprecision(17)<<"{\"family\":\"evaluation\",\"job\":"<<index<<",\"profile\":"<<(o.profile?"true":"false");
 if(o.profile){report<<",\"m0\":";write_vector(basis[0]);report<<",\"mp\":";write_vector(basis[1]);report<<",\"mm\":";write_vector(basis[2]);}
 else{report<<",\"moments\":";write_vector(basis[0]);}
 report<<",\"macro_output\":";write_vector(simulation::last_macro_output);report<<"}\n";
 return bad?2:0;
}
int run_batch(int argc,char**argv){try{
 if((argc!=3&&argc!=5)||std::string(argv[1])!="--batch")throw std::runtime_error("Usage: native --batch jobs.txt [--simulation-workers N]");
 int workers=4;if(argc==5){if(std::string(argv[3])!="--simulation-workers")throw std::runtime_error("Unknown option");workers=std::stoi(argv[4]);if(workers<1||workers>64)throw std::runtime_error("Simulation workers must be in 1..64");}
 auto start=std::chrono::steady_clock::now();std::ifstream file(argv[2]);if(!file)throw std::runtime_error("Cannot open jobs file");std::string line;std::vector<std::unique_ptr<OwnedJob>> owners;std::vector<BatchJob> jobs;
 while(std::getline(file,line)){if(line.empty()||line[0]=='#')continue;std::istringstream in(line);std::vector<std::string> text={"evaluator"};std::string word;while(in>>std::quoted(word))text.push_back(word);std::vector<char*> args;for(auto& x:text)args.push_back(x.data());if(owners.size()>=64)throw std::runtime_error("Maximum batch size is 64");owners.emplace_back(new OwnedJob(parse_options(args.size(),args.data())));jobs.push_back(owners.back()->j);}
 int P=jobs.size();if(!P||P>64)throw std::runtime_error("Batch must contain 1..64 cases");DeviceBuffer<BatchJob> dj(P);upload(dj,jobs);auto* jp=dj.get();CUDA_CHECK(cudaDeviceSynchronize());auto solve_start=std::chrono::steady_clock::now();
 DeviceBuffer<int> hw_batch_counts(P);cudaEvent_t hs,he;cudaEventCreate(&hs);cudaEventCreate(&he);double hw_total=0;
 auto refresh=[&](){jobs=download(dj);};
 for(int outer=1;outer<=50;++outer){
  begin_outer<<<P,1>>>(jp,outer);batch_prices<<<dim3(2198,P),256>>>(jp);
  for(int pol=1;pol<=20;++pol){
   batch_dividend<<<dim3(2198,P),256>>>(jp);clear_policy<<<dim3(2198,P),256>>>(jp);batch_reach<<<dim3(2198,P),256>>>(jp);
   cudaEventRecord(hs);batch_hw_flags<<<dim3(88,P),256>>>(jp);
   for(int c=0;c<P;++c)if(pol==1?jobs[c].active:jobs[c].policy_active)owners[c]->howard_select();
   batch_hw_count<<<P,1>>>(jp,hw_batch_counts.get());auto hc=download(hw_batch_counts);int hm=*std::max_element(hc.begin(),hc.end());
   batch_hw_inverse<<<dim3(88,P),256>>>(jp);batch_hw_plan<<<dim3((hm+127)/128,P),128>>>(jp);
   for(int it=0;it<99;++it)batch_hw_step<<<dim3((hm+127)/128,P),128>>>(jp,it);
   batch_hw_scatter<<<dim3((hm+127)/128,P),128>>>(jp);
   batch_default<<<dim3(3,P),256>>>(jp,99);batch_howard<<<dim3(2198,P),256>>>(jp,99);cudaEventRecord(he);
   batch_expected<<<dim3(2198,P),256>>>(jp);batch_shortlist<<<dim3(18750,P),256>>>(jp);batch_action<<<dim3(2198,P),256>>>(jp);
   batch_error<<<dim3(2198,P),256>>>(jp);decide_policy<<<P,1>>>(jp);accept_policy<<<dim3(2198,P),256>>>(jp);refresh();float hms;cudaEventElapsedTime(&hms,hs,he);hw_total+=hms;bool any=false;for(auto& j:jobs)any|=j.policy_active;if(!any)break;
  }
  for(int c=0;c<P;++c)if(jobs[c].active)owners[c]->adjacency();begin_dist<<<P,1>>>(jp);
  for(int it=1;it<=1500;++it){
   batch_compact<<<dim3((22500*32+255)/256,P),256>>>(jp);batch_dist<<<dim3(2198,P),256>>>(jp);batch_pack<<<dim3(2198,P),256>>>(jp);
   batch_prefix<<<dim3(1,P),256>>>(jp,false);batch_scatter<<<dim3(2198,P),256>>>(jp,false);batch_sum<<<dim3(1,P),256>>>(jp,false);
   decide_dist<<<P,1>>>(jp);batch_normalize<<<dim3(2198,P),256>>>(jp);refresh();bool any=false;for(auto& j:jobs)any|=j.dist_active;if(!any)break;
  }
  batch_pack_default<<<dim3(2198,P),256>>>(jp);batch_prefix<<<dim3(1,P),256>>>(jp,true);batch_scatter<<<dim3(2198,P),256>>>(jp,true);batch_sum<<<dim3(1,P),256>>>(jp,true);
  decide_outer<<<P,1>>>(jp);accept_defaults<<<dim3(2198,P),256>>>(jp);refresh();bool any=false;for(auto& j:jobs)any|=j.active;if(!any)break;
 }
 CUDA_CHECK(cudaDeviceSynchronize());auto solved=std::chrono::steady_clock::now();batch_profit<<<dim3(2198,P),256>>>(jp);std::cerr<<"HOWARD_BATCH_MS "<<hw_total<<"\n";CUDA_CHECK(cudaDeviceSynchronize());auto finish_start=std::chrono::steady_clock::now();
 std::vector<std::ostringstream> reports(P),errors(P);std::vector<int> statuses(P,0);std::vector<std::exception_ptr> exceptions(P);
 std::vector<double> finish_seconds(P),start_seconds(P);std::vector<int> worker_ids(P);std::atomic<int> next_job{0};std::vector<std::thread> threads;
 workers=std::min(workers,P);
 for(int worker=0;worker<workers;++worker)threads.emplace_back([&,worker](){
  for(;;){int c=next_job.fetch_add(1);if(c>=P)break;auto t=std::chrono::steady_clock::now();start_seconds[c]=std::chrono::duration<double>(t-finish_start).count();worker_ids[c]=worker;
   try{CUDA_CHECK(cudaSetDevice(0));statuses[c]=finish_job(*owners[c],jobs[c],c,reports[c],errors[c]);}catch(...){exceptions[c]=std::current_exception();}
   finish_seconds[c]=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
  }
 });
 for(auto& thread:threads)thread.join();
 double finish_wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-finish_start).count();int status=0;
 for(int c=0;c<P;++c){if(exceptions[c])std::rethrow_exception(exceptions[c]);status|=statuses[c];std::cout<<reports[c].str();std::cerr<<errors[c].str();}
 std::cout<<std::setprecision(17)<<"{\"family\":\"parallel_candidate_simulation\",\"workers\":"<<workers<<",\"finish_wall_seconds\":"<<finish_wall<<",\"candidate_finish_seconds\":[";
 for(int c=0;c<P;++c){if(c)std::cout<<",";std::cout<<finish_seconds[c];}std::cout<<"],\"candidate_start_seconds\":[";
 for(int c=0;c<P;++c){if(c)std::cout<<",";std::cout<<start_seconds[c];}std::cout<<"],\"worker_ids\":[";
 for(int c=0;c<P;++c){if(c)std::cout<<",";std::cout<<worker_ids[c];}std::cout<<"]}\n";
 std::cout<<std::setprecision(17)<<"{\"family\":\"native_candidate_axis_batch\",\"jobs\":"<<P<<",\"policy_seconds\":"<<std::chrono::duration<double>(solved-solve_start).count()<<",\"total_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"status\":"<<status<<"}\n";CUDA_CHECK(cudaEventDestroy(hs));CUDA_CHECK(cudaEventDestroy(he));return status;
 }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}

// One CUDA context serves successive optimizer generations. Each job still
// initializes and solves independently; only allocations/context are reused.
int main(int argc,char** argv) {
 if(argc<2 || std::string(argv[1])!="--serve")return run_batch(argc,argv);
 int workers=4;
 if(argc==4 && std::string(argv[2])=="--simulation-workers")workers=std::stoi(argv[3]);
 else if(argc!=2){std::cerr<<"Usage: --serve [--simulation-workers N]\n";return 1;}
 if(workers<1||workers>64)return 1;
 if(cudaFree(nullptr)!=cudaSuccess)return 1;
 std::cout<<"{\"family\":\"ready\",\"protocol\":1}\n"<<std::flush;
 std::string jobs;int sequence=0;
 while(std::getline(std::cin,jobs)) {
   if(jobs=="QUIT")break;
   std::string count=std::to_string(workers);
   std::vector<std::string> text={argv[0],"--batch",jobs,"--simulation-workers",count};
   std::vector<char*> args;for(auto& t:text)args.push_back(t.data());
   int status=run_batch(args.size(),args.data());
   std::cout<<"{\"family\":\"batch_complete\",\"sequence\":"<<sequence++<<",\"status\":"<<status<<"}\n"<<std::flush;
   if(status)return status;
 }
 return 0;
}
