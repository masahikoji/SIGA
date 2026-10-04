// Original equal-weight absolute-range Pocock-Simon allocation.
// Output: int16 product sums [replicate, pair(01/02), profile].
// Independent streams per replicate and per sequence; no outcomes generated.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#else
#include <chrono>
static void omp_set_num_threads(int) {}
static double omp_get_wtime() {
 return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
#endif
static uint64_t mix64(uint64_t z) {
 z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
 z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
 return z^(z>>31);
}
static double unif(std::mt19937_64 &g) { return (g()>>11)*0x1.0p-53; }
static int sgn(int x) { return (x>0)-(x<0); }
static int assign(std::array<int,11>& state, int F, int s, double p, double u) {
 int b=sgn(state[0]);
 for(int j=0;j<F;++j) b+=sgn(state[1+2*j+((s>>j)&1)]);
 const double pr=(b<0?p:(b>0?1-p:0.5));
 const int z=(u<pr?1:-1);
 state[0]+=z;
 for(int j=0;j<F;++j) state[1+2*j+((s>>j)&1)]+=z;
 return z;
}
static int reference(std::array<int,11>& st,int F,int s,double p,double u){
 int pos=std::abs(st[0]+1),neg=std::abs(st[0]-1);
 for(int j=0;j<F;++j){int x=st[1+2*j+((s>>j)&1)];pos+=std::abs(x+1);neg+=std::abs(x-1);}
 int z=(u<(pos<neg?p:pos>neg?1-p:0.5)?1:-1);
 st[0]+=z;for(int j=0;j<F;++j)st[1+2*j+((s>>j)&1)]+=z;
 return z;
}
static void selftest(){
 std::mt19937_64 g(938475U);long checked=0;
 for(int F: {2,5})for(double p:{0.5,0.8,0.95})for(int b=0;b<50;++b){
  std::array<int,11> a{},r{};
  for(int i=0;i<2000;++i){int s=g()%(1<<F);double u=unif(g);
   if(assign(a,F,s,p,u)!=reference(r,F,s,p,u)||a!=r)throw std::runtime_error("rule mismatch");++checked;
  }
 }
 std::cout<<"checked_assignments="<<checked<<"\n";
}
int main(int argc,char**argv){try{
 if(argc==2&&std::string(argv[1])=="--self-test"){selftest();return 0;}
 if(argc!=8)throw std::runtime_error("usage: F n B p master_seed threads output.bin");
 int F=std::stoi(argv[1]),n=std::stoi(argv[2]),B=std::stoi(argv[3]);
 double p=std::stod(argv[4]);uint64_t seed=std::stoull(argv[5]);int threads=std::stoi(argv[6]);
 if((F!=2&&F!=5)||n<1||n>30000||B<2||!std::isfinite(p)||p<0.5||p>=1||threads<1)throw std::runtime_error("invalid input");
 const int J=1<<F;std::array<double,5> probs{{.5,.4,.3,.2,.1}};
 std::vector<double> cdf(J);double acc=0;
 for(int s=0;s<J;++s){double q=1;for(int j=0;j<F;++j)q*=((s>>j)&1)?probs[j]:1-probs[j];acc+=q;cdf[s]=acc;}cdf.back()=1;
 std::vector<int16_t> out(static_cast<size_t>(B)*2*J);
 omp_set_num_threads(threads);double start=omp_get_wtime();
 #pragma omp parallel for schedule(static)
 for(int r=0;r<B;++r){
  const uint64_t base=mix64(seed+UINT64_C(0x9e3779b97f4a7c15)*(static_cast<uint64_t>(r)+1));
  std::mt19937_64 gp(mix64(base+11)),g0(mix64(base+21)),g1(mix64(base+31)),g2(mix64(base+41));
  std::array<int,11> st0{},st1{},st2{};std::array<int,32> u{},v{};
  for(int i=0;i<n;++i){
   int s=std::upper_bound(cdf.begin(),cdf.end(),unif(gp))-cdf.begin();
   int z0=assign(st0,F,s,p,unif(g0)),z1=assign(st1,F,s,p,unif(g1)),z2=assign(st2,F,s,p,unif(g2));
   u[s]+=z0*z1;v[s]+=z0*z2;
  }
  auto off=static_cast<size_t>(r)*2*J;for(int s=0;s<J;++s){out[off+s]=static_cast<int16_t>(u[s]);out[off+J+s]=static_cast<int16_t>(v[s]);}
 }
 std::ofstream f(argv[7],std::ios::binary);if(!f)throw std::runtime_error("cannot open output");
 f.write(reinterpret_cast<char*>(out.data()),out.size()*sizeof(int16_t));if(!f)throw std::runtime_error("write failed");
 std::cout<<"F="<<F<<" n="<<n<<" B="<<B<<" p="<<p<<" master_seed="<<seed<<" threads="<<threads<<" seconds="<<omp_get_wtime()-start<<"\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
