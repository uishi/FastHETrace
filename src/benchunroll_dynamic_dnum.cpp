// Define PROFILE to enable TIC-TOC timing measurements
#define PROFILE

#include "palisade.h"

#include <cassert>
#include <fstream>
#include <omp.h>
using namespace std;

using namespace lbcrypto;
using PolyType = DCRTPoly;
using Ctxt = Ciphertext<PolyType>;
using EvkAut = shared_ptr<map<usint, LPEvalKey<PolyType>>>;
using FHEContext = CryptoContext<PolyType>;

#include "stat.hpp"
#include "mem_usage.h"
#include "trace_keygen.hpp"
#include "automo_opt.hpp"
#include "hoist.hpp"

double EvalTraceUnrollHKOnline(
 const size_t M,
 const CryptoContext<DCRTPoly> cc,
 const LPKeyPair<DCRTPoly>& keys,
 const uint32_t   L,
 const size_t num_unroll,
 const size_t num_mod_reduce,
 vector<vector<usint>>& auto_indices_parallel,
 map<usint, std::pair< std::vector<DCRTPoly>, std::vector<DCRTPoly> >>& inv_evks
)
{
  // Input
  vector<complex<double>> x = { 0, 0, 0, 0, 0, 0, 0, 1 };
  Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x);

  cout << "Input x: " << ptxt << endl;
  auto c = cc->Encrypt(keys.publicKey, ptxt);
  vector<complex<double>> x1(cc->GetRingDimension()/2, 1);
	auto c1 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(x1));
	for (size_t i = 0; i < num_mod_reduce; ++i)
	{
   c = cc->EvalMult(c, c1);
	 c = cc->Rescale(c);
	}
  cout << "# Dropped = " << c->GetLevel()  << " # MODs Available = " << (L+1) - c->GetLevel()<< endl;

  // First, we perform 7 regular (non-hoisted) rotations
  // and measure the runtime.
  TimeVar t;
  TIC(t);

  for (size_t numitr = 0; numitr < num_unroll; ++numitr)
	{
    size_t num_par = auto_indices_parallel[numitr].size();
    // Hoist out expensive
    vector< shared_ptr<vector<DCRTPoly>> > c_rot_hoisted(num_par);
  
    // Hoisting
    auto cPrecomp = cc->EvalFastRotationPrecompute(c);

    auto c_tmp = c->Clone();
    #pragma omp parallel for
    for (size_t j = 0; j < num_par; ++j)
      c_rot_hoisted[j] = EvalIPAndAutomorph(c, auto_indices_parallel[numitr][j], cPrecomp, inv_evks);

    // Post-Hoisting
    shared_ptr<vector<DCRTPoly>> merge_sum =  EvalAddLazyComponents(c_rot_hoisted);

    // ModDown and addit to the original
    c += EvalModDownIntermediate(c_tmp, merge_sum);
	}
  double timeHoisting = TOC(t);

  Plaintext result;
  cc->Decrypt(keys.secretKey, c, &result);
  result->SetLength(5);
  cout << "\t Result = " << result << endl;
  cout << "Timing = " << timeHoisting << " ms...";

	auto m = cc->GetCyclotomicOrder();

	size_t num_ans = 0;
	int ans = 1;
	size_t num_expected_ans = 1UL << M;
  if (M == size_t(log2(m/2)))
	{
	  ans = 2;
	  num_expected_ans = (1UL << (M-1));
	}

  // See all the slot values
	for (size_t i = 0; i <m/4; ++i)
	{
	  if (round(result->GetCKKSPackedValue()[i].real()) == ans)
		  num_ans += 1;
	}

  if (num_ans == num_expected_ans)
    cout << " correct :) #1 = "  << num_ans  << " Expected = " << num_expected_ans << endl;
  else
    cout << "      wrong #1 = " << num_ans << " Expected = " << num_expected_ans <<endl;

  return timeHoisting;
}

double EvalTraceUnrollHKMult(
 const size_t M,
 const CryptoContext<DCRTPoly> cc,
 const LPKeyPair<DCRTPoly>& keys,
 const uint32_t   L,
 const size_t num_unroll,
 const size_t num_mod_reduce,
 vector<vector<usint>>& auto_indices_parallel,
 map<usint, std::pair< std::vector<DCRTPoly>, std::vector<DCRTPoly> >>& inv_evks
)
{
  // Input
  vector<complex<double>> x = { 0, 0, 0, 0, 0, 0, 0, 1 };
  Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x);

  cout << "Input x: " << ptxt << endl;
  auto c = cc->Encrypt(keys.publicKey, ptxt);
  vector<complex<double>> x1(cc->GetRingDimension()/2, 1);
	auto c1 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(x1));
	for (size_t i = 0; i < num_mod_reduce; ++i)
	{
   c = cc->EvalMult(c, c1);
	 c = cc->Rescale(c);
	}
  cout << "# Dropped = " << c->GetLevel()  << " # MODs Available = " << (L+1) - c->GetLevel()<< endl;

  TimeVar t;
  TIC(t);

  for (size_t numitr = 0; numitr < num_unroll; ++numitr)
      c += EvalHoistedAutomorph(c, auto_indices_parallel[numitr], inv_evks);
  double timeHoisting = TOC(t);

  Plaintext result;
  cc->Decrypt(keys.secretKey, c, &result);
  result->SetLength(5);
  cout << "\t Result = " << result << endl;
  cout << "Timing = " << timeHoisting << " ms...";

	auto m = cc->GetCyclotomicOrder();

	size_t num_ans = 0;
	int ans = 1;
	size_t num_expected_ans = 1UL << M;
  if (M == size_t(log2(m/2)))
	{
	  ans = 2;
	  num_expected_ans = (1UL << (M-1));
	}

  // See all the slot values
	for (size_t i = 0; i <m/4; ++i)
	{
	  if (round(result->GetCKKSPackedValue()[i].real()) == ans)
		  num_ans += 1;
	}

  if (num_ans == num_expected_ans)
    cout << " correct :) #1 = "  << num_ans  << " Expected = " << num_expected_ans << endl;
  else
    cout << "      wrong #1 = " << num_ans << " Expected = " << num_expected_ans <<endl;

  return timeHoisting;
}

double EvalTraceUnrollSingleTh(
 const size_t M,
 const CryptoContext<DCRTPoly> cc,
 const LPKeyPair<DCRTPoly>& keys,
 const uint32_t   L,
 const size_t num_unroll,
 const size_t num_mod_reduce,
 vector<vector<usint>>& auto_indices_parallel,
 map<usint, std::pair< std::vector<DCRTPoly>, std::vector<DCRTPoly> >>& inv_evks
)
{
  // Input
  vector<complex<double>> x = { 0, 0, 0, 0, 0, 0, 0, 1 };
  Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x);

  cout << "Input x: " << ptxt << endl;
  auto c = cc->Encrypt(keys.publicKey, ptxt);
  vector<complex<double>> x1(cc->GetRingDimension()/2, 1);
	auto c1 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(x1));
	for (size_t i = 0; i < num_mod_reduce; ++i)
	{
   c = cc->EvalMult(c, c1);
	 c = cc->Rescale(c);
	}
  cout << "# Dropped = " << c->GetLevel()  << " # MODs Available = " << (L+1) - c->GetLevel()<< endl;


  // First, we perform 7 regular (non-hoisted) rotations
  // and measure the runtime.
  TimeVar t;
  TIC(t);

  for (size_t numitr = 0; numitr < num_unroll; ++numitr)
    c += EvalHoistedAutomorphHKSingleTh(c, auto_indices_parallel[numitr], inv_evks);

  double timeHoisting = TOC(t);

  Plaintext result;
  cc->Decrypt(keys.secretKey, c, &result);
  result->SetLength(5);
  cout << "\t Result = " << result << endl;
  cout << "Timing = " << timeHoisting << " ms...";

	auto m = cc->GetCyclotomicOrder();

	size_t num_ans = 0;
	int ans = 1;
	size_t num_expected_ans = 1UL << M;
  if (M == size_t(log2(m/2)))
	{
	  ans = 2;
	  num_expected_ans = (1UL << (M-1));
	}

  // See all the slot values
	for (size_t i = 0; i <m/4; ++i)
	{
	  if (round(result->GetCKKSPackedValue()[i].real()) == ans)
		  num_ans += 1;
	}

  if (num_ans == num_expected_ans)
    cout << " correct :) #1 = "  << num_ans  << " Expected = " << num_expected_ans << endl;
  else
    cout << "      wrong #1 = " << num_ans << " Expected = " << num_expected_ans <<endl;

  return timeHoisting;
}



template<class T>
void ShowVec (const std::vector<T>& v) {
  std::size_t e = v.size();
	std::cout << "[";
	for (std::size_t i = 0, e = v.size(); i < e-1; ++i) {
	  std::cout << v[i] << ",";
	}
	std::cout << v[e-1] << "]" << std::endl;
}



int main(int argc, char* argv[])
{
#ifdef USE_FAST_AUTOMORPH
  cout << "Fast Automorph" << endl;
#endif

  int num_thread = 1;
  // For some reason, the number of threads is shown only in the parallel section
  #pragma omp parallel
  {
    #pragma omp single
    num_thread = omp_get_num_threads();
  }
  // try to unroll evenly distributed
	// In case of evalsum, logN-1/#Unroll
  if (argc != 9)
  {
     cerr << "Usage: ./xxx logN M #Unroll depth ell scale dnum  numexp" << endl;
		 cerr << "  # Argc (Expected) = 9" << endl;
		 cerr << "  # Argc (Actual)   = "<< argc << endl;
     return 1;
  }

  size_t      logN  = stoi(argv[1]);
  size_t         M  = stoi(argv[2]);
  size_t num_unroll = stoi(argv[3]);
  uint32_t   depth  = stoi(argv[4]);
  uint32_t     ell  = stoi(argv[5]);
  uint32_t   scale  = stoi(argv[6]);
  uint32_t    dnum  = stoi(argv[7]);
  size_t    numExp  = stoi(argv[8]);

  assert(logN >= M);
  assert(M >= num_unroll);
  assert(depth >= ell);

  size_t num_mod_reduce = depth - ell;

  SecurityLevel securityLevel = HEStd_128_classic;
  RescalingTechnique rsTech = APPROXRESCALE;
  KeySwitchTechnique ksTech = HYBRID;

  CryptoContext<DCRTPoly> cc;
  cc =
      CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
         depth,
         scale,
         1 << (logN - 1),
         securityLevel,
         1 << logN,
         rsTech,
         ksTech,
         dnum
  );

  cc->Enable(ENCRYPTION);
  cc->Enable(SHE);
  cc->Enable(LEVELEDSHE);
  vector<int> indices;
	LPKeyPair<DCRTPoly> keys = cc->KeyGen();
	string outfile;

  // evk
  map<usint, std::pair< std::vector<DCRTPoly>, std::vector<DCRTPoly> >> inv_evks;
  // ctxt
  vector<complex<double>> x = {0};
  Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x);
  auto cipher = cc->Encrypt(keys.publicKey, ptxt);
 
  
  vector<vector<usint>> auto_indices_par;

  auto time_offline = EvalSumKeyGenUnroll(cc, keys, M, num_unroll, auto_indices_par);
  cc->EvalMultKeyGen(keys.secretKey);

  for(auto vec: auto_indices_par)
    for (auto id: vec)
      PreComputeEvk(cc, id, cipher, inv_evks);

	{
    string off_outfile = "/tmp/result/LweRlwePalisadeKeyGen/evalsum_AutoOpt_unroll_div_optimize_HK_logN" + to_string(logN)
      + "_M"  + to_string(M)
      + "_depth"  + to_string(depth)
      + "_numunroll"   + to_string(num_unroll)
      + "_scale"  + to_string(scale)
      + "_numth"  + to_string(num_thread)
      + "_dnum"   + to_string(dnum)
      + "_numexp" + to_string(numExp)
      + ".txt";
		ofstream offofs(off_outfile);
		offofs << time_offline << "," << ((double)(getPeakRSS()) / (1 << 20)) << endl;
	}
	int numkeys = 0;
	for (auto key_vec:auto_indices_par)
	{
	  cout << "["<< key_vec.size() <<",";
		numkeys += key_vec.size();
	}
	cout << "]";
	cout << "# keys =" << numkeys << endl;

	vector<double> timings;
	if (num_thread == 1)
	{
    for (size_t i = 0; i < numExp; ++i)
      timings.push_back(EvalTraceUnrollSingleTh(M, cc, keys, depth, num_unroll, num_mod_reduce, auto_indices_par, inv_evks));
	}
	else
	{
    for (size_t i = 0; i < numExp; ++i)
      timings.push_back(EvalTraceUnrollHKMult(M, cc, keys, depth, num_unroll, num_mod_reduce, auto_indices_par, inv_evks));
	}


  
	int num_special_mod = (depth + 1)/dnum;

  string outprefix = "";
  outfile = "/tmp/result/TraceRuntime/" + outprefix+ "dynamic_dnum_unroll_HK_logN" + to_string(logN)
    + "_M"           + to_string(M)
    + "_L"           + to_string(depth)
    + "_ell"         + to_string(ell)
    + "_numunroll"   + to_string(num_unroll)
    + "_scale"       + to_string(scale)
    + "_numth"       + to_string(num_thread)
    + "_k"           + to_string(num_special_mod)
    + "_numexp"      + to_string(numExp)
    + ".txt";
  

 StatDropFirst s_unroll(timings, { (double)(getPeakRSS()) / (1 << 20) });
 s_unroll.Write(cout);
 ofstream ofs(outfile);
 s_unroll.Write(ofs);
 ofs.close();
 
	return 0;
}
