#pragma once

#include "palisade.h"


void AutoIdComposeTensor(
 const uint32_t M,
 const vector<usint>& indices_src1,
 const vector<usint>& indices_src2,
 vector<usint>& indices_dest
) {

  indices_dest.clear();
  uint32_t n1 = indices_src1.size();
  uint32_t n2 = indices_src2.size();

  for (uint32_t i = 0; i < n1; ++i)
  {
    for (uint32_t j = 0; j < n2; ++j)
    {
      usint id = indices_src1[i] * indices_src2[j];
      id %= M;
      indices_dest.push_back(id);
    }
  }
}


void MergeAutoIds(
 const uint32_t M,
 const vector<vector<usint>>& base_indices,
 const uint32_t beg,
 const uint32_t end,
 vector<usint>& auto_indices_each_unroll)
{
  uint32_t num_pack = end - beg + 1;
  uint32_t num_node_bin_tree = 2 * num_pack - 1;

  vector<vector<usint>> merge_ids(num_node_bin_tree);

  for (uint32_t i = 0; i < num_pack; ++i)
    merge_ids[i] = base_indices[i+beg];
  
  
  uint32_t lim = num_pack - 1;
  for (uint32_t i = 0, j = 0; j < lim; i=i+2, ++j)
  {
//    cour << "j=" << j << endl;
//    cour << "ACcess to " << 2*j  << " " << 2*j + 1 << " " << num_pack + j << endl;
    AutoIdComposeTensor(M, merge_ids[2*j], merge_ids[2*j + 1], merge_ids[num_pack + j]);
  }

  auto_indices_each_unroll =  merge_ids.back();
}


template<class T>
void TraceKeyGen(
  lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& cc,
  const size_t numitr,
  T& keys,
  vector<usint>& index_list,
  bool is_evalsums = false)
{
  uint32_t N = cc->GetRingDimension();

  size_t diff = static_cast<size_t>(log2(N)) - numitr; 
  if (is_evalsums)
	{
    uint32_t M = 2 * N;
    uint32_t g = 5;
    index_list.push_back(g);
    for (usint i = N >> 1; i != (1UL << diff); i = i >> 1)
		{
		  g *= g;
			g %= M;
      index_list.push_back(g);
		}
		int32_t gsign = NativeInteger(5).ModInverse(M).ConvertToInt();
		g = gsign;
    index_list.push_back(g);
    auto eval_automorph_ks = cc->EvalAutomorphismKeyGen(keys.secretKey, index_list);
    cc->InsertEvalAutomorphismKey(eval_automorph_ks);
	  return;
	}

  index_list.push_back(N + 1); // from the current ring to the next
// N/2 + 1, N/4 + 1, ..., N/(N/2) + 1
  for (usint i = N >> 1; i != (1UL << diff); i = i >> 1) {
    index_list.push_back(i + 1); // ring recursively...
  }
  //index_list.push_back(2 * N - 1); // for negating things
  //EvkAut eval_automorph_ks;
  auto eval_automorph_ks = cc->EvalAutomorphismKeyGen(keys.secretKey, index_list);
  cc->InsertEvalAutomorphismKey(eval_automorph_ks);
}

// Pure Unroll Keygen
template<class T>
double EvalSumKeyGenUnroll(
 CryptoContext<DCRTPoly>& cc,
 T& keys,
 const size_t num_itr,
 const size_t h, // num unroll
 vector<vector<usint>>& vec_auto_indices_par,
 bool is_trace = true
)
{
  TimeVar v;
  uint32_t N = cc->GetRingDimension();
  uint32_t M = N << 1;
  const size_t logN = log2(N);
  const size_t num_loop = num_itr;
  // number of times to do pack base case
  size_t unroll_base = num_loop / h;
  size_t unroll_rem  = unroll_base + 1;
  size_t num_repeat_rem  = num_loop % h;
  size_t num_repeat_base = h - num_repeat_rem;
  cout << " #loop (" << num_loop << ") is divided into " << h
       << " chunks with base chunk size " 
       << unroll_base << " (repeating " << num_repeat_base << " times) and "
       << unroll_rem  << " (repeating " << num_repeat_rem << " times)" << endl;

  usint diff = logN - num_loop;
  vector<vector<usint>> base_indices;
  if (is_trace)
  {
    base_indices.push_back({1, N + 1}); // from the current ring to the next
    // N/2 + 1, N/4 + 1, ..., N/(N/2) + 1
    for (usint i = N >> 1; i != (1UL << diff); i = i >> 1)
      base_indices.push_back({1, i + 1}); // ring recursively...
  }
  else
  {
    usint g = 5;
    base_indices.push_back({1, g});
    for (usint i = N >> 1; i != (1UL << diff); i = i >> 1)
    {
      g *= g;
      g %= M;
      base_indices.push_back({1, g});
    }
  }

  vec_auto_indices_par.resize(h);
  uint32_t beg = 0;
  uint32_t end = 0;
  for (size_t i = 0; i < num_repeat_base; ++i)
  {
    beg = i * unroll_base;
    end = beg + (unroll_base - 1);
    MergeAutoIds(M, base_indices, beg, end, vec_auto_indices_par[i]);
  }

  // Ceil case
  uint32_t base = end + 1;
  for (size_t i = 0; i < num_repeat_rem; ++i)
  {
    beg = base + i * unroll_rem;
    end = beg + (unroll_rem - 1);
    MergeAutoIds(M, base_indices, beg, end, vec_auto_indices_par[num_repeat_base + i]);
  }

  cout << "Indices (including 1):[\n";
  for (auto vec: vec_auto_indices_par)
  {
    cout << "  [";
    for (auto a: vec)
      cout << a << " ";
    cout << "] #keys = " << vec.size() <<  endl;
  }

  // Drop "1" as it does nothing on automorpshim
  for (size_t j = 0; j != h; ++j)
    vec_auto_indices_par[j].erase(vec_auto_indices_par[j].begin());

  vector<usint> auto_all;
  cout << "\nIndices (w/o 1):[\n";
  for (auto vec: vec_auto_indices_par)
  {
    cout << "  [";
    for (auto a: vec)
    {
      cout << a << " ";
      auto_all.push_back(a);
    }

    cout << "] #keys = " << vec.size() <<  endl;
  }

  auto eval_automorph_ks = cc->EvalAutomorphismKeyGen(keys.secretKey, auto_all);
  cc->InsertEvalAutomorphismKey(eval_automorph_ks);
  return TOC(v);
}


void GenZmstar(uint32_t logN, vector<uint32_t>& odds, vector<int>& is_consumed)
{
  uint32_t N = 1UL << logN;
	uint32_t M = 2 * N;
	for (uint32_t i = 0; i < M; ++i)
	{
	  if (i % 2 == 1)
		{
		  odds.push_back(i);
			is_consumed[i] = 1;
		// init as false; true when actually used
		}
	}
	return;
}

size_t GenUnroll(
 const uint32_t logN,
 const uint32_t M,
 const size_t num_unroll,
 vector<vector<uint32_t>>& vec_auto_indices
)
{
  size_t num_loop = logN ;
  size_t base_chunk = num_loop / num_unroll;
  size_t rem = num_loop % num_unroll;
  //cout << " #loop (" << num_loop << ") for Z_{N/2} is divided into "
  //     << num_unroll << " chunks with base chunk size " << base_chunk
  //     << " and " << rem << " chunks are of size " << base_chunk + 1;

  vector<uint32_t> chunk_sizes(num_unroll, base_chunk);

  int cnt = 0;
  while (rem != 0)
  {
    ++chunk_sizes[cnt];
    --rem;
    ++cnt;
  }

  cout << endl;
  for (auto chunk: chunk_sizes)
     cout << chunk << " ";
  cout << endl;
  uint32_t mm = M -1; 

  vector<uint32_t> auto_all;
  uint32_t gen = 5;

  for (size_t j = 0; j != num_unroll; ++j)
  {
    vector<uint32_t> auto_indices{1, gen};
    // NOTE 1 never gets registered
    // Unrolling
    // Add gen^2 * (1, gen)
    for (size_t i = 1; i != chunk_sizes[j]; ++i)
    {
      //cout << "INNER LOOP i = " << i << endl;
      // power of gen
      uint32_t new_id = gen * gen;
      //cout << " Multiplying thins by gen = " << new_id  << endl;
      new_id &= mm;
      const vector<uint32_t> prev_indices = auto_indices;
      // Composition of the automorphisms
      for (auto prev_index: prev_indices)
      {
        auto id =  (prev_index * new_id);
        //cout << "->" << prev_index << "*" << new_id << "=" << id;
        id &= mm;
        //cout << "  mod " << M << "=" << id << endl;
        auto_indices.push_back(id);
      }
      gen = new_id;
    }
    // next base gen
    gen *= gen;
    gen &= M-1;

    vec_auto_indices.emplace_back(auto_indices);
  }
  // Drop "1" as it does nothing on automorpshim
  for (size_t j = 0; j != num_unroll; ++j)
    vec_auto_indices[j].erase(vec_auto_indices[j].begin());

  size_t numkeys = 0;
	for(auto v: vec_auto_indices)
	  numkeys += v.size();
	return numkeys;
}

double TraceSeqHKOffline(
  const size_t logN,
  const size_t num_itr,
  const uint32_t multDepth,
  const uint32_t scaleFactorBits,
  const uint32_t dnum,
  CryptoContext<DCRTPoly>& cc,
  LPKeyPair<DCRTPoly>& keys,
  vector<usint>& auto_indices_seq,
  bool is_evalsum = false
 )
{
  TimeVar t;
  TIC(t);
  uint32_t batchSize = 8;
  SecurityLevel securityLevel = HEStd_128_classic;
  RescalingTechnique rsTech = APPROXRESCALE;
  KeySwitchTechnique ksTech = HYBRID;

   cc =
      CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
         multDepth,
         scaleFactorBits,
         batchSize,
         securityLevel,
         1 << logN,
         rsTech,
         ksTech,
         dnum
   );

  uint32_t N = cc->GetRingDimension();  cout << "CKKS scheme is using ring dimension " << N << " logN=" << std::log2(N) << endl << endl;

  cc->Enable(ENCRYPTION);
  cc->Enable(LEVELEDSHE);
  cc->Enable(SHE);

  keys = cc->KeyGen();
  vector<int> indices;

  cout << "Trace keygen..." << endl;
  TraceKeyGen(cc, num_itr, keys, auto_indices_seq, is_evalsum);
  cout << "done.." << endl;

  return TOC(t);
}
