#pragma once

#ifdef USE_FAST_AUTOMORPH
void GenAutomorphTable(const usint n, const usint k, std::vector<usint> &perm) {
  // usint n = GetRingDimension();
  usint m = n << 1;
  usint logn = log2(n);
  usint logm = logn + 1;
  perm.resize(n);
  for (usint j = 1; j < m; j += 2) {
    usint idx = (j * k) - (((j * k) >> logm) << logm);
    usint jrev = ReverseBits(j >> 1, logn);
    usint idxrev = ReverseBits(idx >> 1, logn);
    // result.m_values->operator[](jrev) = GetValues().operator[](idxrev);
    perm[idxrev] = jrev;
  }
}
#endif

void PreComputeEvk(
    const CryptoContext<DCRTPoly> &cc, const usint autoIndex,
    ConstCiphertext<DCRTPoly> ciphertext,
    map<usint, std::pair<std::vector<DCRTPoly>, std::vector<DCRTPoly>>> &evks) {
  uint32_t n = cc->GetRingDimension();
  uint32_t m = n << 1;

  // 1. Compute autoIndex^{-1} mod m
  NativeInteger index_int(to_string(autoIndex));
  auto inv_autoid = index_int.ModInverse(m).ConvertToInt();
  //  cout << "Autoid = ," << autoIndex <<  "InvAuto = " << inv_autoid << " mod
  //  " << m << endl;

  // 2. Retrieve the automorphism key that corresponds to the auto index.
  auto evalKey = ciphertext->GetCryptoContext()
                     ->GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag())
                     .find(autoIndex)
                     ->second;

  // Get the parts of the automorphism key
  std::vector<DCRTPoly> b = evalKey->GetBVector();
  std::vector<DCRTPoly> a = evalKey->GetAVector();

  // 3. Automorphism by its inverse
  for (size_t i = 0; i < b.size(); i++)
    b[i] = b[i].AutomorphismTransform(inv_autoid);

  for (size_t i = 0; i < a.size(); i++)
    a[i] = a[i].AutomorphismTransform(inv_autoid);

  evks[autoIndex] = make_pair(b, a);
}
