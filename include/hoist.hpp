#pragma once

Ciphertext<DCRTPoly> EvalHoistedAutomorphHKSingleTh(
    ConstCiphertext<DCRTPoly> ciphertext, // original ciphertext to  rotate
		const vector<usint>& autoids, 
    map<usint, std::pair< std::vector<DCRTPoly>, std::vector<DCRTPoly> >>& inv_evks
)
{

  Ciphertext<DCRTPoly> newCiphertext = ciphertext->CloneEmpty();

  const shared_ptr<vector<DCRTPoly>> expandedCiphertext = ciphertext->GetCryptoContext()->EvalFastRotationPrecompute(ciphertext);

  // 1. Retrieve the automorphism key that corresponds to the auto index.
  std::vector<DCRTPoly> c(2);

  c[0] = ciphertext->GetElements()[0].Clone();
  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
	      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(ciphertext->GetCryptoParameters()); 

  const shared_ptr<typename DCRTPoly::Params> paramsQ = ciphertext->GetElements()[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP = cryptoParamsLWE->GetAuxElementParams();
  size_t cipherTowers = paramsQ->GetParams().size();
  size_t towersToSkip = cryptoParamsLWE->GetElementParams()->GetParams().size() - cipherTowers;


  size_t n = autoids.size();
  DCRTPoly c_sum(paramsQ, Format::EVALUATION, true);
  DCRTPoly cTilda0_sum((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);
  DCRTPoly cTilda1_sum((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);
  for (size_t ii = 0; ii < n; ++ii)
	{
    const usint autoIndex = autoids[ii];

    auto evalKey = ciphertext->GetCryptoContext()->
      GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag()).find(autoIndex)->second;
#ifdef GET_RAM
    const std::vector<DCRTPoly> &b = evalKey->GetBVector();
    const std::vector<DCRTPoly> &a = evalKey->GetAVector();
#else
    const std::vector<DCRTPoly> &b = inv_evks[autoIndex].first;
    const std::vector<DCRTPoly> &a = inv_evks[autoIndex].second;
#endif

    DCRTPoly ct0;
    DCRTPoly ct1;

  // 1. Retrieve the automorphism key that corresponds to the auto index.

  // Inner Prod with Evk
    // Now, we are preparing the resulting storage mod PQ (filled with 0)
    DCRTPoly cTilda0((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);
    DCRTPoly cTilda1((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);
  
    for (uint32_t j=0; j<expandedCiphertext->size(); j++) {
      for (usint i=0; i< (*expandedCiphertext)[j].GetNumOfElements(); i++) {
        usint idx = ( i < cipherTowers ) ? i : i + towersToSkip;
        cTilda0.SetElementAtIndex(i, cTilda0.GetElementAtIndex(i) + (*expandedCiphertext)[j].GetElementAtIndex(i) * b[j].GetElementAtIndex(idx));
  			cTilda1.SetElementAtIndex(i, cTilda1.GetElementAtIndex(i) + (*expandedCiphertext)[j].GetElementAtIndex(i) * a[j].GetElementAtIndex(idx));    
  		 }
  	}

#ifdef TRACE_DEBUG
    // We see that the last k RNS components are relative to P
		uint32_t num_digits = expandedCiphertext->size();
	  cout << " k                             = "  << paramsP->GetParams().size() << endl;
	  cout << " cipherTows (ell+1)            = "  << cipherTowers << endl;
	  cout << " k + ell + 1                   = "  << cipherTowers + paramsP->GetParams().size() << endl;
	  cout << " TowsToSkip (#Drop = L - ell ) = "  << towersToSkip << endl;
	  cout << " # Digtit                      = "  << num_digits << endl;
	  cout << " # RSN on EVK                  = "  << towersToSkip + cipherTowers + paramsP->GetParams().size() << endl;
    for (uint32_t j=0; j< num_digits; j++) {
	    cout << "\tDigit No: "  << j<< endl;
			usint num_rns_on_digit = (*expandedCiphertext)[j].GetNumOfElements();
	    cout << "\t\t #RNS Components on digit: "  << num_rns_on_digit  << endl;
      for (usint i=0; i< num_rns_on_digit ; i++) {
        usint idx = ( i < cipherTowers ) ? i : i + towersToSkip;
	      cout << "\t\t EVk RNS-id to refer     : "  << idx << endl;
			}
		}
#endif

#ifdef USE_FAST_AUTOMORPH
    vector<usint> perm;
    GenAutomorphTable(cTilda0.GetRingDimension(), autoIndex, perm);
    c_sum       += std::move(ciphertext->GetElements()[0].Permute(perm));
    cTilda0_sum += std::move(cTilda0.Permute(perm));
    cTilda1_sum += std::move(cTilda1.Permute(perm));
#else
    c_sum       += std::move(ciphertext->GetElements()[0].AutomorphismTransform(autoIndex));
    cTilda0_sum += std::move(cTilda0.AutomorphismTransform(autoIndex));
    cTilda1_sum += std::move(cTilda1.AutomorphismTransform(autoIndex));
#endif
  }

  cTilda0_sum.SetFormat(Format::COEFFICIENT);
  cTilda1_sum.SetFormat(Format::COEFFICIENT);

  DCRTPoly cHat0 = cTilda0_sum.ApproxModDown(
      paramsQ,
      paramsP,
      cryptoParamsLWE->GetPInvModQTable(),        // P^{-1} mod q_i
      cryptoParamsLWE->GetPInvModQPreconTable(),  // Barrett Const to multiply with P^{-1} mod q_i
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  DCRTPoly cHat1 = cTilda1_sum.ApproxModDown(
      paramsQ,
      paramsP,
      cryptoParamsLWE->GetPInvModQTable(),
      cryptoParamsLWE->GetPInvModQPreconTable(),
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  cHat0.SetFormat(Format::EVALUATION);
  cHat1.SetFormat(Format::EVALUATION);

  DCRTPoly ct0;
  DCRTPoly ct1;
  ct0 = c_sum + cHat0;
  ct1 = cHat1;


  newCiphertext->SetElements({ ct0, ct1 });
  newCiphertext->SetDepth(ciphertext->GetDepth());
  newCiphertext->SetLevel(ciphertext->GetLevel());
  newCiphertext->SetScalingFactor(ciphertext->GetScalingFactor());

  return newCiphertext;
}

Ciphertext<DCRTPoly> EvalHoistedAutomorph(
    ConstCiphertext<DCRTPoly> ciphertext, // original ciphertext to  rotate
		const vector<usint>& autoids, 
    map<usint, std::pair< std::vector<DCRTPoly>, std::vector<DCRTPoly> >>& inv_evks
)
{

  Ciphertext<DCRTPoly> newCiphertext = ciphertext->CloneEmpty();

  const shared_ptr<vector<DCRTPoly>> expandedCiphertext = ciphertext->GetCryptoContext()->EvalFastRotationPrecompute(ciphertext);

  // 1. Retrieve the automorphism key that corresponds to the auto index.
  std::vector<DCRTPoly> c(2);

  c[0] = ciphertext->GetElements()[0].Clone();
  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
	      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(ciphertext->GetCryptoParameters()); 

  const shared_ptr<typename DCRTPoly::Params> paramsQ = ciphertext->GetElements()[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP = cryptoParamsLWE->GetAuxElementParams();
  const shared_ptr<typename DCRTPoly::Params> paramsPQ = (*expandedCiphertext)[0].GetParams();
  size_t cipherTowers = paramsQ->GetParams().size();
  size_t towersToSkip = cryptoParamsLWE->GetElementParams()->GetParams().size() - cipherTowers;


  size_t n = autoids.size();
  vector<DCRTPoly>   c_vec    (n);
  vector<DCRTPoly> cTilda0_vec(n);
  vector<DCRTPoly> cTilda1_vec(n);
	#pragma omp parallel for
  for (size_t ii = 0; ii < n; ++ii)
	{
    const usint autoIndex = autoids[ii];

    auto evalKey = ciphertext->GetCryptoContext()->
      GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag()).find(autoIndex)->second;
#ifdef GET_RAM
    const std::vector<DCRTPoly> &b = evalKey->GetBVector();
    const std::vector<DCRTPoly> &a = evalKey->GetAVector();
#else
    const std::vector<DCRTPoly> &b = inv_evks[autoIndex].first;
    const std::vector<DCRTPoly> &a = inv_evks[autoIndex].second;
#endif

    DCRTPoly ct0;
    DCRTPoly ct1;

  // 1. Retrieve the automorphism key that corresponds to the auto index.

  // Inner Prod with Evk
    // Now, we are preparing the resulting storage mod PQ (filled with 0)
    DCRTPoly cTilda0(paramsPQ, Format::EVALUATION, true);
    DCRTPoly cTilda1(paramsPQ, Format::EVALUATION, true);
  
    for (uint32_t j=0; j<expandedCiphertext->size(); j++) {
      for (usint i=0; i< (*expandedCiphertext)[j].GetNumOfElements(); i++) {
        usint idx = ( i < cipherTowers ) ? i : i + towersToSkip;
        cTilda0.SetElementAtIndex(i, cTilda0.GetElementAtIndex(i) + (*expandedCiphertext)[j].GetElementAtIndex(i) * b[j].GetElementAtIndex(idx));
  			cTilda1.SetElementAtIndex(i, cTilda1.GetElementAtIndex(i) + (*expandedCiphertext)[j].GetElementAtIndex(i) * a[j].GetElementAtIndex(idx));    
  		 }
  	}

#ifdef USE_FAST_AUTOMORPH
    vector<usint> perm;
    GenAutomorphTable(cTilda0.GetRingDimension(), autoIndex, perm);
    c_vec      [ii] = std::move(ciphertext->GetElements()[0].Permute(perm));
    cTilda0_vec[ii] = std::move(cTilda0.Permute(perm));
    cTilda1_vec[ii] = std::move(cTilda1.Permute(perm));
#else
    c_vec      [ii] = std::move(ciphertext->GetElements()[0].AutomorphismTransform(autoIndex));
    cTilda0_vec[ii] = std::move(cTilda0.AutomorphismTransform(autoIndex));
    cTilda1_vec[ii] = std::move(cTilda1.AutomorphismTransform(autoIndex));
#endif
  }
  DCRTPoly       c_sum(paramsQ, Format::EVALUATION, true);
  DCRTPoly cTilda0_sum(paramsPQ, Format::EVALUATION, true);
  DCRTPoly cTilda1_sum(paramsPQ, Format::EVALUATION, true);

#if 0
  #pragma omp parallel for reduction(+:c_sum)
  for (size_t ii = 0; ii < n; ++ii)
	   c_sum += c_vec[ii];

  #pragma omp parallel for reduction(+:cTilda0_sum, cTilda1_sum)
  for (size_t ii = 0; ii < n; ++ii)
	{
	   cTilda0_sum += cTilda0_vec[ii];
	   cTilda1_sum += cTilda1_vec[ii];
  }
#else
  for (size_t ii = 0; ii < n; ++ii)
	{
	   c_sum += c_vec[ii];
	   cTilda0_sum += cTilda0_vec[ii];
	   cTilda1_sum += cTilda1_vec[ii];
  }
#endif


  cTilda0_sum.SetFormat(Format::COEFFICIENT);
  cTilda1_sum.SetFormat(Format::COEFFICIENT);

  DCRTPoly cHat0 = cTilda0_sum.ApproxModDown(
      paramsQ,
      paramsP,
      cryptoParamsLWE->GetPInvModQTable(),        // P^{-1} mod q_i
      cryptoParamsLWE->GetPInvModQPreconTable(),  // Barrett Const to multiply with P^{-1} mod q_i
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  DCRTPoly cHat1 = cTilda1_sum.ApproxModDown(
      paramsQ,
      paramsP,
      cryptoParamsLWE->GetPInvModQTable(),
      cryptoParamsLWE->GetPInvModQPreconTable(),
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  cHat0.SetFormat(Format::EVALUATION);
  cHat1.SetFormat(Format::EVALUATION);

  DCRTPoly ct0;
  DCRTPoly ct1;
  ct0 = c_sum + cHat0;
  ct1 = cHat1;


  newCiphertext->SetElements({ ct0, ct1 });
  newCiphertext->SetDepth(ciphertext->GetDepth());
  newCiphertext->SetLevel(ciphertext->GetLevel());
  newCiphertext->SetScalingFactor(ciphertext->GetScalingFactor());

  return newCiphertext;
}

inline shared_ptr<vector<DCRTPoly>> EvalIPAndAutomorph(
    ConstCiphertext<DCRTPoly> ciphertext, // original ciphertext to do rotate
    const usint autoIndex,
    const shared_ptr<vector<DCRTPoly>> expandedCiphertext,
    map<usint, std::pair< std::vector<DCRTPoly>, std::vector<DCRTPoly> >>& inv_evks
    // Decomposed Ciphertext mod PQ where each DCRTPoly corresponds to a decomposed digit    // to rotate
)
{
  // 1. Retrieve the automorphism key that corresponds to the auto index.
  auto evalKey = ciphertext->GetCryptoContext()->
      GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag()).      find(autoIndex)->second;
  std::vector<DCRTPoly> c(2);

  c[0] = ciphertext->GetElements()[0].Clone();

  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(evalKey->GetCryptoParameters());

  //Ciphertext<DCRTPoly> newCiphertext = ciphertext->CloneEmpty();

  // Eval Keys
  // Get the parts of the automorphism key
  std::vector<DCRTPoly> b = inv_evks[autoIndex].first;
  std::vector<DCRTPoly> a = inv_evks[autoIndex].second;

  DCRTPoly ct0;
  DCRTPoly ct1;

  const shared_ptr<typename DCRTPoly::Params> paramsQ = c[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP = cryptoParamsLWE->GetAuxElementParams();

  size_t cipherTowers = c[0].GetParams()->GetParams().size();
  size_t towersToSkip = cryptoParamsLWE->GetElementParams()->GetParams().size() - cipherTowers;

  // 1. Retrieve the automorphism key that corresponds to the auto index.

// Inner Prod with Evk
  // Now, we are preparing the resulting storage mod PQ (filled with 0)
  DCRTPoly cTilda0((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);
  DCRTPoly cTilda1((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);

  for (uint32_t j=0; j<expandedCiphertext->size(); j++) {
    for (usint i=0; i< (*expandedCiphertext)[j].GetNumOfElements(); i++) {
      usint idx = ( i < cipherTowers ) ? i : i + towersToSkip;
      cTilda0.SetElementAtIndex(i, cTilda0.GetElementAtIndex(i) + (*expandedCiphertext)[j].GetElementAtIndex(i) * b[j].GetElementAtIndex(idx));
			cTilda1.SetElementAtIndex(i, cTilda1.GetElementAtIndex(i) + (*expandedCiphertext)[j].GetElementAtIndex(i) * a[j].GetElementAtIndex(idx));    
		 }
	}

#ifdef USE_FAST_AUTOMORPH
  vector<usint> perm;
  GenAutomorphTable(cTilda0.GetRingDimension(), autoIndex, perm);
  c[0]    = ciphertext->GetElements()[0].Permute(perm);
  cTilda0 = cTilda0.Permute(perm);
  cTilda1 = cTilda1.Permute(perm);
#else
  c[0]    = ciphertext->GetElements()[0].AutomorphismTransform(autoIndex);
  cTilda0 = cTilda0.AutomorphismTransform(autoIndex);
  cTilda1 = cTilda1.AutomorphismTransform(autoIndex);
#endif


//  vector<DCRTPoly> result{c[0], cTilda0, cTilda1};
  shared_ptr<vector<DCRTPoly>> resultPtr = make_shared<vector<DCRTPoly>>(
	 vector<DCRTPoly>{c[0], cTilda0, cTilda1});

  return resultPtr;
}

inline shared_ptr<vector<DCRTPoly>> EvalFastRotationPrecomputeHybridOpt(ConstCiphertext<DCRTPoly> ciphertext)
{
  Ciphertext<DCRTPoly> newCiphertext = ciphertext->CloneEmpty();

  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(ciphertext->GetCryptoParameters());

  const std::vector<DCRTPoly> &c = ciphertext->GetElements();

  DCRTPoly ct0;
  DCRTPoly ct1;

  const shared_ptr<typename DCRTPoly::Params> paramsQ = c[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP = cryptoParamsLWE->GetAuxElementParams();

  size_t cipherTowers = c[0].GetParams()->GetParams().size();

  DCRTPoly cTmp = c[1].Clone();

  uint32_t l = cipherTowers - 1;
  uint32_t alpha = cryptoParamsLWE->GetNumberOfTowersPerDigit();
  uint32_t beta = ceil(((double)(l+1))/alpha); // The number of digits of the current ciphertext
  if (beta > cryptoParamsLWE->GetNumberOfQPartitions())
    beta = cryptoParamsLWE->GetNumberOfQPartitions();

  vector<DCRTPoly> digitsCTmp(beta);

  // Digit decomposition
  // Zero-padding and split
  uint32_t numTowersLastDigit = cryptoParamsLWE->GetQPartition(beta-1)->GetParams().size();
  for (uint32_t j=0; j<beta; j++) {
    cout << "****Decomposed Index ="<< j << endl;
    if (j == beta-1) {
      auto part = cryptoParamsLWE->GetQPartition(j);
      part->GetParams();

      numTowersLastDigit = cipherTowers - alpha*j;

      vector<NativeInteger> moduli(numTowersLastDigit);
      vector<NativeInteger> roots(numTowersLastDigit);

      for (uint32_t i=0; i<numTowersLastDigit; i++) {
        moduli[i] = part->GetParams()[i]->GetModulus();
        roots[i] = part->GetParams()[i]->GetRootOfUnity();
      }

      auto params = DCRTPoly::Params(part->GetCyclotomicOrder(),
                    moduli, roots, {}, {}, 0);

      digitsCTmp[j] = DCRTPoly(std::make_shared<typename DCRTPoly::Params>(params), EVALUATION, true);

    } else
      digitsCTmp[j] = DCRTPoly(cryptoParamsLWE->GetQPartition(j), Format::EVALUATION, true);

    cout << "Mod Size " << j << " th Decomp:" << log2(digitsCTmp[j].GetModulus().ConvertToDouble()) <<endl;

    uint32_t iters = (j == beta-1) ? numTowersLastDigit : alpha;
    for (uint32_t i=0; i<iters; i++) {
    cout << " beta="<< beta << ","
             << "#itr=" << iters << ", "
             << " j*alpha + i " << j * alpha + i << ", "
             << " l " << l<< ", "
             << " RNS Capacity " << l+1
             << endl;
      if (j*alpha + i <= l) {
        auto tmp = cTmp.GetElementAtIndex(j*alpha+i);
        digitsCTmp[j].SetElementAtIndex(i, tmp);
    cout << "  Set Decomp[" << j << "].Poly[" << i << "]"
         << "=" << "tmp.Poly[" <<  j*alpha + i << "] where it has full of RNS"<< endl;
      }
    }
  }

  // RNS decompose
  // j: decomposed stuff
  for (uint32_t j=0; j<beta; j++) {
    for (uint32_t i=0; i<alpha; i++) {
      if (j*alpha + i <= l) {
        auto tmp = digitsCTmp[j].GetElementAtIndex(i).Times(cryptoParamsLWE->GetQHatInvModqTable()[j][j*alpha+i]);
        digitsCTmp[j].SetElementAtIndex(i, tmp);
    cout << " ==> Update " << j << "-th decomp at " << i << "-th poly" << " using " <<  j*alpha + i<< endl;
      } else
    cout << " ==> Zerois " << j << "-th decomp at " << i << "-th poly"<< endl;
    }
  }
  for (size_t i = 0; i < l; ++i)
      cout <<  "\t\t" << digitsCTmp[0].GetNumOfElements()<< endl;

  vector<DCRTPoly> pPartExtC(digitsCTmp.size());
  vector<DCRTPoly> expandedC(digitsCTmp.size());
  for (uint32_t j=0; j<digitsCTmp.size(); j++) {
    auto tmpDigit = digitsCTmp[j].Clone();

    tmpDigit.SetFormat(Format::COEFFICIENT);

    const shared_ptr<typename DCRTPoly::Params> params = cryptoParamsLWE->GetComplementaryPartition(cipherTowers-1, j);

    pPartExtC[j] = tmpDigit.ApproxSwitchCRTBasis(cryptoParamsLWE->GetQPartition(j), params, //paramsP,
        cryptoParamsLWE->GetPartitionQHatInvModQTable(j)[digitsCTmp[j].GetNumOfElements()-1],
        cryptoParamsLWE->GetPartitionQHatInvModQPreconTable(j)[digitsCTmp[j].GetNumOfElements()-1],
        cryptoParamsLWE->GetPartitionQHatModPTable(cipherTowers-1)[j],
        cryptoParamsLWE->GetPartitionPrecon(cipherTowers-1)[j]);

    pPartExtC[j].SetFormat(Format::EVALUATION);
    cout << " # DCRT RNS Ring Element of PPart" <<     pPartExtC[j].GetNumOfElements() << endl;
    // Mod PQ
    expandedC[j] = DCRTPoly(cTmp.GetExtendedCRTBasis(paramsP), Format::EVALUATION, true);
    cout << "1. Setting Expanded Ctxt Mod Q part for " << j <<"th Decom " << "\n";
    // NOTE we are working on j-th decomposed stuff
    for (usint i=0; i<cipherTowers; i++) {
      if (i/alpha == j)
      {
        cout << " i = " << i << " set from DecomposeModQ Stuff on poly index i mod alpha " << i%alpha << endl;
        expandedC[j].SetElementAtIndex(i, digitsCTmp[j].GetElementAtIndex(i % alpha));
    }
      else {
        if (i/alpha < j) {
          cout << " i = " << i << " set from ExtendedToPQ_iHat Stuff on poly index i " << i << endl;
          expandedC[j].SetElementAtIndex(i, pPartExtC[j].GetElementAtIndex(i));
        } else {
          cout << " i = " << i << " set from ExtendedToPQ_iHat Stuff on poly index i-alpha" << i-alpha << endl;
          expandedC[j].SetElementAtIndex(i, pPartExtC[j].GetElementAtIndex(i - alpha));
        }
      }
    }

    cout << "2. Setting Expanded Ctxt Mod P part for " << j <<"th Decom " << "\n";
    for (usint i=0; i<paramsP->GetParams().size(); i++) {
      expandedC[j].SetElementAtIndex(i+cipherTowers, pPartExtC[j].GetElementAtIndex(i + params->GetParams().size() - paramsP->GetParams().size()));
          cout << " i = " << i+cipherTowers << " set from ExtendedToPQ_iHat Stuff on poly index " << i + params->GetParams().size() - paramsP->GetParams().size() << endl;
    }
  }

  // for (size_t dcrt_id = 0; dcrt_id < expandedC.size(); ++dcrt_id)
  // {
    // cout << "[" << dcrt_id << "]" << endl;
    // for (usint polyid = 0; polyid < expandedC[dcrt_id].GetNumOfElements(); ++polyid)
      // cout << "\t\t" << expandedC[dcrt_id].GetElementAtIndex(polyid)[0] << endl;
  // }

  vector<DCRTPoly> result(expandedC.size());
  for (uint32_t i=0; i<expandedC.size(); i++) {
    result[i] = expandedC[i];
  }

  shared_ptr<vector<DCRTPoly>> resultPtr = make_shared<vector<DCRTPoly>>(result);

  return resultPtr;
}


// Each pointer refers to (c(X^i), IP.b, IP.a) \in R_q \times (R_PQ)^2
inline shared_ptr<vector<DCRTPoly>> EvalAddLazyComponents(
  const vector<shared_ptr<vector<DCRTPoly>>>& lazy_components
)
{
  size_t n = lazy_components.size();
  vector<DCRTPoly> result{
    (*lazy_components[0])[0],
    (*lazy_components[0])[1],
    (*lazy_components[0])[2]
  };
  //cout << result[0].IsEmpty() << endl;
  //cout << result[1].IsEmpty() << endl;
  //cout << result[2].IsEmpty() << endl;

  for (size_t i = 1; i < n; ++i)
  {
     result[0] = result[0] + (*lazy_components[i])[0];
     result[1] = result[1] + (*lazy_components[i])[1];
     result[2] = result[2] + (*lazy_components[i])[2];
  }
  //cout << std::log2(result[0].GetModulus().ConvertToDouble()) << endl;
  //cout << std::log2(result[1].GetModulus().ConvertToDouble()) << endl;
  //cout << std::log2(result[2].GetModulus().ConvertToDouble()) << endl;
  shared_ptr<vector<DCRTPoly>> resultPtr = make_shared<vector<DCRTPoly>>(result);
  return resultPtr;
}

// TODO Given summed ctxt mod Q and mod PQ
inline Ciphertext<DCRTPoly> EvalModDownIntermediate(
  ConstCiphertext<DCRTPoly> ciphertext,
  const shared_ptr<vector<DCRTPoly>> expandedCiphertext_sum // b part summed mod Q, and Decomposed Ciphertext mod PQ summed
)
{

  // Resulting Storate
  Ciphertext<DCRTPoly> newCiphertext = ciphertext->CloneEmpty();

  // Preparation for ModDown
  (*expandedCiphertext_sum)[1].SetFormat(Format::COEFFICIENT);
  (*expandedCiphertext_sum)[2].SetFormat(Format::COEFFICIENT);

  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(ciphertext->GetCryptoParameters());
  const std::vector<DCRTPoly> &cc = ciphertext->GetElements();
  const shared_ptr<typename DCRTPoly::Params> paramsQ = cc[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP = cryptoParamsLWE->GetAuxElementParams();

  // TODO
  // ModDown (May need modifications)
  DCRTPoly cHat0 = (*expandedCiphertext_sum)[1].ApproxModDown(
      paramsQ,
      paramsP,
      cryptoParamsLWE->GetPInvModQTable(),        // P^{-1} mod q_i
      cryptoParamsLWE->GetPInvModQPreconTable(),  // Barrett Const to multiply with P^{-1} mod q_i
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

// TODO
// ModDown (May need modifications)
  DCRTPoly cHat1 = (*expandedCiphertext_sum)[2].ApproxModDown(
      paramsQ,
      paramsP,
      cryptoParamsLWE->GetPInvModQTable(),
      cryptoParamsLWE->GetPInvModQPreconTable(),
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

// TODO If this is the intermediate step, skip this transform
// But this lazy ntt should be done later
  cHat0.SetFormat(Format::EVALUATION);
  cHat1.SetFormat(Format::EVALUATION);

  DCRTPoly ct0;
  DCRTPoly ct1;
// NOTE cHat0, b-part of innerprodcucted result which is mod-downed
// (c: original automorphed ctxt mod Q)
  ct0 = (*expandedCiphertext_sum)[0] + cHat0;
// NOTE cHat1, a-part of innerprodcucted result which is mod-downed
  ct1 = cHat1;

  newCiphertext->SetElements({ ct0, ct1 });
  newCiphertext->SetDepth(ciphertext->GetDepth());
  newCiphertext->SetLevel(ciphertext->GetLevel());
  newCiphertext->SetScalingFactor(ciphertext->GetScalingFactor());

  return newCiphertext;
}

#if 0
inline Ciphertext<DCRTPoly> EvalFastAutomorphHybrid(
    ConstCiphertext<DCRTPoly> ciphertext,
    const usint autoIndex,
    const shared_ptr<vector<DCRTPoly>> expandedCiphertext
)
{

  // Retrieve the automorphism key that corresponds to the auto index.
  auto evalKey = ciphertext->GetCryptoContext()->
      GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag()).
      find(autoIndex)->second;

  // Apply the automorphism to the first component of the ciphertext.
  DCRTPoly psi_c0(ciphertext->GetElements()[0].AutomorphismTransform(autoIndex));

  std::vector<DCRTPoly> c(2);

  c[0] = psi_c0;
  c[1] = ciphertext->GetElements()[1].Clone();

  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(evalKey->GetCryptoParameters());

  Ciphertext<DCRTPoly> newCiphertext = ciphertext->CloneEmpty();

  std::vector<DCRTPoly> b = evalKey->GetBVector();
  std::vector<DCRTPoly> a = evalKey->GetAVector();

  DCRTPoly ct0;
  DCRTPoly ct1;

  const shared_ptr<typename DCRTPoly::Params> paramsQ = c[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP = cryptoParamsLWE->GetAuxElementParams();

  size_t cipherTowers = c[0].GetParams()->GetParams().size();
  size_t towersToSkip = cryptoParamsLWE->GetElementParams()->GetParams().size() - cipherTowers;

  DCRTPoly cTmp, cOrig;

  cOrig = c[1];
  cTmp = c[1].Clone();

  DCRTPoly cTilda0((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);
  DCRTPoly cTilda1((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);

  for (uint32_t j=0; j<expandedCiphertext->size(); j++) {
    DCRTPoly expandedC((*expandedCiphertext)[j].AutomorphismTransform(autoIndex));

    for (usint i=0; i<expandedC.GetNumOfElements(); i++) {
      usint idx = ( i < cipherTowers ) ? i : i + towersToSkip;
      cTilda0.SetElementAtIndex(i, cTilda0.GetElementAtIndex(i) + expandedC.GetElementAtIndex(i) * b[j].GetElementAtIndex(idx));
      cTilda1.SetElementAtIndex(i, cTilda1.GetElementAtIndex(i) + expandedC.GetElementAtIndex(i) * a[j].GetElementAtIndex(idx));
    }
  }

  cTilda0.SetFormat(Format::COEFFICIENT);
  cTilda1.SetFormat(Format::COEFFICIENT);

  DCRTPoly cHat0 = cTilda0.ApproxModDown(paramsQ, paramsP,
      cryptoParamsLWE->GetPInvModQTable(),
      cryptoParamsLWE->GetPInvModQPreconTable(),
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  DCRTPoly cHat1 = cTilda1.ApproxModDown(paramsQ, paramsP,
      cryptoParamsLWE->GetPInvModQTable(),
      cryptoParamsLWE->GetPInvModQPreconTable(),
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  cHat0.SetFormat(Format::EVALUATION);
  cHat1.SetFormat(Format::EVALUATION);

  ct0 = c[0] + cHat0;
  ct1 = cHat1;

  newCiphertext->SetElements({ ct0, ct1 });

  newCiphertext->SetDepth(ciphertext->GetDepth());
  newCiphertext->SetLevel(ciphertext->GetLevel());
  newCiphertext->SetScalingFactor(ciphertext->GetScalingFactor());

  return newCiphertext;
}

// Almost logically done
// TODO: check what's really needed and unneeded
//
//
inline shared_ptr<vector<DCRTPoly>> EvalFastAutomorphHybridLazy(
    ConstCiphertext<DCRTPoly> ciphertext, // original ciphertext to do rotate
    const usint autoIndex,
    const shared_ptr<vector<DCRTPoly>> expandedCiphertext
    // Decomposed Ciphertext mod PQ where each DCRTPoly corresponds to a decomposed digit
    // to rotate
)
{
  // 1. Retrieve the automorphism key that corresponds to the auto index.
  auto evalKey = ciphertext->GetCryptoContext()->
      GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag()).
      find(autoIndex)->second;

  // 2. Apply the automorphism to the first component of the ciphertext.
  DCRTPoly psi_c0(ciphertext->GetElements()[0].AutomorphismTransform(autoIndex));


  std::vector<DCRTPoly> c(2);

  c[0] = psi_c0;
  // TODO seems unused anywhere
  c[1] = ciphertext->GetElements()[1].Clone();

  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(evalKey->GetCryptoParameters());

  //Ciphertext<DCRTPoly> newCiphertext = ciphertext->CloneEmpty();

  // Eval Keys
  std::vector<DCRTPoly> b = evalKey->GetBVector();
  std::vector<DCRTPoly> a = evalKey->GetAVector();

  DCRTPoly ct0;
  DCRTPoly ct1;

  const shared_ptr<typename DCRTPoly::Params> paramsQ = c[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP = cryptoParamsLWE->GetAuxElementParams();

  size_t cipherTowers = c[0].GetParams()->GetParams().size();
  size_t towersToSkip = cryptoParamsLWE->GetElementParams()->GetParams().size() - cipherTowers;

// Inner Prod with Evk
  // Now, we are preparing the resulting storage mod PQ (filled with 0??) Perhaps
  DCRTPoly cTilda0((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);
  DCRTPoly cTilda1((*expandedCiphertext)[0].GetParams(), Format::EVALUATION, true);

  // Perform automorphisms on extended ctxt
  for (uint32_t j=0; j<expandedCiphertext->size(); j++) {
    DCRTPoly expandedC((*expandedCiphertext)[j].AutomorphismTransform(autoIndex));

    // And then perform innerprod b/w evk
    for (usint i=0; i<expandedC.GetNumOfElements(); i++) {
      usint idx = ( i < cipherTowers ) ? i : i + towersToSkip;
      cTilda0.SetElementAtIndex(i, cTilda0.GetElementAtIndex(i) + expandedC.GetElementAtIndex(i) * b[j].GetElementAtIndex(idx));
      cTilda1.SetElementAtIndex(i, cTilda1.GetElementAtIndex(i) + expandedC.GetElementAtIndex(i) * a[j].GetElementAtIndex(idx));
    }
  }

  vector<DCRTPoly> result{psi_c0, cTilda0, cTilda1};
  //cout << psi_c0 << endl;

  shared_ptr<vector<DCRTPoly>> resultPtr = make_shared<vector<DCRTPoly>>(result);
  //cout << " " << result[0].GetElementAtIndex(0)[0]  << << endl;
  //cout << " " << result[1].GetElementAtIndex(0)[0]  << << endl;
  //cout << " " << result[2].GetElementAtIndex(0)[0]  << << endl;

  return resultPtr;
}
#endif
