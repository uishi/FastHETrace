from enum import IntEnum, auto

class KSMethod(IntEnum):
    ORGHK20   = auto()
    PALIHK20  = auto()
    LATGOHK20 = auto()

def invNTT(logN, nummod, numringelem):
    N = 1 << logN
    return N * logN * nummod * numringelem

def NTT(logN, nummod, numringelem):
    N = 1 << logN
    return N * logN * nummod * numringelem

def approx_base_conv(logN, nummod_src, nummod_dst, numringelem):
    N = 1 << logN
    base_conv_per_int = nummod_src * (1 + nummod_dst)
    return  base_conv_per_int * numringelem * N

def approx_base_conv_HK20(logN, nummod_src, nummod_dst, numringelem):
    N = 1 << logN
    base_conv_per_int = nummod_src * (nummod_dst - nummod_src)
    return  base_conv_per_int * numringelem * N

def base_conv(logN, nummod_src, nummod_dst, numringelem, m_ks: KSMethod):
    bconv_cost = 0
    if   m_ks == KSMethod.ORGHK20:
        bconv_cost  = approx_base_conv_HK20(logN, nummod_src, nummod_dst, numringelem)
    elif m_ks == KSMethod.PALIHK20 or m_ks == KSMethod.LATGOHK20:
        bconv_cost  = approx_base_conv(logN, nummod_src, nummod_dst, numringelem)
    else:
        raise Exception("KS Method invalid BASECONV")

    # This should not happen
    if bconv_cost == 0:
        raise Exception("Bconv_Cost must not be zero")

    return bconv_cost

def cost_decomp(logN, num_ctxt_mod):
    N = 1 << logN
    return N * num_ctxt_mod

def cost_modup(logN, ell, k, num_digit, m_ks: KSMethod):
    # For each digit, (ell + 1) RNS components are produced from k RNS components.
    bconv_cost  = base_conv(logN, k, ell + 1, num_digit, m_ks)
    if bconv_cost == 0:
        raise Exception("KS Method invalid")

    return bconv_cost

# Multipication over R_{PQ} and R_{PQ} twice
def cost_innerprod(logN, ell, k, d):
    N = 1 << logN
    return   2 * d * (k + ell + 1) * N

# ModDown on R_{PQ}^{numringelem}
def cost_moddown(logN, ell, k, numringelem, m_ks:KSMethod):
    N = 1 << logN
    bconv_cost  =  0
    if m_ks == KSMethod.ORGHK20:
        raise Exception("Unknown")
    else:
        bconv_cost = base_conv(logN, k, ell + 1, 1, m_ks)

    if bconv_cost == 0:
        raise Exception("BaseConv cost should not be 0")

    # Mult Pinv mod qi for each i = 0 to ell
    mult_by_pinv  =  (ell + 1) * N
    return  numringelem * (bconv_cost + mult_by_pinv)
