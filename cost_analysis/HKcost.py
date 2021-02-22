import math
from cost_basic import *

# Computes unit costs for C_{H1}, C_{H2} (per-unrolled iteration), and C_{IP} (on one evaluation key and a decomposed ring vector).
# Both (input, output) are in coeff representation.
def compute_HKcost_cc(logN, ell, k, d, m_ks:KSMethod):
    c_hoist1 = 0
    c_hoist2 = 0
    # LATTIGO made this step free
    # 1) RNS-decomposition
    if m_ks != KSMethod.LATGOHK20:
        c_hoist1  += cost_decomp(logN, ell + 1)
    # 2) ModUp
    c_hoist1  += cost_modup(logN, ell, k, d, m_ks)
    # 3) NTTs on extended digits (k + ell + 1 RNS components)
    NTT_on_all_digit  = NTT(logN, (k + ell + 1), d)
    c_hoist1 += NTT_on_all_digit
    # 4) IPs
    c_ip      = cost_innerprod(logN, ell, k, d)

    # 5) INTTs to do moddown (baseconv)
    # INTTs on all RNS components ensure ModDown output is in coeff rep
    iNTTonRpqSquare = invNTT(logN, ell + k + 1, 2)
    c_hoist2  += iNTTonRpqSquare

    # 6) ModDown on RqP^2 (Base Conv and drop)
    c_hoist2  += cost_moddown(logN, ell, k, 2, m_ks)
    return c_hoist1, c_ip, c_hoist2

# Computes unit costs for C_{H1}, C_{H2} (per-unrolled iteration), and C_{IP} (on one evaluation key and a decomposed ring vector).
# Both (input, output) are in evaluation representation.
def compute_HKcost_ee(logN, ell, k, d, m_ks:KSMethod):
    c_hoist1 = 0
    c_hoist2 = 0

    # C_{H1}
    # 1) INTT on c1 for base_conv
    iNTTonRq = invNTT(logN, ell + 1, 1)
    c_hoist1 += iNTTonRq
    # LATTIGO made this step free
    if m_ks != KSMethod.LATGOHK20:
        # 2) Split and zero-padding (BaseDecomp)
        c_hoist1 += cost_decomp(logN, ell) # 0 maybe
    # 3) ModUp
    c_hoist1 += cost_modup(logN, ell, k, d, m_ks)
    # 4) NTTs on extended digits (NOTE we only need ell+1 RNSs as k RNSs are already there)
    NTT_on_all_digit  = NTT(logN, (ell + 1), d)
    c_hoist1 += NTT_on_all_digit

    # 5) IP (C_{IP})
    c_ip     = cost_innerprod(logN, ell, k, d)
    # 6) Sum: FREE

    # C_{H2}
    # Less # NTTs
    if m_ks == KSMethod.LATGOHK20:
        # 7) INTTs on P part to do compute [x]_p mod Ql (baseconv)
        iNTTonRpSquare = invNTT(logN, k, 2)
        c_hoist2  += iNTTonRpSquare
        # 8) (After Baseconv) NTTs on [x]_p mod Ql to do rounding
        NTTonRqSquare = NTT(logN, ell +  1, 2)
        c_hoist2  += NTTonRqSquare
        # 9) BaseConv and MultByPinv and Rounding to complete ModDown  on 2 ring elements
        c_hoist2  += cost_moddown(logN, ell, k, 2, m_ks)

    # More # NTTs
    elif m_ks == KSMethod.PALIHK20:
        # 7) INTTs on input ctxt
        iNTTonRpSquare = invNTT(logN, ell + 1 + k, 2)
        c_hoist2  += iNTTonRpSquare
        # 8) BaseConv and MultByPinv and Rounding to complete ModDown  on 2 ring elements
        c_hoist2  += cost_moddown(logN, ell, k, 2, m_ks)
        # 9) Get result back to NTT (Extra NTT)
        NTTonRqSquare = NTT(logN, ell +  1, 2)
        c_hoist2  += NTTonRqSquare
    else:
        raise Exception("ModDown is currently only for LATTIGO and PALISADE")

    return c_hoist1, c_ip, c_hoist2

# Computes unit costs for C_{H1}, C_{H2} (per-unrolled iteration), and C_{IP} (inner-product on one evaluation key and a decomposed ring array) by paying a special attention on ModUp Phase.
# Both (input, output) are in evaluation representation.
def compute_HKcost_ee_advanced(logN, ell, k, beta, m_ks:KSMethod):
    c_hoist1 = 0
    c_hoist2 = 0
    # C_{H1}
    # 1) INTT on c1 for base_conv
    iNTTonRq  = invNTT(logN, ell + 1, 1)
    c_hoist1 += iNTTonRq
    # 2) RNS-decomposition
    # If not lattigo, constant multiplications are needed
    if m_ks != KSMethod.LATGOHK20:
        c_hoist1 += cost_decomp(logN, ell)

    if beta < 1:
        raise Exception("beta cannot be less than 1")
    # 3) ModUp (PALISADE)
    # NTT is performed on dest_mods after modup.
    if (ell + 1) % k != 0:
        # d - 1 digits follows the easy case, but the last digit does not.
        # simple modup on beta-1 digits
        if beta > 1:
            c_hoist1     += cost_modup(logN, ell, k, beta - 1, m_ks)
            c_hoist1     += NTT(logN, (ell + 1), beta - 1)

        # Modup on the lsat digit
        num_src_mod_last  = (ell + 1) % k
        num_dest_mod_last = (ell + 1 + k) - num_src_mod_last 
        c_hoist1         += base_conv(logN, num_src_mod_last, num_dest_mod_last, 1, m_ks)
        c_hoist1         += NTT(logN, num_dest_mod_last, 1)
    else:
        # 3) simple
        c_hoist1 += cost_modup(logN, ell, k, beta, m_ks)
        # 4) NTTs on extended digits (NOTE we only need ell+1 RNS components as k RNS components are already available in eval. rep..)
        c_hoist1 += NTT(logN, (ell + 1), beta)
    # In either way, the #NTTs are beta * (ell + 1 + k)

    # 5) IPs (C_{IP})
    c_ip     = cost_innerprod(logN, ell, k, beta)
    # 6) Sum: FREE

    # C_{H2}
    if m_ks == KSMethod.LATGOHK20:
        # 7) INTTs on P part to do compute [x]_p mod Ql (baseconv)
        iNTTonRpSquare = invNTT(logN, k, 2)
        c_hoist2      += iNTTonRpSquare
        # 8) (After Baseconv) NTTs on [x]_p mod Ql to do rounding
        NTTonRqSquare  = NTT(logN, ell +  1, 2)
        c_hoist2      += NTTonRqSquare
        # 9) BaseConv and MultByPinv and Rounding to complete ModDown  on 2 ring elements
        c_hoist2      += cost_moddown(logN, ell, k, 2, m_ks)

    elif m_ks == KSMethod.PALIHK20:
        # 7) INTTs on input ctxt
        iNTTonRpSquare = invNTT(logN, ell + 1 + k, 2)
        c_hoist2      += iNTTonRpSquare
        # 8) BaseConv and MultByPinv and Rounding to complete ModDown  on 2 ring elements
        c_hoist2      += cost_moddown(logN, ell, k, 2, m_ks)
        # 9) Get result back to NTT (Extra NTT)
        NTTonRqSquare  = NTT(logN, ell +  1, 2)
        c_hoist2      += NTTonRqSquare
    else:
        raise Exception("ModDown is currently only for LATTIGO and PALISADE")

    return c_hoist1, c_ip, c_hoist2

# Computes overall costs where both (input, output) are in coefficient representation
def unrolled_cc_cost(h, total_num_ks, logN, ell, k, d, m_ks:KSMethod):
    c_hoist1, c_ip, c_hoist2   = compute_HKcost_cc(logN, ell, k, d, m_ks)
    hoist1_cost                = h * c_hoist1
    ip_cost                    = total_num_ks * c_ip
    hoist2_cost                = h * c_hoist2
    return [hoist1_cost, ip_cost, hoist2_cost]

# Computes overall costs where both (input, output) are in evaluation representation
def unrolled_ee_cost(h, total_num_ks, logN, ell, k, d, m_ks:KSMethod):
    c_hoist1, c_ip, c_hoist2   = compute_HKcost_ee(logN, ell, k, d, m_ks)
    hoist1_cost                = h * c_hoist1
    ip_cost                    = total_num_ks * c_ip
    hoist2_cost                = h * c_hoist2
    return [hoist1_cost, ip_cost, hoist2_cost]

# Computes overall costs where both (input, output) are in evaluation representation, with a special attention on the last-digit for modup
def unrolled_ee_cost_advanced(h, total_num_ks, logN, ell, k, d, m_ks:KSMethod):
    c_hoist1, c_ip, c_hoist2   = compute_HKcost_ee_advanced(logN, ell, k, d, m_ks)
    hoist1_cost                = h * c_hoist1
    ip_cost                    = total_num_ks * c_ip
    hoist2_cost                = h * c_hoist2
    return [hoist1_cost, ip_cost, hoist2_cost]
