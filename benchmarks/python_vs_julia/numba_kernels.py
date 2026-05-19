"""Numba @njit translations of SigmaTau's mtotdev/htotdev greenhall kernels.

Faithful 0-based ports of `_mtotdev_greenhall` and `_htotdev_greenhall` from
`src/stab/core/total.jl`. Single-threaded — `parallel=False` — to make the
comparison apples-to-apples vs single-threaded allantools.
"""
import numpy as np
from numba import njit


@njit(cache=True, fastmath=False)
def mtotdev_numba(x: np.ndarray, m_values: np.ndarray, tau0: float) -> np.ndarray:
    N = x.shape[0]
    K = m_values.shape[0]
    devs = np.empty(K, dtype=np.float64)

    max_m = 0
    for k in range(K):
        if m_values[k] > max_m:
            max_m = m_values[k]
    ext = np.empty(3 * 3 * max(max_m, 1), dtype=np.float64)

    for k in range(K):
        m = m_values[k]
        nsubs = N - 3 * m + 1
        if nsubs < 1:
            devs[k] = np.nan
            continue

        seg_len = 3 * m
        half_n = seg_len / 2.0
        outer_sum = 0.0

        # n is 0-based subsequence start: x window is x[n .. n+seg_len-1]
        for n in range(nsubs):
            if m == 1:
                slope = (x[n + 2] - x[n]) / (2.0 * tau0)
            else:
                hi = int(np.floor(half_n))
                s1 = 0.0
                for i in range(hi):
                    s1 += x[n + i]
                s1 /= hi
                s2 = 0.0
                for i in range(hi, seg_len):
                    s2 += x[n + i]
                s2 /= (seg_len - hi)
                slope = (s2 - s1) / (half_n * tau0)

            # per-window detrended values + time-reversed mirror, tripled
            for j in range(seg_len):
                val = x[n + j] - slope * tau0 * j
                rev_val = x[n + seg_len - 1 - j] - slope * tau0 * (seg_len - 1 - j)
                ext[j] = rev_val
                ext[seg_len + j] = val
                ext[2 * seg_len + j] = rev_val

            # sliding triple-m-window 2nd-difference
            a1 = 0.0
            a2 = 0.0
            a3 = 0.0
            for i in range(m):
                a1 += ext[i]
                a2 += ext[i + m]
                a3 += ext[i + 2 * m]

            d2 = (a3 - 2.0 * a2 + a1) / m
            block_sum = d2 * d2

            for j in range(6 * m - 1):
                a1 += ext[j + m] - ext[j]
                a2 += ext[j + 2 * m] - ext[j + m]
                a3 += ext[j + 3 * m] - ext[j + 2 * m]
                d2 = (a3 - 2.0 * a2 + a1) / m
                block_sum += d2 * d2

            outer_sum += block_sum / (6.0 * m)

        devs[k] = np.sqrt(outer_sum / (2.0 * m * m * tau0 * tau0 * nsubs))

    return devs


@njit(cache=True, fastmath=False)
def htotdev_numba(x: np.ndarray, m_values: np.ndarray, tau0: float) -> np.ndarray:
    N = x.shape[0]
    K = m_values.shape[0]
    devs = np.empty(K, dtype=np.float64)

    Ny = N - 1
    y = np.empty(Ny, dtype=np.float64)
    for i in range(Ny):
        y[i] = (x[i + 1] - x[i]) / tau0

    max_m = 0
    for k in range(K):
        if m_values[k] > max_m:
            max_m = m_values[k]
    ext = np.empty(3 * 3 * max(max_m, 1), dtype=np.float64)

    for k in range(K):
        m = m_values[k]
        if m == 1:
            L = N - 3
            if L <= 0:
                devs[k] = np.nan
                continue
            s = 0.0
            for i in range(L):
                d3 = x[i + 3] - 3.0 * x[i + 2] + 3.0 * x[i + 1] - x[i]
                s += d3 * d3
            devs[k] = np.sqrt(s / (6.0 * L * tau0 * tau0))
            continue

        n_iter = Ny - 3 * m + 1
        if n_iter < 1:
            devs[k] = np.nan
            continue

        seg_len = 3 * m
        # The julia kernel's outer iterator `i` is 0-based and accesses
        # y[i+j] with j 1-based, i.e. y at offsets i..i+seg_len-1 in 0-based.
        # Translate the entire block to 0-based offsets directly.
        hi_lo = seg_len // 2            # was floor(seg_len/2)
        lo_lo = (seg_len + 1) // 2      # was ceil(seg_len/2) -> python "+1)//2"
        # julia: hi=floor(L/2); s1 over j=1..hi -> 0-based: y[i+0..i+hi-1]
        #         lo_start=ceil(L/2)+1; s2 over j=lo_start..L -> 0-based: y[i+lo_lo .. i+L-1]
        mid = float(seg_len // 2)
        dev_sum = 0.0

        for i_start in range(n_iter):
            s1 = 0.0
            for j in range(hi_lo):
                s1 += y[i_start + j]
            m1 = s1 / hi_lo

            s2 = 0.0
            for j in range(lo_lo, seg_len):
                s2 += y[i_start + j]
            m2 = s2 / (seg_len - lo_lo)

            if seg_len % 2 == 1:
                slope = (m2 - m1) / (0.5 * (seg_len - 1) + 1.0)
            else:
                slope = (m2 - m1) / (0.5 * seg_len)

            # julia: val = y[i+j] - slope*(j-1-mid) for j 1-based 1..L
            #   0-based j' = j-1: val_j' = y[i + j'] - slope*(j' - mid)
            # julia: rev = y[i + L - j + 1] - slope*(L - j - mid)
            #   0-based j': rev_j' = y[i + L - 1 - j'] - slope*(L - 1 - j' - mid)
            for j in range(seg_len):
                val = y[i_start + j] - slope * (j - mid)
                rev = y[i_start + seg_len - 1 - j] - slope * (seg_len - 1 - j - mid)
                ext[j] = rev
                ext[seg_len + j] = val
                ext[2 * seg_len + j] = rev

            a1 = 0.0
            a2 = 0.0
            a3 = 0.0
            for j in range(m):
                a1 += ext[j]
                a2 += ext[j + m]
                a3 += ext[j + 2 * m]

            d3 = (a3 - 2.0 * a2 + a1) / m
            sq = d3 * d3

            for j in range(6 * m - 1):
                a1 += ext[j + m] - ext[j]
                a2 += ext[j + 2 * m] - ext[j + m]
                a3 += ext[j + 3 * m] - ext[j + 2 * m]
                d3 = (a3 - 2.0 * a2 + a1) / m
                sq += d3 * d3

            dev_sum += sq / (6.0 * m)

        devs[k] = np.sqrt(dev_sum / (6.0 * n_iter))

    return devs
