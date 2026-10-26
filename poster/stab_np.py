"""Minimal numpy port of the overlapping ADEV / MDEV / HDEV / TOTDEV kernels,
used only to draw poster figures when Julia is unavailable. Verified against
the Stable32 fixture values in docs/src/validation.md (see check_port.py)."""
import numpy as np

def read_phase(path, header=10):
    with open(path) as f:
        lines = f.readlines()[header:]
    return np.array([float(l) for l in lines if l.strip()])

def adev(x, m, tau0=1.0):
    x = np.asarray(x); n = len(x)
    d = x[2*m:] - 2*x[m:n-m] + x[:n-2*m]
    return np.sqrt(np.mean(d**2) / (2*(m*tau0)**2)), len(d)

def hdev(x, m, tau0=1.0):
    x = np.asarray(x); n = len(x)
    d = x[3*m:] - 3*x[2*m:n-m] + 3*x[m:n-2*m] - x[:n-3*m]
    return np.sqrt(np.mean(d**2) / (6*(m*tau0)**2)), len(d)

def mdev(x, m, tau0=1.0):
    x = np.asarray(x); n = len(x)
    d = x[2*m:] - 2*x[m:n-m] + x[:n-2*m]          # second differences, length n-2m
    c = np.cumsum(np.concatenate(([0.0], d)))
    s = c[m:] - c[:-m]                            # moving sums of m consecutive 2nd differences
    return np.sqrt(np.mean(s**2) / (2*m**2*(m*tau0)**2)), len(s)

def tdev(x, m, tau0=1.0):
    md, k = mdev(x, m, tau0)
    return md * m*tau0 / np.sqrt(3), k

def totdev(x, m, tau0=1.0):
    """Raw TOTDEV (SP1065 eq. 25): reflect the mean-removed phase at both ends."""
    x = np.asarray(x, dtype=float); n = len(x)
    xs = x - np.mean(x)          # SP1065 uses the phase after removing the linear slope? no: raw kernel uses x as-is in allantools (no detrend)
    xs = x
    ext = np.concatenate((2*xs[0] - xs[1:n-1][::-1], xs, 2*xs[-1] - xs[1:n-1][::-1]))
    off = n - 2                  # index of xs[0] in ext
    tot = 0.0
    for i in range(1, n-1):      # i = 1 .. n-2 (0-based), SP1065 i = 2 .. N-1
        j = off + i
        tot += (ext[j-m] - 2*ext[j] + ext[j+m])**2
    return np.sqrt(tot / (2*(m*tau0)**2*(n-2))), n-2
