import sys; sys.path.insert(0, '/home/user/SigmaTau.jl/poster')
import numpy as np, stab_np as S
x = S.read_phase('/home/user/SigmaTau.jl/test/fixtures/validation/stable32gen.DAT')
print('N =', len(x))
ref = {}  # Stable32 sigma from docs/src/validation.md
ref['adev'] = {1:1.0097e+00,2:5.0444e-01,4:2.5021e-01,8:1.3063e-01,16:6.2750e-02,32:3.1804e-02,64:1.6085e-02,128:8.1063e-03,256:4.5004e-03,512:3.1642e-03,1024:3.1373e-03,2048:4.6131e-03}
ref['mdev'] = {1:1.0097e+00,2:3.5650e-01,4:1.2611e-01,8:4.6674e-02,16:1.5818e-02,32:6.0448e-03,64:2.5945e-03,128:1.7347e-03,256:1.8372e-03,512:2.2115e-03,1024:2.7671e-03,2048:5.2906e-03}
ref['hdev'] = {1:1.0646e+00,2:5.3205e-01,4:2.6192e-01,8:1.3854e-01,16:6.6018e-02,32:3.3500e-02,64:1.6951e-02,128:8.4500e-03,256:4.5672e-03,512:2.9510e-03,1024:2.4608e-03,2048:2.8169e-03}
ref['totdev'] = {1:1.0097e+00,2:5.0443e-01,4:2.5026e-01,8:1.3060e-01,16:6.2782e-02,32:3.1842e-02,64:1.6090e-02,128:8.1074e-03,256:4.4756e-03,512:3.0558e-03,1024:2.8694e-03,2048:3.4123e-03,4096:4.0832e-03}  # 512: raw value
for name, fn in [('adev',S.adev),('mdev',S.mdev),('hdev',S.hdev),('totdev',S.totdev)]:
    worst = 0
    for m, r in ref[name].items():
        v, k = fn(x, m)
        rel = abs(v-r)/r; worst = max(worst, rel)
    print(f'{name}: worst rel diff vs Stable32 = {worst:.2e}')
