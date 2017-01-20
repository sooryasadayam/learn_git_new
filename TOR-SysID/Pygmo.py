from scipy.integrate import odeint
from scipy.optimize import differential_evolution as d_e
import matplotlib.pyplot as plt
#dS_dt(S_vec, t, v, w, eps, zhi, sig, mu, lam, kap)
import numpy as np
P_vec = np.array([8.64925618e-05,   4.11695577e+00,   5.02092999e-01,
         1.48110096e+00,   8.30433556e+01,   2.47557115e-03,
         5.57898571e-01,   9.96937216e-01])
t=np.linspace(0,20,201)
def dS_dt(S_vec, t, v, w, eps, zhi, sig, mu, lam, kap):
    dx_dt = ((1 - S_vec[0]) * S_vec[3] - mu * S_vec[0]) / v
    dm_dt = (1 - S_vec[0] - S_vec[1]) / w
    dp_dt = S_vec[1] - S_vec[2] - sig * (S_vec[2] - lam * S_vec[3])
    dz_dt = (S_vec[2] - (lam + kap) * S_vec[3] - zhi * ((1 - S_vec[0]) * S_vec[3] - sig * S_vec[0])) / eps
    return np.array([dx_dt, dm_dt, dp_dt, dz_dt])


S_in = np.array([0, 0, 0, 0])
Solution = odeint(dS_dt, S_in, t, args=tuple(P_vec))
P_sat = Solution[len(t) - 1, 3]
a = np.max(np.where(Solution[:, 3] <= (P_sat * 0.5)))
b = np.min(np.where(Solution[:, 3] >= (P_sat * 0.5)))
c = np.max(np.where(Solution[:, 3] <= (P_sat * 0.99)))
d = np.min(np.where(Solution[:, 3] >= (P_sat * 0.99)))
if (abs(Solution[a,3] - P_sat * 0.5) > abs(Solution[b, 3] - P_sat * 0.5)):
    t50 = b
    print t[t50]
else:
    t50 = a
    print t[50]
if (abs(Solution[c, 3] - P_sat * 0.99) > abs(Solution[d, 3] - P_sat * 0.99)):
    t99 = d
    print t[t99]
else:
    t99 = c
    print t[t99]
plt.plot(t, Solution[:, 3], t[t50], Solution[t50, 3], 'bs', t[t99], Solution[t99,3], 'rs')
plt.show()