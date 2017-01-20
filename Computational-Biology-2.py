from scipy.integrate import odeint
from scipy.optimize import differential_evolution as d_e
import matplotlib.pyplot as plt
#dS_dt(S_vec, t, v, w, eps, zhi, sig, mu, lam, kap)
import numpy as np
# S_vec=(np.array([X,M,P,Z]))
# P_vec=np.array([v,w,eps,zhi,sig,mu,lam,kap])
P_vec = np.array([10e-3, 10e-1, 10e-2, 10, 0.01, 10e-2, 0.1, 10e-6])
T = np.linspace(0,200,10000)


def dS_dt(S_vec, t, v, w, eps, zhi, sig, mu, lam, kap):
    dx_dt = ((1-S_vec[0])*S_vec[3]**15-mu*S_vec[0])/v
    dm_dt = (1-S_vec[0]-S_vec[1])/w
    dp_dt = S_vec[1]-S_vec[2]-sig*(S_vec[2] - lam*S_vec[3])
    dz_dt = (S_vec[2]-(lam + kap)*S_vec[3]-zhi*((1 - S_vec[0])*S_vec[3]**15 - mu*S_vec[0]))/eps
    return np.array([dx_dt, dm_dt, dp_dt, dz_dt])

S_in = np.array([0, 0, 0, 0])
Solution = odeint(dS_dt, S_in, T, args=tuple(P_vec))
D_Solution = (Solution[1:len(T)-1,3]-Solution[0:len(T)-2,3])/(T[1]-T[0])
P_sat = Solution[len(T) - 1, 3]
a = np.max(np.where(Solution[:, 3] <= (P_sat * 0.5)))
b = np.min(np.where(Solution[:, 3] >= (P_sat * 0.5)))
c = np.max(np.where(Solution[:, 3] <= (P_sat * 0.99)))
d = np.min(np.where(Solution[:, 3] >= (P_sat * 0.99)))
if (abs(Solution[a,3] - P_sat * 0.5) > abs(Solution[b,3] - P_sat * 0.5)):
    t50 = b
else:
    t50 = a
if (abs(Solution[c, 3] - P_sat * 0.99) > abs(Solution[d, 3] - P_sat * 0.99)):
    t99 = d
else:
    t99 = c
#plt.plot(T, Solution[:, 3], T[t50], Solution[t50, 3], 'bs', T[t99], Solution[t99, 3], 'rs')
#plt.show()
plt.plot(Solution[1:len(T)-1,3],D_Solution)
plt.show()
print len(Solution[0:len(T)-2,3])
print len(D_Solution)
