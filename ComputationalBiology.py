from lmfit import minimize, Parameters, Parameter, report_fit
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def f(xs, t, ps):
    try:
        n = ps['a'].value
        phi = ps['b'].value
        mu = ps['c'].value
    except:
        a, b, c = ps
    x = xs[0]
    p = xs[1]
    dx_dt = (1-x)*(1-phi*x-p)/n-mu*x/n
    dp_dt = phi*x
    return np.array([dx_dt,dp_dt])

def g(t, x0, ps):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(f, x0, t, args=(ps,))
    return x
def residual(ps, ts, data):
    x0 = ps['x0'].value
    p0 = ps['p0'].value
    X0=np.array([x0,p0])
    model = g(ts, X0, ps)
    return (model - data).ravel()
n = 2.0
phi = 0.5
mu = 0.8
true_params = [n,phi,mu]
X0 = np.array([10.0,0])

t = np.linspace(0, 10, 10)
data = g(t, X0, true_params)
data += np.random.normal(size=data.shape)

# set parameters incluing bounds
params = Parameters()
params.add('x0', value=float(data[0]), min=0, max=100)
params.add('n', value= 1.0, min=0, max=10)
params.add('phi', value= 1.0, min=0, max=10)
params.add('mu',value=1.0,min=0,max=10)
# fit model and find predicted values
result = minimize(residual, params, args=(t, data), method='leastsq')
final = data + result.residual.reshape(data.shape)

# plot data and fitted curves
plt.plot(t, data, 'o')
plt.plot(t, final, '--', linewidth=2, c='blue')

# display fitted statistics
report_fit(result)
