import numpy as np
from diffusion import Diffusion
from constant_drift_constant_diffusion import Constant_Drift_Constant_Diffusion


"""
We want to evaluate V_0 = E[g(X_T)]. 
Here, g(x) is the call option payoff max(x-k,0).
"""

#np.random.seed(0)

k=0
def g(x):
    return np.maximum(x-k, 0)
    
"""
This case demonstrates Standard Brownian Motion dX_t = dW_t, X_0=0
where W_t is standard Brownian Motion.
"""
obj = Diffusion(num_of_paths=100, num_of_steps=100, T=1)
obj.simulate_W()
print obj.get_V0(g, obj.W)
#obj.plot_paths(obj.W)

"""
This case demonstrates the constant diffusion case dX_t = mu*dt + sigma*dW_t, X_0=0
where W_t is standard Brownian Motion.
"""
obj1 = Constant_Drift_Constant_Diffusion(num_of_paths=100, num_of_steps=100, T=1, mu=20, sigma=10)
obj1.simulate_process()
print obj1.get_V0(g, obj1.X)
#obj1.plot_paths(obj1.X)
