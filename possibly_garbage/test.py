from __future__ import  division
import numpy as np
from matplotlib import pyplot as plt

from  ars import ARS


######################################
# Example 1: sample 10000 values 
# from the normal distribution N(2,3)
######################################
T = 1.0
def f(x, mu=1, sigma=0.8):
    """ 
    Log Normal distribution 
    """
    #return -1/(2*sigma**2)*(x-mu)**2
    return -sigma*x*(0.5*x - mu) - x**2/(2.0*T)

def fprima(x, mu=1, sigma=0.8):
    """
    Derivative of Log Normal distribution
    """
    #return -1/sigma**2*(x-mu)
    return -sigma*(x - mu) - x/T
    
def h(x, mu=1, sigma=0.8):
    denom = np.sqrt(2*np.pi*T/(1+sigma*T)) * np.exp(T*mu**2*sigma**2/(2+2*sigma*T))
    return np.exp(-sigma*x*(0.5*x - mu) - x**2/(2.0*T)) / denom

x = np.linspace(-100,100,100)
ars = ARS(f, fprima, xi = [-4,1,40], mu=1, sigma=0.8)
samples = ars.draw(10000)
print "mean", np.mean(samples), T*1*0.8/(1+0.8*T)
print "std", np.std(samples), np.sqrt(T/(1+0.8*T))
plt.hist(samples, bins=100, normed=True)
x_dist = np.arange(-4,4,0.1)
plt.plot(x_dist, [h(y) for y in x_dist], color='red')
plt.show()
