import numpy as np
from numpy import exp, sin, cos, sqrt
from math import factorial
import matplotlib.pyplot as plt
from scipy.stats import norm, t

class Beskos_Roberts:
    def __init__(self, alpha, alpha_deriv, alpha_integral, k1, k2):
        self.alpha = alpha
        self.alpha_deriv = alpha_deriv
        self.A = alpha_integral
        self.k1 = k1
        self.T = 1.0 / (k2 - k1)

    def phi(self, u):
        return 0.5 * ((self.alpha(u))**2 + self.alpha_deriv(u)) - self.k1

    def h(self, u):
        # The value of c only matters if we want int_{-inf}^{inf} h(u)du = 1.
        # Since h(u) is only used to draw random samples, the value
        # of c does not matter, only the ratio of h(u1) to h(u2) where 
        # u1,u2 are real numbers
        c = 1 
        return exp(self.A(u) - u**2/(2*self.T)) / c


    def random_h(self):
        # generate a random variable X distributed according to h
        # this just works for mean reverting for now
        return np.random.normal(T*1*0.8/(1+0.8*T), np.sqrt(T/(1+0.8*T)))

    def is_even(self, x):
        return not x % 2

    def generate_skeleton_path(self):
        while True:
            # Step 1
            omega_0 = 0.0
            #omega_T = 0.0 # fix this. Should be distributed according to P(omega_T=x) = h(x)
            omega_T = self.random_h() # only works for mean reverting for now
            omega = {0.0: omega_0, self.T: omega_T}
            V = []

            # Step 2
            U = np.random.uniform(0, 1)
            #print "U", U
            i = 0
            
            fill_another_step = True
            while fill_another_step:
                # Step 3
                V_i, W = np.random.uniform(0, [self.T, 1.0/self.T])
                i += 1
                #print "\ni", i
                #print "V_i =", V_i
                #print "W =", W
                

                # Step 4
                V_neg = max([0] + [e for e in V if e < V_i]) # the greatest V less than V_i
                #print "V_neg =", V_neg
                V_pos = min([T] + [e for e in V if e > V_i]) # the smallest V greater than V_i
                #print "V_pos =", V_pos
                mean = omega[V_neg] + (omega[V_pos] - omega[V_neg]) / (V_pos - V_neg) * (V_i - V_neg) # linear interpolation
                std = sqrt((V_pos - V_i) * (V_i - V_neg) / (V_pos - V_neg))
                omega[V_i] = np.random.normal(mean, std)
                #print "mean, std = ", mean, std
                #print "phi =", phi(omega[V_i])
                V.append(V_i)
                #print "V", V

                # Step 5
                if self.phi(omega[V_i]) < W or U > 1.0/factorial(i):
                    fill_another_step = False
                    if self.is_even(i):
                        I = 0 # why do we need I?
                        #print "not done 1 - reject path"
                    else:
                        I = 1 # why do we need I?
                        #print "done"
                        # Step 6
                        return omega
                # else go back to to step 3
                #print "omega = ", omega

def alpha(u):
    # radians or degrees?
    return sin(u)

def alpha_deriv(u):
    # radians or degrees?
    return cos(u)

def A(u):
    # A(u) = int_0^u alpha(y)dy ~ int alpha(y)dy
    return -cos(u)

k1 = -1.0/2
k2 = 5.0/8
T = 1.0 / (k2 - k1)

process = Beskos_Roberts(alpha, alpha_deriv, A, k1, k2)


omega = process.generate_skeleton_path()
print omega
print len(omega)
plt.plot(sorted(omega), [omega[key] for key in sorted(omega)])
plt.show()

