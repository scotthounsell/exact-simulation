import numpy as np
import matplotlib.pyplot as plt
from constant_drift_constant_diffusion import Constant_Drift_Constant_Diffusion

beta = 1 # not yet optimal

class Poisson_Process():
    def __init__(self, beta, T):
        self.waiting_times = []
        total = 0
        while total <= T:
            y = np.random.exponential(beta)
            total += y
            self.waiting_times.append(y)
        self.arrival_times = np.cumsum(self.waiting_times)


class Test(Constant_Drift_Constant_Diffusion):
    def __init__(self, num_of_paths, num_of_steps, T, X0=0, mu=0, sigma=1):
        super(Test, self).__init__(num_of_paths, num_of_steps, T, X0, mu, sigma)

    def generate_Y(self):
        poisson = Poisson_Process(beta, self.T)
        poisson.arrival_times = np.insert(poisson.arrival_times,0,0)
        print poisson.arrival_times#[1:-1]
        X = np.zeros(len(poisson.arrival_times))
        poisson.waiting_times[-1] = self.T - poisson.arrival_times[-2]

        X[0] = self.X0
        for i, t in enumerate(poisson.arrival_times[1:]):
            X[i+1] = X[i] + self.mu*poisson.waiting_times[i] + self.sigma*np.random.normal(0, np.sqrt(poisson.waiting_times[i]))
        print X
        plt.plot(poisson.arrival_times, X)

    def plot_path(self):
        pass

for i in xrange(1):
    obj = Test(1,1,9000, 2)
    obj.generate_Y()
plt.show()
