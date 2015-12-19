import numpy as np

class Poisson_Process():
    """
    Generates a Poisson process with a truncated last step.
    """
    def __init__(self, beta, T):
        """
        Constructor. Sets two arrays:
        self.dt is the waiting times between arrivals (exponential random variables). 
            The last waiting time is shrunk so that its arrival coincides with T.
            A 0 is inserted in the beginning so that the indices match the reference paper.
        self.t is the arrival times. Note that a 0 is inserted to the front of the array. Also, the last arrival time will always be T.
        
        Args:
            beta (float): Parameter for exponential random variable.
            T (float): Expiry. The Poisson process is generated up to this point.
        """
        self.dt = [0]
        
        # Generate exponential random variables until their sum adds to more than T
        total = 0
        while total <= T:
            random_exponential = np.random.exponential(beta)
            total += random_exponential
            self.dt.append(random_exponential)
            
        # Change the last exponential so that the sum of all exponentials equals T
        self.dt[-1] -= total - T

        # self.t is the cumulative sum of the exponentials
        self.t = np.cumsum(self.dt)


# An example to show how the class works. Keep this commented out unless testing the class.
#~ test = Poisson_Process(1,3)
#~ print test.dt
#~ print test.t
