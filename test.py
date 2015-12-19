from numpy import sqrt, exp, log
import numpy as np

X0 = 1
sigma0 = X0
beta = 0.2
T = 1.0

def mu(y):
    x = convert_y_to_x(y)
    return 0.2*X0*(sqrt(x)-1) - 0.25*X0

def convert_y_to_x(y):
    return X0*exp(0.5*y)


def convert_x_to_y(x):
    return 2*log(x/X0)

def temp_mu(x):
    return 0.1*x*(sqrt(x)-1)

dt = [0, 0.7971959, 0.2028041]
N_T = 1
dW_1 = [0, 0.523418, -0.4584004]
#~ dW_2 = [0, -0.8515401, -0.1628777]

Yhat_1 = mu(0)*dt[1] + X0*dW_1[1]
print Yhat_1
Yhat_2 = Yhat_1 + mu(Yhat_1)*dt[2] + X0*dW_1[2]
print Yhat_2

Xhat_1 = convert_y_to_x(Yhat_1)
print Xhat_1
Xhat_2 = convert_y_to_x(Yhat_2)
print Xhat_2

dW = np.array([[0,  0.8252239, -0.43812191], 
               [0,  0.9794792, -0.07006420], 
               [0,  0.3490172,  0.04498612], 
               [0, -1.5740653, -0.74912387]])
               
Xhat = np.array([[1, 1.3674702, 1.074641],
                 [1, 1.4771139, 1.396648],
                 [1, 1.0777323, 1.075489],
                 [1, 0.4120205, 0.274211]])


diff = temp_mu(Xhat[:,1]) - temp_mu(Xhat[:,0])
product = sum(diff*dW[:,2])/(dt[2]*sigma0)
print 'product', product


Xhat_2 = [1.9124475, 0.8854799, 1.2303173, 1.0732895]

#~ print 'product', (mu(Xhat_1) - mu(1))*dW_1[2]/(dt[2]*sigma0)
#~ product = (temp_mu(Xhat_1) - temp_mu(1))*dW_1[2]/(dt[2]*sigma0)
#~ print 'correct product', product
#~ 
g_multiple = exp(beta*T) * beta**(-N_T) * product
print 'g_multiple', g_multiple
g_multiple = exp(beta*T) * beta**(-0)

strikes = np.arange(0.6, 1.6, 0.1)
print strikes
for k in strikes:
    print k, (max(sum(Xhat_2)/4.0 - k,0) - 0)*g_multiple
    






