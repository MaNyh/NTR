import numpy as np
import matplotlib.pyplot as plt

def test_transientcheck():
    return 0

def transientcheck():
    return 0

def signalcreator():
    return 0


datalen = 1000
amp = 10
oscilation_duration = 10
date_rng = np.arange(0,datalen)

time = np.arange(0, datalen)/10


def tanh(time):
    stationary = np.random.rand()*2-1
    zaehl = 1#np.random.rand()
    nenn = np.e**(2*time*np.random.rand()) - 1
    ans = stationary + zaehl/nenn
    ans /= max(abs(ans))
    return ans

sinus = np.sin(time) * np.e**(-time*oscilation_duration**-1)*amp
tanh = tanh(time)
rausch = np.random.rand(datalen)

signal= sinus+tanh+rausch
#print(signal)
fig, axs = plt.subplots(2, 1)

axs[1].plot(time,signal)
axs[0].plot(time,tanh+sinus)

def mean(series):
    sum = np.sum(series)
    return sum/len(series)



plt.show()
