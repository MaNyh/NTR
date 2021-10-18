import numpy as np
import matplotlib.pyplot as plt

"""
this module is supposed to return a timesstamp from a time-series, that equals the time when the transient signal ends

Ries 2018
https://link.springer.com/content/pdf/10.1007/s00162-018-0474-0.pdf

numerical solutions have a transient behaviour at the initial process
it is assumed that this initial transient can be reproduced by a sine, tanh and a noise-function
with these functions given, we can analytically define where the transient process ends
"""

class signal_generator:
    # resolution
    datalen = 10000 # length of timeseries
    # value between [0;1]  defines when a signal is stationary
    transientlimit = 0.95

    def __init__(self):
        #some kind of factor for the duration of an abating signal
        self.sin_lasting = np.random.randint(0, 100)
        self.tanh_lasting = np.random.randint(0, 100)

        self.sin_stationary_ts = -self.sin_lasting * (1 + np.log(1 - self.transientlimit))
        self.tanh_stationary_ts = np.arctanh(self.transientlimit) * self.tanh_lasting

        # as defined in Ries 2018, the signal must be at least two times as long as the transient process
        if self.sin_stationary_ts > self.tanh_stationary_ts:
            self.time = 2*self.sin_stationary_ts
        else:
            self.time = 2*self.tanh_stationary_ts

        #defining the used timesteps (unnesessary, the signal-length could be normed to 1!)
        self.timesteps = np.arange(0, self.time, self.datalen ** -1)

        #damping the sine with euler
        abate = np.e ** (-self.timesteps * self.sin_lasting ** -1)
        self.sin_abate = abate / max(abate)
        #randomising the sign of tanh
        self.tanh_sign = np.random.choice([-1, 1])

        #values for the 'stationarity'. this can be defined because the function is analytical. but is this correct?
        self.sin_stationary = np.argmax(self.sin_abate < (1 - self.transientlimit))
        self.tanh_stationarity = np.argmax( np.tanh(self.timesteps * self.tanh_lasting ** -1) > self.transientlimit)

        print()
        print("time when stationary")
        sin_statts = int(self.sin_stationary_ts / self.time * len(self.timesteps))
        tanh_statts = int(self.tanh_stationary_ts / self.time * len(self.timesteps))
        print("sin : ", self.sin_stationary_ts , " timestep ", sin_statts)
        print("tanh : ", self.tanh_stationary_ts , " timestep ", tanh_statts)
        self.stat_time = self.tanh_stationary_ts if self.tanh_stationary_ts > self.sin_stationary_ts else self.sin_stationary_ts
        print("signal stationary at", self.stat_time)

        self.stationarity_relt = round(int(self.stat_time/self.time*len(self.timesteps))/len(self.timesteps),1)
        self.stationarity_time = self.stationarity_relt * self.time

        print("this equals a timestep of ", int(self.stat_time / self.time * len(self.timesteps)), " from ",
              len(self.timesteps), " or ", self.stationarity_relt)


    def tanh_signal(self):
        ans = np.tanh(self.timesteps * self.tanh_lasting ** -1)
        return ans * self.tanh_sign

    def sin_signal(self):
        sinus = np.sin(self.timesteps*np.random.rand()) * self.sin_abate
        return sinus

    def noise_signal(self):
        mu, sigma = 0, np.random.rand() # mean and standard deviation
        s = np.random.normal(mu, sigma, size=len(self.timesteps))
        if max(s)>1:
            s /= max(s)
        return s

    def generate(self):
        sinus = self.sin_signal()
        tanh = self.tanh_signal()
        rausch = self.noise_signal()

        sin_stats = (-1 + self.sin_abate) * -1
        tanh_stats = np.sinh(self.timesteps * self.tanh_lasting ** -1) / np.cosh(self.timesteps * self.tanh_lasting ** -1)

        signal = sinus + tanh + rausch

        return sinus, tanh, rausch, signal, sin_stats , tanh_stats


    def plot(self,sinus, tanh, rausch, signal, stat_sin , stat_tanh):
        fig, axs = plt.subplots(6, 1)

        axs[0].plot(np.arange(0, self.time, self.datalen ** -1), sinus, color="orange", label="abating sine signal")
        axs[0].axvline(self.sin_stationary_ts)
        axs[1].plot(np.arange(0, self.time, self.datalen ** -1), stat_sin, color="orange", label="stationarity sine signal")
        axs[1].axvline(self.sin_stationary_ts)
        axs[1].fill_between(np.arange(0, self.time, self.datalen ** -1),stat_sin, color="orange")
        axs[2].plot(np.arange(0, self.time, self.datalen ** -1), tanh, color="blue", label="tanh signal")
        axs[2].axvline(self.tanh_stationary_ts)
        axs[3].plot(np.arange(0, self.time, self.datalen ** -1), stat_tanh, color="blue", label="stationarity tanh signal")
        axs[3].axvline(self.tanh_stationary_ts)
        axs[3].fill_between(np.arange(0, self.time, self.datalen ** -1),stat_tanh, color="blue")
        axs[4].plot(np.arange(0, self.time, self.datalen ** -1), rausch, color="black", label="noise")
        axs[5].plot(np.arange(0, self.time, self.datalen ** -1), signal, color="red", label="signal")
        axs[5].axvline(self.stationarity_time)

        for a in axs:
            a.legend(loc="upper right")

        plt.show()


def test_transientcheck(verbose=True):
    siggen=signal_generator()

    sinus, tanh, rausch, signal, stat_sin , stat_tanh = siggen.generate()
    if verbose:
        siggen.plot(sinus, tanh, rausch, signal, stat_sin, stat_tanh)

    timestationarity_ries = transientcheck(signal)
    print("transientcheck-rückgabe (ries2018): ",timestationarity_ries)
    assert siggen.stat_time < timestationarity_ries, "von transientcheck zurückgegebene zeit der stationarität zu klein"
    return 0


def transientcheck(signal):
    """
    Hello Berk
    This function needs to analyze a time-series.
    It should return a value for the time, that equals the
    :param signal: timeseries
    :return: time_of_stationarity
    """
    return 0
