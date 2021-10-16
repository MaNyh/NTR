import numpy as np
import matplotlib.pyplot as plt


class signal_generator:
    datalen = 10000 # length of timeseries
    transientlimit = 0.95

    def __init__(self):
        self.sin_lasting = np.random.randint(0, 100)
        self.tanh_lasting = np.random.randint(0, 100)
        self.sin_stationary_ts = -self.sin_lasting * (1 + np.log(1 - self.transientlimit))
        self.tanh_stationary_ts = np.arctanh(self.transientlimit) * self.tanh_lasting

        if self.sin_stationary_ts > self.tanh_stationary_ts:
            self.time = 2*self.sin_stationary_ts
        else:
            self.time = 2*self.tanh_stationary_ts

        self.timesteps = np.arange(0, self.time, self.datalen ** -1)

        self.sin_abate = np.e ** (-self.timesteps * self.sin_lasting ** -1) / np.e
        self.tanh_sign = np.random.choice([-1, 1])
        self.sin_stationary = np.argmax(self.sin_abate < (1 - self.transientlimit))
        self.tanh_stationarity = np.argmax( np.tanh(self.timesteps * self.tanh_lasting ** -1) > self.transientlimit)

        print()
        print("time when stationary")
        ts = self.sin_stationary_ts
        ss = self.tanh_stationary_ts
        print("sin : ", ts)
        print("tanh : ", ss)
        stat_time = ts if ts > ss else ss
        print("signal stationary at", stat_time)
        print("this equals a timestep of ", int(stat_time/self.time*len(self.timesteps)), " from ", len(self.timesteps) , " or " , round(int(stat_time/self.time*len(self.timesteps))/len(self.timesteps),1))


    def tanh_signal(self):
        ans = np.tanh(self.timesteps * self.tanh_lasting ** -1)
        return ans * self.tanh_sign

    def sin_signal(self):
        sinus = np.sin(self.timesteps) * self.sin_abate
        return sinus

    def noise_signal(self):
        mu, sigma = 0, 1 # mean and standard deviation
        s = np.random.normal(mu, sigma, size=len(self.timesteps))
        s /= max(s)
        return s

    def generate(self):
        sinus = self.sin_signal()
        tanh = self.tanh_signal()
        rausch = self.noise_signal()

        sin_stats = (-1 + self.sin_abate) * -1
        tanh_stats = np.sinh(self.timesteps * self.tanh_lasting ** -1) / np.cosh(self.timesteps * self.tanh_lasting ** -1)
        if sin_stats[-1] < self.transientlimit:
            print("sin still transient ", str(sin_stats[-1]))
        if tanh_stats[-1] < self.transientlimit:
            print("tanh still transient ", str(tanh_stats[-1]))

        signal = sinus + tanh + rausch

        return sinus, tanh, rausch, signal, sin_stats , tanh_stats


    def plot(self,sinus, tanh, rausch, signal, stat_sin , stat_tanh):
        fig, axs = plt.subplots(6, 1)

        axs[0].plot(np.arange(0, self.time, self.datalen ** -1), sinus, color="orange", label="abating sine signal")
        axs[1].plot(np.arange(0, self.time, self.datalen ** -1), stat_sin, color="orange", label="stationarity sine signal")
        axs[2].plot(np.arange(0, self.time, self.datalen ** -1), tanh, color="blue", label="tanh signal")
        axs[3].plot(np.arange(0, self.time, self.datalen ** -1), stat_tanh, color="blue", label="stationarity tanh signal")
        axs[4].plot(np.arange(0, self.time, self.datalen ** -1), rausch, color="black", label="noise")
        axs[5].plot(np.arange(0, self.time, self.datalen ** -1), signal, color="red",label="signal")

        axs[0].legend(loc="upper right")
        axs[1].legend(loc="upper right")
        axs[2].legend(loc="upper right")
        axs[3].legend(loc="upper right")
        axs[4].legend(loc="upper right")
        axs[5].legend(loc="upper right")

        plt.show()


def test_transientcheck():
    siggen=signal_generator()

    sinus, tanh, rausch, signal, stat_sin , stat_tanh = siggen.generate()
    siggen.plot(sinus, tanh, rausch, signal, stat_sin , stat_tanh)

    return 0


def transientcheck():
    """
    this method must check when a steady-state is
    https://link.springer.com/content/pdf/10.1007/s00162-018-0474-0.pdf
    """
    return 0



