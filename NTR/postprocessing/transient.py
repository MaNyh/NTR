import numpy as np
import matplotlib.pyplot as plt


class signal_generator:
    datalen = 10000 # length of timeseries
    transientlimit = 0.95

    def __init__(self):
        self.sin_lasting = np.random.randint(0, 100)
        self.tanh_lasting = np.random.randint(0, 100)
        self.time = 100

        self.turn = np.random.choice([-1, 1])

        self.timesteps = np.arange(0, self.time, self.datalen ** -1)

        self.sinabate =np.e ** (-self.timesteps * self.sin_lasting ** -1)
        self.sinabate /= max(self.sinabate)

        self.sin_stationary = np.argmax(self.sinabate<(1-self.transientlimit))
        self.tanh_stationarity = np.argmax(np.sinh(self.timesteps * self.tanh_lasting ** -1) / np.cosh(self.timesteps * self.tanh_lasting ** -1)>self.transientlimit)

    def tanh_signal(self):
        ans = np.sinh(self.timesteps * self.tanh_lasting ** -1) / np.cosh(self.timesteps * self.tanh_lasting ** -1)
        return ans * self.turn

    def sin_signal(self):
        sinus = np.sin(self.timesteps) * self.sinabate
        return sinus

    def noise_signal(self):
        return np.random.rand(len(self.timesteps))*np.random.rand()

    def generate(self):
        sinus = self.sin_signal()
        tanh = self.tanh_signal()
        rausch = self.noise_signal()

        sin_stats = (-1 + self.sinabate)*-1
        tanh_stats = np.sinh(self.timesteps * self.tanh_lasting ** -1) / np.cosh(self.timesteps * self.tanh_lasting ** -1)
        if sin_stats[-1] < self.transientlimit:
            print("sin still transient ", str(sin_stats[-1]))
        if tanh_stats[-1] < self.transientlimit:
            print("tanh still transient ", str(tanh_stats[-1]))

        signal = sinus + tanh + rausch

        return sinus, tanh, rausch, signal, sin_stats , tanh_stats

    def analytical_stat_solution(self):
        return 0

    def plot(self,sinus, tanh, rausch, signal, stat_sin , stat_tanh):
        fig, axs = plt.subplots(6, 1)

        axs[0].plot(np.arange(0, self.time, self.datalen ** -1), sinus)
        axs[1].plot(np.arange(0, self.time, self.datalen ** -1), stat_sin)
        axs[2].plot(np.arange(0, self.time, self.datalen ** -1), tanh)
        axs[3].plot(np.arange(0, self.time, self.datalen ** -1), stat_tanh)
        axs[4].plot(np.arange(0, self.time, self.datalen ** -1), rausch)
        axs[5].plot(np.arange(0, self.time, self.datalen ** -1), signal)

        plt.show()


def test_transientcheck():
    siggen=signal_generator()
    infval = siggen.turn

    sin_stationary_ts = -siggen.sin_lasting**-1*np.log(1-siggen.transientlimit)
    print(sin_stationary_ts)
    tanh_stationary_ts = np.arctanh(siggen.transientlimit)*siggen.tanh_lasting**-1
    print(tanh_stationary_ts)

    print(siggen.timesteps)
    print(siggen.time)

    sinus, tanh, rausch, signal, stat_sin , stat_tanh = siggen.generate()
    siggen.plot(sinus, tanh, rausch, signal, stat_sin , stat_tanh)

    return 0


def transientcheck():
    """
    this method must check when a steady-state is
    https://link.springer.com/content/pdf/10.1007/s00162-018-0474-0.pdf
    """
    return 0



