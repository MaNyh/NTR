import numpy as np
import matplotlib.pyplot as plt


class signal_generator:
    def __init__(self):
        self.datalen = 10000
        self.date_rng = np.arange(0, self.datalen)

        self.sin_lasting = np.random.randint(0, 100)
        self.tanh_lasting = np.random.randint(0, 100)
        self.time = 100

        self.timesteps = np.arange(0, self.time, self.datalen ** -1)

    def tanh_signal(self):
        """
        :param time:
        :return:
        """

        turn = np.random.choice([-1, 1])
        ans = np.sinh(self.timesteps * self.tanh_lasting ** -1) / np.cosh(self.timesteps * self.tanh_lasting ** -1)

        stationarity = ans
        return ans * turn, stationarity

    def sin_signal(self, oscilation_duration):
        abate = np.e ** (-self.timesteps * oscilation_duration ** -1)
        abate /= max(abate)
        sinus = np.sin(self.timesteps) * abate
        return sinus, 1 - abate

    def noise(self):
        return np.random.rand(len(self.timesteps))

    def generate(self):
        sinus, stat_sin = self.sin_signal(self.sin_lasting)
        tanh, stat_tanh = self.tanh_signal()
        rausch = self.noise()

        if stat_sin[-1] < 0.95:
            print("sin still transient ", str(stat_sin[-1]))
        if stat_tanh[-1] < 0.95:
            print("tanh still transient ", str(stat_tanh[-1]))

        signal = sinus + tanh + rausch

        return sinus, tanh, rausch, signal, stat_sin , stat_tanh

    def plot(self):
        fig, axs = plt.subplots(6, 1)
        sinus, tanh, rausch, signal, stat_sin , stat_tanh = self.generate()
        axs[0].plot(np.arange(0, self.time, self.datalen ** -1), sinus)
        axs[1].plot(np.arange(0, self.time, self.datalen ** -1), stat_sin)
        axs[2].plot(np.arange(0, self.time, self.datalen ** -1), tanh)
        axs[3].plot(np.arange(0, self.time, self.datalen ** -1), stat_tanh)
        axs[4].plot(np.arange(0, self.time, self.datalen ** -1), rausch)
        axs[5].plot(np.arange(0, self.time, self.datalen ** -1), signal)

        plt.show()


def test_transientcheck():
    siggen = signal_generator()
    siggen.plot()
    return 0


def transientcheck():
    return 0


def signalcreator():
    return 0


transientcheck()
