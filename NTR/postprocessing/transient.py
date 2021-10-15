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

        self.transientlimit = 0.95

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
        return np.random.rand(len(self.timesteps))*np.random.rand()

    def generate(self):
        sinus, stat_sin = self.sin_signal(self.sin_lasting)
        tanh, stat_tanh = self.tanh_signal()
        rausch = self.noise()

        if stat_sin[-1] < self.transientlimit:
            print("sin still transient ", str(stat_sin[-1]))
        if stat_tanh[-1] < self.transientlimit:
            print("tanh still transient ", str(stat_tanh[-1]))


        signal = sinus + tanh + rausch

        return sinus, tanh, rausch, signal, stat_sin , stat_tanh

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
    sinus, tanh, rausch, signal, stat_sin , stat_tanh = siggen.generate()
    siggen.plot(sinus, tanh, rausch, signal, stat_sin , stat_tanh)
    return 0


def transientcheck():
    """
    https://link.springer.com/content/pdf/10.1007/s00162-018-0474-0.pdf
    :return:
    """
    return 0



