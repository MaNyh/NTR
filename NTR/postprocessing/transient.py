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

from scipy.integrate import simps
import scipy.fftpack as fftpack


class signal_generator:
    """
    this is a signal-generator

    only parameters for
    """
    # resolution
    timeresolution = 1000  # resolution of times

    transientlimit = 0.95

    def __init__(self):
        # some kind of factor for the duration of an abating signal
        self.tanh_lasting = np.random.randint(0, 100)

        self.sin_lasting = np.random.randint(0, 100)
        self.sin_omega = np.random.rand()

        self.tanh_stationary_ts = np.arctanh(self.transientlimit) * self.tanh_lasting
        # self.tanh_sign = np.random.choice([-1, 1])

        self.sin_stationary_ts = -self.sin_lasting * (1 + np.log(1 - self.transientlimit))

        # as defined in Ries 2018, the signal must be at least two times as long as the transient process
        # todo current calculation of self.time is not according ries2018
        if self.sin_stationary_ts > self.tanh_stationary_ts:
            self.time = 2 * self.sin_stationary_ts
        else:
            self.time = 2 * self.tanh_stationary_ts
        # todo we need to correct the wrong calculated self.time. this is not nice
        self.time *= 2

        # defining the used timesteps (unnesessary, the signal-length could be normed to 1!)
        self.timesteps = np.arange(0, self.time, self.timeresolution ** -1)

        # damping the sine with euler
        abate = np.e ** (-self.timesteps * self.sin_lasting ** -1)
        self.sin_abate = abate / max(abate)

        # values for the 'stationarity'. this can be defined because the function is analytical. but is this correct?
        self.sin_stationary = np.argmax(self.sin_abate < (1 - self.transientlimit))
        self.tanh_stationarity = np.argmax(np.tanh(self.timesteps * self.tanh_lasting ** -1) > self.transientlimit)

        sin_freq = self.sin_omega ** -1
        sig_time = self.time

        sin_timescales = sig_time * sin_freq

        print()
        print("time when stationary")
        sin_statts = int(self.sin_stationary_ts / self.time * len(self.timesteps))
        tanh_statts = int(self.tanh_stationary_ts / self.time * len(self.timesteps))
        print("sin : ", self.sin_stationary_ts, " timestep ", sin_statts)
        print("sin frequency ", sin_freq, " ; timescales in whole dataset ", round(sin_timescales, 2))
        print("tanh : ", self.tanh_stationary_ts, " timestep ", tanh_statts)
        self.stat_time = self.tanh_stationary_ts if self.tanh_stationary_ts > self.sin_stationary_ts else self.sin_stationary_ts
        print("signal stationary at", self.stat_time)
        self.stationarity_relt = round(int(self.stat_time / self.time * len(self.timesteps)) / len(self.timesteps), 1)
        self.stationarity_time = self.stationarity_relt * self.time
        print("this equals a timestep of ", int(self.stat_time / self.time * len(self.timesteps)), " from ",
              len(self.timesteps), " or ", self.stationarity_relt)
        print("time emulated signal ", self.time)

    def tanh_signal(self):
        ans = np.tanh(self.timesteps * self.tanh_lasting ** -1)
        return ans  # * self.tanh_sign

    def sin_signal(self):
        sinus = np.sin(self.timesteps * self.sin_omega) * self.sin_abate
        return sinus * 0.5

    def noise_signal(self):
        mu, sigma = 0, np.random.rand()  # mean and standard deviation
        s = np.random.normal(mu, sigma, size=len(self.timesteps))

        t = self.time
        Dt = t / len(self.timesteps)
        # dt is the timescale of the noisy signal --> emulated length scale!
        dt = t / 1000
        timescale = int(dt / Dt)
        weights = np.repeat(1.0, timescale) / timescale
        out = np.convolve(s, weights, 'same')
        # todo: when a sine is applied, the integral scale fits. but the signal is not that noisy anymore. whats the best solution?
        out /= max(out)
        out += np.sin(self.timesteps * dt ** -1) * 0.5
        out /= max(out)*0.5
        return out

    def generate(self):
        sinus = self.sin_signal()
        tanh = self.tanh_signal()
        rausch = self.noise_signal()

        sin_stats = (-1 + self.sin_abate) * -1
        tanh_stats = np.sinh(self.timesteps * self.tanh_lasting ** -1) / np.cosh(
            self.timesteps * self.tanh_lasting ** -1)

        signal = sinus + tanh + rausch

        return sinus, tanh, rausch, signal, sin_stats, tanh_stats

    def plot(self, sinus, tanh, rausch, signal, stat_sin, stat_tanh):
        fig, axs = plt.subplots(6, 1)

        axs[0].plot(np.arange(0, self.time, self.timeresolution ** -1), sinus, color="orange",
                    label="abating sine signal")
        axs[0].axvline(self.sin_stationary_ts)
        axs[1].plot(np.arange(0, self.time, self.timeresolution ** -1), stat_sin, color="orange",
                    label="stationarity sine signal")
        axs[1].axvline(self.sin_stationary_ts)
        axs[1].fill_between(np.arange(0, self.time, self.timeresolution ** -1), stat_sin, color="orange")
        axs[2].plot(np.arange(0, self.time, self.timeresolution ** -1), tanh, color="blue", label="tanh signal")
        axs[2].axvline(self.tanh_stationary_ts)
        axs[3].plot(np.arange(0, self.time, self.timeresolution ** -1), stat_tanh, color="blue",
                    label="stationarity tanh signal")
        axs[3].axvline(self.tanh_stationary_ts)
        axs[3].fill_between(np.arange(0, self.time, self.timeresolution ** -1), stat_tanh, color="blue")
        axs[4].plot(np.arange(0, self.time, self.timeresolution ** -1), rausch, color="black", label="noise")
        axs[5].plot(np.arange(0, self.time, self.timeresolution ** -1), signal, color="red", label="signal")
        axs[5].axvline(self.stationarity_time)

        for a in axs:
            a.legend(loc="upper right")

        plt.show()


def test_transientcheck(verbose=True):
    sig_gen = signal_generator()

    sinus, tanh, rausch, signal, stat_sin, stat_tanh = sig_gen.generate()

    if verbose:
        sig_gen.plot(sinus, tanh, rausch, signal, stat_sin, stat_tanh)

    timestationarity_ries = transientcheck(signal, sig_gen.timesteps)
    print("transientcheck-rückgabe (ries2018): ", timestationarity_ries)
    assert sig_gen.stat_time < timestationarity_ries, "von transientcheck zurückgegebene zeit der stationarität zu klein"
    return 0


def autocorr(x):
    norm = np.sum(np.array(x) ** 2)
    result = np.correlate(np.array(x), np.array(x), 'full') / norm
    return result[int(len(result) / 2):]


def transientcheck(signal, timesteps):
    """

    :param signal: timeseries
    :return: time_of_stationarity
    """

    second_half_id = int(len(signal) / 2)
    second_half_of_signal = np.copy(signal[second_half_id:])
    second_half_mean = np.mean(second_half_of_signal)
    second_half_of_signal -= second_half_mean
    second_half_timesteps = np.copy(timesteps[second_half_id:])
    autocorrelated = autocorr(second_half_of_signal)

    # we are integrating from zero to zero-crossing in the autocorrelation, we need the time to begin with zeros
    # probably the used datasample is not beginning with 0. therefore:
    timesteps -= timesteps[0]

    if len(zero_crossings(autocorrelated)) > 0:
        acorr_zero_crossings = zero_crossings(autocorrelated)[0]
    else:
        print("no zero crossing found, using first minimal value (possibly last timestep). check data quality!")
        acorr_zero_crossings = np.where(autocorrelated == min(autocorrelated))[0][0]

    integral_scale = simps(autocorrelated[:acorr_zero_crossings], timesteps[:acorr_zero_crossings])

    integrals_window = 30
    time_window = integrals_window * integral_scale
    windows = []

    window_upperlimit = time_window
    windows_rms = []
    windows_mean = []
    window_signal = []
    window_std = []
    for signal_at_time, time in zip(signal, timesteps):
        if time >= window_upperlimit:
            window_upperlimit += time_window

            windows.append(np.array(window_signal))
            window_std.append(np.std(window_signal))
            windows_mean.append(np.mean(window_signal))
            windows_rms.append(np.sqrt(np.sum(windows[-1] ** 2) / len(windows[-1])))

            window_signal = []
        else:
            window_signal.append(signal_at_time)

    #signal_rms = np.sqrt(np.sum(second_half_of_signal ** 2) / len(second_half_of_signal))
    eps_time_mean = np.std(second_half_of_signal) / second_half_mean * (
                2 * integral_scale / (second_half_timesteps[-1] - second_half_timesteps[0])) ** .5
    eps_time_rms = (integral_scale / (second_half_timesteps[-1] - second_half_timesteps[0])) ** .5

    confidence_mean_high = second_half_mean * (1 + 1.96 * eps_time_mean)
    confidence_mean_low = second_half_mean * (1 - 1.96 * eps_time_mean)
    confidence_rms_high = np.mean(windows_rms) * (1 + 1.96 * eps_time_rms)
    confidence_rms_low = np.mean(windows_rms) * (1 - 1.96 * eps_time_rms)

    fig, axs = plt.subplots(2, 1)
    axs[0].plot(windows_mean)
    axs[0].hlines(confidence_mean_high, xmin=0, xmax=len(windows_mean))
    axs[0].hlines(confidence_mean_low, xmin=0, xmax=len(windows_mean))
    axs[1].plot(windows_rms)
    axs[1].hlines(confidence_rms_high, xmin=0, xmax=len(windows_mean))
    axs[1].hlines(confidence_rms_low, xmin=0, xmax=len(windows_mean))
    plt.show()

    return 0


def zero_crossings(data_series):
    zcs = np.where(np.diff(np.sign(data_series)))[0]
    return zcs


def dominant_freq(signal, timesteps):
    times = timesteps
    length = int(len(times) / 2)

    forcez = signal
    Nev = 1
    N = len(signal)
    t = np.linspace(times[length], times[-1], len(signal))
    forcezint = np.interp(t, times, forcez)

    fourier = fftpack.fft(forcezint[Nev - 1:N - 1])
    frequencies = fftpack.fftfreq(forcezint[Nev - 1:N - 1].size, d=t[1] - t[0])
    freq = frequencies[np.argmax(np.abs(fourier))]
    return freq
