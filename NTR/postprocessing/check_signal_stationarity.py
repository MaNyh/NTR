import numpy as np
import matplotlib.pyplot as plt

from NTR.postprocessing.integralscales_from_signal import integralscales_from_timeseries

"""
this module is supposed to return a timesstamp from a time-series, that equals the time when the transient signal ends

Ries 2018
https://link.springer.com/content/pdf/10.1007/s00162-018-0474-0.pdf

numerical solutions have a transient behaviour at the initial process
it is assumed that this initial transient can be reproduced by a sine, tanh and a noise-function
with these functions given, we can analytically define where the transient process ends
"""


class signal_generator:
    """
    this is a signal-generator

    only parameters for
    """
    # resolution
    timeresolution = 100  # resolution of times

    transientlimit = 0.95

    def __init__(self):
        # some kind of factor for the duration of an abating signal
        self.tanh_lasting = np.random.randint(0, 100)
        self.sin_lasting = np.random.randint(0, 100)

        self.sin_omega = np.random.rand()

        self.tanh_stationary_ts = np.arctanh(self.transientlimit) * self.tanh_lasting
        self.sin_stationary_ts = -self.sin_lasting * (1 + np.log(1 - self.transientlimit))

        # as defined in Ries 2018, the signal must be at least two times as long as the transient process
        # todo current calculation of self.time is not according ries2018
        if self.sin_stationary_ts > self.tanh_stationary_ts:
            self.time = 2 * self.sin_stationary_ts
        else:
            self.time = 2 * self.tanh_stationary_ts

        # time equals approx 300 timescales, adjust to approx 1000
        self.time *= 4
        # defining the used timesteps (unnesessary, the signal-length could be normed to 1!)
        self.timesteps = np.arange(0, self.time, self.timeresolution ** -1)

        # damping the sine with euler
        abate = np.e ** (-self.timesteps * self.sin_lasting ** -1)
        self.sin_abate = abate / max(abate)

        # values for the 'stationarity'. this can be defined because the function is analytical. but is this correct?
        self.sin_stationary = np.argmax(self.sin_abate < (1 - self.transientlimit))
        self.tanh_stationarity = np.argmax(np.tanh(self.timesteps * self.tanh_lasting ** -1) > self.transientlimit)


    def tanh_signal(self):
        ans = np.tanh(self.timesteps * self.tanh_lasting ** -1)
        return ans  # * self.tanh_sign

    def sin_signal(self):
        sinus = np.sin(self.timesteps * self.sin_omega) * self.sin_abate
        return sinus #* 0.5

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
        out += np.sin(self.timesteps * dt ** -1) * 0.25
        out /= max(out) * 2
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
        axs[3].fill_between(np.arange(0, self.time, self.timeresolution ** -1), stat_tanh, color="blue")
        axs[4].plot(np.arange(0, self.time, self.timeresolution ** -1), rausch, color="black", label="noise")
        axs[5].plot(np.arange(0, self.time, self.timeresolution ** -1), signal, color="red", label="signal")

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
    assert False, "Funktion nicht vollständig"
    return 0


def transientcheck(signal, timesteps):
    """

    :param signal: timeseries
    :return: time_of_stationarity
    """

    second_half_id = int(len(signal) / 2)
    second_half_of_signal = np.copy(signal[second_half_id:])

    second_half_mean = np.mean(second_half_of_signal)
    second_half_of_signal_fluctations = second_half_of_signal-second_half_mean
    second_half_timesteps = np.copy(timesteps[second_half_id:])

    integral_time_scale, integral_length_scale = integralscales_from_timeseries(second_half_mean, second_half_of_signal_fluctations, timesteps)

    integrals_window = 30
    time_window = integrals_window * integral_time_scale
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

    #todo: find a suitable name for the following var
    eps_timescale = (integral_time_scale / (second_half_timesteps[-1] - second_half_timesteps[0]))

    eps_time_mean = np.std(second_half_of_signal_fluctations) / second_half_mean * (2 *eps_timescale) ** .5
    eps_time_rms = eps_timescale ** .5

    no_windows_mean = len(windows_mean)
    no_windows_rms = len(windows_rms)

    confidence_mean_high = second_half_mean * (1 + 1.96 * eps_time_mean)
    confidence_mean_low = second_half_mean * (1 - 1.96 * eps_time_mean)
    confidence_rms_high = np.mean(windows_rms) * (1 + 1.96 * eps_time_rms)
    confidence_rms_low = np.mean(windows_rms) * (1 - 1.96 * eps_time_rms)

    mean_in_ci_range = [True if (confidence_mean_high >= i >= confidence_mean_low) else False for i in windows_mean]
    rms_in_ci_range = [True if (confidence_rms_high >= i >= confidence_rms_low) else False for i in windows_rms]

    mean_stationarity_fraction = np.array(
        [sum(mean_in_ci_range[pt:]) / (no_windows_mean - pt) for pt in range(len(mean_in_ci_range))])
    rms_stationarity_fraction = np.array(
        [sum(rms_in_ci_range[pt:]) / (no_windows_rms - pt) for pt in range(len(rms_in_ci_range))])

    statwindow_mean = np.where(mean_stationarity_fraction > 0.95)[0][0]
    statwindow_rms = np.where(rms_stationarity_fraction > 0.95)[0][0]

    fig, axs = plt.subplots(2, 1)

    axs[0].plot(windows_mean)
    axs[0].hlines(confidence_mean_high, xmin=0, xmax=no_windows_mean)
    axs[0].hlines(confidence_mean_low, xmin=0, xmax=no_windows_mean)
    axs[1].plot(windows_rms)
    axs[1].hlines(confidence_rms_high, xmin=0, xmax=no_windows_rms)
    axs[1].hlines(confidence_rms_low, xmin=0, xmax=no_windows_rms)
    plt.show()

    return 0


