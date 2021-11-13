import numpy as np
from scipy.integrate import simps

from NTR.utils.mathfunctions import autocorr, zero_crossings


def integralscales_from_timeseries(mean, fluctations, timesteps):

    autocorrelated = autocorr(fluctations)
    # we are integrating from zero to zero-crossing in the autocorrelation, we need the time to begin with zeros
    # probably the used datasample is not beginning with 0. therefore:
    timesteps -= timesteps[0]
    if len(zero_crossings(autocorrelated)) > 0:
        acorr_zero_crossings = zero_crossings(autocorrelated)[0]
    else:
        print("no zero crossing found, using first minimal value (possibly last timestep). check data quality!")
        acorr_zero_crossings = np.where(autocorrelated == min(autocorrelated))[0][0]

    integral_time_scale = simps(autocorrelated[:acorr_zero_crossings], timesteps[:acorr_zero_crossings])
    integral_length_scale = integral_time_scale * mean

    return integral_time_scale, integral_length_scale
