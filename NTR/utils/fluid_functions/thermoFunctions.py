

def Sutherland_Law(T, As=1.458e-06, Ts=110.4):
    """
    :param T:
    :param As:
    :param Ts:
    :return:
    """
    nu = As * (T) ** (0.5) / (1 + Ts / T)
    return nu


def TByIdealGasLaw(rho, p, R):
    T = p / (rho * R)
    return T

