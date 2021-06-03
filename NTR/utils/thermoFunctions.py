def Sutherland_Law(T, As=1.458e-06, Ts=110.4):
    nu = As * (T) ** (0.5) / (1 + Ts / T)
    return nu


def TByIdealGasLaw(rho, p, R=287.058):
    T = p / (rho * R)
    return T
