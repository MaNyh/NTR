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

class idealgas:
    R = 8.31446261815324
    def __init__(self,**args):
        self.cp = args.get("cp")
        self.cv = args.get("cv")
        self.Rs = args.get("Rs")
        self.kappa = args.get("kappa")
        self.M = args.get("M")
        self.As = args.get("As")

        if self.cp and self.kappa:
            self.cv =  self.cp * self.kappa
        elif self.cp and self.cv:
            self.kappa = self.cp/self.cv
        elif self.cv and self.kappa:
            self.cp = self.cv / self.kappa

        if self.M:
            self.Rs = self.R/self.M
        elif self.Rs:
            self.M = self.R /self.Rs

    def checksettings(self):

        return 0
