class FluidCoeffs:
    # https://www.chemie.de/lexikon/Universelle_Gaskonstante.html
    R = 8.314472

    def __init__(self):
        self.kappa = None
        self.R_L = None
        self.mu = None
        self.M = None
        self.cp = None
        self.p_k = None
        self.As = None
        self.l = None
        self.Ts = None

    def set_kappa(self, kappa):
        """
        :param kappa:heat capacity ratio [ ]
        :return: None
        saves value to object (self)
        """
        assert type(kappa) == float, "kappa needs to be a float"
        self.kappa = kappa

    def set_R_L(self, R_L):
        """
        :param R_L:specific gas constant [J / ( kg K )]
        :return: None
        saves value to object (self)
        """
        assert type(R_L) == float, "R_L needs to be a float"
        self.R_L = R_L

    def set_mu(self, mu):
        """
        :param mu:dynamic viscosity [Pa s]
        :return: None
        saves value to object (self)
        """
        assert type(mu) == float, "mu needs to be a float"
        self.mu = mu

    def set_M(self, M):
        """
        :param M:molar weight [kg/mol]
        :return: None
        saves value to object (self)
        """
        assert type(M) == float, "M needs to be a float"
        self.M = M
        self.R_L = self.R / self.M

    def set_cp(self, cp):
        """
        :param cp:specific heat [J/(kg·K)]
        :return: None
        saves value to object (self)
        """
        assert type(cp) == float, "cp needs to be a float"
        self.cp = cp

    def set_p_k(self, p_k):
        """
        :param p_k:Kinematic pressure (p/rho)[m2/s2]
        :return: None
        saves value to object (self)
        """
        assert type(p_k) == float, "p_k needs to be a float"
        self.p_k = p_k

    def set_As(self, As):
        """
        :param As:??? [ ]
        :return: None
        saves value to object (self)
        """
        assert type(As) == float, "As needs to be a float"
        self.As = As

    def set_l(self, l):
        """
        :param l:??? [ ]
        :return: None
        saves value to object (self)
        """
        assert type(l) == float, "l needs to be a float"
        self.l = l

    def set_Ts(self, Ts):
        """
        :param Ts:??? [ ]
        :return: None
        saves value to object (self)
        """
        assert type(Ts) == float, "Ts needs to be a float"
        self.Ts = Ts
