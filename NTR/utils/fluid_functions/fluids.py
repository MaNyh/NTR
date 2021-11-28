class idealgas:
    R = 8.31446261815324

    def __init__(self, **args):
        self.cp = args.get("cp")
        self.cv = args.get("cv")
        self.Rs = args.get("Rs")
        self.kappa = args.get("kappa")
        self.M = args.get("M")

        if self.cp and self.kappa:
            self.cv = self.cp * self.kappa
        elif self.cp and self.cv:
            self.kappa = self.cv/self.cp
        elif self.cv and self.kappa:
            self.cp = self.cv / self.kappa

        if self.M:
            self.Rs = self.R/self.M
        elif self.Rs:
            self.M = self.R /self.Rs

        self.checksettings()

    def checksettings(self):
        assert self.cv == self.cp * self.kappa
        assert self.kappa == self.cv/self.cp
        assert self.cp == self.cv / self.kappa

        assert self.Rs == self.R/self.M
        assert self.M == self.R /self.Rs

