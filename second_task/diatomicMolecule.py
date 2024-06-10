import numpy as np

planc = 6.582e-16 # [eV * s]

class DiatomicMol:

    def __init__(self, D, B, omega, ga, gb, gab, ma, mb):
        
        self.D = D
        self.B = B
        
        self.omega = omega

        self.ga = ga
        self.gb = gb

        self.gab = gab

        self.ma = ma * 1.672e-24
        self.mb = mb * 1.672e-24 # [gramm]

        self.ma *= 6.242e11
        self.mb *= 6.242e11 # [Ev]

    def ffunc(self, T, P):

        f = self.ga * self.gb / self.gab *                      \
            np.power(self.ma * self.mb / (self.ma + self.mb) *  \
            T / 2 / np.pi / planc / planc, 1.5)                 \
            * (1 - np.exp(-self.omega * 1.2398*10e-4 / T))      \
            * self.B / P                                        \
            * np.exp(-self.D / T)   

        return f

    def dissociationDegree(self, T, P):

        f = self.ffunc(T, P)

        return np.sqrt(f / (f + 1))

    def atomConcentration(self, T, P, alpha):

        res = alpha / (alpha + 1) * P / T

        return res
    
    def moleculConcentration(self, T, P, alpha):

        res = (1 - alpha) / (alpha + 1) * P / T

        return res
 
