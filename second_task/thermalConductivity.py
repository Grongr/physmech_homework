import numpy             as np
import matplotlib.pyplot as plt

from diatomicMolecule import DiatomicMol as DM

Boltzman = 8.6173e-5

class ThermalConductivity(DM):

    def temperatureVelocity(self, T, P):

        self.vt_a  = np.sqrt(3 * T / self.ma)
        self.vt_b  = np.sqrt(3 * T / self.mb)
        self.vt_ab = np.sqrt(3 * T / (self.ma + self.mb))


    def calculateConcentrations(self, T, P):

        h = T[1] - T[0]

        self.alpha  = self.dissociationDegree(T, P)
        self.dalpha = np.diff(self.alpha, n=1) / h

        self.atom_conc  = self.atomConcentration(T, P, self.alpha)
        self.mol_conc   = self.moleculConcentration(T, P, self.alpha)

        self.conc = 2 * self.atom_conc + self.mol_conc

        self.datom_conc = np.diff(self.atom_conc, n=1) / h * self.dalpha
        self.dmol_conc  = np.diff(self.mol_conc,  n=1) / h * self.dalpha

        Oa = self.datom_conc[0]
        Om = self.dmol_conc[0]
        self.datom_conc = np.array( [Oa] + list(self.datom_conc) )
        self.dmol_conc  = np.array( [Om] + list(self.dmol_conc) )

    def getConcentrations(self):
        return self.atom_conc, self.mol_conc, self.conc
    
    def freePath(self):
        sig = 3e-15
        # self.la = 1 / (2 * self.atom_conc + self.mol_conc) / sig
        # self.lm = 1 / (2 * self.atom_conc + self.mol_conc) / sig
        self.la = 1 / (self.atom_conc + self.mol_conc) / sig
        self.lm = 1 / (2 * self.atom_conc) / sig

    def energyPerParticle(self, T, P):
        self.ea = self.eb = 5 / 2 * T
        self.em = 7 / 2 * T - self.D

    def thermalConductivity(self, T, P):

        sig = 3e-15

        self.calculateConcentrations(T, P)
        self.temperatureVelocity(T, P)
        self.freePath()
        self.energyPerParticle(T, P)

        # self.lam_tr  = self.atom_conc * 5 / 2 * self.la
        # self.lam_tr *= (self.vt_a + self.vt_b) / 3
        # self.lam_tr += self.mol_conc * 7 / 2 * self.lm / 3 * self.vt_ab

        self.lam_tr  = 0
        self.lam_tr += 5/6 * self.atom_conc * self.la * self.vt_a
        self.lam_tr += 5/6 * self.atom_conc * self.la * self.vt_b
        self.lam_tr += 7/6 * self.mol_conc  * self.lm * self.vt_ab

        Da  = self.la * self.vt_a / 3
        Db  = self.la * self.vt_b / 3
        Dab = self.lm * self.vt_ab / 3

        # self.lam_chem  = 0
        # self.lam_chem += self.ea * Da * self.datom_conc
        # self.lam_chem += self.eb * Db * self.datom_conc
        # self.lam_chem += self.em * Dab * self.dmol_conc
        # self.lam_chem = self.em * Dab * self.dmol_conc

        self.lam_chem  = 1
        self.lam_chem *= -((7/2) * T - self.D)
        self.lam_chem *= 1 / 3 / sig
        self.lam_chem *= np.sqrt(3 * T / (self.ma + self.mb))
        self.lam_chem *= self.alpha * (1 - self.alpha**2) / (1 + self.alpha**2)
        self.lam_chem *= self.D / (T**2)

        self.lam = self.lam_tr + self.lam_chem
        return self.lam


if __name__ == "__main__":

    T = np.linspace(1e3, 1e4, 10000) * Boltzman

    pressure = list( np.linspace(0.1, 1, 7) * 1e5)
    # pressure = list( np.linspace(0.1, 1, 7) * 1e5 * 6.242e12)

    thCon = ThermalConductivity(D = 4.5723, B = 10.59341, omega = 2990.946,
                                ga = 6, gb = 2, gab = 1,
                                ma = 35.45, mb = 1.008)

    fig, ax = plt.subplots(3)
    for p in pressure:
        thCon.thermalConductivity(T, p)
        ax[0].plot(T / Boltzman, thCon.lam_tr)
    # ax[0].plot(T / Boltzman, 6.76e25 * np.sqrt(T) + 1, color='red')
    ax[0].set_ylabel('Transfer thermal conductivity')
    ax[0].grid()

    for p in pressure:
        thCon.thermalConductivity(T, p)
        ax[1].plot(T / Boltzman, thCon.lam_chem)
    ax[1].set_ylabel('Chemical thermal conductivity')
    ax[1].grid()

    for p in pressure:
        lam = thCon.thermalConductivity(T, p)
        ax[2].plot(T / Boltzman, lam)
    # ax[2].plot(T / Boltzman, 6.76e25 * np.sqrt(T), color='red')
    ax[2].set_xlabel('Temperature')
    ax[2].set_ylabel('Thermal conductivity')
    ax[2].grid()

    plt.show()

    # for p in pressure:
    #     thCon.thermalConductivity(T, p)
    #     plt.plot(T / Boltzman, thCon.dalpha)

    # plt.show()
