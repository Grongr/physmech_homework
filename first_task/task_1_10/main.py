import diatomicMolecule as dm

import matplotlib.pyplot as plt
import numpy as np

def graphDissocDegree(T, alphas, pres):

    plt.grid()
    plt.xlabel("T, K.")
    plt.ylabel("Степень диссоциации.")

    for i in range(len(pres)):
        plt.plot(T, alphas[i], label=str(pres[i]) + "атм")

    plt.legend()
    plt.show()

def calculateDissocDegree(T, pres, diss):

    alphas = []

    for p in pres:
        alphas.append(diss.dissociationDegree(T, p*1e6))

    return alphas

def graphAtomConcentration(T, concs, pres):

    plt.grid()
    plt.xlabel("T, K.")
    plt.ylabel("Концентрация атомов.")

    for i in range(len(pres)):
        plt.plot(T, concs[i], label=str(pres[i]) + "атм")

    plt.legend()
    plt.show()

def calculateAtomConcentration(T, pres, diss, alphas):

    concs = []

    for i in range(len(pres)):
        concs.append(diss.atomConcentration(T, pres[i]*1e6, alphas[i]))

    return concs

def graphMoleculConcentration(T, mconcs, pres):

    plt.grid()
    plt.xlabel("T, K.")
    plt.ylabel("Концентрация молекул.")

    for i in range(len(pres)):
        plt.plot(T, mconcs[i], label=str(pres[i]) + "атм")

    plt.legend()
    plt.show()

def calculateMoleculeConcentration(T, pres, diss, alphas):

    mconcs = []

    for i in range(len(pres)):
        mconcs.append(diss.moleculConcentration(T, pres[i]*1e6, alphas[i]))

    return mconcs



if __name__ == "__main__":

    # HCl
    dissN = dm.DiatomicMol(D = 4.5723, B = 10.59341, omega = 2990.946,
                           ga = 6, gb = 2, gab = 1,
                           ma = 35.45, mb = 1.008)

    T = np.linspace(200, 15000, 1000) * 8.6173e-5
    pres = list(np.linspace(0.1, 2, 10))

    alphas = calculateDissocDegree(T, pres, dissN)
    concs  = calculateAtomConcentration(T, pres, dissN, alphas)
    mconcs = calculateMoleculeConcentration(T, pres, dissN, alphas)

    graphDissocDegree(T / 8.6173e-5, alphas, pres)
    graphAtomConcentration(T / 8.6173e-5, concs, pres)
    graphMoleculConcentration(T / 8.6173e-5, mconcs, pres)

