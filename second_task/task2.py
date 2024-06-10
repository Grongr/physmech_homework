import numpy as np
import matplotlib.pyplot as plt


# Functions for conversion to other units of measurement
def evToKelvin(T_ev):
    return T_ev * 11606


def kelvinToEv(T_k):
    return T_k / 11606


def smToEv(smRev):
    return smRev / 8065


def atmToEv(p):
    return p * 101325 * 6.242e12


def gramToEvs2_sm2(m):
    return m * 6.242e11


# Parameters of the problem to be solved
m_p = 1.67262129e-24

M_H = 1.008
M_F = 18.998
M_HF = 20.006

m_H = m_p * M_H  # g
m_F = m_p * M_F  # g
m_HF = m_p * M_HF  # g

g_H = 2
g_F = 6
g_HF = 8

D_HF = 5.8  # eV
w_HF = 4139  # sm^(-1)
B_HF = 20.96  # sm^(-1)

h_ev = 6.528e-16
sigma = 3e-15  # sm^2


def R(T_K, p_atm):
    T_ev = kelvinToEv(T_K)
    p = atmToEv(p_atm)
    multiplier1 = g_H * g_F / g_HF
    multiplier2 = (gramToEvs2_sm2(m_H * m_F / m_HF) * T_ev / (2. * np.pi * h_ev * h_ev)) ** (3. / 2.)
    multiplier3 = 1. - np.exp(-smToEv(w_HF) / T_ev)
    multiplier4 = smToEv(B_HF) / p
    multiplier5 = np.exp(-D_HF / T_ev)
    return multiplier1 * multiplier2 * multiplier3 * multiplier4 * multiplier5


def calcDegOfDissociation(R):
    return np.sqrt(R / (1 + R))


p_atm = np.linspace(0.1, 1., 9)
T_K = np.linspace(1000, 10000, 1000)

# degreeOfDissociation = calcDegOfDissociation(R(T_K, 0.1))
# ratioOfLambda = (D_HF / kelvinToEv(T_K)) ** 2 * degreeOfDissociation * (1 - degreeOfDissociation) / (1 + degreeOfDissociation)
# plt.plot(T_K, ratioOfLambda)

alpha = calcDegOfDissociation(R(T_K, 0.1))
lambda_tr_H = (5 / 6) / sigma * np.sqrt(3 * kelvinToEv(T_K) / gramToEvs2_sm2(m_H))
lambda_tr_F = (5 / 6) / sigma * np.sqrt(3 * kelvinToEv(T_K) / gramToEvs2_sm2(m_F))
lambda_tr_HF = (7 / 6) / sigma * np.sqrt(3 * kelvinToEv(T_K) / gramToEvs2_sm2(m_HF))
plt.plot(T_K, lambda_tr_H + lambda_tr_F + lambda_tr_HF)

multiplier1 = -((7 / 2) * kelvinToEv(T_K) - D_HF)
multiplier2 = 1 / (3 * sigma)
multiplier3 = np.sqrt(3 * kelvinToEv(T_K) / gramToEvs2_sm2(m_HF))
multiplier4 = alpha * (1 - alpha ** 2) / (1 + alpha ** 2)
multiplier5 = D_HF / (kelvinToEv(T_K) ** 2)
lambda_ch = multiplier1 * multiplier2 * multiplier3 * multiplier4 * multiplier5
# plt.plot(T_K, lambda_ch)
# plt.plot(T_K, lambda_tr_H + lambda_tr_F + lambda_tr_HF + lambda_ch)
plt.show()
