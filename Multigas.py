import hapi
from scipy.optimize import curve_fit
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def exponential_comb_iso1(x, y1, y2):
    return np.exp((- x[0] * y1 - x[1] * y2) * l)


def lin_comb(x, y1, y2, y3, y4, y5, y6, y7, y8):
    return x[0] * y1 + x[1] * y2 + x[2] * y3 + x[3] * y4 + x[4] * y5 + x[5] * y6 + x[6] * y7 + x[7] * y8


#hapi.getHelp(hapi.ISO_ID)
hapi.db_begin('Data')
#hapi.fetch_by_ids('Multigas', [1, 2, 3, 4, 7, 8, 9, 10], 3583.7, 3604.7)
#hapi.fetch_by_ids('Multigas', [1, 7], 3583.7, 3604.7)
nu_0 = 3562
nu_end = 3564
hapi.fetch_by_ids('Check', [1, 7], nu_0 - 10, nu_end + 10)
#hapi.select('Multigas')

molec, iso, nuij, Sref, E, Slf, Air, nair, Delta = hapi.getColumns('Check', ['molec_id',
            'local_iso_id', 'nu', 'Sw', 'elower', 'gamma_self', 'gamma_air', 'n_air', 'delta_air'])
molec_num = 2
iso_num = 1
T = 300  # Заданная температура
Tref = 296  # Опорная температура в HITRAN
c = 29979245800  # Скорость света
c2 = 1.4387769  # Константа ch/k
Na = 6.02214129e23  # Число Авогадро
kb = 1.380649e-16  # Постоянная Больцмана
pref = 1  # Атмосферное давление
ps = np.zeros((molec_num, iso_num))  # Парциальное давление метана
p = 1  # Заданное давление
#los = 2.6867811e19
l = 140
Tr = 273
concentration = [0.008, 0.00042]
for j in range(iso_num):
    ps[0, j] = p * concentration[0] * hapi.abundance(1, j + 1)
    ps[1, j] = p * concentration[1] * hapi.abundance(2, j + 1)

Qref = np.zeros((molec_num, iso_num))
Q = np.zeros((molec_num, iso_num))
M = np.zeros((molec_num, iso_num))
for i in range(iso_num):
    Qref[0, i] = hapi.partitionSum(1, i + 1, Tref)
    Qref[1, i] = hapi.partitionSum(2, i + 1, Tref)
    Q[0, i] = hapi.partitionSum(1, i + 1, T)
    Q[1, i] = hapi.partitionSum(2, i + 1, T)
    M[0, i] = hapi.molecularMass(1, i + 1)
    M[1, i] = hapi.molecularMass(2, i + 1)

#hitran = pd.read_excel("Transmittance spectrums/Mixture.xlsx")
hitran = pd.read_excel("3562-3564/3562.xlsx")
Trans = np.array(hitran.Trans)
nu = [nu_0 + i * 0.001 for i in range((nu_end - nu_0) * 1000 + 2)]
khapi = np.zeros((molec_num, iso_num, len(nu)))
thapi = np.zeros(len(nu))

for i in range(len(nuij)):
    AlphaD = nuij[i] / c * math.sqrt(2 * Na * kb * T * 0.693 / M[molec[i] - 1, iso[i] - 1])
    Gamma = ((Tref / T) ** nair[i]) * (Air[i] * (p - ps[molec[i] - 1, iso[i] - 1]) + Slf[i] * ps[molec[i] - 1, iso[i] - 1])
    S = Sref[i] * Qref[molec[i] - 1, iso[i] - 1] / Q[molec[i] - 1, iso[i] - 1] * \
        math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (1 - math.exp(-c2 * nuij[i] / T)) \
        / (1 - math.exp(-c2 * nuij[i] / Tref))
    khapi[molec[i] - 1, iso[i] - 1, :] += S * hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, p * Delta[i], nu)
    #khapi[molec[i] - 1, iso[i] - 1, :] += S * hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, 0, nu)

x = [khapi[0, 0, :], khapi[1, 0, :]]
#x = [khapi[0, 0, :], khapi[0, 1, :], khapi[0, 2, :], khapi[0, 3, :], khapi[1, 0, :], khapi[1, 1, :], khapi[1, 2, :], khapi[1, 3, :]]
#Trans_linear = - np.log(Trans) / l
#los, cov = curve_fit(lin_comb, x, Trans_linear, p0=1e13 * np.ones(8), bounds=[1e10, 1e22])
#print(los)
los, cov = curve_fit(exponential_comb_iso1, x, Trans, p0=1e13 * np.ones(2), bounds=[1e8, 1e22])
print(los)

for i in range(iso_num):
    thapi += los[i] * khapi[0, i, :] + los[i + iso_num] * khapi[1, i, :]
thapi = np.exp(- l * thapi)

plt.plot(nu, thapi, nu, Trans)
plt.show()

plt.plot(nu, (thapi - Trans) / Trans)
plt.show()
