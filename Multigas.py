import hapi
from scipy.optimize import curve_fit
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#hapi.getHelp(hapi.ISO_ID)
hapi.db_begin('Data')
hapi.fetch_by_ids('Mixture', [1, 2, 3, 4, 7, 8, 9, 10], 3593.6, 3594.8)
#hapi.fetch_by_ids('CO2', [7, 8, 9, 10], 3593.6, 3594.8)
#hapi.select('Mixture')

molec, iso, nuij, Sref, E, Slf, Air, nair, Delta = hapi.getColumns('Mixture', ['molec_id',
            'local_iso_id', 'nu', 'Sw', 'elower', 'gamma_self', 'gamma_air', 'n_air', 'delta_air'])
Qref = np.zeros((2, 4))
Q = np.zeros((2, 4))
M = np.zeros((2, 4))
for i in range(4):
    Qref[0, i] = hapi.partitionSum(1, i + 1, 296)
    Qref[1, i] = hapi.partitionSum(2, i + 1, 296)
    Q[0, i] = hapi.partitionSum(1, i + 1, 300)
    Q[1, i] = hapi.partitionSum(2, i + 1, 300)
    M[0, i] = hapi.molecularMass(1, i + 1)
    M[1, i] = hapi.molecularMass(2, i + 1)

T = 300  # Заданная температура
Tref = 296  # Опорная температура в HITRAN
c = 29979245800  # Скорость света
c2 = 1.4387769  # Константа ch/k
Na = 6.02214129e23  # Число Авогадро
kb = 1.380649e-16  # Постоянная Больцмана
pref = 1  # Атмосферное давление
ps = 1  # Парциальное давление метана
p = 1  # Заданное давление
#los = 2.6867811e19
l = 140
Tr = 273

#hitran = pd.read_excel("Lorentz on the Web.xlsx")
#Trans = np.array(hitran.Trans[9000:11001])
nu = [3593.6 + i * 0.001 for i in range(1201)]
khapi = np.zeros((2, 4, len(nu)))
thapi = np.zeros(len(nu))

for i in range(len(nuij)):
    AlphaD = nuij[i] / c * math.sqrt(2 * Na * kb * T * 0.693 / M[molec[i] - 1, iso[i] - 1])
    Gamma = ((Tref / T) ** nair[i]) * (Air[i] * (p - ps) + Slf[i] * ps)
    S = Sref[i] * Qref[molec[i] - 1, iso[i] - 1] / Q[molec[i] - 1, iso[i] - 1] * \
        math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (1 - math.exp(-c2 * nuij[i] / T)) \
        / (1 - math.exp(-c2 * nuij[i] / Tref))
    khapi[molec[i] - 1, iso[i] - 1, :] += S * hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, p * Delta[i], nu)
    #khapi[molec[i] - 1, iso[i] - 1, :] += S * hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, 0, nu)

"""
Trans = - np.log(trans) / l
los,_ = curve_fit(lambda x1, x2, x3, x4, x5, x6, x7, x8, y1, y2, y3, y4, y5, y6, y7, y8: x1 * y1 + 
    x2 * y2 + x3 * y3 + x4 * y4 + x5 * y5 + x6 * y6 + x7 * y7 + x8 * y8, khapi[0, 0, :], khapi[0, 1, :],
    khapi[0, 2, :], khapi[0, 3, :], khapi[1, 0, :], khapi[1, 1, :], khapi[1, 2, :], khapi[1, 3, :],
    Trans, p0=2.4e19, bounds=[2.3e19, 2.5e19])
#print(los)
"""

los = np.array([2.6867811e19, 2.6867811e19, 2.6867811e19, 2.6867811e19, 2.6867811e19, 2.6867811e19, 2.6867811e19, 2.6867811e19])
for i in range(4):
    thapi += los[i] * khapi[0, i, :] + los[i + 4] * khapi[1, i, :]
thapi = np.exp(- l * thapi)
plt.plot(nu, thapi)
plt.show()
