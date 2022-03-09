import hapi
from scipy.optimize import curve_fit
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#hapi.getHelp(hapi.partitionSum)
hapi.db_begin('Data')
hapi.fetch('CH4', 6, 1, 6047, 6067)
nuij, Sref, E, Slf, Air, nair, Delta = hapi.getColumns('CH4', ['nu', 'Sw', 'elower', 'gamma_self',
                                                               'gamma_air', 'n_air', 'delta_air'])
Qref = hapi.partitionSum(6, 1, 296)
Q = hapi.partitionSum(6, 1, 300)
M = hapi.molecularMass(6, 1)

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
l = 5
Tr = 273

hitran = pd.read_excel("Lorentz on the Web.xlsx")
Trans = np.array(hitran.Trans[9000:11001])
nu1 = [6056 + i * 0.001 for i in range(2001)]
khapi = np.zeros(len(nu1))

for i in range(len(nuij)):
    AlphaD = nuij[i] / c * math.sqrt(2 * Na * kb * T * 0.693 / M)
    Gamma = ((Tref / T) ** nair[i]) * Slf[i] * ps
    S = Sref[i] * Qref / Q * math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (
            1 - math.exp(-c2 * nuij[i] / T)) / (1 - math.exp(-c2 * nuij[i] / Tref))
    #khapi += S * hapi.PROFILE_VOIGT(nuij_H2O[i], AlphaD, Gamma, p * Delta[i], nu)
    khapi += S * hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, 0, nu1)

los,_ = curve_fit(lambda q, a: np.exp(- a * q * l), khapi, Trans, p0=2.4e19, bounds=[2.3e19, 2.5e19])
print(los)
thapi = np.exp(- los * l * khapi)
plt.plot(nu1, (thapi - Trans) / thapi)
plt.show()
