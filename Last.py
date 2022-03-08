import hapi
from scipy.optimize import curve_fit
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.special

hapi.db_begin('Data')
hapi.fetch('CH4', 6, 1, 6047, 6067)
#hapi.getHelp(hapi.absorptionCoefficient_Voigt)

nuhapi, coef = hapi.absorptionCoefficient_Voigt(((6, 1),), 'CH4', WavenumberRange=[6047, 6067],
    WavenumberStep=0.001, HITRAN_units=False, GammaL='gamma_self', Environment={'p':1., 'T':300.})

nuhapi, transhapi = hapi.transmittanceSpectrum(nuhapi, coef, Environment={'l':5.})
#plt.plot(nu[9000:11001], transhapi[9000:11001])
#plt.show()

T = 300  # Заданная температура
Tref = 296  # Опорная температура в HITRAN
Qref = 590.52834
Q = 602.86666
c = 29979245800  # Скорость света
c2 = 1.4387769  # Константа ch/k
Na = 6.02214129e23  # Число Авогадро
kb = 1.380649e-16  # Постоянная Больцмана
M = 16  # Молярная масса метана
pref = 1  # Атмосферное давление
ps = 1  # Парциальное давление метана
p = 1  # Заданное давление
#los = 2.6867811e19
l = 5
Tr = 273

hitran = pd.read_excel("Lorentz on the Web.xlsx")
methane = pd.read_excel("Lorentz.xlsx")
Trans = np.array(hitran.Trans[9000:11001])
Delta = methane.Delta  # Сдвиг частоты
nuij = methane.nu   # Волновые частоты переходов
Sref = methane.S  # Интенсивности
E = methane.E  # Энергии нижних состояний переходов
Air = methane.Air  # Полуширина на полувысоте в воздухе
Slf = methane.Self  # Полуширина на полувысоте при парциальном давлении
nair = methane.Nair  # Коэффициент зависимости nair от температуры

z = 20001
x = 4000
nu = [6047 + i * 20 / (z - 1) for i in range(z)]  # Сетка частот
k = np.zeros(x + 1)
khapi = np.zeros(round((z - 1) / 10) + 1)
kS = np.zeros(round((z - 1) / 10) + 1)
#kspecial = np.zeros(z)
fl = np.zeros(z + 1)
fg = np.zeros(z + 1)

for i in range(len(nuij)):
    AlphaD = nuij[i] / c * math.sqrt(2 * Na * kb * T * 0.693 / M)
    sigma = AlphaD / math.sqrt(2 * math.log(2))
    Gamma = ((Tref / T) ** nair[i]) * Slf[i] * ps
    S = Sref[i] * Qref / Q * math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (
            1 - math.exp(-c2 * nuij[i] / T)) / (1 - math.exp(-c2 * nuij[i] / Tref))
    for j in range(z):
        fl[j] = Gamma / math.pi / (Gamma ** 2 + (nu[j] - (nuij[i])) ** 2)
        fg[j] = math.sqrt(0.693 / math.pi / (AlphaD ** 2)) * math.exp(
            -((nu[j] - nuij[i]) ** 2) * 0.693 / (AlphaD ** 2))
    #khapi = hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, p * Delta[i], nu[round(z / 2) - round(z / 20):round(z / 2) + round(z / 20) + 1])
    #kspecial += S * scipy.special.voigt_profile(nu - nuij[i], sigma, Gamma)
    #kS += S * khapi
    n1 = np.argmax(fg)
    if n1 > round((z - 1) / 2) and 2 * n1 < z + round((z - 1 - x) / 2):
        fl2 = fl[round((z - 1 - x) / 2):2 * n1 - round((z - 1 - x) / 2) + 1]
        fg2 = fg[round((z - 1 - x) / 2):2 * n1 - round((z - 1 - x) / 2) + 1]
        k += S * np.convolve(fg2, fl2, mode='same')[:x + 1] * 20 / (z - 1)
    if n1 <= round((z - 1) / 2) and 2 * n1 >= round((z - 1 + x) / 2):
        fl2 = fl[2 * n1 - round((z - 1 + x) / 2):round((z - 1 + x) / 2) + 1]
        fg2 = fg[2 * n1 - round((z - 1 + x) / 2):round((z - 1 + x) / 2) + 1]
        k += S * np.convolve(fg2, fl2, mode='same')[z - 1 - 2 * n1:z + x - 2 * n1] * 20 / (z - 1)

#los,_ = curve_fit(lambda q, a: np.exp(- a * q * l), k[round(x / 2) - round(z / 20):round(x / 2) + round(z / 20) + 1], Trans, p0=2.4e19, bounds=[2.3e19, 2.5e19])
#print(los)
los = 2.44584129e19
t = np.exp(- los * l * k[round(x / 2) - round(z / 20):round(x / 2) + round(z / 20) + 1])
#tS = np.exp(- los * l * kS)
#tspecial = np.exp(- los * l * kspecial[round(z / 2) - 1000:round(z / 2) + 1001])
plt.plot(nu[round((z - 1) / 2) - round(z / 20):round((z - 1) / 2) + round(z / 20) + 1], (t - transhapi[9000:11001]) / t)
#plt.plot(nu[round((z - 1) / 2) - 1000:round((z - 1) / 2) + 1001], t, nu[round((z - 1) / 2) - 1000:round((z - 1) / 2) + 1001], Trans)
plt.show()

#f = hapi.PROFILE_VOIGT(6057, 5, 7, -1, 6057)
#print(f)
#x = np.arange(len(f))
#plt.plot(x, f)
#plt.show()
