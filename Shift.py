from scipy.optimize import curve_fit
import scipy.special
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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

nu = [6056 + i * 0.001 for i in range(2001)]  # Сетка частот
S = np.zeros(len(nuij))
k = np.zeros(2001)
fl = np.zeros(2001)
fg = np.zeros(2001)

for i in range(len(nuij)):
    AlphaD = nuij[i] / c * math.sqrt(2 * Na * kb * T * 0.693 / M)
    sigma = AlphaD / math.sqrt(2 * math.log(2))
    Gamma = ((Tref / T) ** nair[i]) * Slf[i] * ps
    S[i] = Sref[i] * Qref / Q * math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (
            1 - math.exp(-c2 * nuij[i] / T)) / (1 - math.exp(-c2 * nuij[i] / Tref))
    for j in range(2001):
        fl[j] = Gamma / math.pi / (Gamma ** 2 + (nu[j] - (nuij[i] + Delta[i] * p)) ** 2)
        fg[j] = math.sqrt(0.693 / math.pi / (AlphaD ** 2)) * math.exp(
            -((nu[j] - nuij[i]) ** 2) * 0.693 / (AlphaD ** 2))
    for j in range(1001):
        K = 0
        for o in range(1001 + j):
            K += fg[o] * fl[1000 + j - o]
        k[j] += K * S[i] * 0.001
    for j in range(1001, 2001):
        K = 0
        for o in range(j - 1000, 2001):
            K += fg[o] * fl[1000 + j - o]
        k[j] += K * S[i] * 0.001
    #k += S[i - 36] * np.convolve(fg, fl, mode='same') * 0.001
    #k += S[i - 36] * scipy.special.voigt_profile(nu - nuij[i], sigma, Gamma)

los,_ = curve_fit(lambda x, a: np.exp(- a * x * l), k, Trans, p0=3e17, bounds=[3e16, 1e25])
print(los)
t = np.exp(- los * l * k)
#plt.plot(nu[4000:6001], (t - Trans) / t)
plt.plot(nu, t, nu, Trans)
plt.show()
