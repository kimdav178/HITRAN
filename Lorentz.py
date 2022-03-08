import math

from scipy.optimize import curve_fit
import scipy.special
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import copy

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
Delta = methane.Delta[36:174]  # Сдвиг частоты
nuij = methane.nu[36:174] + Delta * p   # Волновые частоты переходов
Sref = methane.S[36:174]  # Интенсивности
E = methane.E[36:174]  # Энергии нижних состояний переходов
Air = methane.Air[36:174]  # Полуширина на полувысоте в воздухе
Slf = methane.Self[36:174]  # Полуширина на полувысоте при парциальном давлении
nair = methane.Nair[36:174]  # Коэффициент зависимости nair от температуры

nu = [6052 + i * 0.001 for i in range(10001)]  # Сетка частот
S = np.zeros((138))
k = np.zeros((10001))

for i in range(36, 174):
    AlphaD = nuij[i] / c * math.sqrt(2 * Na * kb * T * 0.693 / M)
    sigma = AlphaD / math.sqrt(2 * math.log(2))
    Gamma = ((Tref / T) ** nair[i]) * Slf[i] * ps
    S[i - 36] = Sref[i] * Qref / Q * math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (
            1 - math.exp(-c2 * nuij[i] / T)) / (1 - math.exp(-c2 * nuij[i] / Tref))
    k += S[i - 36] * scipy.special.voigt_profile(nu - nuij[i], sigma, Gamma)

los,_ = curve_fit(lambda x, a: np.exp(- a * x * l), k[4000:6001], Trans, p0=2.4e19, bounds=[2e19, 3e19])
print(los)
t = np.exp(- los * l * k[4000:6001])
plt.plot(nu[4000:6001], (t - Trans) / t)
plt.show()


"""
t = np.zeros((20001))
for i in range(20001):
    t[i] = - math.log(T5mb[i]) / k[i] / l / Tr / p * T
print(np.mean(t))

t = np.zeros((2001))
best = []
x = 2.4461e19
min = 1e20

for j in range(20):
    y = 0
    for i in range(2001):
        t[i] = math.exp(- k[4000 + i] * l * x)
        y += (t[i] - Trans[9000 + i]) ** 2
    if y < min:
        min = y
        best = copy.deepcopy(t)
        los = x
    x += 1e14
print(los)
plt.plot(nuij[4000:6001], (best - Trans) / best)
plt.show()
"""
