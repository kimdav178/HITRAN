import math

import copy
import scipy.special
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
ps = 0.005  # Парциальное давление метана
p = 0.005  # Заданное давление 5 мбар
los = 2.6867811e19
l = 5
Tr = 273

hitran = pd.read_excel("Hitran.xlsx")
methane = pd.read_excel("methane.xlsx")
T5mb = hitran.Ts
Ta = hitran.Trans
Delta = methane.Delta  # Сдвиг частоты
nu = methane.nu + (Delta * p)  # Волновые частоты
Sref = methane.S  # Интенсивности
E = methane.E  # Энергии нижних состояний переходов
Air = methane.Air  # Полуширина на полувысоте в воздухе
Slf = methane.Self  # Полуширина на полувысоте при парциальном давлении
nair = methane.Nair  # Коэффициент зависимости nair от температуры

nuij = [6056 + i * 0.001 for i in range(2001)]  # Сетка частот
S = np.zeros((2001))
k = np.zeros((2001))

for i in range(len(nu)):
    AlphaD = nu[i] / c * math.sqrt(2 * Na * kb * T * 0.693 / M)
    sigma = AlphaD / math.sqrt(2 * math.log(2))
    Gamma = ((Tref / T) ** nair[i]) * Slf[i] * ps
    for j in range(2001):
        S[j] = Sref[i] * Qref / Q * math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (
                1 - math.exp(-c2 * nuij[j] / T)) / (1 - math.exp(-c2 * nuij[j] / Tref))
    k += S * scipy.special.voigt_profile(nuij - nu[i], sigma, Gamma)

"""
t = np.zeros((2001))
for i in range(2001):
    t[i] = - math.log(T5mb[i]) / k[i] / l / Tr / p * T
print(np.mean(t))
"""

t = np.zeros((2001))
best = []
x = 1.21e17
min = 1e20
for j in range(2000):
    y = 0
    for i in range(2001):
        t[i] = math.exp(- k[i] * l * x)
        y += (t[i] - T5mb[i]) ** 2
    if y < min:
        min = y
        best = copy.deepcopy(t)
        los = x
    x += 1e12
print(los)
plt.plot(k, (best - T5mb) / best)
plt.show()
