import math

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

methane = pd.read_excel("methane.xlsx")
nu = methane.nu  # Волновые частоты
Sref = methane.S  # Интенсивности
E = methane.E  # Энергии нижних состояний переходов
Air = methane.Air  # Полуширина на полувысоте в воздухе
Slf = methane.Self  # Полуширина на полувысоте при парциальном давлении
Delta = methane.Delta  # Сдвиг частоты
nair = methane.Nair  # Коэффициент зависимости nair от температуры

T = 300  # Заданная температура
Tref = 296  # Опорная температура в HITRAN
Qref = 590.52834
Q = 602.86666
c = 29979245800  # Скорость света
c2 = 1.4387769  # Константа ch/k
Na = 6.02214129e23  # Число Авогадро
k = 1.380649e-16  # Постоянная Больцмана
M = 16  # Молярная масса метана
pref = 1  # Атмосферное давление
ps = 0.005  # Парциальное давление метана
p = 0.005  # Заданное давление 5 мбар
los = 2.6867811e19
l = 5
Tr = 273

nuij = [6056 + i * 0.001 for i in range(2001)]  # Сетка частот
ALphaD = [nuij[j] / c * math.sqrt(2 * Na * k * T * 0.693 / M) for j in range(2001)]
fl = [0 for i in range(2001)]
fg = [0 for i in range(2001)]
S = [0 for i in range(2001)]

for i in range(len(nu)):
    Gamma = ((Tref / T) ** nair[i]) * (Air[i] * (p - ps) + Slf[i] * ps)
    for j in range(2001):
        fl[j] += Gamma / math.pi / (Gamma ** 2 + (nu[i] - (nuij[j] + Delta[i] * p)) ** 2)
        fg[j] += math.sqrt(0.693 / math.pi / (ALphaD[j] ** 2)) * math.exp(-((nu[i] - nuij[j]) ** 2) * 0.693 / (ALphaD[j] ** 2))
        S[j] += Sref[i] * Qref / Q * math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (
                    1 - math.exp(-c2 * nuij[j] / T)) / (
                    1 - math.exp(-c2 * nuij[j] / Tref))
voigt = np.convolve(fg, fl, mode='same')
t = []
for i in range(2001):
    t.append(math.exp(- S[i] * voigt[i] * 0.001 * l * los * Tr * p / T))
"""
for i in range(2001):
    K = 0
    for j in range(i+1):
        K += fg[j] * fl[i - j]
    for j in range(i + 1, 2001):
        K += fg[j] * fl[2001 - j + i]
    #t.append(math.exp(- S[i] * K * 0.001 * l * los * Tr * p / T))
    t.append(K)
"""
plt.plot(nuij, t)
plt.show()
