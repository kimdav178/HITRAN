import math

import matplotlib.pyplot as plt
import pandas as pd

methane = pd.read_excel("methane.xlsx")
nu = methane.nu  # Волновые частоты
Sref = methane.S  # Интенсивности
E = methane.E  # Энергии нижних состояний переходов
Air = methane.Air  # Полуширина на полувысоте в воздухе
Slf = methane.Self  # Полуширина на полувысоте при парциальном давлении
Delta = methane.Delta  # Сдвиг частоты
nair = methane.Nair  # Коэффициент зависимости nair от температуры

nuij = [6056 + i * 0.001 for i in range(2000)]  # Шаг сетки частот
T = 300  # Заданная температура
Tref = 296  # Опорная температура в HITRAN
Qref = 590.52834
Q = 602.86666
c = 29979245800  # Скорость света
c2 = 1.4387769  # Константа ch/k
Na = 6.02214129e23  # Число Авогадро
k = 1.380649e-16  # Постоянная Больцмана
M = 16  # Молярная масса метана
pref = 1000000  # Атмосферное давление
ps = 28.57  # Парциальное давление метана
p = 5000  # Заданное давление 5 мбар

ALphaD = nuij / c * math.sqrt(2 * Na * k * T * 0.693 / M)
Gamma = []
fl = []
fg = []
S = []
K = []
for i in range(len(nu)):
    Gamma.append(((Tref / T) ** nair[i]) * (Air[i] * (p - ps) + Slf[i] * ps))
    for j in range(2000):
        fl.append(Gamma[i] / math.pi / (Gamma[i] ** 2 + (nu[i] - (nuij[j] + Delta[i] * p)) ** 2))
        fg.append(math.sqrt(0.693 / math.pi / (ALphaD ** 2)) * math.exp(-(nu[i] - nuij) ** 2 * 0.693 / (ALphaD ** 2)))
        S.append(
            Sref[i] * Qref / Q * math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (1 - math.exp(-c2 * nuij / T)) / (
                        1 - math.exp(-c2 * nuij / Tref)))
        K.append(S[i] * fl[i])
plt.plot(nu, K)
plt.show()
