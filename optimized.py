from scipy.optimize import curve_fit
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.special

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

hitran = pd.read_excel("Discrete.xlsx")
methane = pd.read_excel("Lorentz.xlsx")
Trans = np.array(hitran.Trans)
Delta = methane.Delta  # Сдвиг частоты
nuij = methane.nu   # Волновые частоты переходов
Sref = methane.S  # Интенсивности
E = methane.E  # Энергии нижних состояний переходов
Air = methane.Air  # Полуширина на полувысоте в воздухе
Slf = methane.Self  # Полуширина на полувысоте при парциальном давлении
nair = methane.Nair  # Коэффициент зависимости nair от температуры

z = 2000001
x = 400000
nu = [6047 + i * 20 / (z - 1) for i in range(z)]  # Сетка частот
k = np.zeros(x + 1)
kspecial = np.zeros(z)
fl = np.zeros(z)
fg = np.zeros(z)

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
    #kspecial += S * scipy.special.voigt_profile(nu - nuij_H2O[i], sigma, Gamma)
    n1 = np.argmax(fg)
    if n1 > round((z - 1) / 2) and 2 * n1 < z + round((z - 1 - x) / 2):
        fl2 = fl[round((z - 1 - x) / 2):2 * n1 - round((z - 1 - x) / 2) + 1]
        fg2 = fg[round((z - 1 - x) / 2):2 * n1 - round((z - 1 - x) / 2) + 1]
        k += S * np.convolve(fg2, fl2, mode='same')[:x + 1] * 20 / (z - 1)
    if n1 <= round((z - 1) / 2) and 2 * n1 >= round((z - 1 + x) / 2):
        fl2 = fl[2 * n1 - round((z - 1 + x) / 2):round((z - 1 + x) / 2) + 1]
        fg2 = fg[2 * n1 - round((z - 1 + x) / 2):round((z - 1 + x) / 2) + 1]
        k += S * np.convolve(fg2, fl2, mode='same')[z - 1 - 2 * n1:z + x - 2 * n1] * 20 / (z - 1)

los,_ = curve_fit(lambda q, a: np.exp(- a * q * l), k[round(x / 2) - round(z / 20):round(x / 2) + round(z / 20) + 1], Trans, p0=2.4e19, bounds=[2.3e19, 2.5e19])
print(los)
t = np.exp(- los * l * k[round(x / 2) - round(z / 20):round(x / 2) + round(z / 20) + 1])
#tspecial = np.exp(- los * l * kspecial[round(z / 2) - 1000:round(z / 2) + 1001])
plt.plot(nu[round((z - 1) / 2) - round(z / 20):round((z - 1) / 2) + round(z / 20) + 1], (t - Trans) / Trans)
#plt.plot(nu[round((z - 1) / 2) - 1000:round((z - 1) / 2) + 1001], t, nu[round((z - 1) / 2) - 1000:round((z - 1) / 2) + 1001], Trans)
plt.show()
