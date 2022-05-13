import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from matplotlib.widgets import Slider
import hapi
from scipy.optimize import curve_fit
import math


def lin_comb(x, y1, y2, y3, y4, y5, y6, y7, y8):
    return x[0] * y1 + x[1] * y2 + x[2] * y3 + x[3] * y4 + x[4] * y5 + x[5] * y6 + x[6] * y7 + x[7] * y8


def baseline(a, b, c):
    return a * nu_period + b + c * (nu_period ** 2)


def normalized(a, b, c):
    return exper_period / (a * nu_period + b + c * (nu_period ** 2))


def nu_adjust(nu_0):
    return np.arange(nu_0, nu_0 - 0.0005 + len(exper_period) * 0.001, 0.001)


nu_step = 0.05
beginning = 23  # Номер первого пика
end = 960
a = beginning
pointmax = [beginning]
valuemax = [0]
fabri = [float(i) for i in open("3562-3564/fp1").read().split()]
exper = [float(i) for i in open("3562-3564/1").read().split()]
hitran = pd.read_excel("3562-3564/3562.xlsx")
Trans = np.array(hitran.Trans)
nu_model = np.array(hitran.nu)

max_fabri = max(fabri)  #Обычно этих строк не должно быть
max_exp = max(exper)
for i in range(len(fabri)):
    fabri[i] = max_fabri - fabri[i]
    exper[i] = max_exp - exper[i]


plt.plot(np.arange(len(fabri)), fabri)
plt.show()
plt.plot(np.arange(len(exper)), exper)
plt.show()

for i in range(beginning, end):
    if ((fabri[i] >= max(fabri[i - 2], fabri[i - 1], fabri[i + 1], fabri[i + 2])) or
        (fabri[i] <= min(fabri[i - 3], fabri[i - 2], fabri[i - 1], fabri[i + 1], fabri[i + 2]))) and i > a + 2:
        l = i - a
        pointmax.append(i)
        valuemax.append(valuemax[len(valuemax) - 1] + nu_step / 2)
        a = i

#for i in range(beginning, a):      Обычно эти строки должны быть
#    exper[i] -= exper[end + 20]
exper1_cut = exper[beginning:a]
interp_grid = np.arange(beginning, pointmax[len(pointmax) - 1])
interp_nu = interp1d(pointmax, valuemax, kind='quadratic')
nu_interp = interp_nu(interp_grid)

plt.plot(np.arange(len(nu_interp)), nu_interp)
plt.show()

# Разворот и интерполяция спектра на периодическую сетку частот
nu_reversed = [nu_interp[len(nu_interp) - i - 1] for i in range(len(nu_interp))]
exper_reversed = [exper1_cut[len(exper1_cut) - i - 1] for i in range(len(exper1_cut))]
nu_period = np.arange(nu_reversed[0], 0, -0.001)
interp_exper = interp1d(nu_reversed, exper_reversed, kind='quadratic')
exper_period = interp_exper(nu_period)


# Подбор базовой линии
def sliders_on_changed(val):
    line.set_ydata(baseline(a_slider.val, b_slider.val, c_slider.val))
    exp.set_xdata(nu_adjust(nu_slider.val))
    exp.set_ydata(normalized(a_slider.val, b_slider.val, c_slider.val))
    fig.canvas.draw_idle()


fig = plt.figure()
ax = fig.add_subplot(121)
bx = fig.add_subplot(122)
fig.subplots_adjust(bottom=0.25)
a_0 = 10000
b_0 = 5000
c_0 = 0
nu_0 = 3562
[line, approximation] = ax.plot(nu_period, baseline(a_0, b_0, c_0), nu_period, exper_period, linewidth=2)
[exp, model] = bx.plot(nu_adjust(nu_0), normalized(a_0, b_0, c_0), nu_model, Trans)
a_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03])
a_slider = Slider(a_slider_ax, 'a', 5000, 20000, valinit=a_0)
b_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03])
b_slider = Slider(b_slider_ax, 'b', 0, 10000, valinit=b_0)
c_slider_ax = fig.add_axes([0.25, 0.05, 0.65, 0.03])
c_slider = Slider(c_slider_ax, 'c', -5000, 5000, valinit=c_0)
nu_slider_ax = fig.add_axes([0.25, 0, 0.65, 0.03])
nu_slider = Slider(nu_slider_ax, 'nu', nu_0 - 1, nu_0 + 1, valinit=nu_0)
nu_slider.on_changed(sliders_on_changed)
a_slider.on_changed(sliders_on_changed)
b_slider.on_changed(sliders_on_changed)
c_slider.on_changed(sliders_on_changed)
plt.show()
print(a_slider.val, b_slider.val, c_slider.val)
for i in range(len(exper_period)):
    exper_period[i] = exper_period[i] / (a_slider.val * nu_period[i] + b_slider.val + c_slider.val * (nu_period[i] ** 2))
print(nu_slider.val)
nu = nu_adjust(nu_slider.val)
nu_end = nu[len(nu) - 1]

#for i in range(len(exper_period)):
#    exper_period[i] = exper_period[i] / (9553 * nu_period[i] + 3621 - 1971 * (nu_period[i] ** 2))


# Начало прямой задачи
#hapi.getHelp(hapi.ISO_ID)
hapi.db_begin('Data')
hapi.fetch_by_ids('Mixture', [1, 2, 3, 4, 7, 8, 9, 10], nu_0, nu_end)
#hapi.select('Mixture')

molec, iso, nuij, Sref, E, Slf, Air, nair, Delta = hapi.getColumns('Mixture', ['molec_id',
            'local_iso_id', 'nu', 'Sw', 'elower', 'gamma_self', 'gamma_air', 'n_air', 'delta_air'])

T = 297  # Заданная температура
Tref = 296  # Опорная температура в HITRAN
c = 29979245800  # Скорость света
c2 = 1.4387769  # Константа ch/k
Na = 6.02214129e23  # Число Авогадро
kb = 1.380649e-16  # Постоянная Больцмана
pref = 1  # Атмосферное давление
ps = np.zeros((2, 4))  # Парциальное давление метана
p = 1  # Заданное давление
l = 100
Tr = 273
concentration = [0.008, 0.00042]
for j in range(4):
    ps[0, j] = p * concentration[0] * hapi.abundance(1, j + 1)
    ps[1, j] = p * concentration[1] * hapi.abundance(2, j + 1)

Qref = np.zeros((2, 4))
Q = np.zeros((2, 4))
M = np.zeros((2, 4))
for i in range(4):
    Qref[0, i] = hapi.partitionSum(1, i + 1, Tref)
    Qref[1, i] = hapi.partitionSum(2, i + 1, Tref)
    Q[0, i] = hapi.partitionSum(1, i + 1, T)
    Q[1, i] = hapi.partitionSum(2, i + 1, T)
    M[0, i] = hapi.molecularMass(1, i + 1)
    M[1, i] = hapi.molecularMass(2, i + 1)

khapi = np.zeros((2, 4, len(nu)))
thapi = np.zeros(len(nu))

for i in range(len(nuij)):
    AlphaD = nuij[i] / c * math.sqrt(2 * Na * kb * T * 0.693 / M[molec[i] - 1, iso[i] - 1])
    Gamma = ((Tref / T) ** nair[i]) * (Air[i] * (p - ps[molec[i] - 1, iso[i] - 1]) + Slf[i] * ps[molec[i] - 1, iso[i] - 1])
    S = Sref[i] * Qref[molec[i] - 1, iso[i] - 1] / Q[molec[i] - 1, iso[i] - 1] * \
        math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (1 - math.exp(-c2 * nuij[i] / T)) \
        / (1 - math.exp(-c2 * nuij[i] / Tref))
    khapi[molec[i] - 1, iso[i] - 1, :] += S * hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, p * Delta[i], nu)
    #khapi[molec[i] - 1, iso[i] - 1, :] += S * hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, 0, nu)

x = [khapi[0, 0, :], khapi[0, 1, :], khapi[0, 2, :], khapi[0, 3, :], khapi[1, 0, :], khapi[1, 1, :], khapi[1, 2, :], khapi[1, 3, :]]
exper_linear = - np.log(exper_period) / l
los, cov = curve_fit(lin_comb, x, exper_linear, p0=1e13 * np.ones(8), bounds=[1e8, 1e22])
print(los)

for i in range(4):
    thapi += los[i] * khapi[0, i, :] + los[i + 4] * khapi[1, i, :]
thapi = np.exp(- l * thapi)
plt.plot(nu, thapi, nu, exper_period)
plt.show()

#plt.plot(nu, (thapi - exper_reversed) / exper_reversed)
#plt.show()
