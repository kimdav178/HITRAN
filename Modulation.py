import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from matplotlib.widgets import Slider
import hapi
from scipy.optimize import curve_fit
import math


def lin_comb_iso4(x, y1, y2, y3, y4, y5, y6, y7, y8):
    return x[0] * y1 + x[1] * y2 + x[2] * y3 + x[3] * y4 + x[4] * y5 + x[5] * y6 + x[6] * y7 + x[7] * y8


def lin_comb_iso1(x, y1, y2):
    return x[0] * y1 + x[1] * y2


def exponential_comb_iso1(x, y1, y2):
    return np.exp((- x[0] * y1 - x[1] * y2) * l)


def baseline(a, b, c):
    return a * nu_period1 + b + c * (nu_period1 ** 2)


def normalized(a, b, c):
    return exper_period1 / (a * nu_period1 + b + c * (nu_period1 ** 2))


def nu_adjust(nu_0):
    return np.arange(nu_0, nu_0 - 0.0005 + len(exper_period1) * 0.001, 0.001)


def search(x, a, b, c):
    return a * x * x + b * x + c


nu_step = 0.048
beginning = 0  # Номер первого пика
beginning1 = 48
beginning2 = 0
beginning3 = 0
beginning4 = 0
end = 960
a = beginning1
pointmax = [beginning1]
valuemax = [0]
fabri = [float(i) for i in open("Modulation/fpmod").read().split()]
exper = [float(i) for i in open("Modulation/1mod").read().split()]
hitran = pd.read_excel("3562-3564/3562.xlsx")
Trans = np.array(hitran.Trans)
nu_model = np.array(hitran.nu)

max_fabri = max(fabri)
max_exp = max(exper)
for i in range(len(fabri)):
    fabri[i] = max_fabri - fabri[i]
    exper[i] = max_exp - exper[i]

#plt.plot(np.arange(len(fabri)), fabri)
#plt.show()
#plt.plot(np.arange(len(exper)), exper)
#plt.show()

exper_cut1 = exper[beginning:end]
fabri_cut = fabri[beginning:end]
points1 = []
points2 = []
points3 = []
points4 = []
exper1 = []
exper2 = []
exper3 = []
exper4 = []
fabri1 = []
fabri2 = []
fabri3 = []
fabri4 = []
i = 0
while i < len(exper_cut1) - 3:
    exper1.append(exper_cut1[i])
    fabri1.append(fabri_cut[i])
    points1.append(i)
    exper2.append(exper_cut1[i + 1])
    fabri2.append(fabri_cut[i + 1])
    points2.append(i + 1)
    exper3.append(exper_cut1[i + 2])
    fabri3.append(fabri_cut[i + 2])
    points3.append(i + 2)
    exper4.append(exper_cut1[i + 3])
    fabri4.append(fabri_cut[i + 3])
    points4.append(i + 3)
    i += 4

interp_grid1 = np.arange(0, end - 4)
interp_signal1 = interp1d(points1, exper1, kind='quadratic')
signal_interp1 = interp_signal1(interp_grid1)
interp_fabri1 = interp1d(points1, fabri1, kind='quadratic')
fabri_interp1 = interp_fabri1(interp_grid1)
interp_grid2 = np.arange(1, end - 3)
interp_signal2 = interp1d(points2, exper2, kind='quadratic')
signal_interp2 = interp_signal2(interp_grid2)
interp_fabri2 = interp1d(points2, fabri2, kind='quadratic')
fabri_interp2 = interp_fabri2(interp_grid2)
interp_grid3 = np.arange(2, end - 2)
interp_signal3 = interp1d(points3, exper3, kind='quadratic')
signal_interp3 = interp_signal3(interp_grid3)
interp_fabri3 = interp1d(points3, fabri3, kind='quadratic')
fabri_interp3 = interp_fabri3(interp_grid3)
interp_grid4 = np.arange(3, end - 1)
interp_signal4 = interp1d(points4, exper4, kind='quadratic')
signal_interp4 = interp_signal4(interp_grid4)
interp_fabri4 = interp1d(points4, fabri4, kind='quadratic')
fabri_interp4 = interp_fabri4(interp_grid4)

plt.plot(interp_grid1, fabri_interp1)
plt.show()
#plt.plot(interp_grid2, fabri_interp2)
#plt.show()
#plt.plot(interp_grid3, fabri_interp3)
#plt.show()
#plt.plot(interp_grid4, fabri_interp4)
#plt.show()

for i in range(beginning1, end - 6):
    if ((fabri_interp1[i] >= max(fabri_interp1[i - 2], fabri_interp1[i - 1], fabri_interp1[i + 1], fabri_interp1[i + 2])) or
        (fabri_interp1[i] <= min(fabri_interp1[i - 3], fabri_interp1[i - 2], fabri_interp1[i - 1], fabri_interp1[i + 1], fabri_interp1[i + 2]))) and i > a + 2:
        l = i - a
        pointmax.append(i)
        valuemax.append(valuemax[len(valuemax) - 1] + nu_step / 2)
        a = i
#pointmax.append(end - 4)
#valuemax.append(valuemax[len(valuemax) - 1] + nu_step / 2)

exper_cut1 = signal_interp1[beginning1:a]
interp_grid_nu1 = np.arange(beginning1, pointmax[len(pointmax) - 1])
interp_nu1 = interp1d(pointmax, valuemax, kind='quadratic')
nu_interp1 = interp_nu1(interp_grid_nu1)

#plt.plot(np.arange(len(nu_interp1)), nu_interp1)
#plt.show()

# Разворот и интерполяция спектра на периодическую сетку частот
nu_reversed1 = [nu_interp1[len(nu_interp1) - i - 1] for i in range(len(nu_interp1))]
exper_reversed1 = [exper_cut1[len(exper_cut1) - i - 1] for i in range(len(exper_cut1))]
nu_period1 = np.arange(nu_reversed1[0], 0, -0.001)
interp_exper1 = interp1d(nu_reversed1, exper_reversed1, kind='quadratic')
exper_period1 = interp_exper1(nu_period1)


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
a_0 = 5028.72
b_0 = 1301.86
c_0 = -984.24
nu_0 = 3561.36
[line, approximation] = ax.plot(nu_period1, baseline(a_0, b_0, c_0), nu_period1, exper_period1, linewidth=2)
[exp, model] = bx.plot(nu_adjust(nu_0), normalized(a_0, b_0, c_0), nu_model, Trans)
bx.set_ylim(0, 1)
a_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03])
a_slider = Slider(a_slider_ax, 'a', 0, 10000, valinit=a_0)
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
for i in range(len(exper_period1)):
    exper_period1[i] = exper_period1[i] / (a_slider.val * nu_period1[i] + b_slider.val + c_slider.val * (nu_period1[i] ** 2))
print(nu_slider.val)
nu = nu_adjust(nu_slider.val)
nu_end = nu[len(nu) - 1]


# Начало прямой задачи
#hapi.getHelp(hapi.ISO_ID)
hapi.db_begin('Data')
#hapi.fetch_by_ids('Mixture', [1, 2, 3, 4, 7, 8, 9, 10], nu_0 - 10, nu_end + 10)
#hapi.fetch_by_ids('Mixture', [1, 7], nu_0 - 10, nu_end + 10)
#hapi.select('Mixture')

molec, iso, nuij, Sref, E, Slf, Air, nair, Delta = hapi.getColumns('Mixture', ['molec_id',
            'local_iso_id', 'nu', 'Sw', 'elower', 'gamma_self', 'gamma_air', 'n_air', 'delta_air'])

molec_num = 2
iso_num = 1
T = 297  # Заданная температура
Tref = 296  # Опорная температура в HITRAN
c = 29979245800  # Скорость света
c2 = 1.4387769  # Константа ch/k
Na = 6.02214129e23  # Число Авогадро
kb = 1.380649e-16  # Постоянная Больцмана
pref = 1  # Атмосферное давление
ps = np.zeros((molec_num, iso_num))  # Парциальное давление метана
p = 1  # Заданное давление
l = 130
Tr = 273
concentration = [0.008, 0.00042]
for j in range(iso_num):
    ps[0, j] = p * concentration[0] * hapi.abundance(1, j + 1)
    ps[1, j] = p * concentration[1] * hapi.abundance(2, j + 1)

Qref = np.zeros((molec_num, iso_num))
Q = np.zeros((molec_num, iso_num))
M = np.zeros((molec_num, iso_num))
for i in range(iso_num):
    Qref[0, i] = hapi.partitionSum(1, i + 1, Tref)
    Qref[1, i] = hapi.partitionSum(2, i + 1, Tref)
    Q[0, i] = hapi.partitionSum(1, i + 1, T)
    Q[1, i] = hapi.partitionSum(2, i + 1, T)
    M[0, i] = hapi.molecularMass(1, i + 1)
    M[1, i] = hapi.molecularMass(2, i + 1)

khapi = np.zeros((molec_num, iso_num, len(nu)))
thapi = np.zeros(len(nu))

for i in range(len(nuij)):
    AlphaD = nuij[i] / c * math.sqrt(2 * Na * kb * T * 0.693 / M[molec[i] - 1, iso[i] - 1])
    Gamma = ((Tref / T) ** nair[i]) * (Air[i] * (p - ps[molec[i] - 1, iso[i] - 1]) + Slf[i] * ps[molec[i] - 1, iso[i] - 1])
    S = Sref[i] * Qref[molec[i] - 1, iso[i] - 1] / Q[molec[i] - 1, iso[i] - 1] * \
        math.exp(-c2 * E[i] / T) / math.exp(-c2 * E[i] / Tref) * (1 - math.exp(-c2 * nuij[i] / T)) \
        / (1 - math.exp(-c2 * nuij[i] / Tref))
    khapi[molec[i] - 1, iso[i] - 1, :] += S * hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, p * Delta[i], nu)
    #khapi[molec[i] - 1, iso[i] - 1, :] += S * hapi.PROFILE_VOIGT(nuij[i], AlphaD, Gamma, 0, nu)

"""
x = [khapi[0, 0, :], khapi[0, 1, :], khapi[0, 2, :], khapi[0, 3, :], khapi[1, 0, :], khapi[1, 1, :], khapi[1, 2, :], khapi[1, 3, :]]
exper_linear = - np.log(exper_period1) / l
los, cov = curve_fit(lin_comb_iso4, x, exper_linear, p0=1e13 * np.ones(8), bounds=[1e8, 1e22])

x = [khapi[0, 0, :], khapi[1, 0, :]]
#exper_linear = - np.log(exper_period1) / l
#los, cov = curve_fit(lin_comb_iso1, x, exper_linear, p0=1e13 * np.ones(2), bounds=[1e8, 1e22])
los, cov = curve_fit(exponential_comb_iso1, x, exper_period1, p0=1e13 * np.ones(2), bounds=[1e8, 1e22])
print(los)

los[0] = 1.909e17
los[1] = 1.02e16
for i in range(iso_num):
    thapi += los[i] * khapi[0, i, :] + los[i + iso_num] * khapi[1, i, :]
thapi = np.exp(- l * thapi)
plt.plot(nu, thapi, nu_model, Trans)
plt.show()

plt.plot(nu, (thapi - exper_period1) / exper_period1)
plt.show()

plt.plot(np.arange(len(exper_period1)), exper_period1 / thapi)
plt.show()

koefs, covar = curve_fit(search, nu_period1, exper_period1 / thapi, p0=5000 * np.ones(3), bounds=[-1e4, 1e4])
print(koefs)

new_baseline = search(nu_period1, koefs[0], koefs[1], koefs[2])
plt.plot(np.arange(len(exper_period1)), exper_period1 / thapi, np.arange(len(exper_period1)), new_baseline)
plt.show()

for i in range(len(exper_period1)):
    exper_period1[i] = exper_period1[i] / (koefs[1] * nu_period1[i] + koefs[2] + koefs[0] * (nu_period1[i] ** 2))
"""

thapi = np.zeros(len(nu))
cut_len = int(len(nu) - 180)
x = [khapi[0, 0, :cut_len], khapi[1, 0, :cut_len]]
los, cov = curve_fit(exponential_comb_iso1, x, exper_period1[:cut_len], p0=[2e17, 1e16], bounds=[1e15, 1e18])
print(los)

for i in range(iso_num):
    thapi[:cut_len] += los[i] * khapi[0, i, :cut_len] + los[i + iso_num] * khapi[1, i, :cut_len]
thapi = np.exp(- l * thapi[:cut_len])
plt.plot(nu[:cut_len], thapi[:cut_len], nu[:cut_len], exper_period1[:cut_len])
plt.show()

plt.plot(nu[:cut_len], (thapi[:cut_len] - exper_period1[:cut_len]) / exper_period1[:cut_len])
plt.show()