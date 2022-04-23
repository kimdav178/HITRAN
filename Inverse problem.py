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

nu_step = 0.05
beginning = 27
end = 960
a = beginning
pointmax = [beginning]
valuemax = [0]
fabri = [float(i) for i in open("New exp/f1").read().split()]
exper1 = [float(i) for i in open("New exp/1").read().split()]
hitran = pd.read_excel("Transmittance spectrums/Reverse 3568 3571.xlsx")
Trans = np.array(hitran.Trans)

"""
plt.plot(np.arange(len(fabri)), fabri)
plt.show()
"""

for i in range(beginning, end):
    if ((fabri[i] >= max(fabri[i - 2], fabri[i - 1], fabri[i + 1], fabri[i + 2])) or
        (fabri[i] <= min(fabri[i - 2], fabri[i - 1], fabri[i + 1], fabri[i + 2]))) and i > a + 2:
        l = i - a
        pointmax.append(i)
        valuemax.append(valuemax[len(valuemax) - 1] + nu_step / 2)
        a = i

for i in range(beginning, a):
    exper1[i] -= exper1[end + 20]
exper1_cut = exper1[beginning:a]
interp_grid = np.arange(beginning, pointmax[len(pointmax) - 1])
interp_nu = interp1d(pointmax, valuemax, kind=2)
nu_interp = interp_nu(interp_grid)

"""
plt.plot(np.arange(len(nu_interp)), nu_interp)
plt.show()
"""

nu_period = np.arange(0, nu_interp[len(nu_interp) - 1], 0.001)
interp_exper = interp1d(nu_interp, exper1_cut, kind=2)
exper_period = interp_exper(nu_period)


fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.25)

a_0 = 10000
b_0 = 5000
c_0 = 0
[line, approximation] = ax.plot(nu_period, baseline(a_0, b_0, c_0), nu_period, exper_period, linewidth=2)

a_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03])
a_slider = Slider(a_slider_ax, 'a', 5000, 20000, valinit=a_0)
b_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03])
b_slider = Slider(b_slider_ax, 'b', 0, 10000, valinit=b_0)
c_slider_ax = fig.add_axes([0.25, 0.05, 0.65, 0.03])
c_slider = Slider(c_slider_ax, 'c', -2000, 2000, valinit=c_0)

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    line.set_ydata(baseline(a_slider.val, b_slider.val, c_slider.val))
    fig.canvas.draw_idle()
a_slider.on_changed(sliders_on_changed)
b_slider.on_changed(sliders_on_changed)
c_slider.on_changed(sliders_on_changed)

#plt.plot(nu_period, exper_period, nu_period, Trans[200:len(nu_period) + 200])

plt.show()

for i in range(len(exper_period)):
    exper_period[i] = exper_period[i] / (9553 * nu_period[i] + 3621 - 1971 * (nu_period[i] ** 2))

"""
print(a_slider.val, b_slider.val, c_slider.val)
for i in range(len(exper_period)):
    exper_period[i] = exper_period[i] / (a_slider.val * nu_period[i] + b_slider.val + c_slider.val * (nu_period[i] ** 2))
"""

nu_0 = 3568.945
nu = [nu_0 + i * 0.001 for i in range(len(exper_period))]
nu_end = nu[len(nu) - 1]
exper_reversed = [exper_period[len(exper_period) - i - 1] for i in range(len(exper_period))]


#hapi.getHelp(hapi.ISO_ID)
hapi.db_begin('Data')
#hapi.fetch_by_ids('Mixture', [1, 2, 3, 4, 7, 8, 9, 10], nu_0, nu_end)
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
#los = 2.6867811e19
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
exper_linear = - np.log(exper_reversed) / l
los, cov = curve_fit(lin_comb, x, exper_linear, p0=1e13 * np.ones(8), bounds=[1e10, 1e22])
print(los)

for i in range(4):
    thapi += los[i] * khapi[0, i, :] + los[i + 4] * khapi[1, i, :]
thapi = np.exp(- l * thapi)
plt.plot(nu, thapi, nu, exper_reversed)
plt.show()

#plt.plot(nu[945:945 + len(exper_reversed)], exper_reversed, nu, Trans)
plt.plot(nu, (thapi - exper_reversed) / exper_reversed)
plt.show()
