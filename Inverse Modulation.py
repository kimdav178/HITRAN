import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from matplotlib.widgets import Slider


def baseline(a, b, c):
    return a * nu_period + b + c * (nu_period ** 2)


def normalized(a, b, c):
    return exper_cut1 / (a * nu_period + b + c * (nu_period ** 2))


def lin_comb_iso1(x, y1, y2):
    return x[0] * y1 + x[1] * y2


def exponential_comb_iso1(x, y1, y2):
    return np.exp((- x[0] * y1 - x[1] * y2) * l)


def nu_adjust(nu_0):
    return np.arange(nu_0, nu_0 - 0.0005 + len(exper_cut1) * 0.001, 0.001)


nu_step = 0.0485
beginning = 0  # Номер первого пика
beginning1 = 48  # Номер первого пика для спектра первых гармоник
beginning2 = 0
beginning3 = 0
beginning4 = 0
end = 960
a = beginning1
pointmax = [beginning1]
valuemax = [0]
fabri = [float(i) for i in open("Modulation/fpmod").read().split()]
exper = [float(i) for i in open("Modulation/1mod").read().split()]
hitran = pd.read_excel("3562-3564/Smooth.xlsx")
Trans2 = np.array(hitran.Trans)
Trans = np.gradient(np.gradient(Trans2, 0.01), 0.01)
nu_model = np.array(hitran.nu)
plt.plot(nu_model, Trans)
plt.show()

max_fabri = max(fabri)
max_exp = max(exper)
for i in range(len(fabri)):
    fabri[i] = max_fabri - fabri[i]
    exper[i] = max_exp - exper[i]

plt.plot(np.arange(len(exper)), exper)
plt.show()

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
deriv = []
i = 0
#while i < end - 111:
while i < end - 3:
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
    deriv.append(2 * (exper_cut1[i + 3] + exper_cut1[i + 1] - 2 * exper_cut1[i]) / (exper_cut1[i] + exper_cut1[i + 2]))
    i += 4

interp_grid1 = np.arange(0, end - 4)
#interp_signal1 = interp1d(points1, exper1, kind='quadratic')
#signal_interp1 = interp_signal1(interp_grid1)
interp_fabri1 = interp1d(points1, fabri1, kind='quadratic')
fabri_interp1 = interp_fabri1(interp_grid1)
interp_deriv = interp1d(points1, deriv, kind='quadratic')
deriv_interp = interp_deriv(interp_grid1)

for i in range(beginning1, end - 6):
    if ((fabri_interp1[i] >= max(fabri_interp1[i - 2], fabri_interp1[i - 1], fabri_interp1[i + 1], fabri_interp1[i + 2])) or
        (fabri_interp1[i] <= min(fabri_interp1[i - 3], fabri_interp1[i - 2], fabri_interp1[i - 1], fabri_interp1[i + 1], fabri_interp1[i + 2]))) and i > a + 2:
        l = i - a
        pointmax.append(i)
        valuemax.append(valuemax[len(valuemax) - 1] + nu_step / 2)
        a = i

#exper_cut1 = signal_interp1[beginning1:a]
interp_grid_nu1 = np.arange(beginning1, pointmax[len(pointmax) - 1])
interp_nu1 = interp1d(pointmax, valuemax, kind='quadratic')
nu_interp1 = interp_nu1(interp_grid_nu1)
deriv_cut = deriv_interp[beginning1:a]


# Разворот и интерполяция спектра на периодическую сетку частот
nu_reversed1 = [nu_interp1[len(nu_interp1) - i - 1] for i in range(len(nu_interp1))]
exper_reversed1 = [exper_cut1[len(exper_cut1) - i - 1] for i in range(len(exper_cut1))]
nu_period1 = np.arange(nu_reversed1[0], 0, -0.001)
interp_exper1 = interp1d(nu_reversed1, exper_reversed1, kind='quadratic')
exper_period1 = interp_exper1(nu_period1)
deriv_reversed = [deriv_cut[len(deriv_cut) - i - 1] for i in range(len(deriv_cut))]
interp_deriv = interp1d(nu_reversed1, deriv_reversed, kind='quadratic')
deriv_period = interp_deriv(nu_period1)
#deriv_period = [deriv_period2[len(deriv_period2) - i - 1] for i in range(len(deriv_period2))]


# Подбор частоты


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
[line, approximation] = ax.plot(nu_period, baseline(a_0, b_0, c_0), nu_period, exper_cut1, linewidth=2)
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


def sliders_on_changed(val):
    exp.set_xdata(nu_adjust(nu_slider.val))
    exp.set_ydata(deriv_period)
    fig.canvas.draw_idle()


fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.15)
nu_0 = 3561.36
[exp, model] = ax.plot(nu_adjust(nu_0), deriv_period, nu_model, Trans)
#ax.set_ylim(-1, 1)
nu_slider_ax = fig.add_axes([0.15, 0.05, 0.65, 0.03])
nu_slider = Slider(nu_slider_ax, 'nu', nu_0 - 1, nu_0 + 1, valinit=nu_0)
nu_slider.on_changed(sliders_on_changed)
plt.show()
print(nu_slider.val)
nu = nu_adjust(nu_slider.val)
nu_end = nu[len(nu) - 1]

