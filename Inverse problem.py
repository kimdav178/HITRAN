import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from matplotlib.widgets import Slider


def baseline(a, b, c):
    return a * nu_period + b + c * (nu_period ** 2)

nu_step = 0.05
beginning = 36
end = 960
a = beginning
pointmax = [beginning]
valuemax = [0]
nu1 = [3563 + i * 0.001 for i in range(2001)]
fabri = [float(i) for i in open("Experiments/f1").read().split()]
exper1 = [float(i) for i in open("Experiments/1").read().split()]
hitran = pd.read_excel("Transmittance spectrums/Reverse.xlsx")
Trans = np.array(hitran.Trans)
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

nu_period = np.arange(0, nu_interp[len(nu_interp) - 1], 0.001)
interp_exper = interp1d(nu_interp, exper1_cut, kind=2)
exper_period = interp_exper(nu_period)


fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.25)

a_0 = 10000
b_0 = 5000
c_0 = 1000
[line, approximation] = ax.plot(nu_period, baseline(a_0, b_0, c_0), nu_period, exper_period, linewidth=2)

a_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03])
a_slider = Slider(a_slider_ax, 'a', 5000, 20000, valinit=a_0)
b_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03])
b_slider = Slider(b_slider_ax, 'b', 0, 10000, valinit=b_0)
c_slider_ax = fig.add_axes([0.25, 0.05, 0.65, 0.03])
c_slider = Slider(c_slider_ax, 'c', 0, 2000, valinit=c_0)

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    line.set_ydata(baseline(a_slider.val, b_slider.val, c_slider.val))
    fig.canvas.draw_idle()
a_slider.on_changed(sliders_on_changed)
b_slider.on_changed(sliders_on_changed)
c_slider.on_changed(sliders_on_changed)

#plt.plot(nu_period, exper_period, nu_period, Trans[200:len(nu_period) + 200])

plt.show()
