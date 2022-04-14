import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

nu_step = 0.05
beginning = 36
end = 960
fabri = [float(i) for i in open("Experiments/f1").read().split()]
calibr = np.zeros(end - beginning)
a = beginning
fabrimax = []
pointmax = []
for i in range(beginning, end):
    if ((fabri[i] >= max(fabri[i - 2], fabri[i - 1], fabri[i + 1], fabri[i + 2])) or
        (fabri[i] <= min(fabri[i - 2], fabri[i - 1], fabri[i + 1], fabri[i + 2]))) and i > a + 2:
        l = i - a
        fabrimax.append(fabri[i])
        pointmax.append(i)
        for j in range(a - beginning, i - beginning):
            calibr[j] = nu_step / l / 2
        a = i
for i in range(a, end):
    calibr[i - beginning] = calibr[a - beginning - 1]

"""
nu_max = [calibr[end - beginning - 1] * (end - a)]
for i in range(len(pointmax)):
    nu_max.append(nu_max[i - 1] + 0.05)

interp_function = interp1d(np.arange(end - beginning), nu, kind=2)
nu_period = interp_function(np.arange(0, 2.001, 0.001))
"""

nu1 = [3563 + i * 0.001 for i in range(2001)]
nu = np.zeros(end - beginning)
for i in range(1, end - beginning):
    nu[i] = nu[i - 1] + calibr[end - beginning - i]

hitran = pd.read_excel("Transmittance spectrums/Reverse.xlsx")
Trans = np.array(hitran.Trans)
exper1 = [float(i) for i in open("Experiments/1").read().split()]
exper1_cut = exper1[beginning:end]
exper1_rev = []
for i in range(end - beginning):
    exper1_rev.append((exper1_cut[end - beginning - 1 - i] - exper1[end + 20]) / (19017 - nu[i] * 9474))

nu_period = np.arange(0, nu[end - beginning - 1], 0.001)
interp_exper = interp1d(nu, exper1_rev, kind=2)
exper_period = interp_exper(nu_period)
plt.plot(nu_period, exper_period, nu_period, Trans[200:len(nu_period) + 200])

#t = np.arange(len(fabri))
#plt.plot(t, fabri, pointmax, fabrimax)
#t = np.arange(end - beginning)
#plt.plot(t, calibr)
#plt.plot(t, nu)
#t = np.arange(len(exper1_rev))
#plt.plot(nu, exper1_rev)
plt.show()
