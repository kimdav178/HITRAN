import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

def baseline(a, b):
    return a * t + b

beginning = 36
end = 960
t = np.arange(end - beginning)
exper1 = [float(i) for i in open("Experiments/1").read().split()]
exper1_cut = exper1[beginning:end]
exper1_rev = []
for i in range(end - beginning):
    exper1_rev.append(exper1_cut[end - beginning - 1 - i] - exper1[end + 20])


fig = plt.figure()
ax = fig.add_subplot(111)

# Adjust the subplots region to leave some space for the sliders and buttons
fig.subplots_adjust(left=0.25, bottom=0.25)

a_0 = -10
b_0 = 20000

# Draw the initial plot
# The 'line' variable is used for modifying the line later
[line, approximation] = ax.plot(t, baseline(a_0, b_0), t, exper1_rev, linewidth=2)
#ax.set_xlim([0, 1])
#ax.set_ylim([-10, 10])

# Add two sliders for tweaking the parameters

# Define an axes area and draw a slider in it
a_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03])
a_slider = Slider(a_slider_ax, 'a', -20, 0, valinit=a_0)

# Draw another slider
b_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03])
b_slider = Slider(b_slider_ax, 'b', 18000, 25000, valinit=b_0)

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    line.set_ydata(baseline(a_slider.val, b_slider.val))
    fig.canvas.draw_idle()
a_slider.on_changed(sliders_on_changed)
b_slider.on_changed(sliders_on_changed)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    b_slider.reset()
    a_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

# Add a set of radio buttons for changing color
color_radios_ax = fig.add_axes([0.025, 0.5, 0.15, 0.15])
color_radios = RadioButtons(color_radios_ax, ('red', 'blue', 'green'), active=0)
def color_radios_on_clicked(label):
    line.set_color(label)
    fig.canvas.draw_idle()
color_radios.on_clicked(color_radios_on_clicked)

plt.show()
