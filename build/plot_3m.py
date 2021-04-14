import csv
import matplotlib.pyplot as plt
import numpy as np


def mod(x, y):
    result = []
    for _x, _y in zip(x, y):
        result.append((_x**2 + _y**2) ** 0.5)
    return result


# Specifying formatting for axis
axis_label_font = {'fontname': 'serif','size': 18}
axes_tick_font = {'fontname': 'serif', 'fontsize': 12}
legend_font = {'fontsize': 16}

# Set-up empty dicts for data to be read into
h3, i_h3 = {'n': [], 'time': [], 'real': [], 'imag': []}, {
    'time': [], 'real': [], 'imag': []}

# Open files, write data to respective key, value pair
with open("h3.txt", 'r') as file:
    reader = csv.reader(file)
    for line in reader:
        h3['n'].append(np.float64(line[0]))
        h3['time'].append(np.float64(line[1]))
        h3['real'].append(np.float64(line[2]))
        h3['imag'].append(np.float64(line[3]))
with open("inv_3.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            i_h3[key].append(np.float64(line[key]))


# Plot real & imag part of h3 & i_h3
plt.figure(figsize=(10, 10))
plt.plot(h3['time'], h3['real'], 'r-', linewidth=1.2,
         label=r'h$_3$')
plt.plot(i_h3['time'], i_h3['real'], 'b-', linewidth=1.2,
         label=r'h$_3$`')
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'h$_{3}$(t)', **axis_label_font)
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
plt.xlim(0, 2*np.pi)
plt.ylim(-12.5, 12.5)
plt.legend(**legend_font)
plt.show()

plt.figure(figsize=(10, 10))
plt.plot(h3['time'], h3['imag'], 'r-',
         linewidth=1.2, label=r'h$_3$')
plt.plot(i_h3['time'], i_h3['imag'], 'b-',
         linewidth=1.2, label=r'h$_3$`')
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'h$_{3}$(t)', **axis_label_font)
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
plt.xlim(0, 2*np.pi)
plt.ylim(-12.5, 12.5)
plt.legend(**legend_font)
plt.show()
