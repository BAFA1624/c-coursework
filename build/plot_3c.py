import csv
import matplotlib.pyplot as plt
import numpy as np


def mod(x, y):
    result = []
    for _x, _y in zip(x, y):
        result.append((_x**2 + _y**2) ** 0.5)
    return result


# Specifying formatting for axis
axis_label_font = {'fontname': 'serif',
                   'size': 18}
axes_tick_font = {'fontname': 'serif', 'size': 12}
legend_font = {'fontsize': 16}

h1_prime, h1 = {'time': [], 'real': [], 'imag': []}, {
    'time': [], 'real': [], 'imag': []}
# Open files, write data to respective key, value pair
with open("inv_1.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h1_prime[key].append(np.float64(line[key]))
with open("h1.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h1[key].append(np.float64(line[key]))


# Plot real & imag part of h1
plt.figure(figsize=(10, 10))
plt.plot(h1['time'], h1['real'],
         'r-', linewidth=1.2, label=r'h$_1$')
plt.plot(h1_prime['time'], h1_prime['real'], 'b-',
         linewidth=1.2, label=r'h$_1$`')
# plt.plot(h1['time'], h1['imag'], 'b-',
#         linewidth=1.2, label='imaginary component')
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'h$_{1}$(t)', **axis_label_font)
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
plt.xlim(0, 2*np.pi)
plt.ylim(-2, 2)
plt.legend(**legend_font)
plt.show()

plt.figure(figsize=(10, 10))
plt.plot(h1['time'], h1['imag'],
         'r-', linewidth=1.2, label=r'h$_1$')
plt.plot(h1_prime['time'], h1_prime['imag'], 'b-',
         linewidth=1.2, label=r'h$_1$`')
# plt.plot(h1['time'], h1['imag'], 'b-',
#         linewidth=1.2, label='imaginary component')
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'h$_{1}$(t)', **axis_label_font)
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
plt.xlim(0, 2*np.pi)
plt.ylim(-2, 2)
plt.legend(**legend_font)
plt.show()

