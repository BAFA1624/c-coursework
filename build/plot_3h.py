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
axes_tick_font = {'fontname': 'serif', 'size': 12}
legend_font = {'fontsize': 16}

# Open files, write data to respective key, value pair
h2, h2_prime = {'time': [], 'real': [], 'imag': []}, {
    'time': [], 'real': [], 'imag': []}
# Open files, write data to respective key, value pair
with open("inv_2.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h2_prime[key].append(np.float64(line[key]))
with open("h2.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h2[key].append(np.float64(line[key]))

#
# Plot h2
plt.figure(figsize=(10, 10))
plt.plot(h2['time'], h2['real'], 'r-', 
        linewidth=1.2, label=r'h$_2$')
plt.plot(h2_prime['time'], h2_prime['real'], 'b-', 
        linewidth=1.2, label=r'h$_2$`')
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'h$_{2}$(t)', **axis_label_font)
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
plt.xlim(0, 2*np.pi)
plt.ylim(-20, 140)
plt.legend(**legend_font)
plt.show()

plt.figure(figsize=(10, 10))
plt.plot(h2['time'], h2['imag'], 'r-', 
        linewidth=1.2, label=r'h$_2$')
plt.plot(h2_prime['time'], h2_prime['imag'], 'b-', 
        linewidth=1.2, label=r'h$_2$`')
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'h$_{2}$(t)', **axis_label_font)
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
plt.xlim(0, 2*np.pi)
plt.legend(**legend_font)
plt.show()
