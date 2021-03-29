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

# Set-up empty dicts for data to be read into
h1, h2 = {'time': [], 'real': [], 'imag': []}, {
    'time': [], 'real': [], 'imag': []}

# Open files, write data to respective key, value pair
with open("h1.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h1[key].append(np.float64(line[key]))
with open("h2.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h2[key].append(np.float64(line[key]))

# Find y axis min/max for h1 & h2
h1_r_min, h1_i_min = min(h1['real']), min(h1['imag'])
h2_r_min, h2_i_min = min(h2['real']), min(h2['imag'])

h1_r_max, h1_i_max = max(h1['real']), max(h1['imag'])
h2_r_max, h2_i_max = max(h2['real']), max(h2['imag'])

h1_min, h1_max = min(mod(h1['real'], h1['imag'])), max(
    mod(h1['real'], h1['imag']))
h2_min, m2_max = min(mod(h2['real'], h2['imag'])), max(
    mod(h2['real'], h2['imag']))

# Set p1/p2_ymin/ymax as respective value
p1_ymin = 0
p1_ymax = 0
p2_ymin = 0
p2_ymax = 0

if h1_r_min < h1_i_min:
    p1_ymin = h1_r_min
else:
    p1_ymin = h1_i_min
if h1_r_max > h1_i_max:
    p1_ymax = h1_r_max
else:
    p1_ymax = h1_i_max
if h2_r_min < h2_i_min:
    p2_ymin = h2_r_min
else:
    p2_ymin = h2_i_min
if h2_r_max > h2_i_max:
    p2_ymax = h2_r_max
else:
    p2_ymax = h2_i_max

# Plot real & imag part of h1
plt.figure(figsize=(30, 10))
plt.plot(h1['time'], h1['real'],
         'r-', linewidth=1.2, label='real component')
plt.plot(h1['time'], h1['imag'], 'b-',
         linewidth=1.2, label='imaginary component')
# plt.plot(h1['time'], h1['imag'], 'b-',
#         linewidth=1.2, label='imaginary component')
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'h$_{1}$(t)', **axis_label_font)
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
plt.ylim(min(h1['imag']) - 0.3, max(h1['imag']) + 0.3)
plt.legend()
plt.show()

# Plot real part of h2 (has no imaginary part)
plt.figure(figsize=(10, 10))
plt.plot(h2['time'], mod(h2['real'], h2['imag']), 'r-', linewidth=1.2)
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'h$_{2}$(t)')
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
plt.xlim(0, 2*np.pi)
plt.ylim(p2_ymin, p2_ymax)
plt.show()
