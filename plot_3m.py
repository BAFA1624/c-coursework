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

# Find y axis min/max for h1 & h2
h3_min, h3_max = min(mod(h3['real'], h3['imag'])), max(
    mod(h3['real'], h3['imag']))
ih3_min, ih3_max = min(mod(i_h3['real'], i_h3['imag'])), max(
    mod(i_h3['real'], i_h3['imag']))

h3_r_max, h3_i_max = max(h3['real']), max(h3['imag'])
ih3_r_max, ih3_i_max = max(i_h3['real']), max(i_h3['imag'])

# Set p1/p2_ymin/ymax as respective value
p1_ymin = 0
p1_ymax = 0

if h3_min < ih3_min:
    p1_ymin = h3_min
else:
    p1_ymin = ih3_min
if h3_max > ih3_max:
    p1_ymax = h3_max
else:
    p1_ymax = ih3_max


# Plot real & imag part of h3 & i_h3
plt.figure(figsize=(10, 10))
plt.plot(h3['time'], mod(h3['real'], h3['imag']), 'b-', linewidth=1.2,
         label='Original Signal')
plt.plot(i_h3['time'], mod(i_h3['real'], i_h3['imag']), 'r--', linewidth=1.2,
         label='Transformed Version')
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'|h$_{3}$(t)|', **axis_label_font)
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
plt.xlim(0, 2*np.pi)
plt.ylim(p1_ymin, p1_ymax)
plt.legend()
plt.show()


# plt.figure(figsize=(10, 10))
# plt.plot(i_h3['time'], i_h3['real'], 'b-',
#         linewidth=1.2, label='Transformed real component')
# plt.plot(i_h3['time'], i_h3['imag'], 'r--',
#         linewidth=1.2, label='Transformed imaginary component')
# plt.xlabel('Time / s', **axis_label_font)
# plt.ylabel(r'h$_{3}$`(t)', **axis_label_font)
# plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
#           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
#           **axes_tick_font)
# plt.yticks(**axes_tick_font)
# plt.tick_params(direction='in')
# plt.xlim(0, 2*np.pi)
# plt.ylim(p1_ymin, p1_ymax)
# plt.legend()
# plt.show()
