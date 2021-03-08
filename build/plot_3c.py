import csv
import matplotlib.pyplot as plt
import numpy as np

axis_label_font = {'fontname': 'serif',
                   'size': 18}
axes_tick_font = {'fontname': 'serif', 'size': 12}

h1, h2 = {'time': [], 'h1.r': [], 'h1.i': []}, {
    'time': [], 'h2.r': [], 'h2.i': []}

with open("inv_1.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h1[key].append(np.float64(line[key]))
with open("inv_2.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h2[key].append(np.float64(line[key]))

h1_r_min, h1_i_min = min(h1['h1.r']), min(h1['h1.i'])
h2_r_min, h2_i_min = min(h2['h2.r']), min(h2['h2.i'])

h1_r_max, h1_i_max = max(h1['h1.r']), max(h1['h1.i'])
h2_r_max, h2_i_max = max(h2['h2.r']), max(h2['h2.i'])

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

plt.figure(figsize=(10, 10))
plt.plot(h1['time'], h1['h1.r'], 'r-', linewidth=1.2, label='real component')
plt.plot(h1['time'], h1['h1.i'], 'b-',
         linewidth=1.2, label='imaginary component')
plt.xlabel('Time / s', **axis_label_font)
plt.ylabel(r'h$_{1}$(t)', **axis_label_font)
plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
           ['0', r'$\pi_2$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'],
           **axes_tick_font)
plt.yticks(**axes_tick_font)
plt.tick_params(direction='in')
# plt.ylim(-2, 2)
plt.xlim(0, 2*np.pi)
plt.ylim(p1_ymin, p1_ymax)
plt.legend()
plt.show()

plt.figure(figsize=(10, 10))
plt.plot(h2['time'], h2['h2.r'], 'r-', linewidth=1.2)
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
