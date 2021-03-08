import csv
import matplotlib.pyplot as plt
import numpy as np

axis_label_font = {'fontname': 'serif',
                   'size': 18}
axes_tick_font = {'fontname': 'serif', 'size': 12}

h1, h2 = {'time': [], 'h1.r': [], 'h1.i': []}, {
    'time': [], 'h2.r': [], 'h2.i': []}

with open("inverse.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h1[key].append(np.float64(line[key]))
with open("3_b_2.txt", 'r') as file:
    reader = csv.DictReader(file)
    for line in reader:
        for key in line.keys():
            h2[key].append(np.float64(line[key]))

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
#plt.ylim(-2, 2)
plt.xlim(0, 2*np.pi)
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
plt.ylim(0, 140)
plt.show()
