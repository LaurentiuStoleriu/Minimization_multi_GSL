

import numpy as np
import matplotlib.pyplot as plt

l0 = 0.4
r_mic = 0.25
R_mare = 0.45
A = 0.001
period = 5.0 * np.pi

ini_poz = []
ini_r = []
fin_poz = []
fin_r = []

with open('minim_10.dat', 'r') as fis:
    for line in fis:
        y = line.split()
        ini_poz.append(float(y[1]))
        ini_r.append(float(y[2]))

ini_poz = np.asarray(ini_poz)
ini_r = np.asarray(ini_r)
some_y_ini = np.ones(len(ini_poz)) * (-0.0015)
energ_ini = -A*np.sin(ini_poz * period)


with open('minim_10_sw.dat', 'r') as fis:
    for line in fis:
        y = line.split()
        fin_poz.append(float(y[1]))
        fin_r.append(float(y[2]))

fin_poz = np.asarray(fin_poz)
fin_r = np.asarray(fin_r)
some_y_fin = np.ones(len(fin_poz)) * (+0.0015)
energ_fin = -A*np.sin(fin_poz * period)
        
x_coords = np.linspace(-1,6, 1000)
energ_profile = -A*np.sin(x_coords * period)



fig, ax = plt.subplots()

plt.style.use('ggplot')

ax.plot(x_coords, energ_profile, linewidth=0.5)

ax.scatter(ini_poz, energ_ini, s = 100, color='blue')
ax.scatter(ini_poz, some_y_ini, s = 100, color='blue')
ax.plot(ini_poz, some_y_ini, color='blue')

ax.scatter(fin_poz, energ_fin, s=100, color='red')
ax.scatter(fin_poz, some_y_fin, s = 100, color='red')
ax.plot(fin_poz, some_y_fin, color='red')

plt.show()