import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

r1 = 14.7  # medido
nesp = 11
df = pd.read_csv('ensaio_livre_travado.csv', comment='#')

travado = df.loc[0]
livre = df.loc[1:]
# print(travado)
# print(livre)

# ensaio rotor travado
v1 = travado['Vab']
i1 = travado['Ib'] / np.sqrt(3) / nesp
p1 = (travado['Pal'] + travado['Pbe']) / nesp / 3
q1 = np.sqrt((v1*i1)**2 - p1**2)
print(f'travado: {v1 = :.3f}V, {i1 = :.3f}A, {p1 = :.3f}W, {q1 = :.3f}Var')

req = p1 / i1**2
xeq = q1 / i1**2

# print(f'{req = :.3f} Ohm')
# print(f'{xeq = :.3f} Ohm')

r2 = req - r1
x2 = xeq / 1.68
x1 = xeq - x2



# aproximação das Prot
p = (livre['Pal'] + livre['Pbe']) / nesp
v = livre['Vab']

ks = np.polyfit(v**2, p, 1)

prot = ks[1]
print(f'{prot = :.3f}')

vf = np.linspace(0, 400, 200)
pf = ks[0]* vf**2 + ks[1]

# plt.plot(v, p, marker='*', linestyle='')
# plt.plot(vf, pf)
# plt.ylim(0, 300)
# plt.xlim(0, 400)
# plt.show()

v1 = livre['Vab'].iloc[0]
i1 = livre['Ib'] .iloc[0]/ np.sqrt(3) / nesp
p1 = (livre['Pal'].iloc[0] + livre['Pbe'].iloc[0]) / nesp / 3
q1 = np.sqrt((v1*i1)**2 - p1**2)
print(f'livre: {v1 = :.3f}V, {i1 = :.3f}A, {p1 = :.3f}W, {q1 = :.3f}Var')

z1 = r1 + 1j*x1
i1f = i1 * np.exp(-1j*math.atan2(q1, p1))
e0 = abs(v1 - i1f*z1)
prf = p1 - i1**2*r1 - prot/3
qxm = q1 - i1**2*x1
rf = e0**2 / prf
xm = e0**2 / qxm
print(f'{r1 = :.3f} Ohm')
print(f'{r2 = :.3f} Ohm')
print(f'{x1 = :.3f} Ohm')
print(f'{x2 = :.3f} Ohm')
print(f'{rf = :.3f} Ohm')
print(f'{xm = :.3f} Ohm')

# curva torque
np.seterr('ignore')
v = 380
ns = 1800
nr = np.linspace(-200, ns+100, 500)
s = (ns - nr) / ns
rconv = r2 * (1-s) / s

z0 = 1j*xm*rf / (1j*xm + rf)
z1 = r1 + 1j*x1
z2 = r2 + rconv + 1j*x2

z02 = z0*z2 / (z0 + z2)
zeq = z02 + z1

i1 = v1 / zeq
e0 = i1 * z02
i2 = e0 / z2

pef = 3*abs(i2)**2*r2/s
ws = ns * 2 * np.pi / 60
wr = nr * 2 * np.pi / 60
tind = pef / ws
teixo = tind - prot/ws * abs(wr/ws)
peixo = teixo * wr

carga = pd.read_csv('ensaio_carga.csv', comment='#')
g = 9.8
carga['teixo'] = carga['F']*g * .278

plt.figure(2)
plt.plot(nr, tind)
# plt.plot(nr, teixo)
plt.plot(carga['nr'], carga['teixo'], marker='*', linestyle='')
plt.title('torque')
plt.axhline(color='black')
plt.axvline(color='black')
plt.axvline(ns,color='black')

plt.figure(3)
plt.plot(nr, abs(i1)*np.sqrt(3))
plt.title('In')
plt.axhline(color='black')
plt.axvline(color='black')
plt.axvline(ns,color='black')

plt.figure(4)
plt.plot(nr, peixo)
plt.title('Peixo')
plt.axhline(color='black')
plt.axvline(color='black')
plt.axvline(ns,color='black')

plt.show()


