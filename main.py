import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math

nesp = 11
r1 = 14.7
df = pd.read_csv('ensaio_livre_travado.csv',
                 comment='#')
travado = df.loc[0]


v1 = travado['Vab']
i1 = travado['Ib'] / np.sqrt(3) / nesp
p1 = (travado['Pal']+travado['Pbe']) / 3 / nesp

q1 = np.sqrt((i1*v1)**2 - p1**2)
print(f'travado: {v1 = :.3f} V, {i1 = :.3f} A, {p1 = :.3f} W, {q1 = :.3f} VAr')


req = p1 / i1**2
r2 = req - r1


xeq = q1 / i1**2
x2 = xeq / 1.68  # cat N
x1 = xeq - x2





# determinar Prot
livre = df.loc[1:]
vm = livre['Vab']
pm = (livre['Pal'] + livre['Pbe']) / nesp

ks = np.polyfit(vm**2, pm, 1)

prot = ks[1]
print(f'{prot = :.3} W')

vf = np.linspace(0, 400, 300)
pf = ks[0]*vf**2 + ks[1]


# plt.plot(vm, pm, marker='*', linestyle='')
# plt.plot(vf, pf)
# plt.xlim(0, 400)
# plt.ylim(0, 300)
# plt.xlabel('V')
# plt.ylabel('Perdas')
# plt.show()


# ensaio rotor livre
v1 = livre['Vab'].iloc[0]
i1 = livre['Ib'].iloc[0] / np.sqrt(3) / nesp
p1 = (livre['Pal'].iloc[0]+livre['Pbe'].iloc[0]) / 3 / nesp
q1 = np.sqrt((i1*v1)**2 - p1**2)
print(f'livre: {v1 = :.3f} V, {i1 = :.3f} A, {p1 = :.3f} W, {q1 = :.3f} VAr')

z1 = r1 + 1j*x1
phi = -math.atan2(q1, p1)
i1f = i1 * np.exp(1j*phi)
e0f = v1 - i1f*z1
e0 = abs(e0f)
prf = p1 - i1**2*r1 - prot/3
rf = e0**2 / prf

qxm = q1 - i1**2*x1
xm = e0**2 / qxm

print(f'{r1 = :.3f} Ohm')
print(f'{r2 = :.3f} Ohm')
print(f'{x1 = :.3f} Ohm')
print(f'{x2 = :.3f} Ohm')
print(f'{rf = :.3f} Ohm')
print(f'{xm = :.3f} Ohm')


v = 380
ns = 1800
np.seterr('ignore')
nr = np.linspace(-1000, ns+100, 500)

s = (ns - nr) / ns
rconv = r2 * (1-s)/s


z0 = 1j*xm*rf / (1j*xm +rf)
z1 = r1 + 1j*x1
z2 = r2 + rconv + 1j*x2

z02 = z0*z2 / (z0 + z2)
zeq = z02 + z1

i1 = v1 / zeq

e0 = z02*i1
i2 = e0/z2

ws = ns*2*np.pi/60
pef = 3 * abs(i2)**2 * r2/s
tind = pef / ws

carga = pd.read_csv('ensaio_carga.csv', comment='#')

carga['teixo'] = carga['F'] * 9.8 * .278
carga['Peixo'] = carga['nr'] * 2*np.pi/60 * carga['teixo']


plt.figure(1)
plt.plot(nr, abs(i1)*np.sqrt(3))
plt.ylabel('nr')
plt.ylabel('Iin')
plt.grid()
plt.axvline(color='black')
plt.axhline(color='black')
plt.axvline(ns, color='black')

plt.figure(2)
plt.plot(nr, tind)
plt.plot(carga['nr'], carga['teixo'], marker='*', linestyle='')
plt.ylabel('nr')
plt.ylabel('tind')
plt.grid()
plt.axvline(color='black')
plt.axhline(color='black')
plt.axvline(ns, color='black')

wr = nr * 2*np.pi/60
pconv = tind * wr
plt.figure(3)
plt.plot(nr, pconv)
plt.ylabel('nr')
plt.ylabel('pconv')
plt.grid()
plt.axvline(color='black')
plt.axhline(color='black')
plt.axvline(ns, color='black')



plt.show()

