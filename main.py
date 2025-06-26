import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

R1 = 14.7
nesp = 11
cat_factor = 0.68  # cat N

df = pd.read_csv('ensaio_livre_travado.csv', comment="#")

travado = df.loc[0]
v1 = travado['Vab']
p1 = (travado['Pal'] + travado['Pbe']) / nesp / 3
i1 = travado['Ib'] / nesp / np.sqrt(3)
q1 = np.sqrt((v1 * i1) ** 2 - p1 ** 2)

print(f'travado: {v1 = :.3f}V, {i1 = :.3f}A, {p1 = :.3f}W, {q1 = :.3f}')


Req = p1 / i1**2
Xeq = q1 / i1 ** 2

# print(f'{Req = :.3f}')
# print(f'{Xeq = :.3f}')

R2 = Req - R1


X2 = Xeq / (1+cat_factor)
X1 = Xeq - X2




livre = df.loc[1:]
vm = livre['Vab']
pm = (livre['Pal'] + livre['Pbe'])/ nesp
im = livre['Ib'] / nesp
#
#
ks = np.polyfit(vm**2, pm, 1)
# print(ks)
#
vf = np.linspace(0,390,200)
pf = ks[0]*vf**2 + ks[1]

Prot = ks[1]





v1 = vm.iloc[0]
p1 = pm.iloc[0] / 3
i1 = im.iloc[0] / np.sqrt(3)
q1 = np.sqrt((v1 * i1) ** 2 - p1 ** 2)

print(f'livre: {v1 = :.3f}V, {i1 = :.3f}A, {p1 = :.3f}W, {q1 = :.3f}')

Z1 = np.complex64(R1, X1)
i1f = i1*np.exp(1j*np.atan2(-q1, p1))
e0f = v1 - i1f*Z1
e0 = abs(e0f)

# print(e0)

prf = p1 - i1**2*R1 - Prot/3
Rf = e0**2/prf

qxm = q1 - i1**2*X1
Xm = e0**2/qxm

print(f'{Prot = :.3f}')

print(f'{R1 = :.3f}')
print(f'{R2 = :.3f}')
print(f'{X1 = :.3f}')
print(f'{X2 = :.3f}')
print(f'{Rf = :.3f}')
print(f'{Xm = :.3f}')


#
#
# #
# #
# # plt.plot(vm, pm, marker='*', linestyle='')
# # plt.plot(vf, pf)
# # plt.ylim(0,350)
# # # # plt.xlim(0,250)
# # plt.show()



v1 = 380
ns = 1800

nr = np.linspace(0, ns-0.001, 500)
s = (ns-nr)/ns


rconv = R2 * (1-s) / s
z0 = np.complex64(Rf, Xm)
z1 = np.complex64(R1, X1)
z2 = (np.complex64(0, X2) + rconv)

z20 = z2*z0/(z2 + z0)
zeq = z20 + z1

i1 = v1 / zeq
e0 = i1 * z20
i2 = e0 / z2

pef = 3 * abs(i2)**2 * R2 / s


ws = ns * np.pi * 2 / 60
wr = nr * np.pi * 2 / 60
tind = pef / ws

teixo = tind - (1-s)*Prot/ws

carga = pd.read_csv('ensaio_carga.csv', comment='#')
carga['teixo'] = carga['F'] * .5425/2
print(carga)

plt.plot(nr, tind)
plt.plot(nr, teixo)
plt.plot(carga['nr'], carga['teixo'], marker='*', linestyle='')
plt.show()


# tau = np.zeros_like(nr)
# i1 = np.zeros_like(nr)
# peixo = np.zeros_like(nr)
# pconv = np.zeros_like(nr)




