import csv
import matplotlib.pyplot as plt
import numpy as np

def variables(f):
    lector = csv.reader(f,delimiter=' ')
    T,px,py,pz,p=[],[],[],[],[]
    for fila in lector:
        aux = [float(x) for x in fila if x]
        T.append(aux[0])
        px.append(aux[1])
        py.append(aux[2])
        pz.append(aux[3])
        p.append(aux[4])
    return T,px,py,pz,p

with open("polmacro-0.20-sol","r") as f:
    T1,px1,py1,pz1,p1=variables(f)
with open("polmacro-0.20-esfsol","r") as t:
    T2,px2,py2,pz2,p2=variables(t)
with open("polmacro-0.20-esf","r") as g:
    T3,px3,py3,pz3,p3=variables(g)

with open("polmacro-0.40-sol","r") as h:
    T4,px4,py4,pz4,p4=variables(h)
with open("polmacro-0.40-esfsol","r") as j:
    T5,px5,py5,pz5,p5=variables(j)
with open("polmacro-0.40-esf","r") as k:
    T6,px6,py6,pz6,p6=variables(k)

T1=1.45*np.array(T1)
T2=1.45*np.array(T2)
T3=1.45*np.array(T3)
T4=1.45*np.array(T4)
T5=1.45*np.array(T5)
T6=1.45*np.array(T6)

plt.figure(figsize=(8,6))

plt.xlabel(r"Temperatura $[K]$", size=18)
plt.ylabel(r"Polarización $[\mu C/cm^2]$", size=18)
plt.ylim(-1,20)
plt.xticks(np.arange(0, 901, 100))
plt.minorticks_on()
plt.tick_params(axis='both', which='major', labelsize=18)

plt.plot(T1,p1,label=r"$S0$", marker='o')
plt.plot(T2,p2,label=r"$S50$", marker='o')
plt.plot(T3,p3,label=r"$S100$", marker='o')


plt.grid(True,linestyle='--', alpha=0.5)
plt.legend(fontsize=14)
plt.savefig('Pol-macro-0.20.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()

plt.figure(figsize=(8,6))

plt.xlabel(r"Temperatura $[K]$", size=18)
plt.ylabel(r"Polarización $[\mu C/cm^2]$", size=18)
plt.ylim(-1,15)
plt.xticks(np.arange(0, 901, 100))
plt.minorticks_on()
plt.tick_params(axis='both', which='major', labelsize=18)

plt.plot(T4,p4,label=r"$S0$", marker='o')
plt.plot(T5,p5,label=r"$S50$", marker='o')
plt.plot(T6,p6,label=r"$S100$", marker='o')

plt.grid(True,linestyle='--', alpha=0.5)
plt.legend(fontsize=14)
plt.savefig('Pol-macro-0.40.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()


