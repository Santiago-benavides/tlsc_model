#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 22:21:32 2020

@author: santi
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

#se deben ingresar y:condidiones iniciales. t: tiempo de simulación y elresto
#son condiciones fisicias. se deben definir antes
def pend(y, t, Rl, L, w, Ud, Uq, Uz, Vpd, Vpq, Vpz, C, R0):
    # y tiene las condiciones iniciales del problema
    Id,Iq,Iz,Vdc,Ev=y
    # dydt = tiene las ecuaciones del modelo
    dydt= [-(Rl/L)*Id+w*Iq-(Vdc/(2*L))*Ud+Vpd/L,
           -(Rl/L)*Iq-w*Id-(Vdc/(2*L))*Uq+Vpq/L,
           -(Rl/L)*Iz-(3**0.5/L)*Ev-(Vdc/(2*L))*Uz+Vpz/L,
           +(1/C)*Id*Ud+(1/C)*Iq*Uq+(1/C)*Iz*Uz-Vdc/(C*R0),
           +(3**0.5/(2*C))*Iz-(1/(C*R0))*Ev]
    return dydt


f=50
w=2*np.pi*f
Rl=0.5
L=30.0e-3
C=2200e-6
Vpd=120*2**0.5
Vpq=0
Vpz=0
R0=10000
Vdc=0
# Id=0
# Iq=4.3741

#se las considera variables de entrada
Ud=0.9999
Uq=0
Uz=0

#las condiciones iniciales
y0 = [0,0,0,0,0]
#el rango de tiempo que se requeire la simulación
t = np.linspace(0, 6, 10000)
#se resuelve el sistema
sol = odeint(pend, y0, t, args=(Rl, L, w, Ud, Uq, Uz, Vpd, Vpq, Vpz, C, R0))

#en sol quedan todas las variables del problema resuelto.
plt.figure(1)
plt.plot(t, sol[:, 0], 'C0', label='Id')
plt.plot(t, sol[:, 1], 'C1', label='Iq')
plt.plot(t, sol[:, 2], 'C2', label='Iz')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()

plt.figure(2)
plt.plot(t, sol[:, 3], 'C3', label='Vdc')
plt.plot(t, sol[:, 4], 'C4', label='Ev')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()