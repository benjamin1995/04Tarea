#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import scipy as sp
from planeta import Planeta

'''
PARTE2

En este script se realizara el grafico de energia y trayectoria usando euler.
MIENTRAS MAS CHICO EL DT MAS PRECISO ES PARA LA ENERGIA PERO LA TRAYECTORIA NO

'''


'''
 Condiciones iniciales mas la creacion del objeto de clase Planeta.
 Paso arbitrario.
 Arreglos de ceros a los que se les ira anexando valores dado el metodo de euler.
'''

vy0=0.3
condicion_inicial = sp.array([10, 0, 0, vy0])
p = Planeta(condicion_inicial)
t_final =  6000
pasos = 50000
dt= t_final / (float)(pasos)
x = sp.zeros(pasos)
y = sp.zeros(pasos)
vx = sp.zeros(pasos)
vy = sp.zeros(pasos)
energia=sp.zeros(pasos)

[x[0],y[0],vx[0],vy[0]] = condicion_inicial
energia[0]=p.energia_total()

'''Caso euler explicito con alpha=0 .Se recorrera con un for para ir obteniendo para cada tiempo un vector de 4 componentes [x,y,vx,vy] '''


for i in range (1,pasos):
    resultenergy=p.energia_total()
    p.avanza_euler(dt)
    resultados = p.y_actual
    x[i] = resultados[0]
    y[i] = resultados[1]
    vx[i] = resultados[2]
    vy[i] = resultados[3]
    energia[i]=resultenergy


''' Graficos de los resultados obtenidos con el metodo de euler alpha=0 '''

plt.figure(1)
plt.clf()
plt.subplot(2, 1, 1)
plt.subplots_adjust(hspace=.5)
plt.scatter(x,y,15,c=u'darkred',label = "trayectoria")
plt.title("Solucion euler alfa=0")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.ylim([-100,50])
plt.legend(loc='best',fontsize=7)
t_values = sp.linspace(1,t_final,pasos)
plt.subplot(2, 1, 2)
plt.plot(t_values,energia,label='energia total')
plt.xlabel('Tiempo [seg]')
plt.ylabel('Energia [Joule]')
plt.ylim([-0.2,0.2])
plt.plot(t_values,(1/2.)*(vx**2+vy**2),'g',label='energia cinetica')
plt.plot(t_values,-1/(sp.sqrt(x**2+y**2)),'orange',label='energia potencial')
plt.legend(loc='best',fontsize=6.5)
plt.title("Energia vs tiempo")
plt.savefig('Euler.png')
plt.show()
