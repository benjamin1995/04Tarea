import matplotlib.pyplot as plt
import scipy as sp
from planeta import Planeta
import scipy.stats
import numpy as np


'''
PARTE3
Ahora considero el caso alfa no cero, mi rut es 189563752 y se usara verletvelocity para integrar la ecuacion de movimiento.
Luego se determinara la posicion del perihelio y la velocidad angular de precesion. Finalmente se ploteara el grafico de la energia en funcion
del tiempo y de la trayectoria.

'''

vy0=0.3
condicion_inicial = sp.array([10, 0, 0, vy0])
p = Planeta(condicion_inicial, 10**(-2.375))
t_final =  180*30                      #factor 30 por las orbitas
pasos = 40000
dt= t_final / (float)(pasos)

x = sp.zeros(pasos)
y = sp.zeros(pasos)
vx = sp.zeros(pasos)
vy = sp.zeros(pasos)
r = sp.zeros(pasos)
energia = sp.zeros(pasos)
perihelio = [[], [],[] ]             # se trabajara con dt,dy y dx

''' Condiciones iniciales'''
[x[0],y[0],vx[0],vy[0]] = condicion_inicial
r[0] = sp.sqrt(x[0]**2+y[0]**2)
energia[0] = p.energia_total()

''' Creacion de arreglos de valores para [x,y,vx,vy] con metodo de verlet velocity'''
for i in range (1,pasos):
    p.avanza_verletvelocity(dt)
    resultados = p.y_actual
    x[i] = resultados[0]
    y[i] = resultados[1]
    vx[i] = resultados[2]
    vy[i] = resultados[3]
    energia[i] = p.energia_total()
    r[i] = sp.sqrt(x[i]**2+y[i]**2)

'''Tratamos de buscar el perihelio'''
    a = 0.0009
    valor =10
    if valor-a<r[i-1] and r[i-1]<valor+a:
        perihelio[0].append(p.t_actual)
        perihelio[1].append(x[i-1])
        perihelio[2].append(y[i-1])


'''Obtenemos la velocidad angular de precesion'''

vel_angular = sp.zeros(len(perihelio[0])-1)

phi_anterior =  sp.arctan(perihelio[2][0]/(perihelio[1][0]))   #Variable para calcular el angulo phi por trigonometria
t_anterior = perihelio[0][0]




''' calculo de la velocidad angular como deltaphi sobre deltat'''
for i in range(1,len(perihelio[0])):
    phi = sp.arctan(perihelio[2][i]/(perihelio[1][i]))
    t = perihelio[0][i]
    deltaphi = phi - phi_anterior
    deltat = t-t_anterior

    vel_angular[i-1] = deltaphi/deltat

    t_anterior = t
    phi_anterior = phi

vel_angular = np.round(vel_angular,7)                                      #Redondea al septimo decimal el arreglo de velocidades
vel_precesion = scipy.stats.mode(vel_angular)[0][0]                        #Encuentra la moda en un arreglo
print("Arreglo de velocidad angular de precesion = "+(str)(vel_angular))
print("velocidad angular de precesion = "+(str)(vel_precesion))




''' Grafico de los resultados, en este caso de la trayectoria y la energia con metodo de verlet velocity, ademas de la posicion del perihelio'''

plt.figure(1)
plt.subplot(2, 1, 1)
plt.subplots_adjust(hspace=0.8)

'''opcional:plt.scatter(x,y,c=u'gold',label = "trayectoria")'''

plt.plot(x,y,'gold',label = "trayectoria")
plt.title("Integracion verlet $\\alpha = 10^{-2.375}$")
plt.xlabel("x[m]")
plt.ylabel("y[m]")
plt.scatter(perihelio[1],perihelio[2],26,c=u'red',label="posicion perihelio")
plt.legend(loc="best",fontsize=10)
t_values = sp.linspace(0,t_final,pasos)
plt.subplot(2, 1, 2)
plt.plot(t_values,energia,label='energia total')
plt.title("Energia vs tiempo")
plt.plot(t_values,(1/2.)*(vx**2+vy**2),'g',label='energia cinetica')
plt.plot(t_values,-1/(sp.sqrt(x**2+y**2)),'orange',label='energia potencial')
plt.legend(loc="best",fontsize=7)
plt.xlabel("Tiempo [seg]")
plt.ylabel("Energia [Joule]")
plt.savefig("Precesionperihelio.png")
plt.show()
