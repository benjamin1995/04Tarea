import matplotlib.pyplot as plt
import scipy as sp
from planeta import Planeta

'''
PARTE2
En este script se graficara la energia y trayectoria usando verlet version 2.
caso verlet avanzado con alpha 0 con dt muy pequeno queda muy preciso y perfecto, energia constante
'''

''' Condiciones iniciales, se define objeto de clase Planeta, el paso de tiempo y los arreglos de ceros de las distintas variables'''
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

''' Se define el segundo valor del vector yactual mediante rk4 (en este caso) puesto que verlet2 lo necesita para la iteracion'''
p.avanza_rk4(dt)
resultados = p.y_actual
x[1] = resultados[0]
y[1] = resultados[1]
vx[1] = resultados[2]
vy[1] = resultados[3]

''' Se obtienen los arreglos con todos sus valores mediante el metodo de verlet 2'''
for i in range (2,pasos):
    p.verlet_2(dt,x[i-2],y[i-2])
    resultados = p.y_actual
    x[i] = resultados[0]
    y[i] = resultados[1]
    vx[i] = resultados[2]
    vy[i] = resultados[3]
    energia[i] = p.energia_total()

''' Se grafican los resultados, en este caso energia en funcion del tiempo y la trayectoria para el metodo de verlet 2 con alpha=0'''
plt.figure(2)
plt.subplot(2, 1, 1)
plt.subplots_adjust(hspace=.5)
plt.scatter(x,y,6,c=u'darkred',label = "trayectoria")
plt.title(" Solucion verlet 2 alfa = 0 ")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.ylim([-20,20])
plt.legend(loc='best',fontsize=7)
t_values = sp.linspace(1,t_final,pasos)
plt.subplot(2, 1, 2)
plt.plot(t_values,energia,label='energia total')
plt.title("Energia vs tiempo")
plt.xlabel('Tiempo [seg]')
plt.ylabel('Energia [Joule]')
plt.ylim([-0.2,0.2])
plt.plot(t_values,(1/2.)*(vx**2+vy**2),'g',label='energia cinetica')
plt.plot(t_values,-1/(sp.sqrt(x**2+y**2)),'orange',label='energia potencial')
plt.legend(loc='best',fontsize=7)
plt.savefig('Verlet2.png')
plt.show()
