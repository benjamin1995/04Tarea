import matplotlib.pyplot as plt
import scipy as sp
from planeta import Planeta
import scipy.stats

'''
PARTE2
En este script se realizaran graficos de energia y trayectoria usando verlet velocity
INCREIBLE, NECESITA UN DT MUUUUUUUY CHICO PARA TENER LA TRAYECTORIA DE VERLET2 Y RK4, SIN EMBARGO, EN LA ENERGIA SIGUE OSCILANDO ALREDEDOR DEL -0.05

'''

'''
Caso verlet velocity con alpha=0
Caso interesante voy=0.1, se escapa, graficos de energiat con kine y potencial + trayectoria
Condiciones iniciales, definicion de objeto, paso y arreglos de ceros a llenar
'''

vy0=0.3
condicion_inicial = [10, 0, 0, vy0]
p = Planeta(condicion_inicial)
t_final =  6000
pasos = 50000
dt= t_final / (float)(pasos)
x = sp.zeros(pasos)
y = sp.zeros(pasos)
vx = sp.zeros(pasos)
vy = sp.zeros(pasos)
energia=sp.zeros(pasos)
x[0],y[0],vx[0],vy[0] = condicion_inicial
energia[0]=p.energia_total()

''' Calculo de vectores usando metodo de verlet velocity'''
for i in range (1,pasos):
    resultenergy=p.energia_total()
    p.avanza_verletvelocity(dt)
    resultados = p.y_actual
    x[i] = resultados[0]
    y[i] = resultados[1]
    vx[i] = resultados[2]
    vy[i] = resultados[3]
    energia[i]=resultenergy


''' Grafico de los resultados, en particular de la energia en funcion del tiempo y de la trayectoria con metodo de verletvelocity para alpha=0'''

plt.figure(1)
plt.clf()
plt.subplot(2, 1, 1)
plt.subplots_adjust(hspace=0.6)
plt.scatter(x,y,10,c= u'darkgreen',label='trayectoria')
plt.title("Solucion verlet velocity alfa=0")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.legend(loc="best",fontsize=7)
t_values = sp.linspace(1,t_final,pasos)
plt.subplot(2, 1, 2)

'''opcional: plt.scatter(t_values,energia,label='energia total') '''

plt.plot(t_values,energia,label='energia total')
plt.title("Energia vs tiempo")
plt.xlabel('Tiempo [seg]')
plt.ylabel('Energia [Joule]')
plt.ylim([-0.2,0.2])
plt.plot(t_values,(1/2.)*(vx**2+vy**2),'g',label='energia cinetica')
plt.plot(t_values,-1/(sp.sqrt(x**2+y**2)),'orange',label='energia potencial')
plt.legend(loc="best",fontsize=7)
plt.savefig('Verletvelocity.png')
plt.show()
