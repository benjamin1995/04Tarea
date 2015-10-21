import matplotlib.pyplot as plt
import scipy as sp
from planeta import Planeta

'''
PARTE2

En este script se realizara el grafico de energia y trayectoria para rk4
'''

''' Caso rk4 con alpha=0.Se crea al objeto de la clase Planeta. Se define paso y arreglos de cero listos para llenar '''

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

''' Condiciones iniciales'''

[x[0],y[0],vx[0],vy[0]] = condicion_inicial
energia[0]=p.energia_total()


''' Recorre con el metodo de rk4 el objeto p de la clase Planeta para ir obteniendo el vector [x,y,vx,vy] para cada tiempo dado por el paso '''
for i in range (1,pasos):
    resultenergy=p.energia_total()
    p.avanza_rk4(dt)
    resultados = p.y_actual
    x[i] = resultados[0]
    y[i] = resultados[1]
    vx[i] = resultados[2]
    vy[i] = resultados[3]
    energia[i]=resultenergy


''' graficos de los resultados de metodo rk4 alpha=0'''
plt.figure(1)
plt.clf()
plt.subplot(2, 1, 1)
plt.subplots_adjust(hspace=.5)
plt.scatter(x,y,5,c=u'darkorange',label = "trayectoria")
plt.title("Solucion rk4 alfa=0 ")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.ylim([-15,15])
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
plt.legend(loc='best',fontsize=6.3)
plt.savefig('Rk4.png')
plt.show()
