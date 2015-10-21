#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy as sp


''' Constantes definidas'''
G=1
M=1
m=1

class Planeta(object):
    '''
    Complete el docstring.
    '''

    def __init__(self, condicion_inicial, alpha=0):
        '''
        __init__ es un método especial que se usa para inicializar las
        instancias de una clase.

        Ej. de uso:
        >> mercurio = Planeta([x0, y0, vx0, vy0])
        >> print(mercurio.alpha)
        >> 0.
        '''
        self.y_actual = condicion_inicial
        self.t_actual = 0.
        self.alpha = alpha

    def ecuacion_de_movimiento(self,datos=sp.array([0,0,0,0])):
        '''
        Implementa la ecuación de movimiento, como sistema de ecuaciónes de
        primer orden.

        '''
        x, y, vx, vy = self.y_actual

        deltax = datos[0]
        deltay = datos[1]
        deltavx = datos[2]
        deltavy = datos[3]

        x += deltax
        y += deltay
        vx += deltavx
        vy += deltavy

        ''' Se definen las funciones de la aceleracion en base a calculo analitico'''
        a=self.alpha
        fx =  -x*((sp.sqrt(x**2 +y**2))**-3) +  (2*a*x*((x**2+y**2)**-2))
        fy =  -y*((sp.sqrt(x**2 +y**2))**-3) +  (2*a*y*((x**2+y**2)**-2))
        return sp.array([vx, vy, fx, fy])

    def avanza_euler(self,dt):
        '''
        Toma la condición actual del planeta y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de Euler explícito. El
        método no retorna nada, pero re-setea los valores de self.y_actual.

        Tiene un error orden 2 y la forma de iterar es ocupando el vector yactual de la ecuacion de movimiento
        con sus 4 componentes y avanzar cada una de ellas multiplicando el vector por el paso y sumandole el yactual.

        Finalmente se actualiza yactual y se avanza en el tiempo.

        '''

        ynext=self.y_actual + dt*(self.ecuacion_de_movimiento())
        self.y_actual=ynext
        self.t_actual+=dt

        pass

    def avanza_rk4(self, dt):
        '''
        Similar a avanza_euler, pero usando Runge-Kutta 4.
        Error de orden 5. Simplemente se usa la formula para cada componente del vector yactual y se obtiene el siguiente.
        Finalmente se actualizan los valores actuales del vector yactual(4 componentes) y se avanza en el tiempo
        '''
        k1=self.ecuacion_de_movimiento()
        k2=self.ecuacion_de_movimiento(k1*dt/2.)
        k3=self.ecuacion_de_movimiento(k2*dt/2.)
        k4=self.ecuacion_de_movimiento(k3*dt)
        ynext=self.y_actual +(1/6.)*(k1+2*k2+2*k3+k4)*dt
        self.y_actual=ynext
        self.t_actual+=dt
        pass

    def avanza_verletvelocity(self, dt):
        '''
        Similar a avanza_euler, pero usando Verlet-velocity.
        '''


        '''Se define una variable para cada valor actual del vector de 4 componentes'''
        xn , yn, vxn, vyn = self.y_actual


        '''calculo de vector posiciones en x y en y'''
        yn1 = sp.array([ xn , yn ]) + sp.array([ vxn , vyn ]) * dt + 1/2. * self.ecuacion_de_movimiento()[2:] * dt**2


        '''Calculo de vector aux de 4 componentes a usar en ecuacion de movimiento para la aceleracion en x y en y '''
        Yn1 = sp.array([yn1[0] , yn1[1], 0, 0])


        ''' Calculo de vector velocidades en x y en y'''
        vn1 = sp.array([ vxn , vyn ]) + 1/2. * (self.ecuacion_de_movimiento()[2:] + self.ecuacion_de_movimiento(Yn1)[2:]) * dt


        '''Nuevo arreglo actual, avanza en el tiempo '''
        actual = sp.array([yn1[0], yn1[1], vn1[0], vn1[1]])
        self.y_actual = actual
        self.t_actual += dt
        pass

    def verlet_2(self,dt,x_p,y_p):
        '''
        En este metodo de verlet alternativo(verlet2) se requieren 2 puntos para iniciar. Se ocupara tambien rk4 para empezar al correrlo

        '''

        ''' Requiere de valores previos y actual para iterar. Se usara rk4 para obtener el actual y asi comenzar'''
        Y_previo = sp.array([x_p,y_p])


        ''' A cada variable se le asigna una componente del vector yactual'''
        x , y, vx, vy = self.y_actual


        ''' Vector que guarda las posiciones'''
        Yn = sp.array([x,y])


        ''' Vector de posicion de dos componentes con error de orden 4'''
        Yn1 = (2 * Yn) - (Y_previo) + dt**2 * self.ecuacion_de_movimiento()[2:]


        ''' Vector de velocidad de 2 componentes con error de orden 2 (mas grande que orden 4),Yprevio es n-1,Yn es n, Yn1 es n+1'''
        Vn = (Yn1 - Y_previo) / (2*dt)


        ''' Nuevo arreglo (actual) '''
        nuevo = sp.array([Yn1[0],Yn1[1],Vn[0], Vn[1]])


        '''Avanza en el tiempo y actualizo yactual(vector de 4 componentes)'''
        self.y_actual = nuevo
        self.t_actual +=dt
        pass

    def energia_total(self):
        '''
        Calcula la energia total del sistema en las condiciones actuales.
        '''

        '''Se define una variable para cada valor actual del vector de 4 componentes '''
        x, y, vx, vy = self.y_actual


        ''' Energia potencial definicion en funcion de x e y'''
        potencial= (-G*m*M)*(sp.sqrt(x**2+y**2)**-1)+ self.alpha*G*M*m/((x**2+y**2))


        ''' Energia cinetica en funcion de las velocidades vx y vy'''
        kine=(vx**2 + vy**2)*m/2.


        ''' Energia total del sistema'''
        energiaT= kine + potencial

        return energiaT
        pass
