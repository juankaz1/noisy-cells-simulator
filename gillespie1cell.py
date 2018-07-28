
# coding: utf-8

# In[52]:


import numpy as np
import matplotlib.pyplot as plt
#from ecuaciones import *


# In[38]:

#acciones de creacion o destruccion de Variables
actions = [[60.0, 0.0], [-1.0, 0.0], [0.0, 62.0], [0.0, -1.0]]

#probabilidades
def darEventos(xs):
    return [20, 1/30.0*xs[0], 25, 1/32.0*xs[1]]



# In[94]:


#Esta función refresca los valores haciendo suceder una instancia escogida al azar entre los sucesos probables
#Teniendo en cuenta sus probabilidades
#ENTRADAS: tiempo antes del paso "t_old" y el arreglo de valores de las variables del sistema antes del paso "xs".
def st(t_old, xs):

    #Se generan las probabilidades de los eventos teniendo en cuenta las ecuaciones
    ev = darEventos(xs)
    #Se crea un arreglo acumulativo de las probabilidades
    cltive=[]
    #Número de eventos
    n = len(ev)

    #Se calcula la constante de normalización de las probabilidades de los eventos
    s = 0.0
    for j in ev:
        s = s+j

    #Se normalizan las probabilidades de los eventos y se genera el arreglo acumulativo de ellas
    for j in range(n):
        if(j==0):
            cltive.append(ev[0]/s)
        else:
            zzz = float(ev[j]/s)
            cltive.append(cltive[j-1] + zzz)

    #Usamos el número aleatorio r1 para generar el tiempo en el que ocurre el próximo evento de interés
    r1 = np.random.random()
    T = (1/s)*np.log(1/r1) + t_old
    #Usamos el número aleatorio ran para generar el próximo evento de interés
    ran = np.random.random()
    #Aquí decidimos cual de los posibles eventos tomará lugar teniendo en cuenta que ninguna variable sea negativa
    for i in range(n):
        if (i > 0):
            if (cltive[i-1]<ran<=cltive[i]):
                if ( any( l < 0 for l in np.array(xs)+np.array(actions[i]))==False):
                    xs = np.array(xs)+np.array(actions[i])
        elif(ran<=cltive[0]):
            if (any( l < 0 for l in np.array(xs)+np.array(actions[i]))==False):
                xs = np.array(xs)+np.array(actions[i])

    #Devolvemos los valores actuales de tiempo y concentración de compuestos (Variables)
    return T, xs

#Falta


def cell(hours, init):

    #Definimos el arreglo de tiempo y generamos una variable de tiempo
    T = []
    variables = []
    t = 0.0
    T.append(0.0)
    variables.append(init)


    #Simulamos la celula usando la función paso (st)
    while t < hours*60.0:

        #Seguimos usando los mismos nombres para las variables, excepto que ya no serán los valores iniciales.
        t, init = st(t, init)
        #Guardamos los valores generados
        T.append(t)
        #pasar a lista: cosa.tolist()
        variables.append(init.tolist())

    #Se retorna el arreglo de tiempo y un arreglo de arreglos "variables".
    return T, variables


# In[109]:


initVals = [0.0, 0.0]
hours = 5

T, vals = cell(hours, initVals)
vals = np.matrix(vals)


# In[110]:


#Aqui visualizamos los resultados
fig = plt.figure(figsize=(15,10))
plt.plot(T,vals[:,0],label='p1',linewidth=2)
plt.plot(T,vals[:,1],label='p2',linewidth=2)
plt.xlabel('Time(mins)' , fontsize = 20)
plt.ylabel('# of molecules', fontsize = 20)
plt.title('Stochastic Simulation for 1 cell', fontsize = 30)
plt.legend()
#plt.savefig('1cell.png')
plt.show()
