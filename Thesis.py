

from numpy import *
import matplotlib.pyplot as plt
import timeit


# In[40]:


#Constantes comunes proteinas

#(d/dt)g = a + b/(1+(g/k)**h) - gamma*g

a = 20       #1/min
b = 250      #1/min
gamma = 1/30 #1/min
k = 5000
h = 1.3


# In[41]:


def st(t_old, g_old):
    g = g_old

    ks=[]        #Constantes de probabilidad
    cltive=[]    #Constantes acumulativas

    # (d/dt)EAi = a + b/(1+(g/k)**h) - gamma*g

    k0 = a                  #Produccion constitutiva de g                          g+= 60
    k1 = b/(1+(g/k)**h)     #Termino de represion de la producion de g             g+= 60
    k2 = gamma*g            #Degradacion de g                                      g+= -1

    ks.append(k0)
    ks.append(k1)
    ks.append(k2)

    #Constante de normalización, suma de todas las probabilidades
    s = 0.0
    for j in ks:
        s = s+j

    #Creamos el arreglo de probabilidades acumulativas para poder definir cada evento
    for j in range(len(ks)):
        if(j==0):
            cltive.append(ks[0]/s)
        else:
            zzz = float(ks[j]/s)
            cltive.append(cltive[j-1] + zzz)


    #Usamos el número aleatorio r1 para generar el tiempo en el que ocurre el próximo evento de interés
    r1 = random.random()
    T = (1/s)*log(1/r1) + t_old

    #Usamos el número aleatorio ran para generar el próximo evento de interés
    ran = random.random()

    #Aquí decidimos cual de los posibles eventos tomará lugar
    if(ran<cltive[0]):
        g+= 60

    elif(cltive[0]<=ran<cltive[1]):
        g+= 60

    elif(cltive[1]<=ran<=cltive[2] and g!=0):
        g+= -1

    #Retornamos los valores actualizados
    return T, g


def cell(hours, dt, g_0):

    #Definimos los arreglos del tiempo y el compuesto como vacios
    T = []
    g = []

    #Agregamos los valores iniciales de cada uno de los compuestos a los arreglos
    t = 0.0

    T.append(0.0)
    g.append(g_0)

    #Simulamos la celula usando la función step (st)
    while t < hours*60.0:

        #Seguimos usando los mismos nombres para las variables, excepto que ya no serán los valores iniciales.
        t, g_0 = st(t, g_0)
        #Guardamos los valores generados
        T.append(t)
        g.append(g_0)

    #Para hacer el promedio entre muchas células tenemos que estandarizar el tiempo
    #Aquí es donde usamos dt para calcular cuantos intervalos se quieren
    T_s = linspace(0, hours*60.0, int(hours*60.0/dt))

    #Inicializamos el arreglo g estandarizado
    g_s = zeros(len(T_s))

    #Aquí estandarizamos para que se pueda hacer un promedio
    for i in range(len(T_s)):
        j = 0
        while(T[j] < T_s[i]):
            j+=1

        g_s[i] = g[j]

    #Finalmente retornamos los arreglos, primero los estandarizados.
    return T_s, g_s, T, g


def sv_cells(N_cells, hours, dt, g_0):

    #Definimos la contidad de puntos en el eje temporal
    n_points = int(hours*60.0/dt)
    #Definimos una matriz para g
    g = zeros((N_cells, n_points))

    #Ahora llenamos las matrices, corriendo N_cells veces
    for i in range(N_cells):

        #Nombramos los resultados y extraemos los arreglos estandarizados
        stuff = cell(hours, dt, g_0)

        T_s = stuff[0]
        g[i,:] = stuff[1]

    #Iniciamos el arreglo del promedio
    ave_g = zeros( n_points)

    #Llenamos los arreglos de los promedios respectivos para cada punto temporal
    for k in range(n_points):
        ave_g[k] = average(g[:,k])

    #Finalmente retornamos los arreglos
    return T_s, ave_g


# In[42]:


hours = 8
dt = 0.5
#70000 weird, d1=0.437
g_0 = 0

T_s, g_s, T, g = cell(hours, dt, g_0)

#Aqui visualizamos los resultados
plt.figure()
plt.plot(T,g,label='g',linewidth=2)
plt.xlabel('Time(mins)' , fontsize = 20)
plt.ylabel('# of molecules', fontsize = 20)
plt.title('Stochastic Simulation for 1 cell', fontsize = 30)
plt.savefig('1cell.png')
plt.show()


# In[43]:


num = 0.0

alpha1 = 0.0
alpha2 = 0.0
alpha3 = 0.0
alpha4 = 0.0

c0 = 1 - 11/(6*h) -  1/(6*h**2)
c1 = 3/h - 5/(2*h**2) + 1/(2*h**3)
c2 = -3/h + 4/(h**2) - 1/(h**3)
c3 = 2/h - 3/(h**2) + 1/(h**3)

for i in range(len(g)):
    if T[i]>100:

        num +=1
        alpha1 += (g[i]/k)**h
        alpha2 += (g[i]/k)**(2*h)
        alpha3 += (g[i]/k)**(3*h)
        alpha4 += (g[i]/k)**(4*h)

alpha1 = alpha1/num
alpha2 = alpha2/num
alpha3 = alpha3/num
alpha4 = alpha4/num

d1 = b*h/(2*k)*(7/4 - alpha1 + alpha2/4) + a*h/k - gamma*h*(c0 + c1*alpha1 + c2*alpha2/2 + c3*alpha3/6 )

aaa = (6/((c3+0.0001)*gamma*h))*((b*h/(2*k)*(7/4 - alpha1 + alpha2/4) + a*h/k - gamma*h*(c0 + c1*alpha1 + c2*alpha2/2 )))

aaaa = (3/(gamma*h*(c3+0.0001)))*( ((b*h**2)/(2*k**2))*(7/4 - alpha1 + alpha2/4) + (b*h/k)*(7/4*alpha1 - alpha2 + alpha3/4) + (k*gamma*(h/k)**2)*(c0 + c1*alpha1 + c2*alpha2/2 + c3*alpha3/6 ) - 2*gamma*h*(c0*alpha1 + c1*alpha2 + c2*alpha3/2) + a*((h/k)**2 + 2*h*alpha1/k))

print( d1, aaa, alpha3, (aaa-alpha3)*100/aaa, (aaaa-alpha4)*100/aaaa)


# In[32]:





# In[16]:


#El número de celulas que se quieren promedia
N_cells = 5
#Cuanto tiempo queremos correr el programa
hours = 8.0
#Ventana de tiempo entre pasos
dt = 0.5
#Valor inicial
g_0 = 70000

T_s, ave_g = sv_cells(N_cells, hours, dt, g_0)


#Graficamos
plt.figure()
plt.plot(T_s,ave_g,label='g')
plt.xlabel('Time(mins)' , fontsize = 20)
plt.ylabel('# of molecules', fontsize = 20)
plt.title('Stochastic Simulation for %d cells' % (N_cells), fontsize = 30)
plt.legend(fontsize = 10)
plt.savefig('N_cells.png')
plt.show()


# In[5]:


import numpy as np
a = [-2,-3,-4,0]
b = [2,3,5,4]
c = np.array(a)+np.array(b)

bo = any(l < 0 for l in c)
bo

print (c, bo)
