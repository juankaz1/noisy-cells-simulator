import matplotlib.pyplot as plt
import numpy as np
x = np.random.randn(1000)
y = np.random.randn(1000) + 5

# normal distribution center at x=0 and y=5
plt.hist2d(x, y, bins=100, normed = True)
plt.xlabel('Time(mins)' , fontsize = 20)
plt.ylabel('Frecuency', fontsize = 20)
plt.title('Probability density for T (' + str(3)+ ' cells)', fontsize = 25)
plt.colorbar()
plt.savefig('Riboswitch/worked.png')
plt.show()

"""

#ERRORES CON ARREGLOS DINAMICOS

#Funcion que genera los histogramas para cada tiempo en cada molecula
def hist(g, dm, i, k):

    global dens

    lim_mols = len(dens[i,k])
    inf = int(g/dm)
    if ((inf +1) > lim_mols):
        faltante_temp = np.zeros(inf+1-lim_mols)
        dens[i,k] = dens[i,k]+faltante_temp
    dens[i,k,inf] +=1



#Funcion que estandariza la funcion densidad
def dens_s(N_cells):

    global dens

    num_vars = len(dens[0])
    num_tiempo = len(dens)
    for k in range(int(num_vars)):
        l_max = max([len(x) for x in dens[:,k]])
        #cosa = max(dens[:,k], key=len)
        #l_max = len(cosa)
        for i in range(int(num_tiempo)):
            a = int(len(dens[i,k]))
            faltante = np.zeros(l_max-a)
            dens[i,k] = dens[i,k]+faltante
    return dens/N_cells





a = [1,2,3,4,5,6,7,8,9]

b = [1,2,3,4,5]
c = [0]


#genero mi funcion densidad con altura 5 para todos los valores de tiempo (eje-x), o bueno, cada valor de indice muestra un histograma.
dens = np.zeros((len(a),5))
dens

#El problema es: tengo una funcion y=f(x), en este caso los valores de y(x) son a[i].
#Mi funcion histograma quiere sumar un "1" a mi densidad de frecuencia "dens" cada vez que la funcion pasa por una de sus divisiones.
#en un punto la funcion supera el numero de divisiones, para lo cual la funcion hist cuenta con el condicional "if".
#Este condicional mira a ver cuantos ceros le faltan a "dens" en ese punto, los agrega (o eso intenta) y luego ahi si suma uno en esa posicion.
def hist(g, i):

    global dens
    dens = np.array(dens)
    lim_mols = len(dens[i])
    inf = g
    #los conteos de indices de las siguientes dos lineas estan bien.
    if ((inf +1) > lim_mols):
        faltante_temp = np.zeros(inf+1-lim_mols, dtype = dens[i].dtype)
        print(dens[i])
        dens[i] = np.append(dens[i],faltante_temp)
        print(dens[i])
        print('did add for '+str(g))
    dens[i,inf] +=1
    print('did it for '+str(g))
#aqui lo hago efectivo. En efecto muestra que agrega el "1" para a[0], a[1], a[2]. Pero luego no, lo cual revela un problema al agrandar el arreglo
# de la posicion i-esima para i >= 3. generar un arreglo de distinto tamaño que los anteriores es el error.
for i in range(len(a)):
    hist(a[i], i)
"""
