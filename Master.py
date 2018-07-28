
import timeit
import numpy as np
import matplotlib.pyplot as plt


#Esta función refresca los valores haciendo suceder una instancia (evento) escogida al azar entre los sucesos probables
#Teniendo en cuenta sus probabilidades
#ENTRADAS: tiempo antes del paso "t_old", el arreglo de valores de las variables del sistema antes del paso "xs",
# funcion "darEventos" que retorna el arreglo de las probabilidades de cada evento dado "xs" (diferente para cada sistema)
# el arreglo "actions" de acciones posibles dado el evento. Es decir, hay una probabilidad "darEventos(xs)[i]" para la accion "actions[i]"
#mas informacion en los codigos de cada sistema
def st(t_old, xs, darEventos, actions):

    #Se generan las probabilidades de los eventos teniendo en cuenta las ecuaciones
    ev = darEventos(xs)
    #Se crea un arreglo acumulativo de las probabilidades
    cltive=[]
    #Número de eventos posibles
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
    #Aquí decidimos cual de los posibles eventos tomará lugar teniendo en cuenta que ninguna variable tome un valor negativo
    #ponemos como condicion no tener valores negativos para ninguna molecula
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



#Esta funcion simula el comportamiento de una celula, guardando su trayectoria en el espacio de moleculas contra tiempo paso a paso
#retorna los arreglos de tiempo y moleculas estandarizados y no estandarizados (para poder comparar con otras celulas)
#ENTRADAS: numero de horas a simular "hours", resolucion minima de tiempo para estandarizar "dt", condiciones iniciales "init" de las moleculas
# funcion "darEventos" que retorna el arreglo de las probabilidades de cada evento dado "xs" (diferente para cada sistema)
# el arreglo "actions" de acciones posibles dado el evento.
def cell(hours, dt, init, darEventos, actions):

    #Definimos el arreglo de tiempo y generamos una variable de tiempo
    #Definimos el arreglo donde ubicaremos la evolucion de las variables
    T = []
    variables = []
    t = 0.0
    #se agrega el tiempo "0"
    T.append(0.0)
    variables.append(init)

    #Simulamos la celula usando la función paso (st)
    while t < hours*60.0:

        #Seguimos usando los mismos nombres para las variables, excepto que ya no serán los valores iniciales.
        t, init = st(t, init, darEventos, actions)
        #Guardamos los valores generados
        T.append(t)
        #pasar a lista: cosa.tolist()
        variables.append(init.tolist())


    #Para hacer el promedio entre muchas células tenemos que estandarizar el tiempo
    #Aquí es donde usamos dt para calcular cuantos intervalos se quieren
    T_s = np.linspace(0, hours*60.0, int(hours*60.0/dt))

    #Inicializamos la matriz de variables a estandarizar en el tiempo.
    vars_s = np.matrix(np.zeros((int(len(T_s)),int(len(init)))))
    variables = np.matrix(variables)

    #Aquí estandarizamos para que se pueda hacer un promedio
    #Se toman ventanas con la resolucion dt, y para cada una de ellas se toma el valor de moleculas en el tiempo justo anterior
    # al limite derecho del intervalo
    j = 0
    for i in range(len(T_s)):
        while(T[j] < T_s[i]):
            j+=1
        for k in range(len(init)):
            vars_s[i,k] = variables[j,k]


    #Se retorna el arreglo de tiempo, un arreglo de arreglos "variables", el arreglo de tiempo estandarizado y
    #una matriz con el valor de cada variable para cada tiempo estandarizado.
    return T, variables, T_s, vars_s


#Esta funcion simula el comportamiento de muchas celulas, retorna el arreglo tiempo estandarizado, el arreglo del comportamiento promedio
# para las diferentes moleculas, y la matriz de evolucion de los compuestos (numero celulas, tiempo, molecula).
#ENTRADAS: numero de celulas a simular "N_cells", horas a simular "hours", resolucion temporal minima "dt", condiciones iniciales "init" de las moleculas
# funcion "darEventos" que retorna el arreglo de las probabilidades de cada evento dado "xs" (diferente para cada sistema),
# el arreglo "actions" de acciones posibles dado el evento.
def sv_cells(N_cells, hours, dt, init, darEventos, actions):

    #Definimos la contidad de puntos en el eje temporal
    n_points = int(hours*60.0/dt)

    #Para hacer el promedio entre muchas células tenemos que estandarizar el tiempo
    #Aquí es donde usamos dt para calcular cuantos intervalos se quieren
    T_s = np.linspace(0, hours*60.0, int(hours*60.0/dt))

    #Definimos un cubo de informacion donde la primera dimension recorre las copias de los procesos,
    #la segunda la evolucion temporal y la tercera las variables.
    Many_vars = np.zeros((N_cells, n_points, len(init)))

    #Ahora llenamos las varias matrices de evolucion de una celula, corriendo N_cells veces
    for i in range(N_cells):

        #Nombramos los resultados y extraemos los arreglos estandarizados
        stuff = cell(hours, dt, init, darEventos, actions)

        Many_vars[i,:,:] = stuff[3]

    #Iniciamos el arreglo del promedio
    ave_variables = np.zeros(( n_points, len(init)))

    #Llenamos los arreglos de los promedios respectivos para cada punto temporal
    ave_variables = np.mean(Many_vars, axis=0)
    #Finalmente retornamos el arreglo de tiempo y la matriz de evolucion de los compuestos
    return T_s, ave_variables, Many_vars

#Esta funcion calcula la densidad de probabilidad para las diferentes variables y retorna un arreglo con cada una de ellas. "dens(info)[i]" es la
# funcion de probabilidad (no normalizada) de la molecula i-esima. La posicion la indica el arreglo inicial y el arreglo de acciones.
#ENTRADAS: cubo de informacion "info"; primera entrada: copia de la celula simulada, segunda entrada: paso temporal, tercera entrada: molecula.
def dens(info):
    #Se inicializa el arreglo que tendra las funciones de probabilidad
    density = []
    #Definicion de variables a partir del
    num_cells = len(info[:,0,0])
    num_tiempos = len(info[0,:,0])
    num_variables = len(info[0,0,:])

    #Se recorre para el numero de variables (o moleculas)
    for k in range(num_variables):
        #Se verifica el mayor numero que alcanza a la variable (mayor numero de moleculas alcanzado)
        lim_max = info[:,:,k].max()

        #se inicializa la matriz de probabilidad
        prob = np.zeros((int(lim_max +1), int(num_tiempos)))
        #Se suma un "1" cada vez que la funcion pasa por un punto de la cuadricula. De ese modo se generan los histogramas para cada tiempo.
        for n in range(num_tiempos):
            for c in range(num_cells):
                molecs_actual = info[c,n,k]
                #aca la ventaja es que el indice de prob[indice,tiempo] es a la vez el numero de moleculas. El valor de prob[m,t] final será
                #el conteo de funciones que pasan por el punto (t,m). Con t = tiempo, m = moleculas.

                #si no se hiciera asi, sino teniendo en cuenta una resolucion dm, necesitariamos otro arreglo para guardar
                #el label de cada entrada, para cuando grafiquemos.
                #como lo que necesitamos hacer con el tiempo.
                prob[int(molecs_actual),n] += 1.0
        density.append(prob)

    #Se retorna un arreglo donde cada entrada es la funcion de probabilidad (no normalizada) de cada variable.
    return density
