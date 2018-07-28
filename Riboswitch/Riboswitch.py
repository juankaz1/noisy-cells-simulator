import matplotlib.pyplot as plt
import Master as ma
import numpy as np

#RIBOSWITCH
#Las siguientes ecuaciones diferenciales dictaminan la evolucion estocastica del sistema en cuestion.

# (d/dt)Ta = KT + mum*R - mup*E*Ta - gammaT*Ta
# (d/dt)P = KP*Ta - gammaP*P
#(d/dt)R = mup*E*Ta - gammaR*R - mum*R

#Definiciones de constantes
KT = 20.0
gammaT = 1/5.0
gammaR = 1/5.0
KP = 1.0
gammaP = 1/30.0
mum = 1/3.0
mup = 1/20000.0
#Molecula de entrada (nota: pensar como acoplar variaciones en el tiempo)
E = 10000.0

#Acciones de creacion o destruccion de Variables. Hay 7 eventos posibles
actions = [[1.0, 0.0, 0.0], [1.0, 0.0, -1.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]

#Probabilidades para este sistema dadas las ecuaciones. Hay 7 eventos posibles.
def darEventos(xs):
    return [KT, mum*xs[2], mup*E*xs[0], gammaT*xs[0], KP*xs[0], gammaP*xs[1], gammaR*xs[2]]


#CAMBIAR AL GUSTO
#Valores iniciales de cada compuesto
initVals = [0.0, 0.0, 0.0]
#horas
hours = 10
#resolucion minima en minutos, NO CAMBIAR POR EL MOMENTO.
#nota: para hacer que esto funcione bien toca usar numeros imaginarios.
dt = 1.0
#Numero de celulas a simular
N_cells = 100000
#nota: pensar como acoplar resolucion minima de moleculas "dm"

#Se genera la simulacion para las condiciones dadas. Definiendo asi arreglos de tiempo y promedios estandarizados.
#Junto con el cubo de informacion "many" que abarca toda la simulacion.
T_s, s_vals, many = ma.sv_cells(N_cells, hours, dt, initVals, darEventos, actions)
s_vals = np.matrix(s_vals)
#Se definen las funciones de probabilidad no normalizadas
densidad = ma.dens(many)


#Aqui visualizamos los resultados

#Comportamiento promedio
fig = plt.figure(figsize=(15,10))
plt.plot(T_s,s_vals[:,0],label='T',linewidth=2)
plt.plot(T_s,s_vals[:,1],label='P',linewidth=2)
plt.plot(T_s,s_vals[:,2],label='R',linewidth=2)
plt.xlabel('Time(mins)' , fontsize = 20)
plt.ylabel('# of molecules', fontsize = 20)
plt.title('Stochastic Simulation for ' + str(N_cells)+ ' cells', fontsize = 30)
plt.legend()
plt.savefig('Riboswitch/'+str(N_cells)+'cells.png')
plt.show()



#Se guarda la informacion sin normalizar en archivos de texto para mirar luego si se desea.
np.savetxt('Riboswitch/RIBOSWITCH__T_s('+ str(N_cells)+'cells).txt', T_s, delimiter=',')
np.savetxt('Riboswitch/RIBOSWITCH__dens_T('+ str(N_cells)+'cells).txt', densidad[0], delimiter=',')
np.savetxt('Riboswitch/RIBOSWITCH__dens_P('+ str(N_cells)+'cells).txt', densidad[1], delimiter=',')
np.savetxt('Riboswitch/RIBOSWITCH__dens_R('+ str(N_cells)+'cells).txt', densidad[2], delimiter=',')


#Generamos las graficas de las distribuciones de probabilidad (normalizadas) para las diferentes variables.
fig = plt.figure(figsize=(9,7))
plt.imshow(np.matrix(densidad[0]/N_cells), origin='lower', aspect='auto', cmap='hot' )
plt.xlabel('Time(mins)' , fontsize = 20)
plt.ylabel('Molecules', fontsize = 20)
plt.title('Probability density for T (' + str(N_cells)+ ' cells)', fontsize = 20)
plt.colorbar(label = 'Normed frecuency')
plt.savefig('Riboswitch/RIBOSWITCH__prob_T.png')
plt.show()

fig = plt.figure(figsize=(9,7))
plt.imshow(np.matrix(densidad[1]/N_cells), origin='lower', aspect='auto', cmap='hot' )
plt.xlabel('Time(mins)' , fontsize = 20)
plt.ylabel('Molecules', fontsize = 20)
plt.title('Probability density for P (' + str(N_cells)+ ' cells)', fontsize = 20)
plt.colorbar(label = 'Normed frecuency')
plt.savefig('Riboswitch/RIBOSWITCH__prob_P.png')
plt.show()

fig = plt.figure(figsize=(9,7))
plt.imshow(np.matrix(densidad[2]/N_cells), origin='lower', aspect='auto', cmap='hot' )
plt.xlabel('Time(mins)' , fontsize = 20)
plt.ylabel('Molecules', fontsize = 20)
plt.title('Probability density for R (' + str(N_cells)+ ' cells)', fontsize = 20)
plt.colorbar(label = 'Normed frecuency')
plt.savefig('Riboswitch/RIBOSWITCH__prob_R.png')
plt.show()
