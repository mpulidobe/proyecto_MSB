'''Cultivo continuo con recirculacion de biomasa biorreactor 20L'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def continuo_jacket(t, Y):
    X, S, P, Tr, Tj, I = Y
    
    X_calc = max(0, X)
    S_calc = max(0, S)
    P_calc = max(0, P)
    
    #Tasa especifica de crecimiento
    miu = (miu_max * S_calc * Kix * np.exp(-P_calc/Kpx))/((Ksx + S_calc)*(Kix + S_calc))
    qs = (qs_max * S_calc * Kis * np.exp(-P_calc/Kps))/((Kss + S_calc)*(Kis + S_calc))
    qp = (qp_max * S_calc * Kip * np.exp(-P_calc/Kpp))/((Ksp + S_calc)*(Kip + S_calc))
    
    #Tasa volumetrica de generacion de energia
    rQ = (Yqs * miu * X_calc) / Yxs
    
    #Controlador
    Error = Tr - T_setpoint
    F_control = F0 + Kp * Error #+ (Kp/Ti) * I
    F = np.clip(F_control, F_min, F_max)
    
    #Evitar windup
    if F_control > F_max and Error > 0:
        dI = 0
    elif F_control < F_min and Error < 0:
        dI = 0
    else:
        dI = Error
        
    #Ecuaciones diferenciales de las variables de estado
    dX = (miu - Kd)* X_calc #Porque se tienen membrana ideal
    dS = (F_feed/V_reactor) * (S_in - S_calc) - qs * X_calc
    dP = alpha * dX + qp * X_calc - (F_feed/V_reactor) * P_calc
    dTr = (F_feed/V_reactor) * (Tr0 - Tr) + (rQ / (rho * Cp)) - (UA * (Tr - Tj)) / (rho * V_reactor * Cp)
    dTj = (F/V_jacket) * (Tj_entrada - Tj) + (UA * (Tr - Tj)) / (rho * V_jacket * Cp)
    return [dX, dS, dP, dTr, dTj, dI]

##Parametros
V_reactor = 20 #[L]
S_in = 34 #[g/L]
F_feed = 10 #[L/h]
Kd = 0.0001 #coeficiente de muerte celular [h^-1]
miu_max = 1.09 #tasa de crecimiento especifica maxima [h^-1]
qs_max = 4.16 #tasa de utilizacion de sustrato especifica maxima [g/g*h]
qp_max = 1.863 #tasa de produccion de lactato especifica maxima [g/g*h]
alpha = 0.017 #constante asociada al crecimiento en Luedeking-Piret [g/g]

#Constantes de limitacion de sustrato
Ksx = 4.229 #Limitacion de sustrato para el crecimiento de la biomasa [g/L]
Kss = 0.15 #Limitacion de sustrato para el consumo de sustrato [g/L]
Ksp = 0.065 #Limitacion de sustrato para la produccion de lactato [g/L]

#Constantes de inhibicion por sustrato
Kix = 394.20 #Inhibicion del sustrato para el crecimiento de la biomasa [g/L]
Kis = 143.391 #Inhibicion del sustrato para el consumo de sustrato [g/L]
Kip = 373.89 #Inhibicion del sustrato para la produccion de lactato [g/L]

#Constantes de inhibicion por producto
Kpx = 5.001 #Inhibicion del producto para el crecimiento de la biomasa [g/L]                 
Kps = 20.07 #Inhibicion del producto para el consumo de sustrato [g/L]
Kpp = 42.83 #Inhibicion del producto para la produccion de lactato [g/L]

#Propiedades fisicas
rho = 1000 #Densidad del medio [g/L]
Cp = 4.182 #Capacidad calorifica del medio [J/g*°C]

#Chaqueta de enfriamiento
V_jacket = 2 #Volumen de la chaqueta [L]
Tj_entrada = 10 #[°C]
UA = 75 * 3600 #[J/h*°C]

#Rendimientos (Ypx, Yps, Yxs)
Yps = 0.72 #Rendimiento de producto[g/g]
Yxs = 0.074 #Rendimiento de biomasa [g/g]
Ypx = Yps/Yxs #Rendimiento de producto [g/g]
Yqs = 3963 #Rendimiento termico [J/g]

#Parametros del controlador PI
T_setpoint = 30 #[°C]
Kp = 0.5 #[L/h*°C]
Ti = 3.66 #[h]
F0 = 0 #[L/h]
F_min = 0 #[L/h]
F_max = 10 #[L/h]

#Condiciones iniciales
X0 = 0.5 #[g/L]
S0 = 0 #[g/L]
P0 = 0 #[g/L]
Tr0 = 25 #[°C]
Tj0 = 25 #[°C]
I0 = 0 #[°C*h]
array_iniciales = np.array([X0, S0, P0, Tr0, Tj0, I0])

#Tiempo de ejecucion
t_start = 0
t_stop = 80
tspan = (t_start, t_stop)
t_array = np.linspace(t_start, t_stop, num=10000)

#Metodo numerico
solucion = solve_ivp(continuo_jacket, tspan, array_iniciales, t_eval = t_array)

Tiempo = solucion.t
Biomasa = solucion.y[0]
Sustrato = solucion.y[1]
Producto = solucion.y[2]
T_reactor = solucion.y[3]
T_jacket = solucion.y[4]
Integral_error = solucion.y[5]
Error = T_reactor - T_setpoint
F_valor = F0 + Kp * Error + (Kp/Ti)*Integral_error
F = np.clip(F_valor, F_min, F_max)

#Grafica 
f1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.plot(Tiempo, Biomasa, label = 'Biomasa', color = 'red')
ax1.plot(Tiempo, Sustrato, label = 'Sustrato', color = 'blue')
ax1.plot(Tiempo, Producto, label = 'Producto', color = 'green')
ax2.plot(Tiempo, T_reactor, label = 'Temperatura reactor', color = 'orange')
ax2.plot(Tiempo, T_jacket, label = 'Temperatura chaqueta', color = 'purple')
ax2.axhline(T_setpoint, color = 'darkred', linestyle='--', label = 'Set point')
ax3.plot(Tiempo, Error, label = 'Error', color = 'gold')
ax4.plot(Tiempo, F, label = 'Flujo real', color = 'magenta')

for i in [ax1, ax2, ax3, ax4]:
    i.grid()
    i.legend()