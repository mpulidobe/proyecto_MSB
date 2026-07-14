'''Cultivo continuo con recirculacion de biomasa biorreactor 20L HCW'''

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
    rQ = Yqs * qs * X_calc
    
    #Controlador
    Error = Tr - T_setpoint
    F_control = F0 + Kp * Error + (Kp/Ti) * I
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
    dTr = (F_feed/V_reactor) * (T_feed - Tr) + (rQ / (rho * Cp)) - (UA * (Tr - Tj)) / (rho * V_reactor * Cp)
    dTj = (F/V_jacket) * (Tj_entrada - Tj) + (UA * (Tr - Tj)) / (rho * V_jacket * Cp)
    return [dX, dS, dP, dTr, dTj, dI]

##Parametros
V_reactor = 20 #[L]
S_in = 10 #[g/L]
F_feed = 1 #[L/h]
T_feed = 25 #[°C]
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
Kp = 9.59 #[L/h*°C]
Ti = 3.66 #[h]
F0 = 5 #[L/h]
F_min = 0 #[L/h]
F_max = 10 #[L/h]

#Condiciones iniciales
X0 = 0.43 #[g/L]
S0 = 33 #[g/L]
P0 = 0 #[g/L]
Tr0 = 25 #[°C]
Tj0 = 25 #[°C]
I0 = 0 #[°C*h]
array_iniciales = np.array([X0, S0, P0, Tr0, Tj0, I0])

#Tiempo de ejecucion
t_start = 0
t_stop = 24
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
Volumen = np.ones_like(Tiempo) * V_reactor

Error = T_reactor - T_setpoint
F_valor = F0 + Kp * Error + (Kp/Ti)*Integral_error
F = np.clip(F_valor, F_min, F_max)

print(f"Volumen del reactor: {V_reactor:.2f} L (constante)")
print(f"Tasa de dilución (D): {F_feed/V_reactor:.3f} h⁻¹")
print(f"Biomasa final: {Biomasa[-1]:.3f} g/L (máxima: {Biomasa.max():.3f} g/L en t={Tiempo[np.argmax(Biomasa)]:.2f} h)")
print(f"Sustrato final: {Sustrato[-1]:.3f} g/L")
print(f"Producto (lactato) final: {Producto[-1]:.3f} g/L")
print(f"Productividad volumétrica final: {(F_feed/V_reactor)*Producto[-1]:.3f} g/L·h")
print(f"Temperatura del reactor: min={T_reactor.min():.2f} °C, max={T_reactor.max():.2f} °C")
print(f"Flujo de refrigerante: min={F.min():.2f} L/h, max={F.max():.2f} L/h")

#Grafica 
f1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.plot(Tiempo, Biomasa, label = 'Biomasa', color = 'red')
ax1.plot(Tiempo, Sustrato, label = 'Sustrato', color = 'blue')
ax1.plot(Tiempo, Producto, label = 'Producto (lactato)', color = 'green')
ax1.set_xlabel('Tiempo [h]'); ax1.set_ylabel('Concentracion [g/L]')

ax2.plot(Tiempo, T_reactor, label = 'Temperatura reactor', color = 'orange')
ax2.plot(Tiempo, T_jacket, label = 'Temperatura chaqueta', color = 'purple')
ax2.axhline(T_setpoint, color = 'darkred', linestyle='--', label = 'Set point')
ax2.set_xlabel('Tiempo [h]'); ax2.set_ylabel('Temperatura [°C]')

ax3.plot(Tiempo, Error, label = 'Error', color = 'gold')
ax3.set_xlabel('Tiempo [h]'); ax3.set_ylabel('Error [°C]')

ax4.plot(Tiempo, F, label = 'Flujo refrigerante real', color = 'magenta')
ax4.plot(Tiempo, Volumen, label='Volumen reactor', color='teal')
ax4.set_xlabel('Tiempo [h]')

for i in [ax1, ax2, ax3, ax4]:
    i.grid()
    i.legend()