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
    
    # Cambio de modo de operación
    if t < 5:
        D = 0          # Batch
    else:
        D = F_feed / V_reactor
    
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
    dS = D * (S_in - S_calc) - qs * X_calc
    dP = alpha * dX + qp * X_calc - D * P_calc
    dTr = D * (T_feed - Tr) + (rQ / (rho * Cp)) - (UA * (Tr - Tj)) / (rho * V_reactor * Cp)
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
Tr0 = 30 #[°C]
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
D_final = F_feed / V_reactor

#Resumen numerico
print(f"Biomasa final: {Biomasa[-1]:.3f} g/L (maxima: {Biomasa.max():.3f} g/L en t={Tiempo[np.argmax(Biomasa)]:.2f} h)")
print(f"Sustrato final: {Sustrato[-1]:.3f} g/L")
print(f"Producto (lactato) final: {Producto[-1]:.3f} g/L")
print(f"Productividad volumétrica final: {D_final*Producto[-1]:.3f} g/L·h")
print(f"Temperatura del reactor: min={T_reactor.min():.2f} °C, max={T_reactor.max():.2f} °C")

#Graficas
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
    
#%%
'''Cultivo continuo con recirculacion de biomasa biorreactor 20L HCW - optimizacion'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import shgo

def objective (x):
    S_in, F_feed = x
    def continuo_jacket(t, Y):
        X, S, P, Tr, Tj, I, P_acumulado = Y
        
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
        dP_acumulado = F_feed * P_calc
        return [dX, dS, dP, dTr, dTj, dI, dP_acumulado]

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
    Tr0 = 30 #[°C]
    Tj0 = 25 #[°C]
    I0 = 0 #[°C*h]
    P_acumulado0 = 0 #[g/L]
    array_iniciales = np.array([X0, S0, P0, Tr0, Tj0, I0, P_acumulado0])

    #Tiempo de ejecucion
    t_start = 0
    t_stop = 24
    tspan = (t_start, t_stop)
    t_array = np.linspace(t_start, t_stop, num=10000)

    #Metodo numerico
    solucion = solve_ivp(continuo_jacket, tspan, array_iniciales, t_eval = t_array)

    if not solucion.success:
        return 0.0 

    # Cálculos finales de masa de producto
    P_reactor_final = solucion.y[2][-1] * V_reactor  # Gramos que se quedaron en el reactor
    P_reactor_inicial = P0 * V_reactor              # Gramos con los que empezamos
    P_cosechado_final = solucion.y[6][-1]            # Gramos que salieron en el flujo continuo
    
    # Masa neta total producida (g)
    masa_total_producida = P_cosechado_final + P_reactor_final - P_reactor_inicial
    
    # Productividad Global Real g/(L·h)
    PROD_global = masa_total_producida / (V_reactor * t_stop)
    
    return -PROD_global

print("Optimizando Sfeed y F_feed...")
bounds = [(10, 300), (0.1, 5.0)]
result = shgo(objective, bounds=bounds)
Sfeed_opt, F_feed_opt = result.x
prod_max = -result.fun

print(f"\n{'='*45}")
print(f"  Sfeed óptimo          : {Sfeed_opt:.4f} g/L")
print(f"  F_feed óptimo         : {F_feed_opt:.4f} L/h")
print(f"  Productividad máxima  : {prod_max:.6f} g/(L·h)")
print(f"{'='*45}\n")
#%%
'''Cultivo continuo con recirculacion de biomasa biorreactor 20L HCW - sintonizacion del controlador'''

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Parámetros del sistema 
V_reactor = 20
S_in, F_feed, T_feed = 10, 0.1, 25.0
T_setpoint = 30.0
rho, Cp, Yqs = 1000.0, 4.182, 3963.0
UA, V_jacket, Tj_entrada = 75*3600, 2.0, 25
F0, F_min, F_max = 5, 0.0, 10.0

# Parámetros cinéticos 
miu_max, qs_max, qp_max = 1.09, 4.16, 1.863
Ksx, Kix, Kpx = 4.229, 394.20, 5.001
Kss, Kis, Kps = 0.15, 143.391, 20.07
Ksp, Kip, Kpp = 0.065, 373.89, 42.83
alpha, Kd = 0.017, 0.0001

# Modelo dinámico para la simulación
def modelo_control(t, Y, Kp, Ti):
    X, S, P, Tr, Tj, I, P_acumulado = Y
    X_calc, S_calc, P_calc = max(0, X), max(0, S), max(0, P)

    # Cinética
    miu = (miu_max * S_calc * Kix) / ((Ksx + S_calc) * (Kix + S_calc)) * np.exp(-P_calc/Kpx)
    qs = (qs_max * S_calc * Kis) / ((Kss + S_calc) * (Kis + S_calc)) * np.exp(-P_calc/Kps)
    qp = (qp_max * S_calc * Kip) / ((Ksp + S_calc) * (Kip + S_calc)) * np.exp(-P_calc/Kpp)
    rQ = Yqs * qs * X_calc # Generación de calor 
    
    D = F_feed/V_reactor if (t>=5) else 0

    # Controlador PI (Variable manipulada: Fc)
    Error = Tr - T_setpoint
    Fc_control = F0 + Kp * Error + (Kp/Ti) * I
    Fc = np.clip(Fc_control, F_min, F_max)
    
    #Ecuaciones diferenciales de las variables de estado
    dX = (miu - Kd)* X_calc #Porque se tienen membrana ideal
    dS = D * (S_in - S_calc) - qs * X_calc
    dP = alpha * dX + qp * X_calc - D * P_calc
    dTr = D * (T_feed - Tr) + (rQ / (rho * Cp)) - (UA * (Tr - Tj)) / (rho * V_reactor * Cp)
    dTj = (F/V_jacket) * (Tj_entrada - Tj) + (UA * (Tr - Tj)) / (rho * V_jacket * Cp)
    dI = 0 if (Fc_control > F_max and Error > 0) or (Fc_control < F_min and Error < 0) else Error
    return [dX, dS, dP, dTr, dTj, dI]

# Función de costo: Minimización del IAE (Error Integral Absoluto) 
def objetivo_iae(parametros):

    Kp_opt, Ti_opt = parametros

    if Kp_opt <= 0 or Ti_opt <= 0:
        return 1e10

    y0 = [0.43, 33, 0, 30, 25, 0]

    t_span = (0,24)
    t_eval = np.linspace(0,24,200)

    sol = solve_ivp(
        modelo_control,
        t_span,
        y0,
        args=(Kp_opt,Ti_opt),
        t_eval=t_eval,
        method='RK45'
    )

    if not sol.success:
        return 1e10

    iae = np.trapezoid(
        np.abs(sol.y[3]-T_setpoint),
        sol.t
    )

    return iae

# Ejecución de la sintonización
print("Sintonizando controlador... Por favor espere.")
inv_inicial = [9.59, 3.66] # Punto de partida literario [1]
resultado = minimize(objetivo_iae, inv_inicial, method='Nelder-Mead', bounds=[(0.1, 50), (0.1, 20)])

Kp_final, Ti_final = resultado.x
print(f"\n--- Sintonización Completada ---")
print(f"Kp óptimo: {Kp_final:.4f} L/h·°C")
print(f"Ti óptimo: {Ti_final:.4f} h")
print(f"IAE Mínimo: {resultado.fun:.4f}")