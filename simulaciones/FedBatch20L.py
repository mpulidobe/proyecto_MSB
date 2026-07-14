'''Cultivo fedbatch biorreactor 20L HCW '''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def fedbatch_jacket(t, Y):
    X, S, P, Tr, Tj, I, V = Y

    X_calc = max(0, X)
    S_calc = max(0, S)
    P_calc = max(0, P)
    V_calc = max(1e-6, V)   

    #Tasa especifica de crecimiento
    miu = (miu_max * S_calc * Kix * np.exp(-P_calc/Kpx))/((Ksx + S_calc)*(Kix + S_calc))
    qs = (qs_max * S_calc * Kis * np.exp(-P_calc/Kps))/((Kss + S_calc)*(Kis + S_calc))
    qp = (qp_max * S_calc * Kip * np.exp(-P_calc/Kpp))/((Ksp + S_calc)*(Kip + S_calc))

    #Tasa volumetrica de generacion de energia metabolica
    rQ = Yqs * qs * X_calc

    #Alimentacion
    if t < 5 :
        F=0 #Fase Inicial Bacth
        dV=0
    else:
        F=F_feed if V_calc < 20 else 0 #Inicio fase fed-bacth
        dV=F
        
    D=F/V_calc #Tasa de dilución   
    
    #Balances de Masa 
    dX = (miu - Kd) * X_calc - D * X_calc
    dS = D * (S_in - S_calc) - qs * X_calc
    dP = (alpha * (miu - Kd) * X_calc + qp * X_calc) - D * P_calc

    #Controlador PI del flujo de refrigerante en la chaqueta
    Error = Tr - T_setpoint
    Fc_control = F0 + Kp * Error + (Kp/Ti) * I
    Fc = np.clip(Fc_control, F_min, F_max) #Caudal de agual de enfriamiento 

    #Anti-windup
    if (Fc_control > F_max and Error > 0) or (Fc_control < F_min and Error < 0):
        dI = 0
    else:
        dI = Error

    #Balance de Energia 
    dTr = (F/V_calc) * (T_feed - Tr) + (rQ / (rho * Cp)) - (UA * (Tr - Tj)) / (rho * V_calc * Cp) # Temperatura en el reactor 
    dTj = (Fc/V_jacket) * (Tj_entrada - Tj) + (UA * (Tr - Tj)) / (rho * V_jacket * Cp) #Temperatura chaqueta 
    return [dX, dS, dP, dTr, dTj, dI, dV]

##Parametros
V0 = 1.5 #Volumen inicial de caldo en el reactor [L]
V_max = 20 #Volumen maximo de operacion del reactor [L]
S_in = 10 #Concentracion de sustrato en la corriente de alimentacion [g/L]
F_feed = 1 #Caudal de alimentacion [L/h]
T_feed = 25 #Temperatura de la corriente de alimentacion [°C]
Kd = 0.0001 #Coeficiente de muerte celular [h^-1]
miu_max = 1.09 #Tasa de crecimiento especifica maxima [h^-1]
qs_max = 4.16 #Tasa de utilizacion de sustrato especifica maxima [g/g*h]
qp_max = 1.863 #Tasa de produccion de lactato especifica maxima [g/g*h]
alpha = 0.017 #Constante asociada al crecimiento en Luedeking-Piret [g/g]


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
Tj_entrada = 10 #Temperatura de entrada a la chaqueta[°C]
UA = 75 * 3600 #[J/h*°C]

#Rendimientos
Yps = 0.72 #Rendimiento de producto [g/g]
Yxs = 0.074 #Rendimiento de biomasa [g/g]
Ypx = Yps/Yxs #Rendimiento de producto [g/g]
Yqs = 3963 #Rendimiento termico [J/g]

#Parametros del controlador PI (flujo de refrigerante en la chaqueta)
T_setpoint = 30 #[°C]
Kp = 9.59 #[L/h*°C]
Ti = 3.66 #[h]
F0 = 1 #[L/h]
F_min = 0 #[L/h]
F_max = 10 #[L/h]

#Condiciones iniciales
X0 = 0.43 #[g/L]
S0 = 33 #[g/L]
P0 = 0 #[g/L]
Tr0 = 30 #[°C]
Tj0 = 25 #[°C]
I0 = 0 #[°C*h]
array_iniciales = np.array([X0, S0, P0, Tr0, Tj0, I0, V0])

#Tiempo de ejecucion
t_start = 0
t_stop = 24
tspan = (t_start, t_stop)
t_array = np.linspace(t_start, t_stop, num=1000)

#Metodo numerico
solucion = solve_ivp(fedbatch_jacket, tspan, array_iniciales, t_eval=t_array, method='LSODA')

Tiempo = solucion.t
Biomasa = solucion.y[0]
Sustrato = solucion.y[1]
Producto = solucion.y[2]
T_reactor = solucion.y[3]
T_jacket = solucion.y[4]
Integral_error = solucion.y[5]
Volumen = solucion.y[6]

Error = T_reactor - T_setpoint
F_valor = F0 + Kp * Error + (Kp/Ti) * Integral_error
F = np.clip(F_valor, F_min, F_max)

#Resumen numerico
t_fin_alimentacion = Tiempo[Volumen >= V_max]
t_fin_alimentacion = t_fin_alimentacion[0] if len(t_fin_alimentacion) > 0 else None

print(f"Volumen final: {Volumen[-1]:.2f} L")
if t_fin_alimentacion is not None:
    print(f"La alimentacion se corta en t = {t_fin_alimentacion:.2f} h (V alcanza V_max)")
else:
    print("La alimentacion nunca alcanza V_max en el tiempo simulado")
print(f"Biomasa final: {Biomasa[-1]:.3f} g/L (maxima: {Biomasa.max():.3f} g/L en t={Tiempo[np.argmax(Biomasa)]:.2f} h)")
print(f"Sustrato final: {Sustrato[-1]:.3f} g/L")
print(f"Producto (lactato) final: {Producto[-1]:.3f} g/L")
print(f"Productividad volumetrica promedio: {Producto[-1]/Tiempo[-1]:.3f} g/L*h")
print(f"Temperatura del reactor: min={T_reactor.min():.2f} °C, max={T_reactor.max():.2f} °C")
print(f"Fraccion de tiempo con flujo de refrigerante saturado (F=F_max): {np.mean(F_valor >= F_max)*100:.1f} %")

#Graficas
f1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 7))
ax1.plot(Tiempo, Biomasa, label='Biomasa', color='red')
ax1.plot(Tiempo, Sustrato, label='Sustrato', color='blue')
ax1.plot(Tiempo, Producto, label='Producto (lactato)', color='green')
ax1.set_xlabel('Tiempo [h]'); ax1.set_ylabel('Concentracion [g/L]')

ax2.plot(Tiempo, T_reactor, label='Temperatura reactor', color='orange')
ax2.plot(Tiempo, T_jacket, label='Temperatura chaqueta', color='purple')
ax2.axhline(T_setpoint, color='darkred', linestyle='--', label='Set point')
ax2.set_xlabel('Tiempo [h]'); ax2.set_ylabel('Temperatura [°C]')

ax3.plot(Tiempo, Error, label='Error', color='gold')
ax3.set_xlabel('Tiempo [h]'); ax3.set_ylabel('Error [°C]')

ax4.plot(Tiempo, F, label='Flujo refrigerante real', color='magenta')
ax4.plot(Tiempo, Volumen, label='Volumen reactor', color='teal')
ax4.set_xlabel('Tiempo [h]')

for i in [ax1, ax2, ax3, ax4]:
    i.grid()
    i.legend()

plt.tight_layout()
plt.show()

# %%
'''Cultivo fedbatch biorreactor 20L HCW - optimizacion'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import shgo

def objective (x):
    def fedbatch_jacket(t, Y):
        X, S, P, Tr, Tj, I, V = Y
    
        X_calc = max(0, X)
        S_calc = max(0, S)
        P_calc = max(0, P)
        V_calc = max(1e-6, V)   
    
        #Tasa especifica de crecimiento
        miu = (miu_max * S_calc * Kix * np.exp(-P_calc/Kpx))/((Ksx + S_calc)*(Kix + S_calc))
        qs = (qs_max * S_calc * Kis * np.exp(-P_calc/Kps))/((Kss + S_calc)*(Kis + S_calc))
        qp = (qp_max * S_calc * Kip * np.exp(-P_calc/Kpp))/((Ksp + S_calc)*(Kip + S_calc))
    
        #Tasa volumetrica de generacion de energia metabolica
        rQ = Yqs * qs * X_calc
    
        #Alimentacion
        if t < 5 :
            F = 0 #Fase Inicial Bacth
            dV = 0
        else:
            F = F_feed  #Inicio fase fed-bacth
            dV = F
            
        D=F/V_calc #Tasa de dilución  
        
        #Balances de Masa 
        dX = (miu - Kd) * X_calc - D * X_calc
        dS = D * (S_in - S_calc) - qs * X_calc
        dP = (alpha * (miu - Kd) * X_calc + qp * X_calc) - D * P_calc
    
        #Controlador PI del flujo de refrigerante en la chaqueta
        Error = Tr - T_setpoint
        Fc_control = F0 + Kp * Error + (Kp/Ti) * I
        Fc = np.clip(Fc_control, F_min, F_max) #Caudal de agual de enfriamiento 
    
        #Anti-windup
        if (Fc_control > F_max and Error > 0) or (Fc_control < F_min and Error < 0):
            dI = 0
        else:
            dI = Error
    
        #Balance de Energia 
        dTr = (F/V_calc) * (T_feed - Tr) + (rQ / (rho * Cp)) - (UA * (Tr - Tj)) / (rho * V_calc * Cp) # Temperatura en el reactor 
        dTj = (Fc/V_jacket) * (Tj_entrada - Tj) + (UA * (Tr - Tj)) / (rho * V_jacket * Cp) #Temperatura chaqueta 
        return [dX, dS, dP, dTr, dTj, dI, dV]
    
    ##Parametros
    V0 = 1.5 #[L] volumen inicial de caldo en el reactor 
    V_max = 20 #[L] volumen maximo de operacion del reactor 
    S_in = x[0] #Concentracion de sustrato en la corriente de alimentacion [g/L]
    F_feed = 1 #Caudal de alimentacion [L/h]
    T_feed = 25 #Temperatura de la corriente de alimentacion [°C]
    Kd = 0.0001 #Coeficiente de muerte celular [h^-1]
    miu_max = 1.09 #Tasa de crecimiento especifica maxima [h^-1]
    qs_max = 4.16 #Tasa de utilizacion de sustrato especifica maxima [g/g*h]
    qp_max = 1.863 #Tasa de produccion de lactato especifica maxima [g/g*h]
    alpha = 0.017 #Constante asociada al crecimiento en Luedeking-Piret [g/g]
    
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
    
    #Rendimientos
    Yps = 0.72 #Rendimiento de producto [g/g]
    Yxs = 0.074 #Rendimiento de biomasa [g/g]
    Ypx = Yps/Yxs #Rendimiento de producto [g/g]
    Yqs = 3963 #Rendimiento termico [J/g]
    
    #Parametros del controlador PI (flujo de refrigerante en la chaqueta)
    T_setpoint = 30 #[°C]
    Kp = 9.59 #[L/h*°C]
    Ti = 3.66 #[h]
    F0 = 1 #[L/h]
    F_min = 0 #[L/h]
    F_max = 10 #[L/h]
    
    #Condiciones iniciales
    X0 = 0.43 #[g/L]
    S0 = 33 #[g/L]
    P0 = 0 #[g/L]
    Tr0 = 30 #[°C]
    Tj0 = 25 #[°C]
    I0 = 0 #[°C*h]
    array_iniciales = np.array([X0, S0, P0, Tr0, Tj0, I0, V0])
    
    #Tiempo de ejecucion
    t_start = 0
    t_stop = 24
    tspan = (t_start, t_stop)
    t_array = np.linspace(t_start, t_stop, num=1000)
    
    #Metodo numerico
    solucion = solve_ivp(fedbatch_jacket, tspan, array_iniciales, t_eval=t_array, method='LSODA')
    
    Tiempo = solucion.t
    Biomasa = solucion.y[0]
    Sustrato = solucion.y[1]
    Producto = solucion.y[2]
    T_reactor = solucion.y[3]
    T_jacket = solucion.y[4]
    Integral_error = solucion.y[5]
    Volumen = solucion.y[6]
    
    # Encontrar el índice donde el volumen es máximo (Punto de parada del Fed-batch)
    i_opt = np.argmax(Volumen) 
    t_final = Tiempo[i_opt]
    P_final = Producto[i_opt]
    V_final = Volumen[i_opt]
    # Productividad Global (Qp) en g/(L·h) 
    # Nota: Se usa la concentración final porque el volumen es el del sistema completo
    if t_final > 0:
        PROD_global = P_final / t_final
    else:
        PROD_global = 0
    return -PROD_global
print("Optimizando Sfeed...")
result = shgo(objective, bounds=[(0, 500)])
Sfeed_opt = result.x[0]
prod_max   = -result.fun

print(f"\n{'='*45}")
print(f"  Sfeed óptimo  : {Sfeed_opt:.4f} g/L")
print(f"  Productividad máxima: {prod_max:.6f} g/(L·h)")
print(f"{'='*45}\n")
# %%
'''Cultivo fedbatch biorreactor 20L HCW - sintonizacion del controlador'''
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Parámetros del sistema 
V0, V_max = 1.5, 20
S_in, F_feed, T_feed = 57.99, 1.0, 25.0
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
    X, S, P, Tr, Tj, I, V = Y
    X_c, S_c, P_c, V_c = max(0, X), max(0, S), max(0, P), max(1e-6, V)

    # Cinética
    miu = (miu_max * S_c * Kix) / ((Ksx + S_c) * (Kix + S_c)) * np.exp(-P_c/Kpx)
    qs = (qs_max * S_c * Kis) / ((Kss + S_c) * (Kis + S_c)) * np.exp(-P_c/Kps)
    qp = (qp_max * S_c * Kip) / ((Ksp + S_c) * (Kip + S_c)) * np.exp(-P_c/Kpp)
    rQ = Yqs * qs * X_c # Generación de calor 

    # Operación Fed-batch
    F = F_feed if (t >= 5 and V_c < V_max) else 0
    D = F / V_c

    # Controlador PI (Variable manipulada: Fc)
    Error = Tr - T_setpoint
    Fc_control = F0 + Kp * Error + (Kp/Ti) * I
    Fc = np.clip(Fc_control, F_min, F_max)

    # Balances de masa, energía e integral
    dX = (miu - Kd) * X_c - D * X_c
    dS = D * (S_in - S_c) - qs * X_c
    dP = (alpha * (miu - Kd) * X_c + qp * X_c) - D * P_c
    dTr = (F/V_c)*(T_feed - Tr) + (rQ/(rho*Cp)) - (UA*(Tr-Tj))/(rho*V_c*Cp) # [6]
    dTj = (Fc/V_jacket)*(Tj_entrada - Tj) + (UA*(Tr-Tj))/(rho*V_jacket*Cp) # [3]
    dI = 0 if (Fc_control > F_max and Error > 0) or (Fc_control < F_min and Error < 0) else Error
    dV = F

    return [dX, dS, dP, dTr, dTj, dI, dV]

# Función de costo: Minimización del IAE (Error Integral Absoluto) 
def objetivo_iae(parametros):

    Kp_opt, Ti_opt = parametros

    if Kp_opt <= 0 or Ti_opt <= 0:
        return 1e10

    y0 = [0.43, 33, 0, 30.5, 25, 0, 1.5]

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

# %%
'''Cultivo fedbatch biorreactor 20L HCW con sfeed optimizado el controlador sintonizado '''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def fedbatch_jacket(t, Y):
    X, S, P, Tr, Tj, I, V = Y

    X_calc = max(0, X)
    S_calc = max(0, S)
    P_calc = max(0, P)
    V_calc = max(1e-6, V)   

    #Tasa especifica de crecimiento
    miu = (miu_max * S_calc * Kix * np.exp(-P_calc/Kpx))/((Ksx + S_calc)*(Kix + S_calc))
    qs = (qs_max * S_calc * Kis * np.exp(-P_calc/Kps))/((Kss + S_calc)*(Kis + S_calc))
    qp = (qp_max * S_calc * Kip * np.exp(-P_calc/Kpp))/((Ksp + S_calc)*(Kip + S_calc))

    #Tasa volumetrica de generacion de energia metabolica
    rQ = Yqs * qs * X_calc

    #Alimentacion
    if t < 5 :
        F=0 #Fase Inicial Bacth
        dV=0
    else:
        F=F_feed if V_calc < 20 else 0 #Inicio fase fed-bacth
        dV=F
        
    D=F/V_calc #Tasa de dilución 
    
    #Balances de Masa 
    dX = (miu - Kd) * X_calc - D * X_calc
    dS = D * (S_in - S_calc) - qs * X_calc
    dP = (alpha * (miu - Kd) * X_calc + qp * X_calc) - D * P_calc

    #Controlador PI del flujo de refrigerante en la chaqueta
    Error = Tr - T_setpoint
    Fc_control = F0 + Kp * Error + (Kp/Ti) * I
    Fc = np.clip(Fc_control, F_min, F_max) #Caudal de agual de enfriamiento 

    #Anti-windup
    if (Fc_control > F_max and Error > 0) or (Fc_control < F_min and Error < 0):
        dI = 0
    else:
        dI =Error

    #Balance de Energia 
    dTr = (F/V_calc) * (T_feed - Tr) + (rQ / (rho * Cp)) - (UA * (Tr - Tj)) / (rho * V_calc * Cp) # Temperatura en el reactor 
    dTj = (Fc/V_jacket) * (Tj_entrada - Tj) + (UA * (Tr - Tj)) / (rho * V_jacket * Cp) #Temperatura chaqueta 
    return [dX, dS, dP, dTr, dTj, dI, dV]

##Parametros
V0 = 1.5          #[L] volumen inicial de caldo en el reactor 
V_max = 20        #[L] volumen maximo de operacion del reactor 
S_in = 57.99 #[g/L]
F_feed = 1 #[L/h]
T_feed = 25 #[°C] temperatura de la corriente de alimentacion 
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
Tj_entrada = 25 #[°C]
UA = 75 * 3600 #[J/h*°C]

#Rendimientos
Yps = 0.72 #Rendimiento de producto [g/g]
Yxs = 0.074 #Rendimiento de biomasa [g/g]
Ypx = Yps/Yxs #Rendimiento de producto [g/g]
Yqs = 3963 #Rendimiento termico [J/g]

#Parametros del controlador PI (flujo de refrigerante en la chaqueta)
T_setpoint = 30 #[°C]
Kp = 47.3736 #[L/h*°C]
Ti = 0.1000 #[h]
F0 = 1 #[L/h]
F_min = 0 #[L/h]
F_max = 10 #[L/h]

#Condiciones iniciales
X0 = 0.43 #[g/L]
S0 = 33 #[g/L]
P0 = 0 #[g/L]
Tr0 = 30 #[°C]
Tj0 = 29 #[°C]
I0 = 0 #[°C*h]
array_iniciales = np.array([X0, S0, P0, Tr0, Tj0, I0, V0])

#Tiempo de ejecucion
t_start = 0
t_stop = 60
tspan = (t_start, t_stop)
t_array = np.linspace(t_start, t_stop, num=1000)

#Metodo numerico
solucion = solve_ivp(fedbatch_jacket, tspan, array_iniciales, t_eval=t_array, method='LSODA')

Tiempo = solucion.t
Biomasa = solucion.y[0]
Sustrato = solucion.y[1]
Producto = solucion.y[2]
T_reactor = solucion.y[3]
T_jacket = solucion.y[4]
Integral_error = solucion.y[5]
Volumen = solucion.y[6]

Error = T_reactor - T_setpoint
F_valor = F0 + Kp * Error + (Kp/Ti) * Integral_error
F = np.clip(F_valor, F_min, F_max)

#--- Resumen numerico util para el informe ---
t_fin_alimentacion = Tiempo[Volumen >= V_max]
t_fin_alimentacion = t_fin_alimentacion[0] if len(t_fin_alimentacion) > 0 else None

# CÁLCULOS ADICIONALES

conversion = (S0 + (S_in * (Volumen[-1] - V0)) / Volumen[-1] - Sustrato[-1]) / \
             (S0 + (S_in * (Volumen[-1] - V0)) / Volumen[-1])

producto_total = Producto[-1] * Volumen[-1]
biomasa_total = Biomasa[-1] * Volumen[-1]
substrato_alimentado = (Volumen[-1] - V0) * S_in

Yps_real = producto_total / substrato_alimentado
Yxs_real = biomasa_total / substrato_alimentado
Pv = producto_total / (Volumen[-1] * Tiempo[-1])

Q = np.trapezoid(UA * (T_reactor - T_jacket), Tiempo)

# RESULTADOS DE OPERACIÓN

print("\n" + "="*60)
print("          RESULTADOS DE LA SIMULACIÓN")
print("="*60)

print("\n--- Operación del reactor ---")
print(f"Volumen inicial                 : {V0:.2f} L")
print(f"Volumen final                   : {Volumen[-1]:.2f} L")

if t_fin_alimentacion is not None:
    print(f"Fin de la alimentación          : {t_fin_alimentacion:.2f} h")
else:
    print("Fin de la alimentación          : No alcanza Vmax")

print(f"Tiempo total de simulación      : {Tiempo[-1]:.2f} h")

# BIOMASA, SUSTRATO Y PRODUCTO

print("\n--- Resultados biológicos ---")
print(f"Biomasa final                   : {Biomasa[-1]:.3f} g/L")
print(f"Biomasa máxima                  : {Biomasa.max():.3f} g/L")
print(f"Tiempo biomasa máxima           : {Tiempo[np.argmax(Biomasa)]:.2f} h")

print(f"Sustrato final                  : {Sustrato[-1]:.3f} g/L")
print(f"Conversión de sustrato          : {conversion*100:.2f} %")

print(f"Producto final                  : {Producto[-1]:.3f} g/L")
print(f"Producto total                  : {producto_total:.2f} g")

print(f"Biomasa total                   : {biomasa_total:.2f} g")
print(f"Sustrato alimentado             : {substrato_alimentado:.2f} g")

# RENDIMIENTOS

print("\n--- Indicadores de desempeño ---")
print(f"Productividad volumétrica       : {Pv:.3f} g/L·h")
print(f"Rendimiento Yps real            : {Yps_real:.3f} g/g")
print(f"Rendimiento Yxs real            : {Yxs_real:.3f} g/g")


# CONTROL DE TEMPERATURA

print("\n--- Desempeño térmico ---")
print(f"Temperatura mínima reactor      : {T_reactor.min():.2f} °C")
print(f"Temperatura máxima reactor      : {T_reactor.max():.2f} °C")
print(f"Energía removida por chaqueta   : {Q:.2f} J")

print("\n--- Parámetros del controlador PI ---")
print(f"Kp                             : {Kp:.2f} L/h·°C")
print(f"Ti                             : {Ti:.2f} h")
print(f"Set point                      : {T_setpoint:.2f} °C")
print(f"Caudal nominal                 : {F0:.2f} L/h")
print(f"Caudal máximo                  : {F_max:.2f} L/h")

#Grafica
f1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 7))
ax1.plot(Tiempo, Biomasa, label='Biomasa', color='red')
ax1.plot(Tiempo, Sustrato, label='Sustrato', color='blue')
ax1.plot(Tiempo, Producto, label='Producto (lactato)', color='green')
ax1.set_xlabel('Tiempo [h]'); ax1.set_ylabel('Concentracion [g/L]')

ax2.plot(Tiempo, T_reactor, label='Temperatura reactor', color='orange')
ax2.plot(Tiempo, T_jacket, label='Temperatura chaqueta', color='purple')
ax2.axhline(T_setpoint, color='darkred', linestyle='--', label='Set point')
ax2.set_xlabel('Tiempo [h]'); ax2.set_ylabel('Temperatura [°C]')

ax3.plot(Tiempo, Error, label='Error', color='gold')
ax3.set_xlabel('Tiempo [h]'); ax3.set_ylabel('Error [°C]')

ax4.plot(Tiempo, F, label='Flujo refrigerante real', color='magenta')
ax4.plot(Tiempo, Volumen, label='Volumen reactor', color='teal')
ax4.set_xlabel('Tiempo [h]')

for i in [ax1, ax2, ax3, ax4]:
    i.grid()
    i.legend()

plt.tight_layout()
plt.show()