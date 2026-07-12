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
    miu = (miu_max * S_calc * Kix * np.exp(-P_calc/Kpx))/(Ksx + S_calc)*(Kix + S_calc)
    qs = (qs_max * S_calc * Kis * np.exp(-P_calc/Kps))/(Kss + S_calc)*(Kis + S_calc)
    qp = (qp_max * S_calc * Kip * np.exp(-P_calc/Kpp))/(Ksp + S_calc)*(Kip + S_calc)
    
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
    dX = miu * X_calc #Porque se tienen membrana ideal
    dS = (F/V_reactor) * (S_in - S_calc) - qs * X
    dP = alpha * dX + qp * X
    dTr = (F/V_reactor)
    dTj = (F/V_jacket)
    return [dX, dS, dP, dTr, dTj, dI]

##Parametros
V_reactor = 20 #[L]
S_in = 10 #[g/L]
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
rho = 1000 #[g/L]
Cp = 1 #[kcal/kg*°C]

#Chaqueta de enfriamiento
V_jacket = 2 #[L]
Tj_entrada = 10
UA = 30 #[kcal/h*°C]

#Parametros del controlador PI
T_setpoint = 30 #[°C]
Kp = 5 #[L/h*°C]
Ti = 10 #[h]
F0 = 0.2 #[L/h]
F_min = 0 #[L/h]
F_max = 10 #[L/h]

#Rendimientos (Ypx, Yps, Yxs)

