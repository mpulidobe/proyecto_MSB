'''Analisis de sensibilidad reactor continuo 20L - Parametro miu_max'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.widgets import Slider, Button

try:
    from scipy.integrate import trapezoid
except ImportError:
    from scipy.integrate import trapz as trapezoid

def continuo_jacket(t, Y, miu_max):
    X, S, P, Tr, Tj, I, P_acumulado = Y

    X_calc = max(0, X)
    S_calc = max(0, S)
    P_calc = max(0, P)  

    #Tasa especifica de crecimiento
    miu = (miu_max * S_calc * Kix * np.exp(-P_calc/Kpx))/((Ksx + S_calc)*(Kix + S_calc))
    qs = (qs_max * S_calc * Kis * np.exp(-P_calc/Kps))/((Kss + S_calc)*(Kis + S_calc))
    qp = (qp_max * S_calc * Kip * np.exp(-P_calc/Kpp))/((Ksp + S_calc)*(Kip + S_calc))

    #Tasa volumetrica de generacion de energia metabolica
    rQ = Yqs * qs * X_calc

    #Cambio de modo de operación
    if t < 5:
        D = 0          # Batch
    else:
        D = F_feed / V_reactor
    
    # Controlador PI (Variable manipulada: Fc)
    Error = Tr - T_setpoint
    Fc_control = F0 + Kp * Error + (Kp/Ti) * I
    Fc = np.clip(Fc_control, F_min, F_max)
    
    #Ecuaciones diferenciales de las variables de estado
    dX = (miu - Kd)* X_calc #Porque se tienen membrana ideal
    dS = D * (S_in - S_calc) - qs * X_calc
    dP = alpha * dX + qp * X_calc - D * P_calc
    dTr = D * (T_feed - Tr) + (rQ / (rho * Cp)) - (UA * (Tr - Tj)) / (rho * V_reactor * Cp)
    dTj = (Fc/V_jacket) * (Tj_entrada - Tj) + (UA * (Tr - Tj)) / (rho * V_jacket * Cp)
    dI = 0 if (Fc_control > F_max and Error > 0) or (Fc_control < F_min and Error < 0) else Error
    dP_acumulado = D * V_reactor * P_calc
    return [dX, dS, dP, dTr, dTj, dI, dP_acumulado]

##Parametros
V_reactor = 20 #[L] volumen maximo de operacion del reactor 
S_in = 80.4029 #[g/L]
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
Tj_entrada = 10 #[°C]
UA = 75 * 3600 #[J/h*°C]

#Rendimientos
Yps = 0.72 #Rendimiento de producto [g/g]
Yxs = 0.074 #Rendimiento de biomasa [g/g]
Ypx = Yps/Yxs #Rendimiento de producto [g/g]
Yqs = 3963 #Rendimiento termico [J/g]

#Parametros del controlador PI (flujo de refrigerante en la chaqueta)
T_setpoint = 30 #[°C]
Kp = 50 #[L/h*°C]
Ti = 0.1000 #[h]
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
t_stop = 60
tspan = (t_start, t_stop)
t_array = np.linspace(t_start, t_stop, num=1000)

#Metodo numerico
solucion = solve_ivp(continuo_jacket, tspan, array_iniciales, t_eval=t_array, args = (miu_max, ), method='LSODA')

Tiempo = solucion.t
Biomasa = solucion.y[0]
Sustrato = solucion.y[1]
Producto = solucion.y[2]
T_reactor = solucion.y[3]
T_jacket = solucion.y[4]
Integral_error = solucion.y[5]
Producto_acumulado = solucion.y[6]
Volumen = np.ones_like(Tiempo) * V_reactor

Error = T_reactor - T_setpoint
F_valor = F0 + Kp * Error + (Kp/Ti) * Integral_error
F = np.clip(F_valor, F_min, F_max)
miu_max_base = miu_max

# Graficas
f1, ax1 = plt.subplots(figsize=(8,5))
line1, = ax1.plot(Tiempo, Biomasa, label='Biomasa', color='red')
line2, = ax1.plot(Tiempo, Sustrato, label='Sustrato', color='blue')
line3, = ax1.plot(Tiempo, Producto, label='Producto (lactato)', color='green')
ax1.set_xlabel('Tiempo [h]'); ax1.set_ylabel('Concentracion [g/L]')
ax1.grid()
ax1.legend()

plt.subplots_adjust(left=0.28)

axcolor = "lightgoldenrodyellow"

ax_mu = plt.axes([0.05, 0.55, 0.16, 0.03],
                 facecolor=axcolor)

s_mu = Slider(
    ax=ax_mu,
    label=r'$\mu_{max}$',
    valmin=0.8*miu_max_base,
    valmax=1.2*miu_max_base,
    valinit=miu_max_base,
    valfmt='%1.3f'
)

def update(val):

    mu = s_mu.val

    sol = solve_ivp(
        continuo_jacket,
        tspan,
        array_iniciales,
        t_eval=t_array,
        args=(mu,),
        method='LSODA'
    )

    line1.set_ydata(sol.y[0])
    line2.set_ydata(sol.y[1])
    line3.set_ydata(sol.y[2])

    ax1.relim()
    ax1.autoscale_view()

    f1.canvas.draw_idle()

s_mu.on_changed(update)

resetax = plt.axes([0.08,0.42,0.12,0.05])

button = Button(resetax,
                'Reset',
                color='lightcoral')

def reset(event):
    s_mu.reset()

button.on_clicked(reset)
plt.show()
