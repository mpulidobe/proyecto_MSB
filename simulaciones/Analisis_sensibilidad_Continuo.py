'''Analisis de sensibilidad reactor continuo 20L - Parametro miu_max'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.widgets import Slider, Button

def continuo_jacket(t, Y, miu_max, qp_max, Kix, Kip):
    X, S, P, Tr, Tj, I = Y

    X_calc = max(0, X)
    S_calc = max(0, S)
    P_calc = max(0, P)

    #Tasa especifica de crecimiento
    miu = (miu_max * S_calc * Kix * np.exp(-P_calc/Kpx))/((Ksx + S_calc)*(Kix + S_calc))
    qs = (qs_max * S_calc * Kis * np.exp(-P_calc/Kps))/((Kss + S_calc)*(Kis + S_calc))
    qp = (qp_max * S_calc * Kip * np.exp(-P_calc/Kpp))/((Ksp + S_calc)*(Kip + S_calc))

    rQ = Yqs * qs * X_calc

    D = F_feed / V_reactor #Tasa de dilucion [h^-1]

    #Controlador PI del flujo de refrigerante en la chaqueta
    Error = Tr - T_setpoint
    Fc_control = F0 + Kp * Error + (Kp/Ti) * I
    Fc = np.clip(Fc_control, F_min, F_max)

    #Anti-windup
    if (Fc_control > F_max and Error > 0) or (Fc_control < F_min and Error < 0):
        dI = 0
    else:
        dI = Error

    #Balances de masa (membrana ideal: retencion biomasa = 1, retencion producto = 0)
    dX = (miu - Kd) * X_calc
    dS = D * (S_in - S_calc) - qs * X_calc
    dP = alpha * (miu - Kd) * X_calc + qp * X_calc - D * P_calc

    #Balance de energia
    dTr = D * (T_feed - Tr) + (rQ / (rho * Cp)) - (UA * (Tr - Tj)) / (rho * V_reactor * Cp)
    dTj = (Fc/V_jacket) * (Tj_entrada - Tj) + (UA * (Tr - Tj)) / (rho * V_jacket * Cp)
    return [dX, dS, dP, dTr, dTj, dI]

##Parametros
V_reactor = 20 #[L] volumen maximo de operacion del reactor 
S_in = 64.4440 #[g/L]
D_op = 1.0800 #[h^-1] optimizado
F_feed = D_op * V_reactor #[L/h]
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
F_max = 20 #[L/h]

#Condiciones iniciales
X0 = 0.43 #[g/L]
S0 = 33 #[g/L]
P0 = 0 #[g/L]
Tr0 = 30 #[°C]
Tj0 = 29 #[°C]
I0 = 0 #[°C*h]
array_iniciales = np.array([X0, S0, P0, Tr0, Tj0, I0])

#Tiempo de ejecucion
t_start = 0
t_stop = 100
tspan = (t_start, t_stop)
t_array = np.linspace(t_start, t_stop, num=5000)

#Metodo numerico
solucion = solve_ivp(continuo_jacket, tspan, array_iniciales, t_eval=t_array, args = (miu_max, qp_max, Kix, Kip), method='LSODA')

Tiempo = solucion.t
Biomasa = solucion.y[0]
Sustrato = solucion.y[1]
Producto = solucion.y[2]
T_reactor = solucion.y[3]
T_jacket = solucion.y[4]
Integral_error = solucion.y[5]
Volumen = np.ones_like(Tiempo) * V_reactor

Error = T_reactor - T_setpoint
F_valor = F0 + Kp * Error + (Kp/Ti) * Integral_error
F = np.clip(F_valor, F_min, F_max)

#Fijar los valores de los parametros en _base
miu_max_base = miu_max
qp_max_base = qp_max
Kix_base = Kix
Kip_base = Kip

# Graficas
f1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 8))
line1_miu, = ax1.plot(Tiempo, Biomasa, label='Biomasa', color='red')
line2_miu, = ax1.plot(Tiempo, Sustrato, label='Sustrato', color='blue')
line3_miu, = ax1.plot(Tiempo, Producto, label='Producto (lactato)', color='green')
ax1.set_title(r'Sensibilidad local del parámetro $mu_{max}$', fontsize=12)
ax1.set_xlabel('Tiempo [h]'); ax1.set_ylabel('Concentracion [g/L]')

line1_qp, = ax2.plot(Tiempo, Biomasa, label='Biomasa', color='red')
line2_qp, = ax2.plot(Tiempo, Sustrato, label='Sustrato', color='blue')
line3_qp, = ax2.plot(Tiempo, Producto, label='Producto (lactato)', color='green')
ax2.set_title(r'Sensibilidad local del parámetro $qp_{max}$', fontsize=12)
ax2.set_xlabel('Tiempo [h]'); ax1.set_ylabel('Concentracion [g/L]')

line1_Kix, = ax3.plot(Tiempo, Biomasa, label='Biomasa', color='red')
line2_Kix, = ax3.plot(Tiempo, Sustrato, label='Sustrato', color='blue')
line3_Kix, = ax3.plot(Tiempo, Producto, label='Producto (lactato)', color='green')
ax3.set_title(r'Sensibilidad local del parámetro $K_{ix}$', fontsize=12)
ax3.set_xlabel('Tiempo [h]'); ax1.set_ylabel('Concentracion [g/L]')

line1_Kip, = ax4.plot(Tiempo, Biomasa, label='Biomasa', color='red')
line2_Kip, = ax4.plot(Tiempo, Sustrato, label='Sustrato', color='blue')
line3_Kip, = ax4.plot(Tiempo, Producto, label='Producto (lactato)', color='green')
ax4.set_title(r'Sensibilidad local del parámetro $K_{ip}$', fontsize=12)
ax4.set_xlabel('Tiempo [h]'); ax1.set_ylabel('Concentracion [g/L]')

for i in [ax1, ax2, ax3, ax4]:
    i.grid()
    i.legend()

plt.subplots_adjust(left=0.32, right=0.98, bottom=0.12, top=0.95, wspace=0.30, hspace=0.30)

axcolor = "lightgoldenrodyellow"

ax_miu  = plt.axes([0.05,0.78,0.20,0.03], facecolor=axcolor)
ax_qp  = plt.axes([0.05,0.66,0.20,0.03], facecolor=axcolor)
ax_Kix = plt.axes([0.05,0.54,0.20,0.03], facecolor=axcolor)
ax_Kip = plt.axes([0.05,0.42,0.20,0.03], facecolor=axcolor)

s_miu = Slider(
    ax=ax_miu,
    label=r'$\mu_{max}$',
    valmin=0.8*miu_max_base,
    valmax=1.2*miu_max_base,
    valinit=miu_max_base,
    valfmt='%1.3f')

s_qp = Slider(
    ax=ax_qp,
    label=r'$q_{p,max}$',
    valmin=0.8*qp_max_base,
    valmax=1.2*qp_max_base,
    valinit=qp_max_base,
    valfmt='%1.3f')

s_Kix = Slider(
    ax=ax_Kix,
    label=r'$K_{ix}$',
    valmin=0.8*Kix_base,
    valmax=1.2*Kix_base,
    valinit=Kix_base,
    valfmt='%1.3f')

s_Kip = Slider(
    ax=ax_Kip,
    label=r'$K_{ip}$',
    valmin=0.8*Kip_base,
    valmax=1.2*Kip_base,
    valinit=Kip_base,
    valfmt='%1.3f')

def update_miu(val):

    miu = s_miu.val

    sol = solve_ivp(
        continuo_jacket,
        tspan,
        array_iniciales,
        t_eval=t_array,
        args=(miu, qp_max_base, Kix_base, Kip_base),
        method="LSODA"
    )

    line1_miu.set_ydata(sol.y[0])
    line2_miu.set_ydata(sol.y[1])
    line3_miu.set_ydata(sol.y[2])

    ax1.relim()
    ax1.autoscale_view()
    f1.canvas.draw_idle()
    
def update_qp(val):

    qp = s_qp.val

    sol = solve_ivp(
        continuo_jacket,
        tspan,
        array_iniciales,
        t_eval=t_array,
        args=(miu_max_base, qp, Kix_base, Kip_base),
        method="LSODA"
    )

    line1_qp.set_ydata(sol.y[0])
    line2_qp.set_ydata(sol.y[1])
    line3_qp.set_ydata(sol.y[2])

    ax1.relim()
    ax1.autoscale_view()
    f1.canvas.draw_idle()

def update_Kix(val):

    Kix = s_Kix.val

    sol = solve_ivp(
        continuo_jacket,
        tspan,
        array_iniciales,
        t_eval=t_array,
        args=(miu_max_base, qp_max_base, Kix, Kip_base),
        method="LSODA"
    )

    line1_Kix.set_ydata(sol.y[0])
    line2_Kix.set_ydata(sol.y[1])
    line3_Kix.set_ydata(sol.y[2])

    ax1.relim()
    ax1.autoscale_view()
    f1.canvas.draw_idle()

def update_Kip(val):

    Kip = s_Kip.val

    sol = solve_ivp(
        continuo_jacket,
        tspan,
        array_iniciales,
        t_eval=t_array,
        args=(miu_max_base, qp_max_base, Kix_base, Kip),
        method="LSODA"
    )

    line1_Kip.set_ydata(sol.y[0])
    line2_Kip.set_ydata(sol.y[1])
    line3_Kip.set_ydata(sol.y[2])

    ax1.relim()
    ax1.autoscale_view()
    f1.canvas.draw_idle()
    
s_miu.on_changed(update_miu)
s_qp.on_changed(update_qp)
s_Kix.on_changed(update_Kix)
s_Kip.on_changed(update_Kip)

resetax = plt.axes([0.08,0.28,0.12,0.05])

button = Button(resetax, "Reset")

def reset(event):
    s_miu.reset()
    s_qp.reset()
    s_Kix.reset()
    s_Kip.reset()

button.on_clicked(reset)

plt.show()
