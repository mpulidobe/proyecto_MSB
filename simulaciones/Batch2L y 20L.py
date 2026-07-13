import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def modelo_batch(t, Y, p):
    X, S, P = Y

    μ = (p['mu_max'] * S * p['Kix']) / ((p['Ksx'] + S) * (p['Kix'] + S)) * np.exp(-P / p['Kpx'])

    dX = (μ - p['Kd']) * X
    dS = -p['qs_max'] * (S * p['Kis']) / ((p['Kss'] + S) * (p['Kis'] + S)) * np.exp(-P / p['Kps']) * X
    dP = p['alpha'] * dX + p['qp_max'] * (S * p['Kip']) / ((p['Ksp'] + S) * (p['Kip'] + S)) * np.exp(-P / p['Kpp']) * X

    return dX, dS, dP

#Parametros por escenario (unidades en g/L)
Escenarios = {
    '2L_M17': {
        'mu_max': 0.34, 'Kix': 114.06, 'Ksx': 0.023, 'Kpx': 18.11, 'Kd': 0.013,
        'qs_max': 4.92, 'Kis': 290.13, 'Kss': 2.73, 'Kps': 10.44,
        'qp_max': 1.027, 'Kip': 135.89, 'Ksp': 0.025, 'Kpp': 5.31, 'alpha': 0.052,
        'X0': 0.5, 'S0': 18.2, 'P0': 0.5# Condiciones iniciales
    },
    '2L_HCW': {
        'mu_max': 0.99, 'Kix': 399.75, 'Ksx': 0.0023, 'Kpx': 3.77, 'Kd': 0.0001,
        'qs_max': 4.99, 'Kis': 399.99, 'Kss': 0.01, 'Kps': 11.16,
        'qp_max': 2.04, 'Kip': 100.0, 'Ksp': 0.01, 'Kpp': 44.99, 'alpha': 0.01,
        'X0': 0.43, 'S0': 45.0, 'P0': 0.0# Condiciones iniciales
    },
    '20L_HCW': {
        'mu_max': 1.09, 'Kix': 394.20, 'Ksx': 4.229, 'Kpx': 5.001, 'Kd': 0.0001,
        'qs_max': 4.16, 'Kis': 143.391, 'Kss': 0.15, 'Kps': 20.07,
        'qp_max': 1.863, 'Kip': 373.89, 'Ksp': 0.065, 'Kpp': 42.83, 'alpha': 0.017,
        'X0': 0.43, 'S0': 33.0, 'P0': 0.0# Condiciones iniciales
    }
}

#Solver
t_span = (0, 24)
t_array = np.linspace(0, 24, 200) 

sol = {}
for k, v in Escenarios.items():
    Y_init = np.array([v['X0'], v['S0'], v['P0']])
    sol[k] = solve_ivp(modelo_batch, t_span, Y_init, t_eval=t_array, args=(v,))

#Graficas
def plot_escala(nombre, keys, colors):
    f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10))
    axes = (ax1, ax2, ax3)
    titulos = ['Biomasa (g/L)', 'Sustrato (g/L)', 'Acido Lactico (g/L)']

    for k, color in zip(keys, colors):
        cx = sol[k].y[0]
        cs = sol[k].y[1]
        cp = sol[k].y[2]
        line1, = ax1.plot(t_array, cx, label=k.replace('_', ' '), color=color)
        line2, = ax2.plot(t_array, cs, label=k.replace('_', ' '), color=color)
        line3, = ax3.plot(t_array, cp, label=k.replace('_', ' '), color=color)

    for i, ax in enumerate(axes):
        ax.set_title(titulos[i])
        ax.set_xlabel('Tiempo (h)')
        ax.set_ylabel(f'Concentración ({unidades[i]})')
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

    for ax in axes:
        ax.legend()

    for ax in axes:
        ax.grid()

    f.suptitle(f'Escala de Biorreactor: {nombre}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()

plot_escala('2 Litros', ['2L_M17', '2L_HCW'], ['red', 'blue'])
plot_escala('20 Litros', ['20L_HCW'], ['green'])
