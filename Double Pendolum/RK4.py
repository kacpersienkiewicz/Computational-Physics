import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def polar_to_rect(Theta_1, Theta_2):
    x_1 = L1 * np.cos(Theta_1)
    y_1 = -L1 * np.sin(Theta_1)
    x_2 = x_1 + L2 * np.cos(Theta_2)
    y_2 = y_1 - L2 * np.sin(Theta_2)
    return x_1, y_1, x_2, y_2

def RK4_Step(variables: np.array):
    Theta_1, Theta_2, Theta_dot_1, Theta_dot_2 = variables
    delta_Theta = Theta_1 - Theta_2

    Theta_double_dot_1 = (-1 * G * (2 + M_sum) * np.sin(Theta_1) - M2 * G * np.sin(Theta_1 - 2 * Theta_2) - 2 * np.sin(delta_Theta) * M2 * (Theta_dot_2**2 * L2 + Theta_dot_1**2 * L1 * np.cos(delta_Theta)))  \
    / (L1 * (2 * M_sum - M2 * np.cos(2 * (delta_Theta))))

    Theta_double_dot_2 = (2 * np.sin(delta_Theta) * (Theta_dot_1**2 * L1 * (M_sum) + G * (M_sum) + np.cos(Theta_1) + Theta_dot_2**2 * L2 * M2 * np.cos(delta_Theta))) \
    / (L2 * (2 * M_sum - M2 * np.cos(2 * (delta_Theta))))

    return np.array([Theta_dot_1, Theta_dot_2, Theta_double_dot_1, Theta_double_dot_2])

def RK4(t, h, variables: np.array):

    k1 = RK4_Step(variables)
    k2 = RK4_Step(variables + k1 * h / 2)
    k3 = RK4_Step(variables + k2 * h / 2)
    k4 = RK4_Step(variables + k3 * h /2)

    return t + h, variables + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

def calc_Energy(variables: np.array):
    Theta_1, Theta_2, Theta_dot_1, Theta_dot_2 = variables
    delta_Theta = Theta_1 - Theta_2
    kinetic =0.5 + M_sum + L1 ** 2 * Theta_dot_1 ** 2 + 0.5 * M2 * L2 ** 2 * Theta_dot_2 ** 2 + M2 * L1 * L2 * Theta_dot_1 * Theta_dot_2 * np.cos(delta_Theta)
    potential = -1 * M_sum * L1 * np.cos(Theta_1) - M2 * G * L2 * np.cos(Theta_2)
    return kinetic, potential, kinetic + potential


def update_frame(frame):
    animated_string_1.set_data([0,Cartesian_corrdinates[frame][0]], [0,Cartesian_corrdinates[frame][1]])

    animated_bob_1.set_data([Cartesian_corrdinates[frame][0]], [Cartesian_corrdinates[frame][1]])

    animated_string_2.set_data([Cartesian_corrdinates[frame][0],Cartesian_corrdinates[frame][2]], [Cartesian_corrdinates[frame][1],Cartesian_corrdinates[frame][3]])

    animated_bob_2.set_data([Cartesian_corrdinates[frame][2]], [Cartesian_corrdinates[frame][3]])

    if frame < history_len:
        trace_bob_2.set_data([cc[2] for cc in Cartesian_corrdinates[:frame]], [cc[3] for cc in Cartesian_corrdinates[:frame]])
    else:
        trace_bob_2.set_data([cc[2] for cc in Cartesian_corrdinates[frame - history_len:frame]], [cc[3] for cc in Cartesian_corrdinates[frame - history_len:frame]])

    kinetic_plot.set_data(times_list[:frame], [e[0] for e in Energy_totals[:frame]])
    potential_plot.set_data(times_list[:frame], [e[1] for e in Energy_totals[:frame]])
    total_energy_plot.set_data(times_list[:frame], [e[2] for e in Energy_totals[:frame]])

    Theta_1_phase_plot.set_data([v[0] for v in variables_list[:frame]], [v[2] for v in variables_list[:frame]])

    Theta_2_phase_plot.set_data([v[1] for v in variables_list[:frame]], [v[3] for v in variables_list[:frame]])

    time_text.set_text('time = %.1fs' % (frame * h))
    
    return (
        animated_string_1,
        animated_bob_1,
        animated_string_2,
        animated_bob_2,
        trace_bob_2,
        kinetic_plot,
        potential_plot,
        total_energy_plot,
        Theta_1_phase_plot,
        Theta_2_phase_plot
    )

G = 9.8 # [g] = m/s^2
L1 = 1.0  #[L] = m
L2 = 1.0  
M1 = 1.0  #[M] = kg
M2 = 1.0
M_sum = M1 + M2

Theta_1 = np.pi / 3  #[Theta] = rad
Theta_2 = - np.pi / 6
Theta_dot_1 = 0 #[Theta dot] = rad/s
Theta_dot_2 = 0

variables = np.array([Theta_1, Theta_2, Theta_dot_1, Theta_dot_2])

t = 0 # [t] = s
h = 0.025 # step-size, [h] = s
steps = 1000
history_len = 100 # number of steps to keep as a trail

Cartesian_corrdinates = []
x_1, y_1, x_2, y_2 = polar_to_rect(Theta_1, Theta_2)
Cartesian_corrdinates.append([x_1, y_1, x_2, y_2])
kinetic, potential, Energy =  calc_Energy(variables)
Energy_totals= [[kinetic, potential, Energy]]
times_list = [t]
variables_list = [variables]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, figsize=(10,10))


animated_string_1, = ax1.plot([],[], color='green')
animated_bob_1, = ax1.plot([],[], marker= 'o', markersize=4, color='black')
animated_string_2, = ax1.plot([],[], color='red')
animated_bob_2, = ax1.plot([],[], marker= 'o', markersize=4, color='black')
trace_bob_2, = ax1.plot([],[], linestyle='-', color='brown')
ax1.set_xlim(-1* (L1 + L2), L1 + L2)
ax1.set_ylim(-1.5 * (L1 + L2), L1 + L2)
ax1.set_title("Double Pendulum Simulation")
time_text = ax1.text(0.05, 0.9, '', transform=ax1.transAxes)

kinetic_plot, = ax2.plot([],[], color='green', label="Kinetic Energy")
potential_plot, = ax2.plot([],[], color='red', label="Potential Energy")
total_energy_plot, = ax2.plot([],[], color='blue', label="Total Energy")
initial_energy = ax2.axhline(y=Energy_totals[0][2], color='black', label="Initial Energy", linestyle='--')
ax2.set_xlim(0, t + h * steps)
ax2.set_ylim(-25, 25)
ax2.legend()
ax2.set_title("Total Energy")

Theta_1_phase_plot, = ax3.plot([],[], color='black')
ax3.set_xlim(-2 * np.pi,  2 * np.pi)
ax3.set_ylim(-2 * np.pi,  2 * np.pi)
ax3.set_title("$\Theta_1$ and $\dot\Theta_1$ Phase Space")

Theta_2_phase_plot, = ax4.plot([],[], color='black')
ax4.set_xlim(-2 * np.pi,  2 * np.pi)
ax4.set_ylim(-2 * np.pi,  2 * np.pi)
ax4.set_title("$\Theta_2$ and $\dot\Theta_2$ Phase Space")

while (t <= h * steps):
    t, variables = RK4(t, h, variables)
    Theta_1, Theta_2, Theta_dot_1, Theta_dot_2 = variables
    x_1, y_1, x_2, y_2 = polar_to_rect(Theta_1, Theta_2)
    Cartesian_corrdinates.append([x_1, y_1, x_2, y_2])
    kinetic, potential, Energy =  calc_Energy(variables)
    Energy_totals.append([kinetic, potential, Energy])
    times_list.append(t)
    variables_list.append(variables)

animation = FuncAnimation(fig=fig, func=update_frame, frames=steps, interval=25)
animation.save("DoublePendulum.gif")
plt.show()