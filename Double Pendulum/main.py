import numpy as np
import matplotlib.pyplot as plt

# Pendulum Lengths (m)
L1, L2 = 1, 1
# Bob masses (kg), and moments of intertia (kg m^2), I = 1/12ml^2 around the center of each rod.
M1, M2 = 1, 1
I1 = 1/12 * M1 * L1
I2 = 1/12 * M2 * L2
# contants, [g] = m/s^-2
G = 9.81
# thetas
theta1, theta2 = 0, 0
theta1dot, theta2dot = 0, 0
# position varaibles for pendulum 1
x1 = 0.5 * L1 * np.sin(theta1)
y1 = -0.5 * L1 * np.cos(theta1)
# position variables for pendulum 2
x2 = L1 * np.sin(theta1) + 0.5 * L2 * np.sin(theta2)
y2 = -1 * L1 * np.cos(theta1) - 0.5 * L2 * np.cos(theta2)
# velocity variables, derivative of x,y wrt time
x1dot = theta1dot * 0.5 * L1 * np.sin(theta1)
y1dot = theta1dot * 0.5 * L1 * np.cos(theta1)
x2dot = theta1dot * L1 * np.sin(theta1) + 0.5 * theta2dot * L2 * np.sin(theta2)
y2dot = -1 * L1 * np.cos(theta1) - 0.5 * L2 * theta2dot * np.cos(theta2)
# L = T - V

def Calc_Energy(y1,y2,x1dot,x2dot,y1dot,y2dot,theta1dot,theta2dot):
    # Energy is the sum of kinetic and potential energy, v =mgh, t = 1/2mv^2 + 1/2Iomega^2
    v =-1*G*(M1*y1 + M2*y2)
    t = 1/2 * (M1 *(x1dot^2+y1dot^2)+ M2*(x2dot^2+y2dot^2)) + 1/2*(I1*theta1dot^2 + I2*theta2dot^2)
    return t + v