import math
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.interpolate import make_interp_spline as mip

resolution = 500
plot = [[[[], []], [[], []]], [[[], []], [[], []]]] #  [0-%, 1-m] [0-top, 1-bottm] [0-0X, 1-OY]
x_ = []
y_ = []

#   ----------  DANE ŹRÓDŁOWE   ----------
NACA  = "NACA63112"
lenght = 1  # [m]
f = 0.29

print("input:\n", NACA, "\tl = ", lenght, "\tf = ", f, "\notuput:")

"""
x_1     -   wartość dla której liczona jest grubość         w postaci ułamka
x_100   -   wartość dla której liczona jest grubość         w postaci %

y_t     -   dystrybucja grubości    -   funkcja
z_t     -   wartość wzrucona przez Y_t

y_c     -   funcka obliczająca współrzędne linii wygięcia
z_c     -   wartość wzrucona przez Y_c

d_d     -   obliczanie pochodnej 
angle   -   kąt 

t       -   maksymalna grobość pfofilu                      w postaci ułamka
m       -   maksymalna strzałka ugięcia                     w postaci ułamka
p       -   położenie maksymalnej strzałki ugięcia          w postaci ułamka
c       -   współczynnik siły nośnej                        w postaci ułamka

"""

naca = list(NACA)

def y_t(x, t):
    if x == 0:  y = 0
    else:       y = t * ((1.4845 * (x**(1/2))) - (0.63 * (x**1)) - (1.758 * (x**2)) + (1.4215 * (x**3)) - (0.5075 * (x**4)))
    return(y/100)

def y_c4(x, m, p):
    if x < p:   y = (m / (p**2)) * (2*p*x - (x**2))
    else:       y = (m / ((1-p) **2)) * ((1 - 2*p) + 2*p*x - x**2)
    return(y)

def y_c5(x, m, p, f, k1, k2, g):
    if g == 0:
        if x < f:  y = (k1/6) * ((x**3) - (3*f*(x**2)) + ((f**2)*(3-f)*x))
        else:      y = ((k1 * (f**3)) /6) * (1-x)
    else:
        if x < f:  y = (k1/6) * (((x-f)**3) - (k2/k1)*x*((1-f)**3) - x*(f**3) + f**3)
        else:      y = (k1/6) * ((k2/k1) * ((x-f)**3) - (k2/k1)*((1-f)**3)*x - x*(f**3) + (f**3))
    return(y)

def arc4(x, m, p):
    if x < p:   y = (2*m / (p**2)) * (p - x)
    else:       y = (2*m / ((1-p) **2)) * (p-x)
    return(np.arctan(y))

def arc5(x, m, p, f, k1, k2, g):
    if g == 0:
        if x < f:  y = (k1/6) * ((3*(x**2)) - (6*f*x) + ((f**2)*(3-f)))
        else:      y = (k1 * (f**3)) /6
    else:
        if x < f:  y = (k1/6) * (3*((x-f)**2) - (k2/k1)*((1-f)**3) - f**3)
        else:      y = (k1/6) * (3*((k2/k1) * ((x-f)**2)) - (k2/k1)*((1-f)**3) - (f**3))
    return(np.arctan(y))

def f_out(d5, d6):
    d = str(d5) + str(d6)

    f_tab = [[10, 0.0580], [20, 0.1260], [30, 0.2025], [40, 0.2900], [50, 0.3910], [21, 0.1300], [31, 0.2170], [41, 0.3180], [51, 0.4410], [00, 0.5]]
    for i in range(len(f_tab)):
        if d == str(f_tab[i][0]): break
    return f_tab[i][1]

if len(naca) == 8:  # PROFILE 4 CYFROWE

    t = (int(naca[6])*10 + int(naca[7]))
    m = int(naca[4])/100
    p = int(naca[5]) *10 /100
    print("t = ",  t, "\tm = ", m, "\tp = ", p)

    for i in range(resolution + 1):
        x = i / resolution

        z_t = y_t(x, t)
        z_c = y_c4(x, m, p)
        angle = arc4(x, m, p)

        x_u = x - z_t * np.sin(angle)
        y_u = z_c + z_t * np.cos(angle)
        x_l = x + z_t * np.sin(angle)
        y_l = z_c - z_t * np.cos(angle)

        plot[0][0][0].append(x_u)
        plot[0][0][1].append(y_u)
        plot[0][1][0].append(x_l)
        plot[0][1][1].append(y_l)

        plot[1][0][0].append(x_u *lenght)
        plot[1][0][1].append(y_u *lenght)
        plot[1][1][0].append(x_l *lenght)
        plot[1][1][1].append(y_l *lenght)

        print(x, "\t", round(z_t, 3), "\t", round(z_c, 3), "\t", round(angle, 3), "\t\t\t\t\t", round(x_u, 3), "\t", round(y_u, 3), "\t", round(x_l, 3), "\t", round(y_l, 3))

elif len(naca) == 9:  # PROFILE 5 CYFROWE

    c = int(naca[4])*(3/20)
    p = int(naca[5])/20
    t = int(naca[7])*10 + int(naca[8])*1
    u = int(naca[6])
    print("c = ", c, "\tp = ", p, "\tt = ", t)

    f = f_out(naca[5], naca[6])
    m = ((3*f - 7*f**2 + 8*f**3 - 4*f**4) / (math.sqrt(f * (1-f)))) - (3/2) * (1-2*f) * ((math.pi/2) - np.arcsin(1 - 2*f))
    k1 = (6*c) / m
    k2 = k1 * ((3*(f-p)**2 - f**3) / ((1-f)**3))
    print("m = ", round(m, 3), "\tf = ", f, "\tk1 = ", round(k1, 3), "\tk2 = ", round(k2, 3), "\tk2/k1 = ", round((k2/k1),5))
    print("\nx  \t\tz_t\t\tz_c\t\tangle")

    for i in range(resolution + 1):
        x = i / resolution

        #if x > p-0.075 and x < p+0.075: continue

        z_t = y_t(x, t)
        z_c = y_c5(x, m, p, f, k1, k2, u)
        angle = arc5(x, m, p, f, k1, k2, u)

        x_u = x - z_t * np.sin(angle)
        y_u = z_c + z_t * np.cos(angle)
        x_l = x + z_t * np.sin(angle)
        y_l = z_c - z_t * np.cos(angle)

        plot[0][0][0].append(x_u)
        plot[0][0][1].append(y_u)
        plot[0][1][0].append(x_l)
        plot[0][1][1].append(y_l)

        plot[1][0][0].append(x_u *lenght)
        plot[1][0][1].append(y_u *lenght)
        plot[1][1][0].append(x_l *lenght)
        plot[1][1][1].append(y_l *lenght)

        print(x, "\t", round(z_t, 3), "\t", round(z_c, 3), "\t", round(angle, 3))

else:
    print("Profile Error")

"""
spline = mip(plot[0][0][0], plot[0][0][1], k=3, )
x_ = np.linspace(0, 1, 500)
y_ = spline(x_)
"""

#   ----------  RYSOWANIE WYKRESÓW   -----------------
fig, ax = plt.subplots(2, 1, figsize=(8, 4))

ax[0].plot(plot[0][0][0], plot[0][0][1])  # top half
ax[0].plot(plot[0][1][0], plot[0][1][1])  # bottom half
ax[0].set_title(NACA + " [%]")
ax[0].axis('equal')
ax[0].grid()

ax[1].plot(plot[1][0][0], plot[1][0][1])    # top half
ax[1].plot(plot[1][1][0], plot[1][1][1])    # bottom half
ax[1].set_title(NACA + " [m]")
ax[1].axis('equal')
ax[1].grid()

plt.subplots_adjust(hspace=0.5)

plt.show()

print_points = 0
if print_points  == 1:
    print("Data point")
    j = len(plot[0][0][1]) -1
    for i in range(len(plot[0][0][1])):
        l = j-i
        print(round(plot[1][0][0][l], 3), "\t", round(plot[1][0][1][l], 3))
    for i in range(len(plot[0][1][1])):
        print(round(plot[1][1][0][i], 3), "\t", round(plot[1][1][1][i], 3))
