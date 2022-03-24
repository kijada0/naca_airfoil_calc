import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline as mip

resolution = 99
plot = [[[[], []], [[], []]], [[[], []], [[], []]]] #  [0-%, 1-m] [0-top, 1-bottm] [0-0X, 1-OY]
x_ = []
y_ = []

#   ----------  DANE ŹRÓDŁOWE   ----------
NACA  = "NACA25112"
lenght = 1.8  # [m]
f = 0.3913


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
p       -   położenie maksymalnej strzałki ugięcia          w postacji ułamka
c       -   współczynnik siły nośnej                        w postacji ułamka

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

def y_c5(x, m, p, f, k1, k2, u):
    if u == 0:
        if x <= p:  y = (k1/6)*((x**3) - 3*f*(x**2) + (f**2)*(3-f)*x)
        else:       y = (k1/6)*((f**3)*(1-x))
    else:
        if x <= p:  y = (k1/6) * (((x-f)**3) - (k2/k1)*x*((1-f)**3) - x*(f**3) + f**3)
        else:       y = (k1/6) * ((k2/k1) * ((x-f)**3) - (k2/k1)*x*((1-f)**3) - x*(f**3) + (f**3))
    return(y)

def arc4(x, m, p):
    if x < p:   y = (2*m / (p**2)) * (p - x)
    else:       y = (2*m / ((1-p) **2)) * (p-x)
    return(np.arctan(y))

def arc5(x, m, p, f, k1, k2, u):
    if u == 0:
        if x <= p:  y = (k1/6)*((3*(f**2))-(6*f*x)-((f**2)*(3-f)))
        else:       y = (k1/6)*(f**3)
    else:
        if x <= p:  y = (k1/6) * (3*((x-f)**2) - (k2/k1)*x*((1-f)**3) - f**3)
        else:       y = (k1/6) * (3*((k2/k1) * ((x-f)**2)) - (k2/k1)*((1-f)**3) - (f**3))
    return(np.arctan(y))

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

        #print(x_1, "\t", round(z_t, 3), "\t", round(z_c, 3), "\t", round(angle, 3), "\t\t\t\t\t", round(x_u, 3), "\t", round(y_u, 3), "\t", round(x_l, 3), "\t", round(y_l, 3))

elif len(naca) == 9:  # PROFILE 5 CYFROWE

    c = int(naca[4])*(3/20)
    p = int(naca[5])/20
    t = int(naca[7])*10 + int(naca[8])*1
    u = int(naca[6])
    print("c = ", c, "\tp = ", p, "\tt = ", t)

    m = ((3*f - 7*f**2 + 8*f**3 - 4*f**4) / (math.sqrt(f * (1-f)))) - (3/2) * (1-2*f) * ((math.pi/2) - np.arcsin(1 - 2*f))
    k1 = (6*c) / m
    k2 = k1 * ((3*(f-p)**2 - f**3) / (1-f)**3)
    print("m = ", round(m, 3), "\tk1 = ", round(k1, 3), "\tk2 = ", round(k2, 3))

    for i in range(resolution + 1):
        x = i / resolution

        if x > p-0.075 and x < p+0.075:
            print("Abord: ", x)
            continue

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


spline = mip(plot[0][0][0], plot[0][0][1], k=3, )
x_ = np.linspace(0, 1, 500)
y_ = spline(x_)


#   ----------  RYSOWANIE WYKRESÓW   -----------------
plt.subplot(2, 1, 1)
#plt.plot(plot[0][0][0], plot[0][0][1])  # top half
plt.plot(x_, y_)
plt.plot(plot[0][1][0], plot[0][1][1])  # bottom half
plt.title(NACA + " [%]")
plt.axis('equal')
plt.grid()


plt.subplot(2, 1, 2)
plt.plot(plot[1][0][0], plot[1][0][1])    # top half
plt.plot(plot[1][1][0], plot[1][1][1])    # bottom half
plt.title(NACA + " [m]")
plt.axis('equal')
plt.grid()


plt.show()