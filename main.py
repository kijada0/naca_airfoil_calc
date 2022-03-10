import math
import matplotlib.pyplot as plt
import numpy as np

resolution = 20
plot = [[[[], []], [[], []]], [[[], []], [[], []]]] #  [0-%, 1-m] [0-top, 1-bottm] [0-0X, 1-OY]

#   ----------  DANE ŹRÓDŁOWE   ----------
NACA  = "NACA2412"
lenght = 2  # [m]
naca = list(NACA)

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

"""

def y_t(x, t):
    if x == 0:  y = 0
    else:       y = t * ((1.4845 * (x**(1/2))) - (0.63 * (x**1)) - (1.758 * (x**2)) + (1.4215 * (x**3)) - (0.5075 * (x**4)))
    return(y/100)

def y_c(x, m, p):
    if x < p:   y = (m / (p**2)) * (2*p*x - x**2)
    else:       y = (m / ((1-p) **2)) * ((1 - 2*p) + 2*p*x - x**2)
    return(y)

def arc(x, m, p):
    if x < p:   y = (2*m / (p**2)) * (p - x)
    else:       y = (2*m / ((1-p) **2)) * (p-x)
    return(np.arctan(y))

if len(naca) == 8:  # PROFILE 4 CYFROWE

    t = (int(naca[6])*10 + int(naca[7]))
    m = int(naca[4])/100
    p = int(naca[5]) *10 /100
    print("t = ",  t, "\tm = ", m, "\tp = ", p)

    for i in range(resolution + 1):
        x_1 = i / resolution
        x_100 = i *100 / resolution

        z_t = y_t(x_1, t)
        z_c = y_c(x_1, m, p)
        angle = arc(x_1, m, p)

        x_u = x_1 - z_t * np.sin(angle)
        y_u = z_c + z_t * np.cos(angle)
        x_l = x_1 + z_t * np.sin(angle)
        y_l = z_c - z_t * np.cos(angle)

        plot[0][0][0].append(x_u)
        plot[0][0][1].append(y_u)
        plot[0][1][0].append(x_1)
        plot[0][1][1].append(y_l)

        plot[1][0][0].append(x_u *lenght)
        plot[1][0][1].append(y_u *lenght)
        plot[1][1][0].append(x_1 *lenght)
        plot[1][1][1].append(y_l *lenght)

        print(x_1, "\t", round(z_t, 2), "\t", round(z_c, 2), "\t", round(angle, 2), "\t\t\t\t\t", round(x_u, 2), "\t", round(y_u, 2), "\t", round(x_l, 2), "\t", round(y_l, 2))

else:
    print("Profile Error")

#   ----------  RYSOWANIE WYKRESÓW   -----------------
plt.subplot(2,1,1)
plt.plot(plot[0][0][0], plot[0][0][1])  # top half
plt.plot(plot[0][1][0], plot[0][1][1])  # bottom half
plt.title(NACA + " [%]")
plt.axis('equal')
plt.grid()

plt.subplot(2,1,2)
plt.plot(plot[1][0][0], plot[1][0][1])    # top half
plt.plot(plot[1][1][0], plot[1][1][1])    # bottom half
plt.title(NACA + " [m]")
plt.axis('equal')
plt.grid()

plt.show()