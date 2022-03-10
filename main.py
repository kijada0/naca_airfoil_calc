import math
import matplotlib.pyplot as plt
import numpy as np

resolution = 20

z = [[[], []], [[], []]]

def zg(x, c):
    if x == 0: return(0)
    else: return(c * ((1.4845 * (x**(1/2))) - (0.63 * (x**1)) - (1.758 * (x**2)) + (1.4215 * (x**3)) - (0.5075 * (x**4))))

custom = False
if custom:
    NACA = input("Wprowadź profil: ")
    lenght = input("Wprowadź dłogość profilu: ")
else:
    NACA  = "NACA0012"
    lenght = 2

naca = list(NACA)
#print(naca)

if len(naca) == 8:
    #   ---------------------------------------------
    #   ----------  PROFILE 4-RO CYFROWE    ---------
    #   ---------------------------------------------
    if int(naca[4]) == 0 and int(naca[5]) == 0:
        #   ---------------------------------------------
        #   ----------  PROFIL SYMETRYCZNY  -------------
        #   ---------------------------------------------

        #   ----------  DANE ŹRÓDŁOWE   -----------------
        c = (int(naca[6])*10 + int(naca[7]))

        #   ----------  DANE ŹRÓDŁOWE   -----------------
        for i in range(resolution+1):
            x = (i*100)/(resolution)
            zgg=
            z[0][0][].append(x)
            z[0][1].append([zg(x/100, c), -1*zg(x/100, c)])
            z[1][0].append((x/100)*lenght)
            z[1][1].append([lenght*zg(x/100, c)/100, -1*lenght*zg(x/100, c)/100])

            print(i, "\t->\t", x, "\t->\t", round(zg(x/100, c), 3), "\t->\t", round(zg(x/100, c), 3))

    else:
        #   ---------------------------------------------
        #   ----------  PROFIL ASYMETRYCZNY -------------
        #   ---------------------------------------------

        f = int(naca[4])/100
        xf = int(naca[5])/10
        c = (int(naca[6])*10 + int(naca[7]))

        for i in range(resolution+1):
            x = (i*100)/(resolution)


            if x/100 < xf:
                #print(xf, "\t",f)
                zf = (f / (xf**2)) * (2 * xf * x/100 - (x/100)**2)
                a = np.arctan(((2*f) / (xf**2)) * (xf-(x/100)))

            else:
                zf = (f / ((1-xf)**2)) * ((1-(2*xf)) + (2*xf*(x/100) - ((x/100)**2)))/2
                a = np.arctan(((2*f) / ((1-xf)**2)) * (xf-(x/100)))

            zgg = zg(x/100, c)

            z[0][0].append(x - (zgg * np.sin(a)))               #oś X 1
            z[0][1].append(zf + (zgg*np.cos(a)))                #oś Y 1
            z[1][0].append(x + (zgg * np.sin(a)))               #oś X 2
            z[1][1].append(zf - (zgg*np.cos(a)))                #oś Y 2

            print(x, "\t", round(zf, 5), "\t", round(a, 3), "\t", round(zgg, 3), "\t", round(x - (zgg * np.sin(a)), 3), "\t", round(zf + (zgg*np.cos(a)),3))
            #x, zf, a, zgg, zge

elif len(naca) == 9:
    print("5 digit profile")
    if int(naca[6]) == 0:
        print("Profile without undercut")
    elif int(naca[6]) == 1:
        print("Profile with undercut")
    else:
        print("Profile Error")

else:
    print("Profile Error")

#print(z)
#print(max(z[1]))

plt.subplot(2,1,1)
plt.plot(z[0][0][0], z[0][0][1])  # top half
plt.plot(z[0][1][0], z[0][1][1])  # bottom half
plt.title(NACA + " [%]")
plt.axis('equal')
plt.grid()

plt.subplot(2,1,2)
plt.plot(z[1][0][0], z[1][0][1])    # top half
plt.plot(z[1][1][0], z[1][1][1])    # bottom half
plt.title(NACA + " [m]")
plt.axis('equal')
plt.grid()

plt.show()
