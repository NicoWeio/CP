import matplotlib.pyplot as plt
import numpy as np
import math


#A1
#a)################################################################################################
ha, dfa = np.genfromtxt("build/A1_a_h.csv", unpack=True, delimiter= ",")
plt.plot(ha, dfa, "x")
plt.plot(ha[np.where(np.abs(dfa-np.cos(1.5))== np.amin(np.abs(dfa-np.cos(1.5))))[0][0]],
         dfa[np.where(np.abs(dfa-np.cos(1.5))== np.amin(np.abs(dfa-np.cos(1.5))))[0][0]], "x",
         label=r"${|f'_{zweipunkt}(1.5,h)}-f'(1.5)|_{min}}$")
plt.xscale("log")
plt.xlabel(r'$h$')
plt.ylabel(r"$f'_{zweipunkt}(x)$")
plt.title(r"$f=\sin(x), x=1.5$")
plt.legend(loc="best")

plt.savefig('build/plot.pdf')
plt.clf()

######

xa, dfax = np.genfromtxt("build/A1_a_x.csv", unpack=True, delimiter= ",")

plt.plot(xa, np.abs(dfax-np.cos(xa)))
plt.xlabel(r'$x$')
plt.ylabel(r"$|f'_{zweipunkt}-f'_{analytisch}|$")
plt.title(r"$f=\sin(x), h=0.01$")

plt.savefig('build/A1_a_x.pdf')
plt.clf()

#b)#########################################################################################
hb, dfb = np.genfromtxt("build/A1_b_h.csv", unpack=True, delimiter= ",")
plt.plot(hb, dfb, "x")
plt.plot(hb[np.where(np.abs(dfb+np.sin(1.5))== np.amin(np.abs(dfb+np.sin(1.5))))[0][0]],
         dfb[np.where(np.abs(dfb+np.sin(1.5))== np.amin(np.abs(dfb+np.sin(1.5))))[0][0]], "x",
         label=r"${|f''_{zweipunkt}(1.5,h)}-f''(1.5)|_{min}}$")
plt.xscale("log")
plt.xlabel(r'$h$')
plt.ylabel(r"$f''_{zweipunkt}$")
plt.title(r"$f''=-\sin(x), x=1.5, h_{f'}=0.01$")
plt.legend(loc="best")

plt.savefig('build/A1_b_h.pdf')
plt.clf()
#print(hb[np.where(np.abs(dfb+np.sin(1.5))== np.amin(np.abs(dfb+np.sin(1.5))))[0][0]])
#print(-np.sin(1.5))

xb, dfbx = np.genfromtxt("build/A1_b_x.csv", unpack=True, delimiter= ",")

plt.plot(xb, np.abs(dfbx + np.sin(xb)))
plt.xlabel(r'$x$')
plt.ylabel(r"$|f''_{zweipunkt}-f''_{analytisch}|$")
plt.title(r"$f=\sin(x), h=0.01, h_{f'}=0.01$")

plt.savefig('build/A1_b_x.pdf')
plt.clf()

#c)######################################################################
hc, dfc = np.genfromtxt("build/A1_c_h.csv", unpack=True, delimiter= ",")
plt.plot(hc, dfc, "x")
plt.plot(hc[np.where(np.abs(dfc-np.cos(1.5))== np.amin(np.abs(dfc-np.cos(1.5))))[0][0]],
         dfc[np.where(np.abs(dfc-np.cos(1.5))== np.amin(np.abs(dfc-np.cos(1.5))))[0][0]], "x",
         label=r"${|f'_{vierpunt}(1.5,h)}-f'(1.5)|_{min}}$")
plt.xscale("log")
plt.xlabel(r'$h$')
plt.ylabel(r"$f'_{vierpunkt}(x)$")
plt.title(r"$f=\sin(x), x=1.5$")
plt.legend(loc="best")

plt.savefig('build/A1_c_h.pdf')
plt.clf()

#print(hc[np.where(np.abs(dfc-np.cos(1.5))== np.amin(np.abs(dfc-np.cos(1.5))))[0][0]])

xc, dfcx = np.genfromtxt("build/A1_c_x.csv", unpack=True, delimiter= ",")

plt.plot(xc, np.abs(dfcx-np.cos(xc)))
plt.xlabel(r'$x$')
plt.ylabel(r"$|f'_{vierpunkt}-f'_{analytisch}|$")
plt.title(r"$f=\sin(x), h=0.01$")

plt.savefig('build/A1_c_x.pdf')
plt.clf()

#d)#######################################################################################
def f2(x):
    if x % np.pi == 0 :
       return NaN
    elif x >= 0:
        return   np.sin(x % np.pi)
    else :
        return - np.sin(x % np.pi)
#def f2(x):
#    if x >= 0:
#        return  + np.sin(x % np.pi)
#    else :
#        return  - np.sin(x % np.pi)

hd2, dfd2= np.genfromtxt("build/A1_d_h2.csv", unpack=True, delimiter= ",")
hd4, dfd4= np.genfromtxt("build/A1_d_h4.csv", unpack=True, delimiter= ",")

fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle(r"$f_2, x=1.5$")

ax1.plot(hd2, dfd2, "x")
ax1.plot(hd2[np.where(np.abs(dfd2-f2(1.5))== np.amin(np.abs(dfd2-f2(1.5))))[0][0]],
         dfd2[np.where(np.abs(dfd2-f2(1.5))== np.amin(np.abs(dfd2-f2(1.5))))[0][0]], "x",label=r"${|f'_{zweipunkt}(1.5,h)}-f'(1.5)|_{min}}$")
ax1.set(xscale="log", xlabel=r'$h$',ylabel= r"$f'_{zweipunkt}(x)$")
ax1.legend(loc="best")

print(hd2[np.where(np.abs(dfd2-f2(1.5))== np.amin(np.abs(dfd2-f2(1.5))))[0][0]])

ax2.plot(hd4, dfd4, "x")
ax2.plot(hd4[np.where(np.abs(dfd4-f2(1.5))== np.amin(np.abs(dfd4-f2(1.5))))[0][0]],
         dfd4[np.where(np.abs(dfd4-f2(1.5))== np.amin(np.abs(dfd4-f2(1.5))))[0][0]], "x",label=r"${|f'_{vierpunkt}(1.5,h)}-f'(1.5)|_{min}}$")
ax2.set(xscale="log", xlabel=r'$h$',ylabel= r"$f'_{vierpunkt}(x)$")
ax2.legend(loc="best")

print(hd4[np.where(np.abs(dfd4-f2(1.5))== np.amin(np.abs(dfd4-f2(1.5))))[0][0]])


plt.savefig('build/A1_d_h.pdf')
plt.clf()

########################

def f22(x):
    y=np.copy(x)
    for i,e in enumerate(x):
     y[i]= f2(e)
    return(y)

xd2, dfdx2= np.genfromtxt("build/A1_d_x2.csv", unpack=True, delimiter= ",")
xd4, dfdx4= np.genfromtxt("build/A1_d_x4.csv", unpack=True, delimiter= ",")

fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle(r"$f_2, x=1.5, h_{zwei}= 0.0001, h_{vier}=0.01$")

ax1.plot(xd2, np.abs(dfdx2-f22(xd2)), label= r"$|f'_{zweipunkt}-f'_{analytisch}|$")
ax1.plot(xd4, np.abs(dfdx4-f22(xd4)),label= r"$|f'_{vierpunkt}-f'_{analytisch}|$")
ax1.set(xlabel=r'$x$', ylabel= r"$f'(x)$")
ax1.legend(loc="best")

#ax2.plot(xd4, np.abs(dfdx4-f22(xd4)))
ax2.plot(xd4, np.abs(dfdx4-dfdx2))
ax2.set(xlabel=r'$x$', ylabel= r"$|f'_{vierpunkt}-f'_{zweipunkt}|$")

print()

plt.savefig('build/A1_d_x.pdf')
plt.clf()
