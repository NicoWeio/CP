import matplotlib.pyplot as plt
import numpy as np


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