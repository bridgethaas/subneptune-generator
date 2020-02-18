import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr

mearth = 5.97e27
msun = 1.9892e33
rearth = 6.371008e8
rsun = 6.9598e10
rfrac = rsun/rearth
mfrac = msun/mearth

plt.figure()
ax = plt.gca()

def formatstring(myfloat):
    return '%.5f'%myfloat

#mpList = [7.0,7.5,8.0]
mpList = list(np.linspace(7.0,8.0,11))
for i, m in enumerate(mpList):
    mpList[i] = formatstring(m)
#print(['%.5f'%item for item in mpList])

#fList = [0.009,0.01,0.011]
fList = list(np.logspace(-4,np.log10(2*(10**-2)),10))
for i, f in enumerate(fList):
    fList[i] = formatstring(f) 

entropyList = [7.28,7.28,7.29,7.29,7.3,7.3,7.3,7.31,7.31,7.32,7.32]
for i, ent in enumerate(entropyList):
    entropyList[i] = formatstring(ent)

#styles = ['-','--','-.']
alphas = [0.4,0.6,1]
colors = ['r','b','g']

def envelope_fraction(history):
    return history.envelope_mass/(history.star_mass*msun)

def calcX(history):
    return history.Hydrogen_Mass / history.envelope_mass

def calcY(history):
    return history.He4_Mass / history.envelope_mass

def calcZ(history):
    return 1 - calcX(history) - calcY(history)

for i, m in enumerate(mpList):
    ent = entropyList[i]

    for j, f in enumerate(fList):
        h = mr.MesaData('data/biggrid/hist_evolve_%s_%s_0.24000_0.02000_0.03392_%s_0.10000.data'%(m,f,ent) ,file_type='log')

        plt.subplot(4, 3, i+1)
        plt.plot(h.star_mass[-1]*mfrac,h.envelope_mass[-1]/(h.star_mass[-1]*msun),marker='.',color='b')
        #ax.plot(h.star_age,envelope_fraction(h))

#ax.set_xlabel('time, Gyr')
ax.set_ylabel('envelope fraction')
#ax.set_yscale('log')

#lines = ax.get_lines()
#ax.legend([lines[k] for k in [2,5,8]],['$M_0 = 7.0 M_E$','$M_0 = 7.5 M_E$','$M_0 = 8.0 M_E$']) 

plt.show()
