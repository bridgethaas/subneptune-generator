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

#mpList = [7.0,7.5,8.0]
mpList = [7.5]
#fList = [0.009,0.01,0.011]
fList = np.logspace(-3,np.log10(2),10)
#entropyList = [7.28,7.3,7.32]
entropyList = [7.3]

#styles = ['-','--','-.']
alphas = [0.4,0.6,1]
colors = ['r','b','g']

for i, m in enumerate(mpList):
    ent = entropyList[i]

    for j, f in enumerate(fList):
        h = mr.MesaData('feb10/hist_evolve_%s_%s_0.24_0.02_0.03392_%s_0.1.data'%(m,f,ent) ,file_type='log')
        #ax.plot(h.star_mass[-1]*mfrac,h.radius[-1]*rfrac,color=colors[i],marker='.')
        ax.plot(h.star_age,h.envelope_mass/(h.star_mass*msun))

ax.set_xlabel('time, Gyr')
ax.set_ylabel('envelope fraction')

lines = ax.get_lines()
ax.legend([lines[k] for k in [2,5,8]],['$M_0 = 7.0 M_E$','$M_0 = 7.5 M_E$','$M_0 = 8.0 M_E$']) 

plt.show()
