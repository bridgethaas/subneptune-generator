import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr

mearth = 5.97e27
msun = 1.9892e33
rearth = 6.371008e8
rsun = 6.9598e10
rfrac = rsun/rearth
mfrac = msun/mearth

fig = plt.figure()

h1 = mr.MesaData('LOGS/hist_evolve_5.77_0.01_0.24_0.02_0.02584_7.23_0.1.data',file_type='log')
plt.plot(h1.star_age,h1.radius*rfrac,label='$f=1\%$')

h2 = mr.MesaData('LOGS/hist_evolve_5.77_0.02_0.24_0.02_0.02584_7.23_0.1.data',file_type='log')
plt.plot(h2.star_age,h2.radius*rfrac,label='$f=2\%$')

h5 = mr.MesaData('LOGS/hist_evolve_5.77_0.05_0.24_0.02_0.02584_7.23_0.1.data',file_type='log')
plt.plot(h5.star_age,h5.radius*rfrac,label='$f=5\%$')

h10 = mr.MesaData('LOGS/hist_evolve_5.77_0.1_0.24_0.02_0.02584_7.23_0.1.data',file_type='log')
plt.plot(h10.star_age,h10.radius*rfrac,label='$f=10\%$')

plt.errorbar(5.e9,2.04,0.06,fmt='k.',capsize=2)

ax = plt.gca()
#ax.set_xscale('log')
plt.gca().set_xlabel('Time (yr)')
#plt.gca().set_ylabel('Mass (M_Earth)')
plt.gca().set_ylabel('Radius (R_Earth)')
plt.title('K2-146b, With Mass Loss')
plt.legend(loc=0)

#print(h1.He4_Mass[0]/(h1.star_mass[0]*msun))
#print(h1.He4_Mass[-1]/(h1.star_mass[-1]*msun))

plt.savefig('r_vs_t.png')
plt.show()
