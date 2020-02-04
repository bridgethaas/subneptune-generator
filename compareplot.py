import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr

mearth = 5.97e27
msun = 1.9892e33
rearth = 6.371008e8
rsun = 6.9598e10
rfrac = rsun/rearth

plt.figure()
ax = plt.gca()

h1no = mr.MesaData('no_mass_loss_jan27_v2/LOGS/evolve_7.5_0.01_0.24_0.02_0.03392_7.3_0.1.mod',file_type='log')
ax.plot(h1no.star_age,h1no.radius*rfrac,'b-',label='$f=1\%$, no mass loss')

h1w = mr.MesaData('with_mass_loss_jan28/LOGS/evolve_7.5_0.01_0.24_0.02_0.03392_7.3_0.1.mod', file_type='log')
ax.plot(h1w.star_age,h1w.radius*rfrac,'b--',label='$f=1\%$, with mass loss')

h2n = mr.MesaData('no_mass_loss_jan27_v2/LOGS/evolve_7.5_0.02_0.24_0.02_0.03392_7.3_0.1.mod', file_type='log')
ax.plot(h2n.star_age,h2n.radius*rfrac,'r-',label='$f=2\%$, no mass loss')

h2w = mr.MesaData('with_mass_loss_jan28/LOGS/evolve_7.5_0.02_0.24_0.02_0.03392_7.3_0.1.mod', file_type='log')
ax.plot(h2w.star_age,h2w.radius*rfrac,'r--',label='$f=2\%$, with mass loss')

h5n = mr.MesaData('no_mass_loss_jan27_v2/LOGS/evolve_7.5_0.05_0.24_0.02_0.03392_7.3_0.1.mod',file_type='log')
ax.plot(h5n.star_age,h5n.radius*rfrac,'k-',label='$f=5\%$, no mass loss')

h5w = mr.MesaData('with_mass_loss_jan28/LOGS/evolve_7.5_0.05_0.24_0.02_0.03392_7.3_0.1.mod', file_type='log')
ax.plot(h5w.star_age,h5w.radius*rfrac,'k--',label='$f=5\%$, with mass loss')  

h10n = mr.MesaData('no_mass_loss_jan27_v2/LOGS/evolve_7.5_0.1_0.24_0.02_0.03392_7.3_0.1.mod',file_type='log')
ax.plot(h10n.star_age,h10n.radius*rfrac,'g-',label='$f=10\%$, no mass loss')

h10w = mr.MesaData('with_mass_loss_jan28/LOGS/evolve_7.5_0.1_0.24_0.02_0.03392_7.3_0.1.mod', file_type='log')
ax.plot(h10w.star_age,h10w.radius*rfrac,'g--',label='$f=10\%$, with mass loss')

plt.legend(loc=0,fontsize='xx-small')
ax.set_xlabel('Time (yr)')
ax.set_ylabel('Radius (R_earth)')
ax.set_xscale('log')
plt.savefig('with_mass_loss_jan28/plots/compare_log.png')

plt.show()
