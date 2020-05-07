import os
import numpy as np
import mesa_reader as mr
from scipy import interpolate
import emcee

mearth = 5.97e27
msun = 1.9892e33
rearth = 6.371008e8
rsun = 6.9598e10
rfrac = rsun/rearth
mfrac = msun/mearth

def formatstring(myfloat):
    return '%.5f'%myfloat

def envelope_fraction(history):
    return history.envelope_mass/(history.star_mass*msun)

class Planet():
    def __repr__(self):
        return self.name
    
    def __init__(
        self,
        name,
        mass, mass_unc,
        radius, radius_unc,
        datadir,
        mpList,
        fList,
        orbitalList,
        entropyList,
        **kwargs
    ):
        self.name = name
        self.mass = mass
        self.mass_unc = mass_unc
        self.radius = radius
        self.radius_unc = radius_unc
        self.datadir = datadir
        self.mpList = mpList
        self.fList = fList
        self.orbitalList = orbitalList
        self.entropyList = entropyList 
        
        if len(kwargs.keys()):
            print(kwargs)
        
        self.fnames = self.format_file_names()
        
        self.grid_masses, self.grid_fs, self.grid_radii, self.grid_ages = self.load_models(self.fnames)

        self.final_masses, self.final_fs, self.final_radii, self.final_ages = self.get_finals()
        
        
    def format_file_names(self,formatter=None):
        if formatter is None:
            formatter = formatstring
        
        fnames = []
        for i, m in enumerate(self.mpList):
            ent = self.entropyList[i]
            
            for k, orbital in enumerate(self.orbitalList):
                for j, f in enumerate(self.fList):
                    fname = self.datadir + '/hist_evolve_%s_%s_0.24000_0.02000_%s_%s_0.10000.data'%(
                        formatter(m),
                        formatter(f),
                        formatter(orbital),
                        formatter(ent)
                    )
                    fnames.append(fname)
                    
        return fnames
    
    
    def load_models(self,fnames,loud=False):
        masses = []
        radii = []
        fs = []
        ages = []
        
        max_len = 0
        
        for i, fname in enumerate(fnames):
            h = mr.MesaData(fname,file_type='log')
            if loud:
                print(fname)
            
            masses.append(h.star_mass*mfrac)
            radii.append(h.radius*rfrac)
            fs.append(envelope_fraction(h))
            ages.append(h.star_age)
            
            if len(h.star_age) > max_len:
                max_len = len(h.star_age)
        
        arrays = [masses, fs, radii, ages]
        square_arrays = [np.zeros((len(fnames),max_len))+np.nan for arr in arrays]
        for i, arr in enumerate(arrays):
            for j, model_arr in enumerate(arr):
                this_len = len(model_arr)
                square_arrays[i][j,:this_len] = model_arr
        
        return square_arrays
    
    def get_finals(self):
        final_indices = np.argmax(np.isnan(self.grid_ages),axis=-1) - 1

        final_masses = self.grid_masses[np.arange(self.grid_masses.shape[0]), final_indices]
        final_fs = self.grid_fs[np.arange(self.grid_fs.shape[0]), final_indices]
        final_radii = self.grid_radii[np.arange(self.grid_radii.shape[0]), final_indices]
        final_ages = self.grid_ages[np.arange(self.grid_ages.shape[0]), final_indices]
        
        return final_masses, final_fs, final_radii, final_ages
    

    def make_interp_functions(self, newage):
        self.intpd_radii = []
        self.intpd_masses = []
        for i in range(0, len(self.grid_ages)):
            if self.final_ages[i] < newage: 
                self.intpd_radii.append(np.nan)
                self.intpd_masses.append(np.nan)
            else:
                rad_age_interp = interpolate.interp1d(self.grid_ages[i], self.grid_radii[i],kind='linear')
                self.intpd_radii.append(rad_age_interp(newage))

                mass_age_interp = interpolate.interp1d(self.grid_ages[i], self.grid_masses[i], kind='linear')
                self.intpd_masses.append(mass_age_interp(newage)) 

        self.intpd_radii = np.array(self.intpd_radii).flatten()
        self.intpd_masses = np.array(self.intpd_masses).flatten()

        grid_points = np.column_stack((self.grid_masses[:,0], self.grid_fs[:,0]))

        self.radius_interp = interpolate.LinearNDInterpolator(grid_points, self.intpd_radii)
        self.mass_interp = interpolate.LinearNDInterpolator(grid_points, self.intpd_masses)

    def run_mcmc(self, newage, nwalkers, burn_in=500, nsteps=10000):

        #call method to define the interpolation functions
        self.make_interp_functions(newage)

        #define (log) likelihood
        def log_likelihood(x, mu_m, sigma_m, mu_r, sigma_r):
            lval = -0.5 * ( (( (self.mass_interp(x[0],x[1]) - mu_m) / sigma_m ) ** 2) +
                            (( (self.radius_interp(x[0],x[1]) - mu_r) / sigma_r ) ** 2) )

            if np.isnan(lval): 
                lval = -np.inf

            return lval

        #hardcode ndim to set initial walker positions properly
        ndim = 2
        p0 = np.random.rand(nwalkers, ndim)

        #fixing initial walker positions
        p0[:,0] = p0[:,0] +  np.min(self.mpList)
        p0[:,1] = (p0[:,1] * 0.01) + np.min(self.fList)

        sampler = emcee.EnsembleSampler( nwalkers, ndim, log_likelihood, 
                                         args=[self.mass,self.mass_unc,self.radius,self.radius_unc],
                                         moves=[emcee.moves.StretchMove(a=4.0)] )

        #run burn-in
        state = sampler.run_mcmc(p0, burn_in)
        sampler.reset()

        #run mcmc
        sampler.run_mcmc(state, nsteps)

        samples = sampler.get_chain(flat=True)
        print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
        print("Autocorrelation time: {0:.2f} steps".format(sampler.get_autocorr_time()[0]))
        
        return samples 


