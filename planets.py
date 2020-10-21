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
        
        self.color = 'C%d'%len(self.stages)

        self.fnames = self.format_file_names()
        
        self.grid_masses, self.grid_radii, self.grid_fs, self.grid_ages, self.grid_luminosities = self.load_models(self.fnames)

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
        luminosities = []
        self.datadicts = []

        max_len = 0
        
        for i, fname in enumerate(fnames):
            try:
                h = mr.MesaData(fname, file_type='log')
            except OSError:
                print(fname, 'does not exist')
                continue
            if loud:
                print(fname)
            
            datadict = {}
            for bulk_name in h.bulk_names:
                history = getattr(h, bulk_name)
                try:
                    len(history)
                except:
                    history = [history]
                datadict[bulk_name] = history
            self.datadicts.append(datadict)

            masses.append(h.star_mass*mfrac)
            radii.append(h.radius*rfrac)
            fs.append(envelope_fraction(h))
            ages.append(h.star_age)
            luminosities.append(np.float64(h.luminosity))
            
            if type(h.star_age) == np.ndarray:
                if len(h.star_age) > max_len:
                    max_len = len(h.star_age)
            else:
                if max_len < 1:
                    max_len = 1
                else:
                    pass
        
        arrays = [masses, radii, fs, ages, luminosities]
        square_arrays = [np.zeros((len(masses),max_len))+np.nan for arr in arrays]
        for i, arr in enumerate(arrays):
            for j, model_arr in enumerate(arr):
                if type(model_arr) == np.ndarray:
                    this_len = len(model_arr)
                    square_arrays[i][j,:this_len] = model_arr
                else:
                    this_len = 1
                    square_arrays[i][j] = model_arr
        
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
        self.intpd_fs = []

        for i in range(0, len(self.grid_ages)):
            if self.final_ages[i] < newage: 
                self.intpd_radii.append(np.nan)
                self.intpd_masses.append(np.nan)
                self.intpd_fs.append(np.nan)
            else:
                rad_age_interp = interpolate.interp1d(self.grid_ages[i], self.grid_radii[i],kind='linear')
                self.intpd_radii.append(rad_age_interp(newage))

                mass_age_interp = interpolate.interp1d(self.grid_ages[i], self.grid_masses[i], kind='linear')
                self.intpd_masses.append(mass_age_interp(newage)) 
                
                f_age_interp = interpolate.interp1d(self.grid_ages[i], self.grid_fs[i],kind='linear')
                self.intpd_fs.append(f_age_interp(newage))

        self.intpd_radii = np.array(self.intpd_radii).flatten()
        self.intpd_masses = np.array(self.intpd_masses).flatten()
        self.intpd_fs = np.array(self.intpd_fs).flatten()

        self.grid_points = np.column_stack((self.grid_masses[:,0], self.grid_fs[:,0]))

        self.radius_interp = interpolate.LinearNDInterpolator(self.grid_points, self.intpd_radii)
        self.mass_interp = interpolate.LinearNDInterpolator(self.grid_points, self.intpd_masses)
        self.f_interp = interpolate.LinearNDInterpolator(self.grid_points, self.intpd_fs)

    def run_mcmc(self, newage, nwalkers, burn_in=500, nsteps=10000, prior='flat'):

        #call method to define the interpolation functions
        self.make_interp_functions(newage)


        #define (log) likelihood
        def log_likelihood(x, mu_m, sigma_m, mu_r, sigma_r):
            lval = -0.5 * ( (( (self.mass_interp(x[0],x[1]) - mu_m) / sigma_m ) ** 2) +
                            (( (self.radius_interp(x[0],x[1]) - mu_r) / sigma_r ) ** 2) )

            if np.isnan(lval): 
                lval = -np.inf

            return lval
        
        #define (log) prior
        def log_prior(x,prior):
            if prior == 'flat':
                return 1.0
            if prior == 'logflat':
                return 1.0/x
            else:
                raise 

        #prob = prior * likelihood
        def log_prob(x, mu_m, sigma_m, mu_r, sigma_r,prior):
            return log_prior(x,prior) + log_likelihood(x, mu_m, sigma_m, mu_r, sigma_r) 


        #hardcode ndim to set initial walker positions properly
        ndim = 2
        p0 = np.random.rand(nwalkers, ndim)

        #fixing initial walker positions
        p0[:,0] = p0[:,0] +  np.min(self.mpList)
        p0[:,1] = (p0[:,1] * 0.01) + np.min(self.fList)

        sampler = emcee.EnsembleSampler( nwalkers, ndim, log_prob, 
                                         args=[self.mass,self.mass_unc,self.radius,self.radius_unc,prior],
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


    def final_fs_for_samples(self, samples):
        samples_final_fs = []
        deltas = []
        for i in range(0,len(samples)):
            mi = samples[i, 0]
            fi = samples[i, 1]
            
            mf = self.mass_interp(mi, fi) 
            ff = 1.0 - ((mi/mf) * (1 - fi))
            samples_final_fs.append(ff) 
            
            delta = ((mf * ff) - (mi * fi)) / (mi * fi) 
            deltas.append(delta)
        return np.array(samples_final_fs), np.array(deltas)

    stages = ['pre_reduce', 'pre_core', 'comp', 'corel', 'reduce', 
              'corem', 'heating', 'remove_heating', 'irrad']
    def get_early_stages(self):
        planets = []
        for stage in self.stages:
            p = EarlyPlanet(str(stage),
                    self.mass, self.mass_unc,
                    self.radius, self.radius_unc,
                    self.datadir,
                    self.mpList,
                    self.fList,
                    self.orbitalList,
                    self.entropyList,
                    stage
            )
            planets.append(p)

        return planets



class EarlyPlanet(Planet):

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
            stage,
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
        self.stage = stage

        if len(kwargs.keys()):
            print(kwargs)

        self.color = 'C%d'%self.stages.index(stage)
    
        self.fnames = self.format_file_names()

        self.grid_masses, self.grid_radii, self.grid_fs, self.grid_ages, self.grid_luminosities = self.load_models(self.fnames)

    def __repr__(self):
        return self.name

    def format_file_names(self, formatter=None):
        if formatter is None:
            formatter = formatstring

        fnames = []

        if self.stage == 'pre_reduce':
            for i, m in enumerate(self.mpList):
                fname = self.datadir + '/hist_pre_reduce_%s.data'%(formatter(m))
                fnames.append(fname)

        elif self.stage == 'pre_core':
            for i, m in enumerate(self.mpList):
                for j, f in enumerate(self.fList):
                    fname = self.datadir + '/hist_pre_core_%s_%s.data'%(
                        formatter(m),
                        formatter(f)
                    )
                    fnames.append(fname)

        elif self.stage == 'comp' or self.stage == 'corel' or self.stage == 'reduce' or self.stage == 'corem':
            for i, m in enumerate(self.mpList):
                for j, f in enumerate(self.fList):
                    fname = self.datadir + '/hist_' + str(self.stage) + '_%s_%s_0.24000_0.02000.data'%(
                        formatter(m),
                        formatter(f)
                    )
                    fnames.append(fname)

        elif (self.stage == 'heating' or self.stage == 'remove_heating' or
              self.stage == 'cooling' or self.stage == 'remove_cooling'):
            for i, m in enumerate(self.mpList):
                ent = self.entropyList[i]

                for j, f in enumerate(self.fList):
                    fname = self.datadir + '/hist_' + str(self.stage) + '_%s_%s_0.24000_0.02000_%s.data'%(
                        formatter(m),
                        formatter(f),
                        formatter(ent)
                    )
                    fnames.append(fname)

        elif self.stage == 'irrad':
            for i, m in enumerate(self.mpList):
                ent = self.entropyList[i]

                for k, orbital in enumerate(self.orbitalList):
                    for j, f in enumerate(self.fList):
                        fname = self.datadir + '/hist_' + str(self.stage) + '_%s_%s_0.24000_0.02000_%s_%s.data'%(
                            formatter(m),
                            formatter(f),
                            formatter(orbital),
                            formatter(ent)
                        )
                        fnames.append(fname)

        else:
            print(str(self.stage) + ' is not a valid stage')

        return fnames
    





