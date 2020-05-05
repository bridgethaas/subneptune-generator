import os
import numpy as np
import mesa_reader as mr

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
    
    
