import mysubsprograms as my

evolve_profile = 'profile_evolve5.80000_0.01000_0.24000_0.02000_0.02600_7.23000_0.10000'
inlist_evolve = 'inlist_evolve_5.80000_0.01000_0.24000_0.02000_0.02600_7.23000_0.10000' 
irrad_mod = 'irrad_5.80000_0.01000_0.24000_0.02000_0.02600_7.23000.mod' 
evolve_mod = 'evolve_5.80000_0.01000_0.24000_0.02000_0.02600_7.23000_0.10000.mod' 
n_frac = 0.1
a = 1.0
ms = 0.331 
orb_sep = 0.026
ec = 1.0e9
column_depth = 3.113836479074946 
flux_dayside = 20764393.675119787 
formation_time = 6e6
teq = 550.0713820446309 
BA = 0.2
escape_regime = 0
diff_sep = 1

my.run_evolve(evolve_profile, inlist_evolve, irrad_mod, evolve_mod,
			n_frac, a, ms, orb_sep, ec, column_depth, flux_dayside,
			formation_time, teq, BA, escape_regime, diff_sep)
