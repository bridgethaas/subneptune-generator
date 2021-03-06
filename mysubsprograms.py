#!/usr/bin/env python

import math
import numpy as np
import os
import shutil
import scipy
from scipy import loadtxt, optimize
import sys
import time
import random
import pandas as pd
from scipy.interpolate import interp1d
from scipy import interpolate


msun = 1.9892e33
rsun = 6.9598e10
rearth = 6.371008e8
mjup = 1.8986e30
rjup = 6.9911e9
mearth = 5.97e27
sigma=5.67e-5
au = 1.496e13

def formatstr(myfloat):
    return "%.5f"%myfloat
def expstr(myfloat):
    mystring = "%.5e"%myfloat
    mystring = mystring.replace('e','d')
    return mystring

def calculate_rho(mp, enFrac):
    observed_Mcore, observed_Rcore = loadtxt('coreMRcomp2_v40_all.txt', unpack=True, skiprows =11, usecols=[0,1])
    core_radius_function = interp1d(observed_Mcore, observed_Rcore, fill_value="extrapolate")

    #In units of Rearth and Mearth
    planet_core_mass = (mp * (1.0 - enFrac))

    planet_core_radius = core_radius_function(planet_core_mass)

    core_mass_cgs = planet_core_mass * mearth
    core_radius_cgs = planet_core_radius * rearth

    core_volume = (4./3.) * (np.pi) * (core_radius_cgs ** 3.0)
    rhocore = core_mass_cgs / float(core_volume)
    
    return (rhocore, planet_core_radius)


def calculate_column_depth(, profile, Teff_star):

    if Teff_star < 3500:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,5,6])
    elif 3500 <= Teff_star < 4500:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,7,8])
    elif 4500 <= Teff_star < 5500:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,9,10])
    elif 5500 <= Teff_star < 6500:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,11,12])
    else:
        T,P, k_r, k_p = loadtxt('OpacityTableSolarMetal.txt', unpack=True, skiprows =38, usecols=[0,1,13,14])

    points = np.column_stack((T, np.log10(P)))
    Opacity_function = interpolate.LinearNDInterpolator(points, k_p)

    #R, T, P are in LOG10
    header = loadtxt(profile,
                        unpack=True,
                        skiprows=5,
                        max_rows=1,
                        dtype='str')
    header = list(header)
    cols = [header.index('zone'),
            header.index('mass'),
            header.index('logR'),
            header.index('logT'),
            header.index('logP') ]
    
    zone, mass, radius, temperature, pressure = loadtxt(profile, unpack=True, skiprows=6, usecols=cols)

    switch_zone = []
    for i in range(len(zone)):
        mass_column_depth = ((mass[0] - mass[i]) * msun) / (4 * np.pi * (((10**radius[i]) * rsun) ** 2))
        opacity_column_depth = (2 / (Opacity_function(, (pressure[i]))))
        #LOG pressure
        switch_zone.append((opacity_column_depth - mass_column_depth, zone[i], mass_column_depth))

    column_depth = abs(switch_zone[0][0])

    if switch_zone[0][0] > 0:
        for i in range(len(switch_zone)):
            if switch_zone[i][0] < 0:
                column_depth = switch_zone[i - 1][2]
                break
    else:
        for i in range(len(switch_zone)):
            if switch_zone[i][0] > 0:
                column_depth = switch_zone[i - 1][2]
                break

    return column_depth


def run_pre_reduce(inlist_pre_reduce, initial_mod, pre_reduce_mod, mp):
    if os.path.isfile(pre_reduce_mod):
        return 0
    start_time = time.time()
    print('begin pre_reduce')
    f = open('inlist_pre_reduce', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + initial_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + pre_reduce_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + pre_reduce_mod.replace(".mod",".data") + '"')
    g = g.replace("<<mp>>",expstr((mp * 5 * mearth / msun)))
    #original factor of 30, changed to 5? so all models converge
    

    h = open(inlist_pre_reduce, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_pre_reduce, "inlist")


    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_pre_core(inlist_pre_core, pre_reduce_mod, pre_core_mod, enFrac, core_mass, rho):
    if os.path.isfile(pre_core_mod):
        return 0
    start_time = time.time()
    print('begin pre_core')
    f = open('inlist_pre_core', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + pre_reduce_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + pre_core_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + pre_core_mod.replace(".mod",".data") + '"')
    g = g.replace("<<core_mass>>", expstr(.1 * core_mass * mearth / msun))
    g = g.replace("<<rho>>", expstr(rho))
    
    h = open(inlist_pre_core, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_pre_core, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_comp(inlist_comp, pre_core_mod, comp_mod, z, y):
    if os.path.isfile(comp_mod):
        return 0
    start_time = time.time()
    print('begin comp')
    f = open('inlist_comp', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + pre_core_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + comp_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + comp_mod.replace(".mod",".data") + '"')
    g = g.replace("<<y>>",formatstr(y))
    g = g.replace("<<z>>",formatstr(z))

    h = open(inlist_comp, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_comp, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_corel(inlist_corel, comp_mod, corel_mod):
    if os.path.isfile(corel_mod):
        return 0
    start_time = time.time()
    print('begin corel')
    f = open('inlist_corel', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + comp_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + corel_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + corel_mod.replace(".mod",".data") + '"')

    h = open(inlist_corel, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_corel, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_reduce(inlist_reduce, corel_mod, reduce_mod, mp):
    if os.path.isfile(reduce_mod):
        return 0
    start_time = time.time()
    print('begin reduce')
    f = open('inlist_reduce', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + corel_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + reduce_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + reduce_mod.replace(".mod",".data") + '"')
    g = g.replace("<<mp>>",expstr((mp * mearth / msun)))
    
    h = open(inlist_reduce, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_reduce, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_corem(inlist_corem, reduce_mod, corem_mod, core_mass, rho):
    if os.path.isfile(corem_mod):
        return 0
    start_time = time.time()
    print('begin corem')
    f = open('inlist_corem', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + reduce_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + corem_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + corem_mod.replace(".mod",".data") + '"')
    g = g.replace("<<core_mass>>", expstr(core_mass * mearth / msun))
    g = g.replace("<<rho>>", expstr(rho))
    
    h = open(inlist_corem, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_corem, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_heating(inlist_heating, corem_mod, heating_mod, heating_profile, entropy, luminosity):
    if os.path.isfile(heating_mod):
        return 0
    start_time = time.time()
    print('begin heating')
    f = open('inlist_heating', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + corem_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + heating_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + heating_mod.replace(".mod",".data") + '"')
    g = g.replace("<<heating_profile>>",'"' + heating_profile + '"')
    g = g.replace("<<entropy>>", formatstr(entropy))
    
    # This is to inflate the planet
    g = g.replace("<<luminosity>>", expstr(luminosity))
    
    h = open(inlist_heating, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_heating, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_cooling(inlist_cooling, corem_mod, cooling_mod, cooling_profile, entropy):
    if os.path.isfile(cooling_mod):
        return 0
    start_time = time.time()
    print('begin cooling')
    f = open('inlist_cooling', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + corem_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + cooling_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + cooling_mod.replace(".mod",".data") + '"')
    g = g.replace("<<cooling_profile>>", '"' + cooling_profile + '"')
    g = g.replace("<<entropy>>", formatstr(entropy))
    
    h = open(inlist_cooling, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_cooling, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_remove_heating(inlist_remove_heating, heating_mod, remove_heating_profile,remove_mod):
    if os.path.isfile(remove_mod):
        return 0
    start_time = time.time()
    print('begin remove_heating')
    f = open('inlist_remove_heating', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + heating_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + remove_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + remove_mod.replace(".mod",".data") + '"')
    g = g.replace("<<remove_heating_profile>>", '"' + remove_heating_profile + '"')
    
    h = open(inlist_remove_heating, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_remove_heating, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_remove_cooling(inlist_remove_cooling, cooling_mod, remove_cooling_profile, remove_mod):
    if os.path.isfile(remove_mod):
        return 0
    start_time = time.time()
    print('begin remove_cooling')
    f = open('inlist_remove_cooling', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + cooling_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + remove_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + remove_mod.replace(".mod",".data") + '"')
    g = g.replace("<<remove_cooling_profile>>", '"' + remove_cooling_profile + '"')
    
    h = open(inlist_remove_cooling, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_remove_cooling, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_irrad(inlist_irrad, remove_mod, irrad_mod, irrad_profile, column_depth, flux_dayside):
    if os.path.isfile(irrad_mod):
        return 0
    start_time = time.time()
    print('begin irrad')
    f = open('inlist_irrad', 'r')
    g = f.read()
    f.close()

    g = g.replace("<<loadfile>>",'"' + remove_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + irrad_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + irrad_mod.replace(".mod",".data") + '"')
    g = g.replace("<<irrad_profile>>", '"' + irrad_profile + '"')
    g = g.replace("<<flux_dayside>>", expstr(flux_dayside))
    g = g.replace("<<column_depth>>", expstr(column_depth))

    h = open(inlist_irrad, 'w')
    h.write(g)
    h.close()
    shutil.copyfile(inlist_irrad, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time


def run_evolve(inlist_evolve, irrad_mod, evolve_mod, evolve_profile,
                n_frac, a, ms, orb_sep, ec, column_depth, flux_dayside,
                formation_time, teq, BA, escape_regime, diff_sep):

    start_time = time.time()
    print('begin evolve')
    f = open('inlist_evolve', 'r')
    g = f.read()
    f.close()

    #File Parameters
    g = g.replace("<<loadfile>>",'"' + irrad_mod + '"')
    g = g.replace("<<smwtfname>>", '"' + evolve_mod + '"')
    g = g.replace("<<hist_smwtfname>>", '"hist_' + evolve_mod.replace(".mod",".data") + '"')
    g = g.replace("<<evolve_profile>>", '"' + evolve_profile + '"')

    #Flux Parameters
    g = g.replace("<<formation_time>>",expstr(formation_time))
    g = g.replace("<<column_depth>>",expstr(column_depth))
    g = g.replace("<<flux_dayside>>", expstr(flux_dayside))

    #x-controls
    g = g.replace("<<n_frac>>", formatstr(n_frac))
    g = g.replace("<<a>>", formatstr(a))
    g = g.replace("<<ms>>", expstr(ms))
    g = g.replace("<<BA>>", formatstr(BA))
    g = g.replace("<<orb_sep>>", formatstr(orb_sep))
    g = g.replace("<<ec>>", expstr(ec))
    g = g.replace("<<teq>>", formatstr(teq))
    g = g.replace("<<escape_regime>>", str(escape_regime))
    g = g.replace("<<diff_sep>>", str(diff_sep))
    
    h = open(inlist_evolve, 'w')
    h.write(g)
    h.close()

    shutil.copyfile(inlist_evolve, "inlist")

    os.system('./star_make_planets')
    run_time = time.time() - start_time
    return run_time
