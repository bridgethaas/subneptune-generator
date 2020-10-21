import os
import numpy as np
import mesa_reader as mr
import matplotlib.pyplot as plt
import matplotlib 
from matplotlib import cm
from matplotlib.colors import Normalize, LogNorm
from scipy import interpolate
from planets_grid.planets import *

matplotlib.rcParams['figure.facecolor'] = 'white'


def plot_early_stages(ax, planet, x, y, planet_index=0, plot_evolve=True, **kwargs):
    #ax - the axes you want to plot on
    #planet - Planet class object we're plotting the data from
    #x, y - quantities on the x and y axes 

    #make list containing the early stages for planet
    ep = planet.get_early_stages()

    if plot_evolve == True:
        ep.append(planet)
        ax.set_xscale('log')

    #offset used when plotting vs timestep number
    offset = 1
    for stage in ep:
        #TODO: handle when grid of early planets
        if hasattr(stage, y):
            ys = getattr(stage, y)[planet_index]
        else:   
            ys = stage.datadicts[planet_index][y]

        if x == 'timestep':
            xs = np.arange(len(ys)) + offset
            offset += len(xs)
        else:
            xs = getattr(stage, x)

        ax.plot(xs, ys, '.', label=stage.name, color=stage.color, **kwargs)
    
        add_to_legend(ax, label=stage.name, c=stage.color) 

    ax.set_xlabel(x)
    ax.set_ylabel(y)



def plot_grid_mr(ax, planet, age, coloring=False, color_by=None, colormap='plasma', time_index=None):
    #makes time interpolation functions for M, R and fills planet.intpd_masses/radii
    planet.make_interp_functions(age)

    if coloring == True:
        absmin = np.nanmin(np.hstack(getattr(planet, color_by)))
        absmax = np.nanmax(np.hstack(getattr(planet, color_by)))
        cmap = plt.get_cmap(colormap)

    for i in range(0, len(planet.intpd_masses)):

        if coloring == True:
            if time_index is not None:
                thiscolor = (getattr(planet, color_by)[i][time_index] - absmin) / (absmax - absmin)
            else:   
                thiscolor = (getattr(planet, color_by)[i] - absmin) / (absmax - absmin)

            plt.scatter(planet.intpd_masses[i], planet.intpd_radii[i], color=cmap(thiscolor), marker='.') 
        else:
            ax.plot(planet.intpd_masses, planet.intpd_radii, '.')
        
    ax.set_xlabel('mass (M$_{\oplus}$)')
    ax.set_ylabel('radius (R$_{\oplus}$)')

    #cax = plt.gcf().add_axes([0.93, 0.13, 0.03, 0.75])
    if coloring==True:
        norm = Normalize(vmin=absmin, vmax=absmax)
        cbar=plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)        
        return cbar


def plot_r_vs_t(ax, planet, init_mass=None, init_f=None, label=None):
    #assumes make_interp_functions has already been run

    if init_mass is not None and init_f is not None:
        model_index = np.where( (np.around(planet.grid_points, 3) == 
                                [init_mass, init_f]).all(axis=1) )
        plt.plot(planet.grid_ages[model_index][0], planet.grid_radii[model_index][0], '-', label=label)
    else: 
        plt.plot(planet.grid_ages[0], planet.grid_radii[0], '-', label=label)

    ax.set_xlabel('planet age (yr)')
    ax.set_ylabel('planet radius (R$_{\oplus}$)') 


def plot_f_vs_t(ax, planet, init_mass, init_f, label=None):
    #assumes make_interp_functions has already been run

    model_index = np.where( (np.around(planet.grid_points, 3) == 
                            [init_mass, init_f]).all(axis=1) )
    
    plt.plot(planet.grid_ages[model_index][0], planet.grid_fs[model_index][0], '-', label=label)

    ax.set_xlabel('planet age (yr)')
    ax.set_ylabel('envelope fraction')


#-------------------------------------------------------------------------
from matplotlib.lines import Line2D
def add_to_legend(
    ax,
    label='',
    shape='line',
    loc=0,
    legend_kwargs=None,
    make_new_legend=False,
    **kwargs):

    legend = ax.get_legend()
    ## add the current legend to the tracked artists
    ##  and then pretend we weren't passed one
    if make_new_legend and legend is not None:
        ax.add_artist(legend)
        legend=None

    if legend is not None:
        lines = legend.get_lines()
        labels = [text.get_text() for text in legend.get_texts()]
        if lines[0].get_linestyle() == 'None':
            lines[0].set_marker('x')
    else:
        lines,labels=[],[]

    if legend_kwargs is None:
        legend_kwargs = {}

    ## make the new line
    if shape == 'line':
        line = Line2D(
        [0],[0],
        **kwargs)
    else:
        raise NotImplementedError

    if label not in labels:
        lines.append(line)
        labels.append(label)

    if loc in legend_kwargs:
        loc = legend_kwargs.pop('loc')
    ax.legend(lines,labels,loc=loc,**legend_kwargs)

    return ax
