import os
import numpy as np
import mesa_reader as mr
import matplotlib.pyplot as plt
import matplotlib
from scipy import interpolate
from planets_grid.planets import *

matplotlib.rcParams['figure.facecolor'] = 'white'


def plot_early_stages(ax, planet, x, y, plot_evolve=True, **kwargs):
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
            ys = getattr(stage, y)[0]
        else:   
            ys = stage.datadicts[0][y]

        if x == 'timestep':
            xs = np.arange(len(ys)) + offset
            offset += len(xs)
        else:
            xs = getattr(stage, x)

        ax.plot(xs, ys, '.', label=stage.name, **kwargs)
    
        add_to_legend(ax, label=stage.name, c=stage.color) 

    ax.set_xlabel(x)
    ax.set_ylabel(y)


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
