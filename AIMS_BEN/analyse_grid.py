#!/usr/bin/env python

"""
A simple utility for analysing binary AIMS grids
"""

import sys
import dill
import numpy as np
import matplotlib.pyplot as plt

iage = 0  

def print_number_of_models_in_tracks(grid):
    for i in range(len(grid.tracks)):
        print("Track %d contains %d models"%(i,len(grid.tracks[i].models)))

def print_number_of_modes_in_models(grid):
    for i in range(len(grid.tracks)):
        for j in range(len(grid.tracks[i].models)):
            print("Track %d %s %d"%(i,grid.tracks[i].models[j].name, \
                                 grid.tracks[i].models[j].modes.shape[0]))

def plot_number_of_modes(grid):
    nmodes = []
    for track in grid.tracks:
        for model in track.models:
            nmodes.append(model.modes.shape[0])
    nmax = max(nmodes)
    x = np.array(range(nmax+1))
    y = np.zeros((nmax+1,),dtype=np.int)
    for i in x: y[i] = nmodes.count(i)
    plt.plot(x,y,"b-")
    plt.plot(x,y,"bo")
    plt.show()

def find_range(grid,i):
    param_min = grid.tracks[0].params[i]
    param_max = grid.tracks[0].params[i]

    for track in grid.tracks:
        if (param_min > track.params[i]): param_min = track.params[i] 
        if (param_max < track.params[i]): param_max = track.params[i] 

    print("Range on parameter %s: %8.5e to %8.5e"%(grid.grid_params[i],param_min,param_max))

def duplicate_ages(track):
    """
    Check to see if you track contains models with duplicate ages.

    :return: ``True`` if there are duplicate age(s)
    :rtype: boolean

    .. warning::
        This method should only be applied after the track has been
        sorted.
    """

    return any(track.models[i].glb[iage] == track.models[i+1].glb[iage] for i in range(len(track.models)-1))

if __name__ == "__main__":

    # check number of arguments
    assert (len(sys.argv) > 1), "Usage: test_grid.py binary_grid_file"

    grid = dill.load(open(sys.argv[1],"rb"))
    ntracks = len(grid.tracks)

    print("Model postfix:    "+str(grid.postfix))
    print("Grid parameters:  "+str(grid.grid_params))
    print("User parameters:  "+str(grid.user_params))
    print("Number of dims.:  "+str(grid.ndim+1))
    print("Number of tracks: "+str(ntracks))
    print("Number of models: "+str(sum([len(track.models) for track in grid.tracks])))
    print("Number of modes:  "+str(sum([model.modes.shape[0] for model in track.models for track in grid.tracks])))

    nmax = max([len(track.models) for track in grid.tracks])
    print("Number of models in smallest track: "+str(min([len(track.models) for track in grid.tracks])))
    print("Number of models in largest track:  "+str(nmax))

    for track in grid.tracks:
        if (duplicate_ages(track)):
            print("ERROR: the track %s = %s"%(track.grid_params,track.params))
            print("       has models with the same age.  Please remove")
            print("       duplicate models.")
            for i in range(len(track.models)):
                print("Model[%d]: %e"%(i,track.models[i].glb[iage]))

    for i in range(len(grid.grid_params)):
        find_range(grid,i)

    hist = np.zeros((nmax+1,))
    for track in grid.tracks: hist[len(track.models)] += 1
    plt.plot(list(range(nmax+1)),hist,"b-")
    plt.show()

    #plot_number_of_modes(grid)
    #print_number_of_models_in_tracks(grid)
    #print_number_of_modes_in_models(grid)
