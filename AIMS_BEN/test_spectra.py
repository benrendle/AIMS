#!/usr/bin/env python

"""
A simple utility for analysing binary AIMS grids
"""

import sys
import dill
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

nprocesses = 4

def count_incomplete_spectra_in_track(track):
    imissing  = 0
    for model in track.models:
        for j in range(model.modes.shape[0]-1):
            if (model.modes['l'][j] != model.modes['l'][j+1]): continue
            if (model.modes['n'][j]+1 != model.modes['n'][j+1]):
                imissing += 1
                break
    return (imissing, len(track.models)-imissing)

def count_incomplete_spectra(grid):
    pool = Pool(processes = nprocesses)
    results = pool.map(count_incomplete_spectra_in_track,grid.tracks)
    pool.close()
    pool.join()

    results_bis = zip(*results)
    print("Number of incomplete spectra: %d"%(sum(results_bis[0])))
    print("Number of complete spectra:   %d"%(sum(results_bis[1])))

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

    if (i < 9):
        print("Range on parameter %s: %8.5e to %8.5e"%(grid.grid_params[i],param_min,param_max))
    else:
        print("Range on parameter %s: %8.5e to %8.5e"%(grid.user_params[i-9],param_min,param_max))

def find_range_bis(grid,i):
    param_min = grid.tracks[0].models[0].glb[i]
    param_max = grid.tracks[0].models[0].glb[i]

    for track in grid.tracks:
        for model in track.models:
            if (param_min > model.glb[i]): param_min = model.glb[i] 
            if (param_max < model.glb[i]): param_max = model.glb[i] 

    print("Range on parameter %d: %8.5e to %8.5e"%(i,param_min,param_max))

if __name__ == "__main__":

    # check number of arguments
    assert (len(sys.argv) > 1), "Usage: test_grid.py binary_grid_file"

    grid = dill.load(open(sys.argv[1],"rb"))
    ntracks = len(grid.tracks)

    print("Grid parameters:  "+str(grid.grid_params))
    print("User parameters:  "+str(grid.user_params))
    print("Number of dims.:  "+str(grid.ndim+1))
    print("Number of tracks: "+str(ntracks))
    print("Number of models: "+str(sum([len(track.models) for track in grid.tracks])))
    print("Number of modes:  "+str(sum([model.modes.shape[0] for model in track.models for track in grid.tracks])))
    count_incomplete_spectra(grid)
