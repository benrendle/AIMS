#!/usr/bin/env python
# $Id: plot_interpolation_test.py
# Author: Daniel R. Reese <dreese@bison.ph.bham.ac.uk>
# Copyright (C) Daniel R. Reese and contributors
# Copyright license: GNU GPL v3.0
#
#   AIMS is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with AIMS.  If not, see <http://www.gnu.org/licenses/>.
#

"""
An interactive utility which plots various forms of interpolation error,
stored in a binary file produced by :py:func:`AIMS.test_interpolation`.
It specifically tests the errors from two types of interpolation:

- age interpolation: this is interpolation along a given evolutionary track
- track interpolation: this is interpolation between different evolutionary
  tracks

This utility allows various types of plots:

- 3D plots of interpolation errors as a function of grid structural parameters
- 2D slices which show interpolation errors as a function of age for a
  given evolutionary track
- interactive plots which allow you to select 2D slices

.. note::
  Interpolation errors for models in a given evolutionary track are typically
  stored in arrays as follows:

  - ``result[model_number,ndim+0]`` = maximum error on the radial modes
  - ``result[model_number,ndim+1]`` = RMS error on the radial modes
  - ``result[model_number,ndim+2]`` = RMS error on the radial modes near
    :math:`\\nu_{\mathrm{max}}`
  - ``result[model_number,ndim+3]`` = maximum error on the non radial modes
  - ``result[model_number,ndim+4]`` = RMS error on the non radial modes
  - ``result[model_number,ndim+5]`` = RMS error on the non radial modes near
    :math:`\\nu_{\mathrm{max}}`
  - ``result[model_number,ndim+6+[0:nglb]]`` = errors on the global parameters

  where:

  - ``result`` = the array which containts the interpolation errors
  - ``model_numer`` = an index which represents the model (not necessarily the
    number of the model along the evolutionary track
  - ``ndim`` = the number of dimensions in the grid (including age)
  - ``nglb`` = the number of global parameters for stellar models in the grid
"""

__docformat__ = 'restructuredtext'

import sys
import dill
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as interp
import colormaps
# import pandas as pd
import sys
# from pandas import DataFrame, read_csv

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

# global variables
ndim = 0
""" number of dimension in grid (including age) """

nglb = 0
""" number of global parameters """

titles = None
""" the grid quantities, which will serve as axis labels """

results_age    = None
"""  list which contains the arrays with the errors from age interpolation """

results_track = None
"""  list which contains the arrays with the errors from track interpolation """

def all_nan(array):
    """
    Test to see if all of the elements of an array are nan's.

    :param array: array in which we're checking to see if all elements are nan's.
    :type array: np.array

    :return: ``True`` if all the elements of ``array`` are nan's, and ``False``
       otherwise.
    :rtype: boolean
    """
    return np.all(np.isnan(array))

def plot3D(results,error_ndx,tpe="max",title=None,truncate=0):
    """
    Create 3D plot showing the error as a function of the two
    first grid parameters.


    :param results: list of 2D arrays which contain various types of errors as
      a function of the model number along a given evolutionary track.
    :param error_ndx: value which specifies the type of error to be plotted.
    :param tpe: specifies how to combine errors along the
      evolutionary track.  Options include:

      - "max": corresponds to taking the maximum value.
      - "avg": takes the mean-square value.

    :param title: the title of the plot
    :param truncate: (default = 0): specifies how many models should be omitted
      on both ends of the track.  This is useful for comparing results from
      tests involing different sizes of increments.

    :type results: list of np.arrays
    :type error_ndx: int
    :type tpe: string
    :type title: string
    :type truncate: int

    .. note::
      See above introductory description for a more detailled description
      of the indices which intervene in the 2D arrays contained in ``results``
      and of the relevant values for ``error_ndx``.
    """

    n = len(results)
    x = []
    y = []
    z = []
    for i in range(n):
        start = truncate
        stop  = results[i].shape[0] - truncate
        if (tpe == "max"):
            value = np.nanmax(results[i][start:stop,ndim+error_ndx])
        elif (tpe == "avg"):
            num = den = 0.0
            for j in range(start,stop):
                if (not math.isnan(results[i][j,ndim+error_ndx])):
                    num += results[i][j,ndim+error_ndx]**2
                    den += 1.0
            if (den > 0.0):
                value = math.sqrt(num/den)
            else:
                value = np.nan
        else:
            print("ERROR: unrecognised type: ",tpe)
            sys.exit(1)
        if (value > 0.0):
            z.append(math.log10(value))
            x.append(results[i][0,0])
            y.append(results[i][0,1])

    x = np.array(x,dtype = np.float64)
    y = np.array(y,dtype = np.float64)
    z = np.array(z,dtype = np.float64)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x,y,z,cmap=cm.jet,linewidth=0.2)
    ax.set_xlabel(titles[0],fontsize=40)
    ax.set_ylabel(titles[1],fontsize=40)
    ax.set_zlabel(r"$\log_{10}$(%s. error)"%(tpe),fontsize=40)
    ax.xaxis.labelpad = 20
    ax.yaxis.labelpad = 20
    ax.zaxis.labelpad = 20
    ax.tick_params(labelsize=30)
    if (title is not None): ax.set_title(title,fontsize=40)

def plot_grid(grid):
    """
    Make an interactive plot of the grid.  Clicking on the blue dots
    will produce a 2D slice showing age interpolation errors for
    the associated evolutionary track.

    :param grid: array containing basic grid parameters (excluding age)
    :type grid: np.array

    .. warning::
      This only works for two-dimensional grids (excluding the age
      dimension).
    """

    # remember ndim includes the age dimension.
    if (ndim != 3):
        print("Only able to plot the tessellation in two dimensions.")
        return

    # find bounds:
    xmin = np.nanmin(grid[:,0])
    xmax = np.nanmax(grid[:,0])
    ymin = np.nanmin(grid[:,1])
    ymax = np.nanmax(grid[:,1])
    dx = xmax-xmin
    dy = ymax/ymin
    xmin -= dx*0.03
    xmax += dx*0.03
    ymin /= dy**0.05
    ymax *= dy**0.05

    fig = plt.figure()
    #plt.semilogy(grid[:,0],grid[:,1],'bo', picker=5)
    plt.plot(grid[:,0],grid[:,1],'bo', picker=5)
    plt.xlim((xmin,xmax))
    plt.ylim((ymin,ymax))
    plt.xlabel(titles[0],fontsize=40)
    plt.ylabel(titles[1],fontsize=40)
    plt.tick_params(labelsize=30)
    fig.canvas.mpl_connect('pick_event', onpick_age)

def onpick_age(event):
    """
    Event catcher for the grid plot (which shows the positions of
    the evolutionary tracks as a function of the grid parameters,
    excluding age).

    Parameters:

    :param event: event caught by the grid plot.
    """
    for pos in event.ind:
        plot_slice_age(pos)

def plot_slice_age(pos):
    """
    Plot age interpolation error as a function of age for a
    given track.

    :param pos: index of the relevant track.
    :type pos: int

    .. note::
      This `pos` index applies to results_age, i.e., it is based
      on the original track indices.
    """

    if (all_nan(results_age[0][pos][:,ndim:ndim+6]) and \
        all_nan(results_age[1][pos][:,ndim:ndim+6])):
        print("Cowardly refusing to plot nan's")
        return

    style = ["b:","b-","bv","r:","r-","rv"]
    labels = ["Max. (incr=1)", "Avg. (incr=1)", r"$\nu_{\mathrm{max}}$ avg (incr=1)",
              "Max. (incr=2)", "Avg. (incr=2)", r"$\nu_{\mathrm{max}}$ avg (incr=2)"]

    xmin = np.nanmin(results_age[0][pos][:,ndim-1])
    xmax = np.nanmax(results_age[1][pos][:,ndim-1])
    xlim = (xmin,xmax+0.5*(xmax-xmin))

    plt.figure()
    plt.subplot(2,1,1)
    if (all_nan(results_age[0][pos][:,ndim:ndim+3]) and \
        all_nan(results_age[1][pos][:,ndim:ndim+3])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for j in range(2):
            for i in range(3):
                plt.plot(results_age[j][pos][:,ndim-1],results_age[j][pos][:,ndim+i],
                     style[i+3*j],label=labels[i+3*j])
        plt.title(r"Error at $(M,Z)=(%f,%f)$"%(results_age[0][pos][0,0],results_age[0][pos][0,1]),fontsize=40)
        plt.ylabel(r"Error(radial)",fontsize=40)
        plt.yscale('log')
        plt.xlim(xlim)
        plt.tick_params(labelsize=30)
        plt.legend(fontsize=25)

    plt.subplot(2,1,2)
    if (all_nan(results_age[0][pos][:,ndim+3:ndim+6]) and \
        all_nan(results_age[1][pos][:,ndim+3:ndim+6])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for j in range(2):
            for i in range(3):
                plt.plot(results_age[j][pos][:,ndim-1],results_age[j][pos][:,ndim+3+i],
                     style[i+3*j],label=labels[i+3*j])
        plt.ylabel(r"Error(non radial)",fontsize=40)
        plt.yscale('log')
        plt.xlabel(titles[ndim-1],fontsize=40)
        plt.xlim(xlim)
        plt.tick_params(labelsize=30)
        plt.legend(fontsize=25)	##################################################################################################################
    plt.show()

def plot_partition_tessellation(grid, ndx1, ndx2, tessellation):
    """
    Make an interactive tessellation plot based on the supplied partition
    on the grid. Clicking on the blue dots will produce a 2D slice showing
    track interpolation errors for the associated evolutionary track.

    :param grid: array containing basic grid parameters (excluding age)
    :param ndx1: list with the indices of the first part of the partition.
    :param ndx2: list with the indices of the second part of the partition.
    :param tessellation: grid tessellation associated with ``ndx2``

    :type grid: np.array
    :type ndx1: list of int
    :type ndx2: list of int

    .. warning::
      This only works for two-dimensional tessellations.
    """

    # remember ndim includes the age dimension which is not included
    # in the tessellation:
    if (ndim != 3):
        print("Only able to plot the tessellation in two dimensions.")
        return

    # find bounds:
    xmin = np.nanmin(grid[:,0])
    xmax = np.nanmax(grid[:,0])
    ymin = np.nanmin(grid[:,1])
    ymax = np.nanmax(grid[:,1])
    dx = xmax-xmin
    dy = ymax/ymin
    xmin -= dx*0.03
    xmax += dx*0.03
    ymin /= dy**0.05
    ymax *= dy**0.05

    fig = plt.figure()
    plt.plot(grid[ndx1,0],grid[ndx1,1],'bo', picker=5)
    plt.plot(grid[ndx2,0],grid[ndx2,1],'yo')
    plt.triplot(grid[ndx2,0],grid[ndx2,1],tessellation.simplices.copy())
    plt.xlim((xmin,xmax))
    plt.ylim((ymin,ymax))
    plt.xlabel(titles[0],fontsize=40)
    plt.ylabel(titles[1],fontsize=40)
    plt.tick_params(labelsize=30)
    fig.canvas.mpl_connect('pick_event', onpick_track)

def onpick_track(event):
    """
    Event catcher for the partition tessellation plot (associated with
    tests of track interpolation).

    :param event: event caught by the partition tessellation plot.
    """
    for pos in event.ind:
        plot_slice_track(pos)

def plot_slice_track(pos):
    """
    Plot track interpolation error as a function of age for a
    given track.

    :param pos: index of the relevant track.
    :type pos: int

    .. note::
      This `pos` index applies to results_track, i.e., it is based
      on the indices deduced from the grid partition.
    """

    if (all_nan(results_track[pos][:,ndim:ndim+6])):
        print("Cowardly refusing to plot nan's")
        return

    style = ["b:","b-","bv"]
    labels = ["Max.", "Avg.", r"$\nu_{\mathrm{max}}$ avg"]
    xmin = np.nanmin(results_track[pos][:,ndim-1])
    xmax = np.nanmax(results_track[pos][:,ndim-1])
    xlim = (xmin,xmax+0.3*(xmax-xmin))

    plt.figure()
    plt.subplot(2,1,1)
    if (all_nan(results_track[pos][:,ndim:ndim+3])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for i in range(3):
            plt.plot(results_track[pos][:,ndim-1],results_track[pos][:,ndim+i],
                     style[i],label=labels[i])
        plt.title(r"Error at $(M,Z)=(%f,%f)$"%(results_track[pos][0,0],results_track[pos][0,1]),fontsize=40)
        plt.ylabel(r"Error(radial)",fontsize=40)
        plt.yscale('log')
        plt.xlim(xlim)
        plt.tick_params(labelsize=30)
        plt.legend()

    plt.subplot(2,1,2)
    if (all_nan(results_track[pos][:,ndim+3:ndim+6])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for i in range(3):
            plt.plot(results_track[pos][:,ndim-1],results_track[pos][:,ndim+3+i],
                     style[i],label=labels[i])
        plt.ylabel(r"Error(non radial)",fontsize=40)
        plt.yscale('log')
        plt.xlabel(titles[ndim-1],fontsize=40)
        plt.xlim(xlim)
        plt.tick_params(labelsize=30)
        plt.legend()
    plt.show()

def surface2D(p,results,error_ndx,tpe="max",title=None,truncate=0):
    """
    Create contour plots to demonstrate more clearly how the
    errors within the grid change as a function of the input
    parameters.

    :param results: list of 2D arrays which contain various types of errors as
      a function of the model number along a given evolutionary track.
    :param error_ndx: value which specifies the type of error to be plotted.
    :param tpe: specifies how to combine errors along the
      evolutionary track.  Options include:

      - "max": corresponds to taking the maximum value.
      - "avg": takes the mean-square value.

    :param title: the title of the plot
    :param truncate: (default = 0): specifies how many models should be omitted
      on both ends of the track.  This is useful for comparing results from
      tests involing different sizes of increments.

    :type results: list of np.arrays
    :type error_ndx: int
    :type tpe: string
    :type title: string
    :type truncate: int

    .. note::
      See above introductory description for a more detailled description
      of the indices which intervene in the 2D arrays contained in ``results``
      and of the relevant values for ``error_ndx``.
    """

    n = len(results)
    x = []
    y = []
    z1 = []
    z2 = []
    for i in range(n):
        start = truncate
        stop  = results[i].shape[0] - truncate
        if (tpe == "max"):
            value = np.nanmax(results[i][start:stop,ndim+error_ndx])
        elif (tpe == "avg"):
            num = den = 0.0
            for j in range(start,stop):
                if (not math.isnan(results[i][j,ndim+error_ndx])):
                    num += results[i][j,ndim+error_ndx]**2
                    den += 1.0
            if (den > 0.0):
                value = math.sqrt(num/den)
            else:
                value = np.nan
        else:
            print("ERROR: unrecognised type: ",tpe)
            sys.exit(1)
        if (value > 0.0):
            z1.append(value)	# math.log10
	    z2.append(math.log10(value))
            x.append(results[i][0,0])
            y.append(results[i][0,1])

    x = np.array(x,dtype = np.float64)
    y = np.array(y,dtype = np.float64)
    z1 = np.array(z1,dtype = np.float64)
    xi, yi = np.linspace(x.min(),x.max(),200), np.linspace(y.min(),y.max(),200)
    xi, yi = np.meshgrid(xi,yi)
    rbf = interp.Rbf(x,y,z1,function='linear')
    zi = rbf(xi,yi)
    j = 0
    kkk = 0
    fig = plt.figure()
    cont = plt.contourf(xi,yi,zi,100,cmap=colormaps.parula)
    k = cont.levels

    kk = np.array(k)

    print(kk)

#    try:
#        j = np.where(kk == -1.3979)[0][0]
#        cont.collections[j].set_color('m')

#    except:
#	pass

#    try:
#	if p == 1:
#	    j = 76
#	    cont.collections[j].set_color('m')
 #   except:
  #      pass

#    try:
 #       if p == 2:
#	    j = 39
#	    cont.collections[j].set_color('m')
 #   except:
#	pass

 #   try:
  #      if p == 3:
#	    j = 15
#	    cont.collections[j].set_color('m')
 #   except:
#	pass

#    try:
#        j = np.where(kk == -0.640797)[0][0]
#        cont.collections[j].set_color('k')
#	cont.collections[j].set_linestyle('dashed')
#    except:
#        pass

#    try:
#        j = np.where(kk == -1.775)[0][0]
#        cont.collections[j].set_color('m')
#	cont.collections[j].set_linestyle('dashed')
#    except:
#        pass

    plt.scatter(x, y, c=z1)
    cb = plt.colorbar()
    plt.xlabel(titles[0],fontsize=40)
    plt.ylabel(titles[1],fontsize=40)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    a = min(z1)
    b = max(z1)
    d = a-b
    f1 = a--1.3979
    g1 = f1/d
    f2 = a--0.640797
    g2 = f2/d
    cax = cb.ax
    cax.hlines(g1,0,1,colors='m',linewidth=2)
    cax.hlines(g2,0,1,colors='k',linewidth=2)
    cax.text(3.5,0.7,r"Percentage Error",rotation=270,fontsize=40)	# %s %(tpe)
    cax.tick_params(labelsize=30)
    if (title is not None): plt.title(title,fontsize=40)

    m = n = 0
    for i in z1:
	if i <= 3:
	    m += 1
	else:
	    n += 1
    pp = float(m)/len(z1)
    print('Pass percentage',p,': ',pp*100,'%')



if __name__ == "__main__":
    """
    Plot interpolation tests.

    .. note::
      The user may want to edit this.
    """

    assert (len(sys.argv) > 1), "Usage: plot_interpolation_test.py data_file"

    filename = sys.argv[1]

    input_data = open(filename,"r")
    [ndim, nglb, titles, grid, ndx1, ndx2, tessellation, results_age1, \
        results_age2, results_track] = dill.load(input_data)
    results_age = [results_age1, results_age2]
    input_data.close()

    #print ndx1

    #print ndx2
    # df = pd.DataFrame(list(results_track[0][:]))
    #print df
    # df1 = pd.DataFrame(list(results_age1[0][:]))
    #print df1

    surface2D(1,results_age1,0,tpe="max",title="Max radial error (nincr = 1)",truncate=1)
    #surface2D(results_age2,0,tpe="max",title="Max radial error (nincr = 2)",truncate=1)
    #surface2D(results_track,0,tpe="max",title="Max radial error (nincr = struct)",truncate=1)

    #surface2D(results_age1,1,tpe="max",title="Avg radial error (nincr = 1)",truncate=1)
    #surface2D(results_age2,1,tpe="max",title="Avg radial error (nincr = 2)",truncate=1)
    #surface2D(results_track,1,tpe="max",title="Avg radial error (nincr = struct)",truncate=1)

    #surface2D(1,results_age1,2,tpe="max",title="numax radial error (nincr = 1)",truncate=1)
    #surface2D(2,results_age2,2,tpe="max",title="numax radial error (nincr = 2)",truncate=1)
    #surface2D(3,results_track,2,tpe="max",title="numax radial error (nincr = struct)",truncate=1)

    #plot3D(results_age1,0,tpe="max",title="Max radial error (nincr = 1)",truncate=1)
    #plot3D(results_age2,0,tpe="max",title="Max radial error (nincr = 2)")
    #plot3D(results_track,0,tpe="max",title="Max radial error (struct)")

    #plot3D(results_age1,1,tpe="avg",title="Avg radial error (nincr = 1)",truncate=1)
    #plot3D(results_age2,1,tpe="avg",title="Avg radial error (nincr = 2)")
    #plot3D(results_track,1,tpe="avg",title="Avg radial error (struct)")

    #plot3D(results_age1,2,tpe="avg",title="numax radial error (nincr = 1)",truncate=1)
    #plot3D(results_age2,2,tpe="avg",title="numax radial error (nincr = 2)")
    #plot3D(results_track,2,tpe="avg",title="numax radial error (struct)")

    plot_grid(grid)
    #plot_partition_tessellation(grid, ndx1, ndx2, tessellation)

    plt.show()
