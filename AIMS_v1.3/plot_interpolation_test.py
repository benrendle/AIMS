#!/usr/bin/env python
# coding: utf-8
# $Id: plot_interpolation_test.py
# Author: Daniel R. Reese <daniel.reese@obspm.fr>
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

.. warning::
   This plot utility only works with 3 dimensional grids (incl. the age dimension).
"""

__docformat__ = 'restructuredtext'

import sys
import dill
import math
import numpy as np
import matplotlib.pyplot  as plt
import matplotlib.lines   as mlines
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

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
      - "avg": takes the root mean-square value.

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
            print("ERROR: unrecognised type: "+str(tpe))
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
    ax.set_xlabel(titles[0])
    ax.set_ylabel(titles[1])
    ax.set_zlabel(r"$\log_{10}$(%s. error)"%(tpe))
    if (title is not None): ax.set_title(title)

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
    plt.xlabel(titles[0])
    plt.ylabel(titles[1])
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
        plt.title(r"Error at $(M,Z)=(%f,%f)$"%(results_age[0][pos][0,0],results_age[0][pos][0,1]))
        plt.ylabel(r"Error(radial)")
        plt.yscale('log')
        plt.xlim(xlim)
        plt.legend()

    plt.subplot(2,1,2)
    if (all_nan(results_age[0][pos][:,ndim+3:ndim+6]) and \
        all_nan(results_age[1][pos][:,ndim+3:ndim+6])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for j in range(2):
            for i in range(3):
                plt.plot(results_age[j][pos][:,ndim-1],results_age[j][pos][:,ndim+3+i],
                     style[i+3*j],label=labels[i+3*j])
        plt.ylabel(r"Error(non radial)")
        plt.yscale('log')
        plt.xlabel(titles[ndim-1])
        plt.xlim(xlim)
        plt.legend()
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
    #plt.semilogy(grid[ndx1,0],grid[ndx1,1],'bo', picker=5)
    #plt.semilogy(grid[ndx2,0],grid[ndx2,1],'yo')
    plt.plot(grid[ndx1,0],grid[ndx1,1],'bo', picker=5)
    plt.plot(grid[ndx2,0],grid[ndx2,1],'yo')
    plt.triplot(grid[ndx2,0],grid[ndx2,1],tessellation.simplices.copy())
    plt.xlim((xmin,xmax))
    plt.ylim((ymin,ymax))
    plt.xlabel(titles[0])
    plt.ylabel(titles[1])
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
        plt.title(r"Error at $(M,Z)=(%f,%f)$"%(results_track[pos][0,0],results_track[pos][0,1]))
        plt.ylabel(r"Error(radial)")
        plt.yscale('log')
        plt.xlim(xlim)
        plt.legend()

    plt.subplot(2,1,2)
    if (all_nan(results_track[pos][:,ndim+3:ndim+6])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for i in range(3):
            plt.plot(results_track[pos][:,ndim-1],results_track[pos][:,ndim+3+i],
                     style[i],label=labels[i])
        plt.ylabel(r"Error(non radial)")
        plt.yscale('log')
        plt.xlabel(titles[ndim-1])
        plt.xlim(xlim)
        plt.legend()
    plt.show()

if __name__ == "__main__":
    """
    Plot interpolation tests.

    .. note::
      The user may want to edit this.
    """

    if (len(sys.argv) < 2):
        print("Usage: plot_interpolation_test.py data_file")
        sys.exit(1)

    filename = sys.argv[1]

    input_data = open(filename,"rb")
    [ndim, nglb, titles, grid, ndx1, ndx2, tessellation, results_age1, \
        results_age2, results_track] = dill.load(input_data)

    if (ndim != 3):
        print("I'm sorry, but I can only plot handle 3D grids.")
        sys.exit(1)

    results_age = [results_age1, results_age2]
    input_data.close()

    plot3D(results_age1,0,tpe="max",title="Max radial error (nincr = 1)",truncate=1)
    plot3D(results_age2,0,tpe="max",title="Max radial error (nincr = 2)")
    plot3D(results_track,0,tpe="max",title="Max radial error (struct)")

    plot3D(results_age1,1,tpe="avg",title="Avg radial error (nincr = 1)",truncate=1)
    plot3D(results_age2,1,tpe="avg",title="Avg radial error (nincr = 2)")
    plot3D(results_track,1,tpe="avg",title="Avg radial error (struct)")

    plot3D(results_age1,2,tpe="avg",title="numax radial error (nincr = 1)",truncate=1)
    plot3D(results_age2,2,tpe="avg",title="numax radial error (nincr = 2)")
    plot3D(results_track,2,tpe="avg",title="numax radial error (struct)")

    plot_grid(grid)
    plot_partition_tessellation(grid, ndx1, ndx2, tessellation)

    plt.show()
