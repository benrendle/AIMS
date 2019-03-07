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
import constants
import pandas as pd
import sys
# from pandas import DataFrame, read_csv

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
plt.rcParams["font.family"] = "serif"

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
    for i in xrange(n):
        start = truncate
        stop  = results[i].shape[0] - truncate
        if (tpe == "max"):
            value = np.nanmax(results[i][start:stop,ndim+error_ndx])
        elif (tpe == "avg"):
            num = den = 0.0
            for j in xrange(start,stop):
                if (not math.isnan(results[i][j,ndim+error_ndx])):
                    num += results[i][j,ndim+error_ndx]**2
                    den += 1.0
            if (den > 0.0):
                value = math.sqrt(num/den)
            else:
                value = np.nan
        else:
            print "ERROR: unrecognised type: ",tpe
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
    ax.set_xlabel(titles[0],fontsize=20)
    ax.set_ylabel(titles[1],fontsize=20)
    ax.set_zlabel(r"$\log_{10}$(%s. error)"%(tpe),fontsize=20)
    ax.xaxis.labelpad = 20
    ax.yaxis.labelpad = 20
    ax.zaxis.labelpad = 20
    ax.tick_params(labelsize=15)
    if (title is not None): ax.set_title(title,fontsize=20)

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
        print "Only able to plot the tessellation in two dimensions."
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
    # plt.semilogy(grid[:,0],10**grid[:,1],'bo', picker=5)
    plt.plot(grid[:,0],grid[:,1],'bo', picker=5)
    plt.xlim((xmin,xmax))
    plt.ylim((ymin,ymax))
    plt.xlabel(titles[0],fontsize=20)
    plt.ylabel(titles[1],fontsize=20)
    plt.tick_params(labelsize=15)
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
        print pos
        print results_age[0][0]
        print results_age[1][0]
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
        print "Cowardly refusing to plot nan's"
        return

    style = ["b:","b-","bv","r:","r-","rv"]
    labels = ["Max. (incr=1)", "Avg. (incr=1)", r"$\nu_{\mathrm{max}}$ avg (incr=1)",
              "Max. (incr=2)", "Avg. (incr=2)", r"$\nu_{\mathrm{max}}$ avg (incr=2)"]

    xmin = np.nanmin(results_age[0][pos][:,ndim-1])
    xmax = np.nanmax(results_age[1][pos][:,ndim-1])
    xlim = (xmin,xmax+0.1*(xmax-xmin))

    plt.figure()
    ax = plt.subplot(2,1,1)
    if (all_nan(results_age[0][pos][:,ndim:ndim+3]) and \
        all_nan(results_age[1][pos][:,ndim:ndim+3])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for j in xrange(2):
            for i in xrange(3):
                plt.plot(results_age[j][pos][:,ndim-1],results_age[j][pos][:,ndim+i],
                     style[i+3*j],label=labels[i+3*j])
        plt.title(r"Error at $(M,Z)=(%f,%f)$"%(results_age[0][pos][0,0],results_age[0][pos][0,1]),fontsize=20)
        plt.ylabel(r"Error(radial)",fontsize=20)
        plt.yscale('log')
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.xlim(xlim)
        plt.tick_params(labelsize=15)
        plt.legend(fontsize=15,loc=8,ncol=3)

    plt.subplot(2,1,2,sharex=ax)
    if (all_nan(results_age[0][pos][:,ndim+3:ndim+6]) and \
        all_nan(results_age[1][pos][:,ndim+3:ndim+6])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for j in xrange(2):
            for i in xrange(3):
                plt.plot(results_age[j][pos][:,ndim-1],results_age[j][pos][:,ndim+3+i],
                     style[i+3*j],label=labels[i+3*j])
        plt.ylabel(r"Error(non radial)",fontsize=20)
        plt.yscale('log')
        plt.xlabel(titles[ndim-1],fontsize=20)
        plt.xlim(xlim)
        plt.tick_params(labelsize=15)
    plt.show()

def plot_hrd():
    """
    Make an HRD showing the distribution of uncertainties as a function of evolution.
    """
    # a = [np.log10(0.01),np.log10(0.0125),np.log10(0.0156),np.log10(0.0194),np.log10(0.0240)] #KIC5786154
    # a = [np.log10(0.0090),np.log10(0.0072),np.log10(0.0057),np.log10(0.0046),np.log10(0.0036)] #KIC8430105
    # a = [np.log10(0.0125),np.log10(0.0100),np.log10(0.0080),np.log10(0.0064),np.log10(0.0051)] #KIC9970396
    # a = [np.log10(0.0156),np.log10(0.0125),np.log10(0.0100),np.log10(0.0080),np.log10(0.0064)] #KIC7037405
    # a = [np.log10(0.0194),np.log10(0.0240),np.log10(0.0297),np.log10(0.0365),np.log10(0.0156)] #KIC8410637
    # a = [np.log10(0.0267),np.log10(0.0297),np.log10(0.0329),np.log10(0.0365),np.log10(0.0404),np.log10(0.0446),np.log10(0.0492)] #NGC6791
    a = [np.log10(0.0090),np.log10(0.0112),np.log10(0.0140),np.log10(0.0174),np.log10(0.0216),np.log10(0.0267),np.log10(0.0329)] #NGC6819
    for k in a:
        idx = []
        fig = plt.figure()
        for i in xrange(len(results_age[0])):
            for j in xrange(len(results_age[0][i][:,0])):
                if results_age[0][i][j,ndim-2] == k:
                    idx = np.append(idx,results_age[0][i][j,:])

        plt.scatter(idx[20::ndim+nglb+6+5],np.log10(idx[21::ndim+nglb+6+5]/3.828e33),c=np.log10(idx[ndim::ndim+nglb+6+5]),label=r'Z = %s '%(10**idx[1]))
        cb = plt.colorbar()
        cax = cb.ax
        cax.text(3.5,0.7,r"$\log_{10}$(Max. error)",rotation=270,fontsize=20)
        cax.tick_params(labelsize=15)
        plt.gca().invert_xaxis()
        plt.legend(loc=3)
        # plt.savefig('/home/bmr135/AIMS/AIMS_BEN/KIC5786154/HRD_KIC5786154_Z'+str(10**idx[1])+'.png')
        # plt.show()

def plot_MXc():
    """
    Make an HRD showing the distribution of uncertainties as a function of evolution.
    """
    a = [np.log10(0.03),np.log10(0.0175),np.log10(0.01),np.log10(0.0057),np.log10(0.0032)]
    for k in a:
        idx = []
        fig = plt.figure()
        for i in xrange(len(results_age[0])):
            for j in xrange(len(results_age[0][i][:,0])):
                if results_age[0][i][j,ndim-2] == k:
                    idx = np.append(idx,results_age[0][i][j,:])
        # print idx[::ndim+nglb+6+4]
        plt.scatter(idx[22::ndim+nglb+6+4],idx[0::ndim+nglb+6+4],c=np.log10(idx[ndim::ndim+nglb+6+4]),label=r'Z = %s '%(10**idx[1]))
        cb = plt.colorbar()
        cax = cb.ax
        cax.text(3.5,0.7,r"$\log_{10}$(Max. error)",rotation=270,fontsize=20)
        cax.tick_params(labelsize=15)
        plt.legend(loc=1)
        plt.xlabel(r'X$_{\rm{c}}$ - Inverted')
        plt.ylabel(r'Mass')
        plt.gca().invert_xaxis()
        # plt.savefig('/media/bmr135/SAMSUNG/AIMS-interp-testing2/MXc_MS_Z'+str(10**idx[1])+'_v3.9.png')
        # plt.show()

def plot_MmHe():
    """
    Make an HRD showing the distribution of uncertainties as a function of evolution.
    """
    # a = [np.log10(0.03),np.log10(0.0175),np.log10(0.01),np.log10(0.0057),np.log10(0.0032)]
    a = [np.log10(0.0090),np.log10(0.0112),np.log10(0.0140),np.log10(0.0174),np.log10(0.0216),np.log10(0.0267),np.log10(0.0329)] #NGC6819
    for k in a:
        idx = []
        fig = plt.figure()
        for i in xrange(len(results_age[0])):
            for j in xrange(len(results_age[0][i][:,0])):
                if results_age[0][i][j,ndim-2] == k:
                    idx = np.append(idx,results_age[0][i][j,:])
        # print idx[::ndim+nglb+6+4]
        plt.scatter(idx[23::ndim+nglb+6+5],idx[0::ndim+nglb+6+5],c=np.log10(idx[ndim::ndim+nglb+6+5]),label=r'Z = %s '%(10**idx[1]))
        cb = plt.colorbar()
        cax = cb.ax
        cax.text(3.5,0.7,r"$\log_{10}$(Max. error)",rotation=270,fontsize=20)
        cax.tick_params(labelsize=15)
        plt.legend(loc=1)
        plt.xlabel(r'Mass He-core')
        plt.ylabel(r'Mass')
        # plt.gca().invert_xaxis()
        # plt.savefig('/media/bmr135/SAMSUNG/AIMS-interp-testing2/MmHe_RGB_Z'+str(10**idx[1])+'_v3.9.png')
        # plt.show()

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
        print "Only able to plot the tessellation in two dimensions."
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
    # plt.plot(grid[ndx2,0],grid[ndx2,1],'yo')
    plt.triplot(grid[ndx2,0],grid[ndx2,1],tessellation.simplices.copy())
    plt.xlim((xmin,xmax))
    plt.ylim((ymin,ymax))
    plt.xlabel(titles[0],fontsize=20)
    plt.ylabel(titles[1],fontsize=20)
    plt.tick_params(labelsize=15)
    # fig.canvas.mpl_connect('pick_event', onpick_track)

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
        print "Cowardly refusing to plot nan's"
        return

    style = ["b:","b-","bv"]
    labels = ["Max.", "Avg.", r"$\nu_{\mathrm{max}}$ avg"]
    xmin = np.nanmin(results_track[pos][:,ndim-1])
    xmax = np.nanmax(results_track[pos][:,ndim-1])
    xlim = (xmin,xmax+0.3*(xmax-xmin))

    plt.figure()
    ax = plt.subplot(2,1,1)
    if (all_nan(results_track[pos][:,ndim:ndim+3])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for i in xrange(3):
            plt.plot(results_track[pos][:,ndim-1],results_track[pos][:,ndim+i],
                     style[i],label=labels[i])
        plt.title(r"Error at $(M,Z)=(%f,%f)$"%(results_track[pos][0,0],results_track[pos][0,1]),fontsize=20)
        plt.ylabel(r"Error(radial)",fontsize=20)
        plt.yscale('log')
        plt.xlim(xlim)
        plt.tick_params(labelsize=15)
        plt.legend()

    plt.subplot(2,1,2)
    if (all_nan(results_track[pos][:,ndim+3:ndim+6])):
        plt.plot([0,1,1,0],[0,1,0,1],"k-")
    else:
        for i in xrange(3):
            plt.plot(results_track[pos][:,ndim-1],results_track[pos][:,ndim+3+i],
                     style[i],label=labels[i])
        plt.ylabel(r"Error(non radial)",fontsize=20)
        plt.yscale('log')
        plt.xlabel(titles[ndim-1],fontsize=20)
        plt.xlim(xlim)
        plt.tick_params(labelsize=15)
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

    error_ndx:
    [0] = maximum error on the radial modes
    [1] = RMS error on the radial modes
    [2] = RMS error on the radial modes near :math:`\\nu_{\mathrm{max}}`
    [3] = maximum error on the non radial modes
    [4] = RMS error on the non radial modes
    [5] = RMS error on the non radial modes near :math:`\\nu_{\mathrm{max}}`
    [6] = Age
    [7] = Mass
    [8] = Teff
    [9] = Z0
    [10] = X0
    [11] = -
    [12] = Xc
    [13] = Period Spacing
    [14] = Reference Frequency
    [15] = Radius
    [16] = Luminosity

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
    for i in xrange(n):
        start = truncate
        stop  = results[i].shape[0] - truncate
        if (tpe == "max"):
            value = np.nanmax(results[i][start:stop,ndim+error_ndx])
        elif (tpe == "avg"):
            num = den = 0.0
            for j in xrange(start,stop):
                if (not math.isnan(results[i][j,ndim+error_ndx])):
                    num += results[i][j,ndim+error_ndx]**2
                    den += 1.0
            if (den > 0.0):
                value = math.sqrt(num/den)
            else:
                value = np.nan
        else:
            print "ERROR: unrecognised type: ",tpe
            sys.exit(1)
        if (value > 0.0):
	    # print value
            z1.append(value)	# math.log10
            z2.append(math.log10(value))
            x.append(results[i][0,0])
            y.append(results[i][0,1])

    x = np.array(x,dtype = np.float64)
    y = np.array(y,dtype = np.float64)
    z1 = np.array(z1,dtype = np.float64)
    z2 = np.array(z2,dtype = np.float64)
    xi, yi = np.linspace(x.min(),x.max(),200), np.linspace(y.min(),y.max(),200)
    xi, yi = np.meshgrid(xi,yi)
    rbf = interp.Rbf(x,y,z2,function='linear')
    zi = rbf(xi,yi)
    j = 0
    kkk = 0
    plt.figure()
    cont = plt.contourf(xi,yi,zi,100)
    plt.close()
    k = cont.levels
    # print k

    fig = plt.figure(figsize=(6,2))
    plt.subplots_adjust(left=0.15, right=0.94, top=0.95, bottom=0.25)
    df = pd.DataFrame()
    df['x'] = x
    df['y'] = y
    df['z2'] = z2
    df = df.sort_values(['z2'])
    # kk = np.array(k)

    # print kk

    # print kk[-18], kk[-22]
    # cont.collections[-18].set_color('m')
    # cont.collections[-22].set_color('k')
    # circle1 = plt.Circle((1.46, np.log10(0.0046)), .0075, color='k', fill=False,linewidth=3)
    # plt.gcf().gca().add_artist(circle1)
    plt.scatter(df['x'], df['y'], c=df['z2'], cmap=colormaps.parula, s=75)
    cb = plt.colorbar()
    # plt.xlabel(titles[0],fontsize=10)
    # plt.xlabel(r'Radius $R/R_{\odot}$',fontsize=10)
    plt.ylabel(titles[1],fontsize=10)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    # plt.xlim(0.75,2.25)
    # plt.ylim(-2.54,-1.46)
    a = min(z2)
    b = max(z2)
    d = a-b
    # print a, b, d
    f1 = a--1.7447
    g1 = f1/d
    f2 = a--1.8538
    g2 = f2/d

    cax = cb.ax
    # cax.hlines(g1,0,1,colors='m',linewidth=2)
    # cax.hlines(g2,0,1,colors='k',linewidth=2)
    #cax.text(3.5,0.7,r"Percentage Error",rotation=270,fontsize=20)	# %s %(tpe)
    cax.text(6.5,0.85,r"$\log_{10}$(%s. error)"%(tpe),rotation=270,fontsize=10)
    cax.tick_params(labelsize=8)
    # if (title is not None): plt.title(title,fontsize=15)
    # fig.savefig('rad_interp.pdf', bbox_inches='tight')
    # plt.show()
    # sys.exit()
    # plt.figure()
    # plt.hist(z2,bins=50)
    # plt.xlabel(r'$\sigma_{\nu_{\rm{i}}}$')

    m = n = 0
    for i in z2:
    	if i <= -1.7447:
    	    m += 1
    	else:
    	    n += 1
    pp = float(m)/len(z2)
    # print 'Pass percentage',p,': ',pp*100,'%'



if __name__ == "__main__":
    """
    Plot interpolation tests.

    .. note::
      The user may want to edit this.
    """


    # input_data = open('/home/bmr135/AIMS/AIMS_BEN/interp_MS_v3.5.8',"r")
    # [ndim, nglb, titles, grid, ndx1, ndx2, tessellation, results_age1, \
    #     results_age2, results_track] = dill.load(input_data)
    # results_age_MS = [results_age1, results_age2]
    # input_data.close()
    # input_data = open('/home/bmr135/AIMS/AIMS_BEN/interp_RGB_v3.7',"r")
    # [ndim, nglb, titles, grid, ndx1, ndx2, tessellation, results_age1, \
    #     results_age2, results_track] = dill.load(input_data)
    # results_age_RGB = [results_age1, results_age2]
    # input_data.close()
    #
    # a = np.log10(0.01)
    # fig = plt.figure()
    # plt.scatter(results_age_MS[0][132][:,20],np.log10( results_age_MS[0][132][:,21]/3.828e33))
    # plt.scatter(results_age_RGB[0][132][:,20],np.log10( results_age_RGB[0][132][:,21]/3.828e33))
    # plt.gca().invert_xaxis()
    # plt.legend(loc=3)
    # # plt.savefig('/media/bmr135/SAMSUNG/AIMS-interp-testing2/HRD_RGB_Z'+str(10**idx[1])+'.png')
    # plt.show()
    # #
    # #
    # #
    # #
    # #
    # sys.exit()

    assert (len(sys.argv) > 1), "Usage: plot_interpolation_test.py data_file"

    filename = sys.argv[1]

    input_data = open(filename,"r")
    [ndim, nglb, titles, grid, ndx1, ndx2, tessellation, results_age1, \
        results_age2, results_track] = dill.load(input_data)
    results_age = [results_age1, results_age2]
    input_data.close()

    # plot_hrd()
    # plot_MXc()
    plot_MmHe()

    # print results_track[]

    # print grid # mass/met
    # print results_age1[0][0]
    #print ndx2
    # df = pd.DataFrame(list(results_track[0][:]))
    #print df
    # df1 = pd.DataFrame(list(results_age1[0][:]))
    #print df1

    # surface2D(1,results_age1,0,tpe="max",title="Max radial error (nincr = 1)",truncate=1)
    # surface2D(1,results_age2,0,tpe="max",title="Max radial error (nincr = 2)",truncate=1)
    # surface2D(1,results_track,0,tpe="max",title="Max. radial error between tracks",truncate=1)
    ''' '''
    # surface2D(1,results_age1,1,tpe="max",title="Avg radial error (nincr = 1)",truncate=1)
    # surface2D(1,results_age2,1,tpe="max",title="Avg radial error (nincr = 2)",truncate=1)
    # surface2D(1,results_track,1,tpe="max",title="Avg radial error (nincr = struct)",truncate=1)

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

    # plot_grid(grid)
    # plot_partition_tessellation(grid, ndx1, ndx2, tessellation)

    plt.show()
