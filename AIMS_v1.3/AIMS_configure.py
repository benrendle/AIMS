#!/usr/bin/env python
# coding: utf-8
# $Id: AIMS_configure.py
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

import math

#########################   Parallelisation   ##############################
# NOTE: this is currently implemented with multiprocessing, which duplicates
#       the memory in each process.  To be more memory efficient, turn off
#       parallelisation using the "parallel" parameter.
nprocesses  = 4      # number of processes (if running in parallel)
parallel    = False  # specifies whether to run in parallel

#########################   EMCEE control parameters   #####################
ntemps      = 10     # number of temperatures
nwalkers    = 250    # number of walkers (this number should be even)
nsteps0     = 200    # number of burn-in steps
nsteps      = 200    # number of steps
thin        = 10     # thinning parameter (1 out of thin steps will be kept ...)
thin_comb   = 100    # thinning parameter for output linear combinations of models
PT          = True   # use parallel tempering?

#########################   Initialisation   ###############################
tight_ball   = True  # initialise with a tight ball around best solution
max_iter     = 1000  # maximum number of iterations to find walker

# Ranges used around tight ball configuration for walkers.
# NOTES:
#   - if ranges are not provided, AIMS will comme up with its own default
#     ranges (typically 1/30th of the total range for the relevant parameter,
#     centred around the value from the best grid model)
#   - if ranges are provided, they will be re-centred around the parameters
#     of the best model in the grid
#   - the ranges on parameters related to surface amplitudes will be reset by AIMS
#   - exact names should be used as keys, since AIMS accesses
#     the relevant distributions by using the name as a key.
#   - it doesn't matter if there are supplementary parameters
#     which don't intervene. AIMS will simply ignore them.

tight_ball_range = {}  # do NOT erase this line
#tight_ball_range["Mass"]     = ("Gaussian", [0.0, 0.10])
#########################   Radial orders   ################################
use_n       = True  # use radial orders when comparing observations with models?
read_n      = True  # read radial orders from input file?
assign_n    = False # use best model to reassign the radial order?
                    # NOTE: this supersedes "read_n"
#########################   Constraints   ##################################
# Determines the type of surface correction to include.  Options include:
#   - None: don't use any surface corrections
#   - "Kjeldsen2008": use surface corrections based on Kjeldsen et al. (2008)
#   - "Kjeldsen2008_scaling": use surface corrections based on Kjeldsen et al. (2008);
#                     The b exponent is determined based on the scaling relation
#                     from Sonoi et al. (2015)
#   - "Kjeldsen2008_2": use surface corrections based on Kjeldsen et al. (2008);
#                     The b exponent is a free parameter
#   - "Ball2014": use one-term surface corrections based on Ball & Gizon (2014)
#   - "Ball2014_2": use two-term surface corrections based on Ball & Gizon (2014)
#   - "Sonoi2015": use the surface corrections based on Sonoi et al. (2015)
#   - "Sonoi2015_scaling": use the surface corrections based on Sonoi et al. (2015);
#                     The beta exponent is determined based on the scaling relation
#                     from Sonoi et al. (2015)
#   - "Sonoi2015_2": use the surface corrections based on Sonoi et al. (2015);
#                     The beta exponent is a free parameter
surface_option = "Kjeldsen2008"
b_Kjeldsen2008 = 4.9  # exponent used in the Kjeldsen et al. surface corrections
beta_Sonoi2015 = 4.0  # exponent used in the Sonoi et al. surface corrections

# Set of seismic constraints to be used. Options include:
#   - "nu": individual frequencies
#   - "r02", "r01", "r10": various frequency ratios
#   - "dnu0": individual large frequency separation using l=0
#   - "avg_dnu": average large frequency separation using all l
#   - "avg_dnu0": average large frequency separation using l=0
# NOTE: combining "nu" with the other constraints leads to a (nearly)
#       singular covariance matrix and is not expected to give good
#       results.
#seismic_constraints = ["r02","r01","r10","avg_dnu0","nu_min0"]
seismic_constraints = ["nu"] 

#########################   Weighting   ########################################
# Determines what type of weighting to apply to seismic and classic contraints.
# Options include:
#    - None: no weighting
#    - "Absolute": absolute weight applied to constraints
#    - "Relative": weights applied after normalising the classic and seismic
#                  constraints to have the same weight.
# NOTE: even with the relative weighting, classic_weight is kept as absolute.
weight_option = "Relative"
seismic_weight = 1.0
classic_weight = 1.0

#########################   Input   ########################################
write_data    = False            # set this to True if you want to write a
                                 # binary grid file
mode_format   = "simple"         # specifies the format of the files with
                                 # the mode frequencies.  Options include:
                                 #   - "simple": the original AIMS format
                                 #   - "CLES": the CLES format (almost the same as "simple")
                                 #   - "MESA": the GYRE format
                                 #   - "agsm": the agsm format from ADIPLS
npositive     = True             # if True, only save modes with n >= 0 in
                                 # binary grid file
cutoff        = 5.0              # remove frequencies above this value times
                                 # the acoustic cutoff-frequency
agsm_cutoff   = False            # if True, only keep frequencies with icase=10010
                                 # (i.e. below the cutoff frequency as determined
                                 # by ADIPLS) in agsm files.  This test is in
                                 # addition to the above user-defined cutoff.
list_grid      = "list_MESA_ms"  # file with list of models and characteristics.
                                 # only used when constructing binary file with
                                 # the model grid (i.e. write_data == True)
grid_params = ('Mass', 'log_Z')  # primary grid parameters (excluding age)
                                 # only used when constructing binary file with
                                 # the model grid (i.e. write_data == True)
                                 # These parameters are used to distinguish
                                 # evolutionary tracks
binary_grid = "data_MESA_ms_log" # binary file with model grid
                                 # this file is written to if write_data == True
                                 # this file is read from if write_data = False
track_threshold = 10             # minimal number of models for a stellar evolutionary
                                 # track.  Tracks with fewer models are removed
#########################   User-defined parameters   ######################
# This variable allows the user to introduce supplementary parameters in
# addition to the parameters hard-coded in to AIMS.  These parameters
# can then be used as output parameters (see output_params) and/or even
# as grid parameters used to define evolutionary tracks (see grid_params).
#
# This variable must be a list (or tuple) of pairs of strings.  The first
# string corresponds to the name of the variable which should be used, for
# instance, in the grid_params and output_params variables.  The second
# string is the fancy latex name for this variable.  Allowance needs to
# be made for a prefix and and postfix (hence the two "%s").  These will
# be replaced by appropriate strings if, for instance, one asks for the
# log of this parameter.

user_params = (('alpha_MLT', 'Mixing length parameter, $%s\\alpha_{\\mathrm{MLT}}%s$'), \
               ('Zs', 'Surface metallicity, $%sZ_s%s$'), \
               ('Xs', 'Surface hydrogen, $%sX_s%s$'),    \
               ('Zc', 'Central metallicity, $%sZ_c%s$'), \
               ('Xc', 'Central hydrogen, $%sX_c%s$'))
#########################   Priors    ######################################
# The priors are given in a similar format as the tight-ball ranges above.
# An important difference is that the relevant probability distributions
# will not be recentred or renormalised (in the case of surface term
# amplitudes).
#
# NOTES:
#   - if a prior on a given parameter is not provided, AIMS will comme
#     up with its own default prior, namely an uniformative prior if
#     the tight_ball option is set to True, and a uniform prior on a
#     relevant range in the opposite case (except for parameters 
#     related to surface corrections, in which case AIMS explicitely asks
#     the user to provide a prior).
#   - exact names should be used as keys, since AIMS accesses
#     the relevant distributions by using the name as a key
#   - it doesn't matter if there are supplementary parameters
#     which don't intervene. AIMS will simply ignore them.

priors = {}                      # The priors will be defined thanks to this 
#priors["Mass"]     = ("Uniform", [0.8, 1.5])
#########################   Interpolation    ###############################
scale_age = True                 # use a scaled age when interpolating
#########################   Interpolation tests    #########################
test_interpolation = False       # decide whether to test the interpolation.
                                 # If True, interpolation tests are carried
                                 # out for the above binary grid, and written
                                 # in binary format to a file which can
                                 # subsequently be analysed using plot_test.py.
interpolation_file = "interpolation_test"  # Name of the file to which to
                                 # write the results from the interpolation
                                 # tests.  This file can be analysed using
                                 # plot_test.py.
#########################   Output   #######################################
# choice of parameters: "Mass", "Radius", "Luminosity", "X", "Y", "Z", "zsx_0",
#                       "Fe_H", "M_H", "Age", "Teff", "Dnu", "numax", "Rho", "g",
#                       "b_Kjeldsen2008", "beta_Sonoi2015"
# Also all quantities provided in the user_params variable are available.
# If the quantities "Zs" and "Xs" are defined via user_params, then "Ys" and
#                       "zsx_s" are also available.
# possible prefixes: "log_", "ln_", "exp_"
# example: "log_g" corresponds to log_{10}(g), where $g$ is the surface gravity
# NOTE: the quantities "b_Kjeldsen2008" and "beta_Sonoi2015" will be the ones
#       based on the scaling relation in Sonoi et al. (2015).  This will differ
#       from the values obtained when using the options surface_option="Kjeldsen2008_2"
#       or "Sonoi2015_2"
output_params = ("Radius","Mass","log_g","Rho","Age","Teff","Xc","Luminosity")
output_dir    = "results"      # name of the root folder with the results
output_osm    = "osm"          # name of the root folder with the OSM files
extended_model  = False        # if True, print all models frequencies
                               # (not just those with observed counterparts)
with_osm        = False        # decide whether to write output files for OSM 
                               # (=Optimal Stellar Model by R. Samadi)
with_combinations = True       # decide whether to write file with model combinations
with_walkers    = True         # decide whether to plot walkers
with_distrib_iter=True         # decide whether to plot distrib_iter
with_echelle    = True         # decide whether to plot echelle diagrams
with_histograms = True         # decide whether to plot histograms
with_triangles  = True         # decide whether to make triangle plots
with_rejected   = True         # decide whether to make triangle plots with accepted/rejected models
plot_extensions = ['png']      # extensions (and formats) for all simple plots
tri_extensions  = ['png']      # extensions (and formats) for triangle plots
                               # supported formats (may depend on backend): eps, jpeg,
                               #   jpg, pdf, png, ps, raw, rgba, svg, svgz, tif, tiff
backend         = 'agg'        # matplotlib backend with which to produce plots.
                               # Options (may differ according to installation):
                               #   'pdf', 'pgf', 'Qt4Agg', 'GTK', 'GTKAgg',
                               #   'ps', 'agg', 'cairo', 'MacOSX', 'GTKCairo',
                               #   'WXAgg', 'template', 'TkAgg', 'GTK3Cairo',
                               #   'GTK3Agg', 'svg', 'WebAgg', 'CocoaAgg',
                               #   'emf', 'gdk', 'WX'
                               # None = use default backend
                               # NOTE: some backends (such as 'TkAgg') do not
                               #       work in batch mode.
