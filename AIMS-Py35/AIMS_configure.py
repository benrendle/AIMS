#!/usr/bin/env python
# $Id: AIMS_configure.py
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

import math

#########################   Parallelisation   ##############################
# NOTE: this is currently implemented with multiprocessing, which duplicates
#       the memory in each process.  To be more memory efficient, turn off
#       parallelisation using the "parallel" parameter.
nprocesses  = 2      # number of processes (if running in parallel)
parallel    = False #$$$True   # specifies whether to run in parallel

#########################   EMCEE control parameters   #####################
ntemps      = 3 #$$$5 # number of temperatures
nwalkers    = 60 #$$$400     # number of walkers (this number should be even)
nsteps0     = 900 #$$$200     # number of burn-in steps
nsteps      = 1000 #$$$4000     # number of steps
thin        = 10     # thinning parameter (1 out of thin steps will be kept ...)
thin_comb   = 100    # thinning parameter for output linear combinations of models
PT          = True   # use parallel tempering?

#########################   Initialisation   ###############################
tight_ball   = True  # initialise with a tight ball around best solution
max_iter     = 10000  # maximum number of iterations to find walker

# Ranges used around tight ball configuration for walkers.
# NOTES:
#   - these ranges will be re-centred around the parameters of the
#     best model in the grid
#   - the ranges on parameters related to surface amplitudes will be reset by AIMS
#   - exact names should be used as keys, since AIMS accesses
#     the relevant distributions by using the name as a key.
#   - it doesn't matter if there are supplementary parameters
#     which don't intervene. AIMS will simply ignore them.

tight_ball_range = {}
tight_ball_range["Mass"]     = ("Gaussian", [0.0, 0.01])	# (29/06/16) edited these to 0.005, 0.002, 0.05, 100, 0.1, 0.1, 0.02, 10
tight_ball_range["Z"]        = ("Gaussian", [0.0, 0.01])	# 0.01 0.002 0.05 100 0.5 0.5 0.02 10
tight_ball_range["log_Z"]    = ("Gaussian", [0.0, 0.05])
tight_ball_range["X"]        = ("Gaussian", [0.0, 0.01])	# 0.01 0.002 0.05 100 0.5 0.5 0.02 10
# tight_ball_range["log_X"]    = ("Gaussian", [0.0, 0.05])
tight_ball_range["Age"]      = ("Gaussian", [0.0, 100.0])
tight_ball_range["numax"]    = ("Gaussian", [0.0, 0.5])
tight_ball_range["Dnu"]      = ("Gaussian", [0.0, 0.5])
tight_ball_range["Radius"]   = ("Gaussian", [0.0, 0.01])
tight_ball_range["Teff"]     = ("Gaussian", [0.0, 10.0])
tight_ball_range["A_surf"]   = ("Gaussian", [0.0, 1.0])  # will be reset by AIMS
tight_ball_range["A3_surf"]  = ("Gaussian", [0.0, 1.0])  # will be reset by AIMS
tight_ball_range["Am1_surf"] = ("Gaussian", [0.0, 1.0])  # will be reset by AIMS
#########################   Radial orders   ################################
use_n       = True  # use radial orders when comparing observations with models?
read_n      = True # read radial orders from input file?
assign_n    = True  # use best model to reassign the radial order?
                    # NOTE: this supersedes "read_n"
#########################   Constraints   ##################################
# Determines the type of surface correction to include.  Options include:
#   - None: don't use any surface corrections
#   - "Kjeldsen2008": use surface corrections based on Kjeldsen et al. (2008)
#   - "Ball2014": use one-term surface corrections based on Ball & Gizon (2014)
#   - "Ball2014_2": use two-term surface corrections based on Ball & Gizon (2014)
surface_option = "Ball2014"
b_Kjeldsen2008 = 4.9  # exponent used in the Kjeldsen et al. surface corrections

# Set of seismic constraints to be used. Options include:
#   - "nu": individual frequencies
#   - "r02", "r01", "r10": various frequency ratios
#   - "dnu0": individual large frequency separation using l=0
#   - "avg_dnu": average large frequency separation using all l
#   - "avg_dnu0": average large frequency separation using l=0
#   - "input_dnu": if using with only dnu and no freqs, use this constraint
# NOTE: combining "nu" with the other constraints leads to a (nearly)
#       singular covariance matrix and is not expected to give good
#       results.
#seismic_constraints = ["r02","r01","r10","avg_dnu0","nu_min0",	"input_dnu","alt_input_dnu"]
seismic_constraints = ['nu']

#########################   Weighting   ########################################
# Determines what type of weighting to apply to seismic and classic contraints.
# Options include:
#    - None: no weighting
#    - "Absolute": absolute weight applied to constraints
#    - "Relative": weights applied after normalising the classic and seismic
#                  constraints to have the same weight.
# NOTE: even with the relative weighting, classic_weight is kept as absolute.
weight_option = "Absolute"
seismic_weight = 1.0
classic_weight = 1.0

#########################   Input   ########################################
write_data    = False            # set this to True if you want to write a
                                 # binary grid file
mode_format   = "simple"         # specifies the format of the files with
                                 # the mode frequencies.  Options include:
                                 #   - "simple": the original AIMS format
                                 #   - "agsm": the agsm format from ADIPLS
npositive     = True             # only save modes with n >= 0 in binary file
cutoff        = 5.0              # remove frequencies above this value times
                                 # the acoustic cutoff-frequency
agsm_cutoff   = False            # if True, only keep frequencies with icase=10010
                                 # (i.e. below the cutoff frequency as determined
                                 # by ADIPLS) in agsm files.  This test is in
                                 # addition to the above user-defined cutoff.
list_grid     = "list_MESA_ms_improved"   # file with list of models and characteristics.
                                 # only used when constructing binary file with
                                 # the model grid (i.e. write_data == True)
grid_params = ("Mass","X")   # primary grid parameters (excluding age)	<--------- Can only be the values used in the file name - the set global parameters of each track.
                                 # only used when constructing binary file with
                                 # the model grid (i.e. write_data == True)
                                 # These parameters are used to distinguish
                                 # evolutionary tracks
binary_grid = "grid_MS_MESA" # binary file with model grid
                                 # this file is written to if write_data == True
                                 # this file is read from if write_data = False
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

#user_params = ()
user_params = (("Xc", r'Central hydrogen, $%sX_c%s$'),)#("DNl1", r'Period Spacing, $%sDNl1%s$'),)#("mHe", r'Helium Mass'),)
#user_params = (("Xc", r'Central hydrogen, $%sX_c%s$'), \
#               ("alpha_MLT", r'Mixing length parameter, $%s\alpha_{\mathrm{MLT}}%s$'), \
#               ("alpha_semi_conv", r'Semiconvection parameter, $%s\alpha_{\mathrm{semi. conv.}}%s$'))
#########################   Priors    ######################################
# The priors are given in a similar format as the tight-ball ranges above.
# An important difference is that the relevant probability distributions
# will not be recentred or renormalised (in the case of surface term
# amplitudes).
#
# NOTES:
#   - exact names should be used as keys, since AIMS accesses
#     the relevant distributions by using the name as a key
#   - it doesn't matter if there are supplementary parameters
#     which don't intervene. AIMS will simply ignore them.

priors = {}                      # The priors will be defined thanks to this
priors["Mass"]     = ("Uniform", [0.9, 1.5])
priors["Z"]        = ("Uniform", [0.0023, 0.0175])
priors["log_Z"]    = ("Uniform", [math.log10(0.0023), math.log10(0.0175)])
priors["X"]        = ("Uniform", [0.68, 0.73])
# priors["log_X"]    = ("Uniform", [math.log10(0.68), math.log10(0.73)])
priors["Age"]      = ("Uniform", [0.0, 2e4])
priors["numax"]    = ("Uniform", [0.0, 5.0e3])
priors["A_surf"]   = ("Uniform", [-1.0, 1.0])  # this is too broad and will be sent by AIMS
priors["A3_surf"]  = ("Uniform", [-1e-12, 1e-12])  # this too broad and should be set experimentally
priors["Am1_surf"] = ("Uniform", [-1e-6, 1e-6])  # this too broad and should be set experimentally
#########################   Interpolation    ###############################
scale_age = True                 # use a scaled age when interpolating
#########################   Interpolation tests    #########################
test_interpolation = False       # decide whether to test the interpolation.
                                 # If True, interpolation tests are carried
                                 # out for the above binary grid, and written
                                 # in binary format to a file which can
                                 # subsequently be analysed using plot_test.py.
interpolation_file = "interp_RGB"  # Name of the file to which to
                                 # write the results from the interpolation
                                 # tests.  This file can be analysed using
                                 # plot_test.py.
#########################   Output   #######################################
# choice of parameters: "Mass", "Radius", "Luminosity", "Z", "X", "Fe_H",
#                       "M_H", "Age", "Teff", "Dnu", "Rho", "g"
# possible prefixes: "log_", "ln_", "exp_"
# example: "log_g" corresponds to log_{10}(g), where $g$ is the surface gravity
output_params = ("Radius","Mass","log_g","Rho","Age","Teff","X","numax","Dnu","Luminosity","Fe_H","M_H")#,"DNl1")
output_dir    = "results"      # name of the root folder with the results
output_osm    = "osm"          # name of the root folder with the OSM files
with_osm       = False         # decide whether to write output files for
                               # OSM (=Optimal Stellar Model by R. Samadi)
with_combinations = True       # decide whether to write file with model combinations
with_walkers    = True         # decide whether to plot walkers
with_echelle    = True         # decide whether to plot echelle diagrams
with_histograms = True         # decide whether to plot histograms
with_triangles  = True         # decide whether to make triangle plots
with_rejected   = True         # decide whether to make triangle plots with accepted/rejected models
plot_extensions = ['png']      # extensions (and formats) for all simple plots
tri_extensions  = ['png']      # extensions (and formats) for triangle plots
# supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff
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
