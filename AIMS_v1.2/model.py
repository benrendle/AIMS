#!/usr/bin/env python
# coding: utf-8
# $Id: model.py
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
A module which contains various classes relevant to the grid of models:

- :py:class:`Model`: a model
- :py:class:`Track`: an evolutionary track
- :py:class:`Model_grid`: a grid of models

These different classes allow the program to store a grid of models and perform
a number of operations, such as:

- retrieving model properties
- interpolate within the grid models
- sort the models within a given evolutionary track
- ...

"""

__docformat__ = 'restructuredtext'

# AIMS configuration:
import AIMS_configure as config

# various packages needed for AIMS
import math
import numpy as np
import random
import sys
from operator import methodcaller
import matplotlib
if (config.backend is not None): matplotlib.use(config.backend)
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

# packages from within AIMS:
import constants
import utilities
import aims_fortran

# Global parameters:
tol     = 1e-10
""" tolerance level for slightly negative interpolation coefficients """

eps     = 1e-6
""" relative tolerance on parameters used for setting up evolutionary tracks """

ntype = np.int16
""" type used for the n values """

ltype = np.int8   # int8 ranges from -128 to 127 
""" type used for the l values """

ftype = np.float64   
"""  type used for the frequencies """

gtype = np.float64
""" type used for grid data """

modetype = [('n',ntype),('l',ltype),('freq',ftype),('inertia',ftype)] 
""" structure for modes """

# quantities related to user-defined parameters:
user_params_index = {}
""" dictionary which will supply the appropriate index for the user-defined parameters"""

user_params_latex = {}
""" dictionary which will supply the appropriate latex name for the user-defined parameters"""

# integer indices for various global quantities which will be stored in a np array:
nglb         = 8 + len(config.user_params)
""" total number of global quantities in a model (see :py:data:`Model.glb`)."""

nlin         = 5 + len(config.user_params)
"""
total number of global quantities which are interpolated in a linear way (see
:py:func:`combine_models`).  These quantities are numbered 0:nlin-1
"""

iage         = 0  
""" index of the parameter corresponding to age in the :py:data:`Model.glb` array """

imass        = 1
""" index of the parameter corresponding to mass in the :py:data:`Model.glb` array """

itemperature = 2
""" index of the parameter corresponding to temperature in the :py:data:`Model.glb` array """

iz0          = 3
"""
index of the parameter corresponding to the initial metallicity
the :py:data:`Model.glb` array
"""

ix0          = 4
"""
index of the parameter corresponding to the initial hydrogen content
in the :py:data:`Model.glb` array
"""

ifreq_ref    = 5 + len(config.user_params)
"""
index of the parameter corresponding to the reference frequency
(used to non-dimensionalise the pulsation frequencies of the model)
in the :py:data:`Model.glb` array
"""

iradius      = 6 + len(config.user_params)
""" index of the parameter corresponding to radius in the :py:data:`Model.glb` array """

iluminosity  = 7 + len(config.user_params)
""" index of the parameter corresponding to luminosity in the :py:data:`Model.glb` array """

def string_to_latex(string,prefix="",postfix=""):
    """
    Return a fancy latex name for an input string.
    
    :param string: string that indicates for which parameter we're seeking a latex name
    :param prefix: optional prefix to add to the string
    :param postfix: optional postfix to add to the string

    :type string: string
    :type prefix: string
    :type postfix: string

    :return: a fancy latex name
    :rtype: string

    .. note::
      This also works for the names of the amplitude parameters for surface corrections.
    """

    if (string.startswith("log_")):
        return string_to_latex(string[4:],prefix+r"\log_{10}\left(",r"\right)"+postfix)
    if (string.startswith("ln_")):
        return string_to_latex(string[3:],prefix+r"\ln\left(",r"\right)"+postfix)
    if (string.startswith("exp_")):
        return string_to_latex(string[4:],prefix+r"\exp\left(",r"\right)"+postfix)
    if (string == "Mass"):       return r'Mass, $%sM/M_{\odot}%s$'%(prefix,postfix)
    if (string == "Radius"):     return r'Radius, $%sR/R_{\odot}%s$'%(prefix,postfix)
    if (string == "Luminosity"): return r'Luminosity, $%sL/L_{\odot}%s$'%(prefix,postfix)
    if (string == "Z"):          return r'Metallicity, $%sZ_0%s$'%(prefix,postfix)
    if (string == "Y"):          return r'Helium content, $%sY_0%s$'%(prefix,postfix)
    if (string == "X"):          return r'Hydrogen content, $%sX_0%s$'%(prefix,postfix)
    if (string == "Ys"):         return r'Helium content, $%sY_s%s$'%(prefix,postfix)
    if (string == "zsx_s"):      return r'Ratio, $%s(Z/X)_s%s$'%(prefix,postfix)
    if (string == "zsx_0"):      return r'Ratio, $%s(Z/X)_0%s$'%(prefix,postfix)
    if (string == "Fe_H"):       return r'Iron content, $%s\mathrm{[Fe/H]}%s$'%(prefix,postfix)
    if (string == "M_H"):        return r'Metal content, $%s\mathrm{[M/H]}%s$'%(prefix,postfix)
    if (string == "Age"):        return r'Age (in Myrs), $%st%s$'%(prefix,postfix)
    if (string == "Teff"):       return r'$%sT_{\mathrm{eff}}%s$ (in K)'%(prefix,postfix)
    if (string == "Dnu"):        return r'Large separation (in $\mu$Hz), $%s\Delta \nu%s$'%(prefix,postfix)
    if (string == "numax"):      return r'$%s\nu_{\mathrm{max}}%s$ (in $\mu$Hz)'%(prefix,postfix)
    if (string == "Rho"):        return r'Density (in g.cm$^{-3}$), $%s\rho%s$'%(prefix,postfix)
    if (string == "g"):          return r'Surface gravity (in cm.s$^{-2}$), $%sg%s$'%(prefix,postfix)
    if (string == "A_surf"):     return r'$%sA_{\mathrm{surf}}%s$'%(prefix,postfix)
    if (string == "A3_surf"):    return r'$%sA_3^{\mathrm{surf}}%s$'%(prefix,postfix)
    if (string == "Am1_surf"):   return r'$%sA_{-1}^{\mathrm{surf}}%s$'%(prefix,postfix)
    if (string == "alpha_surf"): return r'$%s\alpha_{\mathrm{surf}}%s$'%(prefix,postfix)
    if (string == "b_Kjeldsen2008"): return r'$%sb_{\mathrm{Kjeldsen\,et\,al.\,(2008)}}%s$'%(prefix,postfix)
    if (string == "beta_Sonoi2015"): return r'$%s\beta_{\mathrm{Sonoi\,et\,al.\,(2015)}}%s$'%(prefix,postfix)

    try:
        return user_params_latex[string]%(prefix,postfix)
    except KeyError:
        sys.exit("ERROR: unrecognised model quantity: "+string)

def get_surface_parameter_names(surface_option):
    """
    Return the relevant parameter names for a given surface correction option.

    :param surface_option: specifies the type of surface correction.
    :type surface_option: string

    :return: names for the surface parameters
    :rtype: tuple of strings
    """
        
    if (surface_option is None):                   return ()
    if (surface_option == "Kjeldsen2008"):         return ("A_surf",)
    if (surface_option == "Kjeldsen2008_scaling"): return ("A_surf",)
    if (surface_option == "Kjeldsen2008_2"):       return ("A_surf","b_Kjeldsen2008")
    if (surface_option == "Ball2014"):             return ("A3_surf",)
    if (surface_option == "Ball2014_2"):           return ("A3_surf","Am1_surf")
    if (surface_option == "Sonoi2015"):            return ("alpha_surf",)
    if (surface_option == "Sonoi2015_scaling"):    return ("alpha_surf",)
    if (surface_option == "Sonoi2015_2"):          return ("alpha_surf","beta_Sonoi2015")
    sys.exit("ERROR: Unknown surface correction: "+surface_option)

class Model:

    """A class which contains a stellar model, including classical and seismic information."""

    def string_to_param(self,string):
        """
        Return a parameter for an input string.
        
        :param string: string that indicates which parameter we're seeking
        :type string: string

        :return: the value of the parameter
        :rtype: float
        """

        if (string.startswith("log_")): return math.log10(self.string_to_param(string[4:]))
        if (string.startswith("ln_")):  return math.log(self.string_to_param(string[3:]))
        if (string.startswith("exp_")): return math.exp(self.string_to_param(string[4:]))
        if (string == "Mass"):       return self.glb[imass]/constants.solar_mass
        if (string == "Radius"):     return self.glb[iradius]/constants.solar_radius
        if (string == "Luminosity"): return self.glb[iluminosity]/constants.solar_luminosity
        if (string == "Z"):          return self.glb[iz0]
        if (string == "Y"):          return 1.0-self.glb[iz0]-self.glb[ix0]
        if (string == "X"):          return self.glb[ix0]
        if (string == "Ys"):         return 1.0-self.glb[user_params_index["Zs"]]-self.glb[user_params_index["Xs"]]
        if (string == "zsx_s"):      return self.zsx_s
        if (string == "zsx_0"):      return self.zsx_0
        if (string == "Fe_H"):       return self.FeH
        if (string == "M_H"):        return self.MH
        if (string == "Age"):        return self.glb[iage]
        if (string == "Teff"):       return self.glb[itemperature]
        if (string == "Dnu"):        return self.find_large_separation()*self.glb[ifreq_ref]
        if (string == "numax"):      return self.numax
        if (string == "Rho"):        return 3.0*self.glb[imass]/(4.0*math.pi*self.glb[iradius]**3)
        if (string == "g"):          return constants.G*self.glb[imass]/self.glb[iradius]**2
        if (string == "beta_Sonoi2015"): return self.beta_Sonoi2015
        if (string == "b_Kjeldsen2008"): return self.b_Kjeldsen2008

        try:
            return self.glb[user_params_index[string]]
        except KeyError:
            sys.exit("ERROR: unrecognised model quantity: "+string)

    def __init__(self, _glb, _name=None, _modes=None):
        """
        :param _glb: 1D array of global parameters for this model.  Its dimension should
          be greater or equal to :py:data:`nglb`
        :param _name: name of the model (typically the second part of its path)
        :param _modes: list of modes in the form of tuples 
          (n,l,freq,inertia) which will be appended to the set of modes in the model.

        :type _glb: np.array
        :type _name: string
        :type _modes: list of (int, int, float, float)
        """

        # check inputs
        assert (_glb[imass] >= 0.0),        "A star cannot have a negative mass!"
        assert (_glb[iradius] >= 0.0),      "A star cannot have a negative radius!"
        assert (_glb[iluminosity] >= 0.0),  "A star cannot have a negative luminosity!"
        assert (_glb[iz0] >= 0.0),          "A star cannot have a negative metallicity!"
        assert (_glb[ix0] >= 0.0),          "A star cannot have a negative hydrogen abundance!"
        assert (_glb[iage] >= 0.0),         "A star cannot have a negative age!"
        assert (_glb[itemperature] >= 0.0), "A star cannot have a negative temperature!"

        self.name = _name
        """Name of the model, typically the second part of its path"""

        self.glb = _glb
        """Array which will contain various global quantities"""

        self.glb[ifreq_ref] = 5e5*math.sqrt(constants.G*self.glb[imass]/self.glb[iradius]**3)/math.pi
        """Characteristic frequency of the model in :math:`\\mathrm{cyclic \\, \\mu Hz}`"""

        self.modes = np.empty([0],dtype=modetype)
        """array containing the modes (n, l, freq, inertia)"""

        if (_modes is not None): self.append_modes(_modes)

    def read_file(self,filename):
        """
        Read in a set of modes from a file.  This method will either
        call :py:meth:`read_file_simple` or :py:meth:`read_file_agsm`
        according to the value of the ``mode_format`` variable in
        ``AIMS_configure.py``.

        :param filename: name of the file with the modes. The format of this file
                         is decided by the ``mode_format`` variable in
                         ``AIMS_configure.py``.

        :type filename: string

        :return: ``True`` if at least one frequency has been discarded (see note below).
        :rtype: boolean

        .. note::
          At this stage the frequencies should be expressed in :math:`\mu\mathrm{Hz}`.
          They will be non-dimensionalised in :py:func:`read_model_list`.
        """

        if (config.mode_format == "simple"): return self.read_file_simple(filename)
        if (config.mode_format == "agsm"):   return self.read_file_agsm(filename)
        sys.exit("ERROR: unrecognised format \""+config.mode_format+"\".\n" \
                +"       Please choose another value for mode_format in AIMS_configure.py")

    def read_file_simple(self,filename):
        """
        Read in a set of modes from a file.  This uses the "simple" format as
        specified in the ``mode_format`` variable in ``AIMS_configure.py``.

        :param filename: name of the file with the modes.
          The file should contain a one-line header followed by five
          columns which correspond to l, n, frequency, dfreq_var, inertia.

        :type filename: string

        :return: ``True`` if at least one frequency has been discarded (see note below).
        :rtype: boolean

        .. note::
          - The dfreq_var column is discarded.
          - Frequencies above :math:`1.1\\nu_{\\mathrm{cut-off}}` are discarded.
        """

        freqlim = config.cutoff*self.cutoff
        exceed_freqlim = False
        freqfile = open(filename)
        freqfile.readline() # skip head
        mode_temp = [] 
        for line in freqfile:
            line = line.strip()
            columns = line.split()
            n = int(columns[1])
            freq = utilities.to_float(columns[2])
            # remove frequencies above AIMS_configure.cutoff*nu_{cut-off}
            if (freq > freqlim): 
                exceed_freqlim = True
                continue
            if (config.npositive and (n < 0)): continue  # remove g-modes if need be
            mode_temp.append((n,int(columns[0]),freq,utilities.to_float(columns[4])))
        freqfile.close()
        self.modes = np.array(mode_temp,dtype=modetype)

        return exceed_freqlim

    def read_file_agsm(self,filename):
        """
        Read in a set of modes from a file.  This uses the "agsm" format as
        specified in the ``mode_format`` variable in ``AIMS_configure.py``.

        :param filename: name of the file with the modes. This file is a
            binary fortran "agsm" file produced by the ADIPLS code.  See
            instructions to the ADIPLS code for a description of this
            format.

        :type filename: string

        :return: ``True`` if at least one frequency has been discarded (see note below).
        :rtype: boolean
        """

        narr,larr,farr,iarr,nn,exceed_freqlim =  \
            aims_fortran.read_file_agsm(filename,config.npositive,config.agsm_cutoff, \
            config.cutoff*self.cutoff)
        self.modes = np.array(zip(narr[0:nn],larr[0:nn],farr[0:nn],iarr[0:nn]),dtype=modetype)

        return exceed_freqlim

    def write_file_simple(self,filename):
        """
        Write a set of modes into a file using the "simple" format as
        described in :py:meth:`read_file_simple`.

        :param filename: name of the file where the modes should
                         be written.

        :type filename: string

        .. note::
          - Frequencies are non-dimensional and should expressed in muHz
        """

        output = open(filename,"w")
        # write header
        output.write("# %1s %3s %22s %6s %22s\n"%("l","n","nu_theo (muHz)","unused","Inertia"))
        for i in xrange(self.modes.shape[0]):
            output.write("  %1d %3d %22.15e    0.0 %22.15e\n"%( \
                self.modes["l"][i],                        \
                self.modes["n"][i],                        \
                self.modes["freq"][i]*self.glb[ifreq_ref], \
                self.modes["inertia"][i]))
        output.close()

    def append_modes(self,modes):
        """
        Append a list of modes to the model.

        :param modes: list of modes which are in the form of tuples: (n,l,freq,inertia).
        :type modes: list of (int, int, float, float)
        """
        self.modes = np.append(self.modes,np.array(modes,dtype=modetype))
        
    def sort_modes(self):
        """
        Sort the modes by l, then n, then freq.
        """
        # sorts by l, then n, then freq
        ind = np.lexsort((self.modes['freq'], self.modes['n'],self.modes['l']))
        self.modes = np.array([self.modes[i] for i in ind],dtype=modetype)

    def remove_duplicate_modes(self):
        """
        Remove duplicate modes.

        Modes are considered to be duplicate if they have the same l
        and n values (regardless of frequency).

        :return: ``True`` if at least one mode has been removed.
        :rtype: boolean

        .. warning::
          This method assumes the modes are sorted.
        """

        ind = [] 
        for i in xrange(len(self.modes)-1,1,-1):
            if ((self.modes['n'][i] == self.modes['n'][i-1]) and \
                (self.modes['l'][i] == self.modes['l'][i-1])):
                ind.append(i)
        self.modes = np.delete(self.modes,ind)
        return (len(ind) == 0)

    def get_age(self):
        """
        Return age of stellar model.

        This is useful for sorting purposes.

        :return: the age of the model
        :rtype: float
        """
        return self.glb[iage]

    def get_freq(self, surface_option=None, a=[]):
        """
        Obtain model frequencies, with optional frequency corrections.
        
        :param surface_option: specifies the type of surface correction.  Options include:
          
          - ``None``: no corrections are applied
          - ``"Kjeldsen2008"``: apply a correction based on Kjeldsen et al. (2008)
          - ``"Kjeldsen2008_scaling"``: apply a correction based on Kjeldsen et al. (2008).
                                        The exponent is based on a scaling relation from
                                        Sonoi et al. (2015).
          - ``"Kjeldsen2008_2"``: apply a correction based on Kjeldsen et al. (2008).
                                        The exponent is a free parameter.
          - ``"Ball2014"``:     apply a one-term correction based on Ball and Gizon (2014)
          - ``"Ball2014_2"``:   apply a two-term correction based on Ball and Gizon (2014)
          - ``"Sonoi2015"``: apply a correction based on Sonoi et al. (2015)
          - ``"Sonoi2015_scaling"``: apply a correction based on Sonoi et al. (2015)
                                        The exponent is based on a scaling relation from
                                        Sonoi et al. (2015).
          - ``"Sonoi2015_2"``: apply a correction based on Sonoi et al. (2015)
                                        The exponent is a free parameter.

        :param a: amplitude parameters which intervene in the surface correction

        :type surface_option: string
        :type a: array-like

        :return: models frequencies (including surface corrections)
        :rtype: np.array

        .. note::        
          If surface_option==None or a==[], the original frequencies are
          returned (hence modifying them modifies the :py:class:`Model` object).
        """
        
        if (surface_option is None) or (len(a) == 0): return self.modes['freq']
        return self.modes['freq'] + self.get_surface_correction(surface_option, a)

    def get_surface_correction(self, surface_option, a):
        """
        Obtain corrections on model frequencies (these corrections should be
        *added* to the *theorectical* frequencies).
        
        :param surface_option: specifies the type of surface correction.
          Options include:
          
          - ``None``: no corrections are applied
          - ``"Kjeldsen2008"``: apply a correction based on Kjeldsen et al. (2008)
          - ``"Kjeldsen2008_scaling"``: apply a correction based on Kjeldsen et al. (2008).
                                        The exponent is based on a scaling relation from
                                        Sonoi et al. (2015).
          - ``"Kjeldsen2008_2"``: apply a correction based on Kjeldsen et al. (2008).
                                        The exponent is a free parameter.
          - ``"Ball2014"``:     apply a one-term correction based on Ball and Gizon (2014)
          - ``"Ball2014_2"``:   apply a two-term correction based on Ball and Gizon (2014)
          - ``"Sonoi2015"``: apply a correction based on Sonoi et al. (2015)
          - ``"Sonoi2015_scaling"``: apply a correction based on Sonoi et al. (2015)
                                        The exponent is based on a scaling relation from
                                        Sonoi et al. (2015).
          - ``"Sonoi2015_2"``: apply a correction based on Sonoi et al. (2015)
                                        The exponent is a free parameter.

        :param a: parameters which intervene in the surface correction.  According to
          the correction they take on the following meanings:

          - ``"Kjeldsen2008"``:         a[0]*freq**b_Kjeldsen2008
          - ``"Kjeldsen2008_scaling"``: a[0]*freq**b_scaling
          - ``"Kjeldsen2008_2"``:       a[0]*freq**a[1]
          - ``"Ball2014"``:             a[0]*freq**3/I 
          - ``"Ball2014_2"``:           a[0]*freq**3/I + a[1]/(freq*I)
          - ``"Sonoi2015"``:            a[0]*[1 - 1/(1 + (nu/numax)**beta_Sonoi2015)]
          - ``"Sonoi2015_scaling"``:    a[0]*[1 - 1/(1 + (nu/numax)**beta_scaling)]
          - ``"Sonoi2015_2"``:          a[0]*[1 - 1/(1 + (nu/numax)**a[1])]

        :type surface_option: string
        :type a: array-like

        :return: surface corrections on the model frequencies
        :rtype: np.array

        .. note::
          The array operations lead to the creation of a new array with the
          result, which avoids modifications of the original frequencies and inertias.
        """
        
        if (surface_option is None):                   return np.zeros(len(self.modes), dtype=ftype)
        if (surface_option == "Kjeldsen2008"):         return a[0]*self.modes['freq']**config.b_Kjeldsen2008
        if (surface_option == "Kjeldsen2008_scaling"): return a[0]*self.modes['freq']**self.b_Kjeldsen2008
        if (surface_option == "Kjeldsen2008_2"):       return a[0]*self.modes['freq']**a[1]
        if (surface_option == "Ball2014"):             return a[0]*self.modes['freq']**3/self.modes['inertia']
        if (surface_option == "Ball2014_2"):           return a[0]*self.modes['freq']**3/self.modes['inertia'] \
                                                            + a[1]/(self.modes['freq']*self.modes['inertia'])
        if (surface_option == "Sonoi2015"):            return a[0]*(1.0-1.0/(1.0+(self.glb[ifreq_ref]*self.modes['freq'] \
                                                            / self.numax)**config.beta_Sonoi2015))
        if (surface_option == "Sonoi2015_scaling"):    return a[0]*(1.0-1.0/(1.0+(self.glb[ifreq_ref]*self.modes['freq'] \
                                                            / self.numax)**self.beta_Sonoi2015))
        if (surface_option == "Sonoi2015_2"):          return a[0]*(1.0-1.0/(1.0+(self.glb[ifreq_ref]*self.modes['freq'] \
                                                            / self.numax)**a[1]))
        sys.exit("ERROR: Unknown surface correction: "+surface_option)

    @property
    def b_Kjeldsen2008(self):
        """
        Return the exponent for the Kjeldsen et al. (2008) surface correction
        recipe, as calculated based on the Sonoi et al. (2015) scaling relation.

        :return: the Kjeldsen et al. exponent
        :rtype: float
        """
        return 10.0**(-3.16*self.string_to_param("log_Teff") + 0.184*self.string_to_param("log_g")+11.7)
 
    @property
    def beta_Sonoi2015(self):
        """
        Return the exponent for the Sonoi et al. (2015) surface correction
        recipe, as calculated based on the Sonoi et al. (2015) scaling relation.

        :return: the Kjeldsen et al. exponent
        :rtype: float
        """
        return 10.0**(-3.86*self.string_to_param("log_Teff") + 0.235*self.string_to_param("log_g")+14.2)

    def multiply_modes(self,constant):
        """
        Multiply the frequencies by constant.

        :param constant: constant by which the mode frequencies are multiplied
        :type constant: float
        """

        # NOTE: inertias are non-dimensionless, so they don't change.
        # (this has been tested numerically)
        self.modes['freq']*=constant

    def find_mode(self, ntarget, ltarget):
        """
        Find a mode with specific n and l values.

        :param ntarget: target n value
        :param ltarget: target l value
        :type ntarget: int
        :type ltarget: int

        :return: the frequency of the mode
        :rtype: float
        """

        size = len(self.modes)
        for i in xrange(size):
            if ((self.modes['n'][i] == ntarget) and (self.modes['l'][i] == ltarget)):
                return self.modes['freq'][i]
        return np.nan

    def find_mode_range(self):
        """ 
        Find n and l ranges of the modes in the model.

        :return: the n and l ranges of the modes
        :rtype: int, int, int, int
        """

        if (len(self.n) < 1): return -1,-1,-1,-1
        nmin = np.nanmin(self.modes['n'])
        nmax = np.nanmax(self.modes['n'])
        lmin = np.nanmin(self.modes['l'])
        lmax = np.nanmax(self.modes['l'])
        return nmin, nmax, lmin, lmax

    def find_large_separation(self):
        """
        Find large frequency separation using only radial modes.

        :return: the large frequency separation
        :rtype: float
        """

        ind = (self.modes['l'] == 0).flat
        n = self.modes['n'].compress(ind)
        if (len(n) < 2): return np.nan
        freq = self.modes['freq'].compress(ind)
        coeff = np.polyfit(n,freq,1)
        return coeff[0]

    def find_epsilon(self, ltarget):
        """
        Find epsilon, the constant offset in a simplified version of Tassoul's
        asymptotic formula:
        
        :math:`\\nu_n = \\Delta \\nu (n + \\varepsilon)`
        
        :param ltarget: target l value.  Only modes with this l value will be
          used in obtaining epsilon.
        :type ltarget: int

        :return: the constant offset
        :rtype: float
        """

        dnu = self.find_large_separation()
        one = n = nu = 0.0
        for i in xrange(len(self.modes)):
            if (self.modes['l'][i] != ltarget): continue
            one   += 1.0
            n     += self.modes['n'][i]
            nu    += self.modes['freq'][i]
        if (one == 0.0):
            return 0.0
        else:
            return (nu/dnu-n)/one

    @property
    def FeH(self):
        """
        Find [Fe/H] value for model.

        The conversion from (Xs,Zs) to [Fe/H] is performed using the
        following formula:

            :math:`\\mathrm{[Fe/H] = \\frac{[M/H]}{A_{FeH}}  \
                                   = \\frac{1}{A_{FeH}} \\log_{10} \
                                     \\left(\\frac{z/x}{z_{\\odot}/x_{\\odot}} \\right)}`

        :return: the :math:`\\mathrm{[Fe/H]}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """

        return math.log10(self.glb[user_params_index["Zs"]]*constants.solar_x/(self.glb[user_params_index["Xs"]]*constants.solar_z))/constants.A_FeH

    @property
    def MH(self):
        """
        Find [M/H] value for model.

        The conversion from (Xs,Zs) to [M/H] is performed using the
        following formula:

            :math:`\\mathrm{[M/H] = \\log_{10} \\left(\\frac{z/x}{z_{\\odot}/x_{\\odot}} \\right)}`

        :return: the :math:`\\mathrm{[M/H]}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """

        return math.log10(self.glb[user_params_index["Zs"]]*constants.solar_x/(self.glb[user_params_index["Xs"]]*constants.solar_z))

    @property
    def zsx_s(self):
        """
        Find the Zs/Xs value

        :return: the Zs/Xs value
        :rtype: float
        """

        return self.glb[user_params_index["Zs"]]/self.glb[user_params_index["Xs"]]

    @property
    def zsx_0(self):
        """
        Find the Z0/X0 value

        :return: the Z0/X0 value
        :rtype: float
        """
        return self.glb[iz0]/self.glb[ix0]

    @property
    def numax(self):
        """
        Find :math:`\\nu_{\\mathrm{max}}` for model.

        The :math:`\\nu_{\\mathrm{max}}` value is obtained from the
        following scaling relation:

            :math:`\\frac{\\nu_{\\mathrm{max}}}{\\nu_{\\mathrm{max},\odot}} \
                      = \\left(\\frac{M}{M_{\\odot}}\\right)                \
                        \\left(\\frac{R}{R_{\\odot}}\\right)^2              \
                        \\left(\\frac{T_{\\mathrm{eff}}}{T_{\\mathrm{eff},\\odot}}\\right)^{-1/2}`

        :return: the :math:`\\nu_{\\mathrm{max}}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """

        return constants.solar_numax*(self.glb[imass]/constants.solar_mass) \
                                    /((self.glb[iradius]/constants.solar_radius)**2 \
                                    *math.sqrt(self.glb[itemperature]/constants.solar_temperature))

    @property
    def cutoff(self):
        """
        Find :math:`\\nu_{\\mathrm{cut-off}}` for model.

        The :math:`\\nu_{\\mathrm{cut-off}}` value is obtained from the
        following scaling relation:

            :math:`\\frac{\\nu_{\\mathrm{cut-off}}}{\\nu_{\\mathrm{cut-off},\odot}} \
                      = \\left(\\frac{M}{M_{\\odot}}\\right)                        \
                        \\left(\\frac{R}{R_{\\odot}}\\right)^2                      \
                        \\left(\\frac{T_{\\mathrm{eff}}}{T_{\\mathrm{eff},\\odot}}\\right)^{-1/2}`

        :return: the :math:`\\nu_{\\mathrm{cut-off}}` value
        :rtype: float

        .. note::
          The relevant values are given in :py:mod:`constants`
        """

        return constants.solar_cutoff*(self.glb[imass]/constants.solar_mass) \
                                     /((self.glb[iradius]/constants.solar_radius)**2 \
                                     *math.sqrt(self.glb[itemperature]/constants.solar_temperature))

    def freq_sorted(self):
        """
        Check to see if the frequencies are in ascending order for each l value.

        :return: ``True`` if the frequencies are in ascending order.
        :rtype: boolean
        """

        for i in xrange(len(self.modes)-1):
            if (self.modes['l'][i] > 0): continue
            if (self.modes['l'][i] != self.modes['l'][i+1]): continue
            if (self.modes['freq'][i] > self.modes['freq'][i+1]):
                return False
        return True

    def print_me(self):
        """ Print classical and seismic characteristics of the model to standard output."""

        print "----- Model:",self.name," -----"
        print "Mass (in M_sun):              %.5f" % (self.glb[imass]/constants.solar_mass)
        print "Radius (in R_sun):            %.5f" % (self.glb[iradius]/constants.solar_radius)
        print "Reference frequency (in uHz): %.3f" % self.glb[ifreq_ref]
        print "Temperature (in K):           %.1f" % self.glb[itemperature]
        print "Luminosity (in L_sun):        %.3g" % (self.glb[iluminosity]/constants.solar_luminosity)
        print "Age (in Myrs):                %.2f" % self.glb[iage]
        print "Z:                            %.4f" % self.glb[iz0]
        print "X:                            %.4f" % self.glb[ix0]
        for (name, latex_name) in config.user_params:
            print "{0:29} {1:.5e}".format(name,self.glb[user_params_index[name]])
        print "Modes (in muHz):"
        size = self.modes.shape[0]
        for i in xrange(size):
            print "  (n,l,freq,IK) = (%d, %d, %.15f, %.5e)" % \
                   (self.modes['n'][i], self.modes['l'][i],   \
                    self.modes['freq'][i]*self.glb[ifreq_ref],\
                    self.modes['inertia'][i])

class Track:
    """
    An evolutionary track.
    """

    def __init__(self, aModel, grid_params):
        """
        :param aModel: first model to be added to evolutionary track (it does
          not need to be the youngest model in an evolutionary sequence).  This
          Model is used to obtain the relevant parameters for the evolutionary
          track (as given by the :py:data:`grid_params` variable).

        :param grid_params: list of strings which are the names of the
          parameters which describe the evolutionary track.

        :type aModel: :py:class:`Model`
        :type grid_params: list of strings
        """

        self.grid_params = grid_params
        """Names of the parameters used to construct the grid"""

        self.params = map(aModel.string_to_param,self.grid_params)
        """Parameters which characterise this evolutionary track"""

        self.nmodes = len(aModel.modes)
        """Total number pulsation modes from all of the models in this evolutionary track"""

        self.models = [aModel]
        """List of models in this evolutionary track"""

    def append(self,aModel):
        """
        Append a model to the evolutionary track.

        
        :param aModel: model which is being appended to the track
        :type aModel: :py:class:`Model`
        """

        self.models.append(aModel)
        self.nmodes += len(aModel.modes)

    def matches(self, aModel):
        """
        Check to see if a model matches the evolutionary track and can therefore
        be included in the track.

        :param aModel: input model being tested
        :type aModel: :py:class:`Model`

        :return: ``True`` only if the model given as an argument has parameters which
          match those of the evolutionary track.
        :rtype: boolean
        """

        params_bis = map(aModel.string_to_param,self.grid_params)
        for param1, param2 in zip(self.params, params_bis):
            if (abs(param1/param2 - 1.0) > eps): return False
        return True

    def sort(self):
        """Sort models within evolutionary track according to age."""

        self.models.sort(key=methodcaller('get_age'))

    def is_sorted(self):
        """
        Check to see of models are in ascending order according to age.

        :return: ``True`` if the models ar in order of increasing age
        :rtype: boolean
        """

        return all(self.models[i].glb[iage] <= self.models[i+1].glb[iage] for i in xrange(len(self.models)-1))

    def duplicate_ages(self):
        """
        Check to see if you track contains models with duplicate ages.

        :return: ``True`` if there are duplicate age(s)
        :rtype: boolean

        .. warning::
            This method should only be applied after the track has been
            sorted.
        """

        return any(self.models[i].glb[iage] == self.models[i+1].glb[iage] for i in xrange(len(self.models)-1))

    def interpolate_model(self, age):
        """
        Return a model at a given age which is obtained using linear interpolation.

        :param age: age of desired model in :math:`\\mathrm{Myrs}`
        :type age: float

        :return: the interpolated model
        :rtype: :py:class:`Model`

        .. warning::
          This method assumes the track is sorted, since it applies
          a binary search algorithm for increased efficiency.
        """
        
        # easy exit:
        if (age < self.models[0].glb[iage]): return None
        if (age > self.models[-1].glb[iage]): return None
        
        istart = 0
        istop  = len(self.models)-1
        while (istop > istart+1):
            itemp = (istop+istart)/2
            if (age < self.models[itemp].glb[iage]):
                istop = itemp
            else:
                istart = itemp
        mu = (age - self.models[istart].glb[iage]) \
           / (self.models[istop].glb[iage] - self.models[istart].glb[iage])

        return combine_models(self.models[istart],1.0-mu,self.models[istop],mu)

    def find_combination(self, age, coef):
        """
        Return a model combination at a given age which is obtained using linear interpolation.

        :param age: age of desired model in :math:`\\mathrm{Myrs}`
        :param coef: coefficient which multiplies this combination

        :type age: float
        :type coef: float

        :return: pairs composed of an interpolation coefficient and a model name
        :rtype: tuple of (float, string)

        .. warning::
          This method assumes the track is sorted, since it applies
          a binary search algorithm for increased efficiency.
        """
        
        # easy exit:
        if (age < self.models[0].glb[iage]): return None
        if (age > self.models[-1].glb[iage]): return None
        
        istart = 0
        istop  = len(self.models)-1
        while (istop > istart+1):
            itemp = (istop+istart)/2
            if (age < self.models[itemp].glb[iage]):
                istop = itemp
            else:
                istart = itemp
        mu = (age - self.models[istart].glb[iage]) \
           / (self.models[istop].glb[iage] - self.models[istart].glb[iage])
        return ((coef*(1.0-mu), self.models[istart].name), (coef*mu, self.models[istop].name))

    def find_modes(self, ntarget, ltarget):
        """
        Return two lists, one with the ages of the models and the other
        with the mode frequencies corresponding to target n and l values.

        This function is useful for seeing how the frequency of a particular
        mode changes with stellar age.

        :param ntarget: target n value
        :param ltarget: target l value

        :type ntarget: int
        :type ltarget: int

        :return: lists of ages and frequencies
        :rtype: list, list
        """

        ages  = []
        freqs = []
        for model in self.models:
           freq = model.find_mode(ntarget, ltarget)
           if (math.isnan(freq)): continue
           freqs.append(freq)
           ages.append(model.glb[iage])
        return ages, freqs 

    def find_mode_range(self):
        """
        Find n and l ranges of modes in models

        :return: the n and l ranges
        :rtype: int, int, int, int
        """

        if (len(self.models) < 1): return -1,-1,-1,-1
        nmin_total,nmax_total,lmin_total,lmax_total = self.models[0].find_mode_range()
        for model in self.models:
            nmin, nmax, lmin, lmax = model.find_mode_range()
            if (nmin < nmin_total): nmin_total = nmin
            if (nmax > nmax_total): nmax_total = nmax
            if (lmin < lmin_total): lmin_total = lmin
            if (lmax > lmax_total): lmax_total = lmax
        return nmin_total, nmax_total, lmin_total, lmax_total

    def test_interpolation(self, nincr):
        """
        Test accuracy of interpolation along evolutionary track.

        This method removes every other model and retrieves its frequencies
        by interpolation from neighbouring models.  The accuracy of the
        interpolated frequencies and global parameters are tested by carrying
        out comparisons with the original models.

        :param nincr: increment with which to carry out the interpolation.
          By comparing results for different values of ``nincr``, one can
          evaluate how the interpolation error depends on the size of the
          interval over which the interpolation is carried out.
        :type nincr: int

        :return: the interpolation errors
        :rtype: np.array 
        """

        # initialisation
        nmodels = len(self.models)
        ndim = len(self.params)+1
        result = np.zeros((nmodels-2*nincr,ndim+nglb+6),dtype=gtype)

        # loop through all models:
        for i in xrange(nincr,nmodels-nincr):
            # carry out interpolation
            mu = (self.models[i].glb[iage]   - self.models[i-nincr].glb[iage]) \
               / (self.models[i+nincr].glb[iage] - self.models[i-nincr].glb[iage])
            aModel = combine_models(self.models[i-nincr],1.0-mu,self.models[i+nincr],mu)

            result[i-nincr,0:ndim-1] = self.params
            result[i-nincr,ndim-1] = self.models[i].glb[iage]
            result[i-nincr,ndim:ndim+nglb+6] = compare_models(aModel,self.models[i])

        return result

class Model_grid:
    """
    A grid of models.
    """

    def __init__(self):
        self.ndim = 0
        """
        Number of dimensions for the grid (excluding age), as based on the
        :py:data:`Model_grid.grid_params` variable
        """

        self.tracks = []
        """List of evolutionary tracks contained in the grid."""

        self.ndx = []
        """List containing track indices"""

        self.grid = None
        """Array containing the grid parameters for each evolutionary track (excluding age)."""

        self.tessellation = None
        """Object containing the tessellation of the grid used for interpolation."""

        self.grid_params = None
        """
        Set of parameters (excluding age) used to construct the grid and do interpolations.

        .. note::        
          For best interpolation results, these parameters should be comparable.
        """

        self.prefix = None
        """Root folder with grid of models (including final slash)."""

        self.postfix = ".freq" 
        """Last part of the filenames which contain the model frequencies (default = ".freq")."""

        self.user_params = config.user_params
        """
        The set of user parameters involved in the grid.  This is to avoid having a different
        set of user parameters in `AIMS_configure.py`
        """

    def read_model_list(self,filename):
        """
        Read list of models from a file and construct a grid.

        :param filename: name of the file with the list.  The first line
          of this file should contain a prefix which is typically the root
          folder of the grid of models.  This followed by a file with multiple
          columns.  The first 8 contain the following information for each model:

          1. the second part of the path.  When concatenated with the prefix
             on the first line, this should give the full path to the model.
          2. The stellar mass in :math:`\mathrm{g}`
          3. The stellar radius in :math:`\mathrm{cm}`
          4. The stellar luminosity in :math:`\mathrm{g.cm^2.s^{-3}}`
          5. The metallicity
          6. The hydrogen content
          7. The stellar age in :math:`\mathrm{Myrs}`
          8. The effective temperature in :math:`\mathrm{K}`

          The following columns contain the parameters specified in the
          :py:data:`AIMS_configure.user_params` variable.

        :type filename: string
        """

        self.grid_params = config.grid_params

        # set the correct dimension:
        self.ndim = len(self.grid_params)

        # set prefix and postfix:
        listfile = open(filename,"r")
        line = listfile.readline().strip();
        columns = line.split()
        if (len(columns) < 1): sys.exit("Erroneous first line in %s."%(filename))
        self.prefix = columns[0]
        if (len(columns) > 1): self.postfix = columns[1]

        # read models and put them into evolutionary tracks:
        nmodels = 0
        nmodes  = 0
        models_small_spectra = []
        for line in listfile:
            line = line.strip()
            columns = line.split()
            glb = np.empty((nglb,),dtype = gtype)
            glb[imass]        = utilities.to_float(columns[1])
            glb[iradius]      = utilities.to_float(columns[2])
            glb[iluminosity]  = utilities.to_float(columns[3])
            glb[iz0]          = utilities.to_float(columns[4])
            glb[ix0]          = utilities.to_float(columns[5])
            glb[iage]         = utilities.to_float(columns[6])
            glb[itemperature] = utilities.to_float(columns[7])

            i = 8
            for (name, name_latex) in config.user_params:
                glb[user_params_index[name]] = utilities.to_float(columns[i])
                i += 1

            aModel = Model(glb, _name = columns[0])
            exceed_freqlim = aModel.read_file(self.prefix + columns[0] + self.postfix)
            aModel.multiply_modes(1.0/aModel.glb[ifreq_ref])  # make frequencies non-dimensional
            aModel.sort_modes()
            aModel.remove_duplicate_modes()
            for track in self.tracks:
                if (track.matches(aModel)):
                    track.append(aModel)
                    break
            else:
                aTrack = Track(aModel,self.grid_params)
                self.tracks.append(aTrack)
            nmodels += 1
            nmodes  += len(aModel.modes)
            if (not exceed_freqlim):
                models_small_spectra.append(aModel.name)
            print nmodels, nmodes
        listfile.close()

        # right list of models with spectra which are too small in a file:
        output = open("models_small_spectra","w")
        for name in models_small_spectra: output.write(name+"\n")
        output.close()

        # sort tracks:
        for track in self.tracks: track.sort()

        # sanity check:
        for track in self.tracks:
            if track.duplicate_ages():
                print "ERROR: the track ",track.grid_params," = ",track.params
                print "       has models with the same age.  Please remove"
                print "       duplicate models."
                sys.exit(1)

        # update list of indices:
        self.ndx = range(len(self.tracks))
        
        # need to create grid from scratch since tracks have been sorted.
        self.grid = np.asarray([track.params for track in self.tracks])

    def tessellate(self):
        """Apply Delauny triangulation to obtain the grid tessellation."""

        self.tessellation = Delaunay(self.grid)
        
    def plot_tessellation(self):
        """
        Plot the grid tessellation.

        .. warning::
          This only works for two-dimensional tessellations.
        """

        if (self.ndim != 2):
            print "Only able to plot the tessellation in two dimensions."
            return

        # find bounds:
        xmin = np.nanmin(self.grid[:,0])
        xmax = np.nanmax(self.grid[:,0])
        ymin = np.nanmin(self.grid[:,1])
        ymax = np.nanmax(self.grid[:,1])
        dx = xmax-xmin
        dy = ymax-ymin
        xmin -= dx*0.03
        xmax += dx*0.03
        ymin -= dy*0.05
        ymax += dy*0.05

        #plt.semilogy(self.grid[:,0],self.grid[:,1],'o')
        plt.plot(self.grid[:,0],self.grid[:,1],'o')
        plt.triplot(self.grid[:,0],self.grid[:,1],self.tessellation.simplices.copy())
        plt.xlim((xmin,xmax))
        plt.ylim((ymin,ymax))
        plt.xlabel(string_to_latex(self.grid_params[0]))
        plt.ylabel(string_to_latex(self.grid_params[1]))
        plt.savefig("tessellation.eps")

    def test_interpolation(self):
        """
        Test interpolation between different evolutionary tracks in a given grid.

        :return: The following four items are returned:

          - the interpolation errors
          - the first half of the partition (where the interpolation is tested)
          - the second half of the partition (used to carry out the interpolation)
          - the tessellation associated with the second half of the partition

        :rtype: np.array, list, list, tessellation object
        """

        ndx1, ndx2 = self.find_partition()
        tessellation = Delaunay(self.grid[ndx2,:])

        # initialisation
        results = []
        ndim = self.ndim+1

        for j in ndx1:
            nmodels = len(self.tracks[j].models)
            aResult = np.empty((nmodels,ndim+nglb+6),dtype=gtype)
            pt = self.tracks[j].params + [0.0,]

            for i in xrange(nmodels):
                aModel1 = self.tracks[j].models[i]
                pt[-1] = aModel1.glb[iage]
                aModel2 = interpolate_model(self,pt,tessellation,ndx2)
                aResult[i,0:ndim] = pt
                if (aModel2 is None):
                    aResult[i,ndim:ndim+nglb+6] = np.nan
                else:
                    aResult[i,ndim:ndim+nglb+6] = compare_models(aModel1,aModel2)

            results.append(aResult)

        return results, ndx1, ndx2, tessellation

    def find_partition(self):
        """
        Find a partition of the grid for use with :py:meth:`Model_grid.test_interpolation`

        :return: a random partition of [0 ... n-1] into two equal halves, where n is
                 the number of tracks in the grid
        :rtype: two lists of int

        """

        ndx = range(len(self.tracks))
        random.shuffle(ndx)
        nn = len(self.tracks)/2
        return ndx[:nn],ndx[nn:]

    def test_freq(self):
        """
        Test to see if frequencies in all of the models of the grid
        are in ascending order for each l value.

        :return: The following items are returned

          - the effective temperatures of the models with frequencies out of order
          - the luminosities of the models with frequencies out of order
          - the effective temperatures of the models with sorted frequencies
          - the luminosities of the models with sorted frequencies

        :rtype: four lists of floats
        """
        
        Teffs, Lums, Teffs_out, Lums_out = [], [], [], []
        for track in self.tracks:
            for model in track.models:
                if (not model.freq_sorted()):
                    print model.name
                    Teffs_out.append(model.string_to_param("Teff"))
                    Lums_out.append(model.string_to_param("log_Luminosity"))
                else:
                    Teffs.append(model.string_to_param("Teff"))
                    Lums.append(model.string_to_param("log_Luminosity"))
        return Teffs_out, Lums_out, Teffs, Lums

    def find_epsilons(self, ltarget):
        """
        Find epsilon values in models from the grid
        
        :param ltarget: target l value for which epsilons are being obtained
        :type ltarget: int

        :return: the epsilon values
        :rtype: list of floats
        """
        
        epsilons = []
        for track in self.tracks:
            for model in track.models:
                epsilon = model.find_epsilon(ltarget)
                if (epsilon != 0.0): epsilons.append(epsilon)
        return epsilons

def init_user_param_dict():
    """
    Initialise the dictionaries which are related to user-defined parameters.  For
    a given parameter, these dictionaries provide the appropriate index for for
    the :py:data:`Model.glb` array as well as the appropriate latex name.
    """

    i = 5
    for (name, latex_name) in config.user_params:
        user_params_index[name] = i
        user_params_latex[name] = latex_name
        i += 1

def combine_models(model1,coef1,model2,coef2):
    """
    Do linear combination of this model with another.

    This method returns a new model which is the weighted sum
    of two models for the purposes of model interpolation.
    The classical parameters are combined in a self-consistent
    way as are the frequencies.

    :param model1: first model
    :param coef1:  weighting coefficient applied to first model
    :param model2: second model
    :param coef2:  weighting coefficient applied to second model

    :type model1: :py:class:`Model`
    :type coef1:  float
    :type model2: :py:class:`Model`
    :type coef2:  float

    :return: the combined model
    :rtype: :py:class:`Model`

    .. warning::
      One should avoid negative or zero coefficients as
      these could lead to undefined results.
    """
        
    # find global parameters (try to be self-consistent):

    # this first part is simply a linear combination:
    glb = np.empty((nglb,),dtype=gtype)
    glb[0:nlin] = coef1*model1.glb[0:nlin] + coef2*model2.glb[0:nlin]

    # this next part depends on previous results:
    glb[iradius] = (glb[imass]/(coef1*model1.glb[imass]/model1.glb[iradius]**3
                 + coef2*model2.glb[imass]/model2.glb[iradius]**3))**(1.0/3.0)
    cnst1 = model1.glb[iluminosity]/(model1.glb[iradius]**2*model1.glb[itemperature]**4)
    cnst2 = model2.glb[iluminosity]/(model2.glb[iradius]**2*model2.glb[itemperature]**4)
    glb[iluminosity] = (coef1*cnst1 + coef2*cnst2)*glb[iradius]**2*glb[itemperature]**4
    # glb[ifreq_ref] will be correctly defined when the Model() constructor is invoked

    # interpolate spectra:    
    size3 = max(model1.modes.shape[0],model2.modes.shape[0])
    nvalues = np.empty((size3,),dtype=ntype)
    lvalues = np.empty((size3,),dtype=ltype)
    fvalues = np.empty((size3,),dtype=ftype)
    ivalues = np.empty((size3,),dtype=ftype)

    nvalues,lvalues,fvalues,ivalues,n3 = aims_fortran.combine_modes( \
            coef1,model1.modes['n'],model1.modes['l'],model1.modes['freq'],model1.modes['inertia'], \
            coef2,model2.modes['n'],model2.modes['l'],model2.modes['freq'],model2.modes['inertia'], \
            nvalues,lvalues,fvalues,ivalues)

    return Model(glb, _modes=zip(nvalues[0:n3],lvalues[0:n3],fvalues[0:n3],ivalues[0:n3]))

def compare_models(model1,model2):
    """
    Compare two models and find the largest frequency different for
    radial and non-radial modes.

    :param model1: first model
    :param model2: second model

    :type model1: :py:class:`Model`
    :type model2: :py:class:`Model`

    :return: a 1D array to be used in ``plot_test_interpolation.py``
      with the following measurements of the differences between the two models:

      - ``result[0]`` = maximum error on the radial modes
      - ``result[1]`` = RMS error on the radial modes
      - ``result[2]`` = RMS error on the radial modes near
        :math:`\\nu_{\mathrm{max}}`
      - ``result[3]`` = maximum error on the non radial modes
      - ``result[4]`` = RMS error on the non radial modes
      - ``result[5]`` = RMS error on the non radial modes near
        :math:`\\nu_{\mathrm{max}}`
      - ``result[6+[0:nglb]]`` = errors on the global parameters

    :rtype: np.array
    """

    # initialisation:
    n_radial = 0
    n_radial_numax = 0
    n_non_radial = 0
    n_non_radial_numax = 0
    result = np.zeros((6+nglb,),dtype=gtype)
    # define frequency interval around numax:
    numax = 0.5*(model1.numax/model1.glb[ifreq_ref] \
          +      model2.numax/model2.glb[ifreq_ref])
    a = 0.8*numax
    b = 1.2*numax

    # compare frequency spectra:
    size1 = len(model1.modes)
    size2 = len(model2.modes)
    i1 = i2 = 0
    while((i1 < size1) and (i2 < size2)):
        if (model1.modes['l'][i1] < model2.modes['l'][i2]): i1+=1; continue
        if (model1.modes['l'][i1] > model2.modes['l'][i2]): i2+=1; continue
        if (model1.modes['n'][i1] < model2.modes['n'][i2]): i1+=1; continue
        if (model1.modes['n'][i1] > model2.modes['n'][i2]): i2+=1; continue

        # now the two modes have the same n and l values:
        diff = abs(model1.modes['freq'][i1] - model2.modes['freq'][i2])
        avg_freq =(model1.modes['freq'][i1] + model2.modes['freq'][i2])/2.0
        if (model1.modes['l'][i1] == 0):
            if (result[0] < diff): result[0] = diff 
            diff *= diff  # square diff
            result[1] += diff
            n_radial += 1
            # in python, this is called an interval comparison:
            if (a <= avg_freq <= b):
                result[2] += diff 
                n_radial_numax += 1
        else:
            if (result[3] < diff): result[3] = diff 
            diff *= diff  # square diff
            result[4] += diff
            n_non_radial += 1
            if (a <= avg_freq <= b):
                result[5] += diff 
                n_non_radial_numax += 1
        i1+=1
        i2+=1

    # avoid divisions by zero:
    if (n_radial > 0):
        result[1] = math.sqrt(result[1]/float(n_radial))
    else:
        result[1] = np.nan

    if (n_radial_numax > 0):
        result[2] = math.sqrt(result[2]/float(n_radial_numax))
    else:
        result[2] = np.nan

    if (n_non_radial > 0):
        result[4] = math.sqrt(result[4]/float(n_non_radial))
    else:
        result[4] = np.nan

    if (n_non_radial_numax > 0):
        result[5] = math.sqrt(result[5]/float(n_non_radial_numax))
    else:
        result[5] = np.nan

    # absolute differences on global parameters:
    result[6:6+nglb] = np.absolute(model1.glb - model2.glb)

    return result

def find_interpolation_coefficients(grid,pt,tessellation,ndx):
    """
    Find interpolation weights from the corresponding simplex.

    Linear interpolation weights are obtained with the simplex
    by finding the barycentric coordinates of the point given
    by ``pt``.

    :param grid: grid of models in which we're carrying out the
      interpolation
    :param pt: set of parameters used for finding the
      interpolation weights.  The first part contains the grid
      parameters (relevant to this interpolation), whereas
      the last element is the age (not used here).  If the
      provided set of parameters lies outside the grid, then
      ``None`` is returned instead of an interpolated model.
    :param tessellation: tessellation with which to carry out the
      interpolation.
    :param ndx: indices of the grid points associated with the
      tessellation

    :type grid: :py:class:`Model_grid`
    :type pt: array-like
    :type ndx: list of int

    :return: lists of interpolation coefficients and tracks
    :rtype: list of floats, list of :py:class:`Track`
    """
    
    if (pt is None): return None, None
    if (tessellation is None): return None, None
    pt1 = np.asarray(pt[0:-1],dtype=gtype)
    val = tessellation.find_simplex(pt1.reshape((1,grid.ndim)))[0]

    # see if point is outside tessellation
    if (val == -1): return None, None
    mat = tessellation.transform[val]

    # make sure the transformation matrix is defined:
    if (math.isnan(np.sum(mat))): return None, None

    b = mat[:grid.ndim].dot(pt1-mat[grid.ndim])
    coefs = np.r_[b, 1.0-b.sum()]
    ind = tessellation.simplices[val]

    # check to make sure you're not outside the grid:
    for coef in coefs:
        if (coef < -tol): return None, None

    # produce results, filtering out zero elements:
    coefs_out = []
    tracks = []
    for coef,i in zip(coefs,ind):
        # remove negative coefficients to avoid problems.
        if (coef > 0.0):
            coefs_out.append(coef)
            tracks.append(grid.tracks[ndx[i]])

    return coefs_out, tracks

def find_ages(coefs, tracks, age):
    """
    Find ages to which each track needs to be interpolated for a specified
    age.  The global variable ``scale_age`` decides between the following
    two options:
    
    1. ``scale_age`` = ``False``: each track is simply interpolated to ``age``.
    2. ``scale_age`` = ``True``: the age of each model along each evolutionary
       track, including the interpolated track, is linearly mapped onto the
       interval [0,1].  A dimensionless parameter ``eta`` is obtained by
       interpolating ``age`` onto the interval [0,1], using the linear
       transformation associated with the interpolated track.  Using the
       parameter eta, a corresponding age is obtained along each track.

    .. figure:: ./figures/age_interpolation.*
      :figclass: align-center
 
      This diagram illustrates both types of age interpolation and shows
      the advantages of selecting ``scale_age`` = ``True``.


    :param coefs: interpolation coefficients used to weight each track.
    :param tracks:  evolutionary tracks involved in the interpolation.
    :param age: target age for the output interpolated model.

    :type coefs: list of floats
    :type tracks:  list of :py:class:`Track`
    :type age: float

    :return: the relevant age for each track
    :rtype: list of floats

    .. note::
      - the interpolation coefficients should add up to 1.0
      - there should be as many tracks as interpolation coefficients.
    """

    assert (len(coefs) == len(tracks)), "Mismatch between len(coefs) and len(tracks)"
    
    if (not config.scale_age): return [age]*len(coefs)
    
    age_s = 0.0
    age_f = 0.0
    for coef,track in zip(coefs,tracks):
        age_s += coef*track.models[0].glb[iage]
        age_f += coef*track.models[-1].glb[iage]
        
    eta = (age-age_s)/(age_f-age_s)
    
    # check to see if the age lies within the interpolated track:
    if (eta < 0.0): return None
    if (eta > 1.0): return None
    
    ages = []
    for coef,track in zip(coefs,tracks):
        ages.append((1.0-eta)*track.models[0].glb[iage] + eta*track.models[-1].glb[iage])

    return ages

def interpolate_model(grid,pt,tessellation,ndx):
    """
    Interpolate model in grid using provided parameters.

    The interpolation is carried out in two steps.  First, linear
    interpolation according to age is carried out on each node of
    the simplex containing the set of parameters.  This interpolation
    is done using the :py:class:`Track.interpolate_model` method.
    Then, linear interpolation is carried out within the simplex.
    This achieved by finding the barycentric coordinates of the
    model (i.e. the weights), before combining the age-interpolated
    models form the nodes using the :py:class:`combine_models` method.
    In this manner, the weights are only calculated once, thereby
    increasing computational efficiency.

    :param grid: grid of models in which we're carrying out the
      interpolation
    :param pt: set of parameters used for the interpolation.
      The first part contains the grid parameters, whereas
      the last element is the age.  If the provided set
      of parameters lies outside the grid, then ``None``
      is returned instead of an interpolated model.
    :param tessellation: tessellation with which to carry out the
      interpolation.
    :param ndx: indices of the grid points associated with the
      tessellation

    :type grid: :py:class:`Model_grid`
    :type pt: array-like
    :type ndx: list of int

    :return: the interpolated model
    :rtype: :py:class:`Model`
    """

    # find simplex interpolation coefficients
    coefs,tracks = find_interpolation_coefficients(grid,pt,tessellation,ndx)
    if (coefs is None): return None

    # find ages:
    ages = find_ages(coefs,tracks,pt[-1])
    if (ages is None): return None

    n = len(tracks)

    # treat the case where there is only 1 model:
    if (n == 1):
        if (abs(coefs[0]-1.0) > eps):
            print "WARNING: erroneous interpolation coefficient: ",coefs[0]
        return tracks[0].interpolate_model(ages[0])

    # treat the case where there are at least 2 models:
    aModel1 = tracks[0].interpolate_model(ages[0])
    if (aModel1 is None): return None
    aModel2 = tracks[1].interpolate_model(ages[1])
    if (aModel2 is None): return None
    aModel1 = combine_models(aModel1,coefs[0],aModel2,coefs[1])
    for i in xrange(2,n):
        aModel2 = tracks[i].interpolate_model(ages[i])
        if (aModel2 is None): return None
        aModel1 = combine_models(aModel1,1.0,aModel2,coefs[i]);
    return aModel1

def find_combination(grid,pt):
    """
    Find linear combination of models which corresponds to interpolating
    the model based on the provided parameters.

    The interpolation is carried out using the same procedure as in
    :py:func:`interpolate_model`.

    :param grid: grid of models in which we're carrying out the
      interpolation
    :param pt: set of parameters used for the interpolation.
      The first part contains the grid parameters, whereas
      the last element is the age.  If the provided set
      of parameters lies outside the grid, then ``None``
      is returned instead of an interpolated model.

    :type grid: :py:class:`Model_grid`
    :type pt: array-like

    :return: pairs of coefficients and model names
    :rtype: tuple of (float,string)
    """

    # find simplex interpolation coefficients
    coefs,tracks = find_interpolation_coefficients(grid,pt,grid.tessellation,grid.ndx)
    if (coefs is None): return None

    # find ages:
    ages = find_ages(coefs,tracks,pt[-1])
    if (ages is None): return None

    n = len(tracks)

    # combine multiple models:
    results = ()
    for coef,track,age in zip(coefs,tracks,ages):
        if (coef < 0.0): return None # make sure we're not outside the grid
        result = track.find_combination(age,coef)
        if (result is None): return None
        results += result
    return results

