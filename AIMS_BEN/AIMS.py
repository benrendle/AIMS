#!/usr/bin/env python
# coding: utf-8
# $Id: AIMS.py
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
A module which contains the main program for AIMS as well as various classes
which intervene when calculating the priors and likelihood function:

- :py:class:`Distribution`: a class which represents a probability distribution
- :py:class:`Prior_list`: a class with a list of priors
- :py:class:`Mode`: a class used to represent observed modes
- :py:class:`Combination`: a class used to represent frequency combinations
- :py:class:`Likelihood`: a class used to represent the likelihood function
- :py:class:`Probability`: a class which groups the priors and likelihood function together

This module relies on the `emcee <http://dan.iel.fm/emcee/current/>`_ package to apply
an MCMC algorithm which will return a representative sample of models for a given set
of seismic an classic constraints.

.. warning::
  In various places in this module, for instance in the :py:class:`Prior_list`
  and :py:class:`Likelihood` classes, various methods return what is described
  as a :math:`\chi^2` value.  Technically, these are not :math:`\chi^2` values,
  but rather :math:`-\chi^2/2`, i.e. the argument of the exponential function
  which intervenes in the Gaussian probability distribution.
"""

__docformat__ = 'restructuredtext'
__version__ = "1.3.0"

# AIMS configuration option
import AIMS_configure as config # user-defined configuration parameters

import dill
import os
import sys
import time
import math
import matplotlib
if (config.backend is not None): matplotlib.use(config.backend)
import matplotlib.pyplot  as plt
import matplotlib.patches as mpatches
import matplotlib.lines   as mlines
import numpy as np
from scipy.stats import truncnorm
import emcee
import ptemcee
import corner
from lxml import etree
from operator import attrgetter
from multiprocessing import Pool
from tqdm import tqdm

# import modules from the AIMS package
import model
import constants
import utilities
import aims_fortran
plt.rcParams["font.family"] = "serif"
# parameters associated with the grid
grid             = None
""" grid of models """

grid_params_MCMC = ()
""" parameters used in the MCMC run (excluding surface correction parameters) """

grid_params_MCMC_with_surf = ()
""" parameters used in the MCMC run (including surface correction parameters) """

ndims            = 0
""" number of dimensions for MCMC parameters (includes :py:data:`nsurf`)  """

nsurf            = 0
""" number of surface term parameters """

# parameters associated with probabilities
prob      = None
"""
:py:class:`Probability` type object that represents the probability function
which includes the likelihood and priors
"""

tight_ball_distributions = None
""" :py:class:`Prior_list` type object with the distributions for the initial tight ball """

log0      = -1e300
""" a large negative value used to represent ln(0) """

threshold = -1e295
""" threshold for "accepted" models.  Needs to be greater than :py:const:`log0` """

# variables associated with the best model from a scan of the entire grid:
best_grid_model   = None
""" best model from a scan of the entire grid """

best_grid_params  = None
""" parameters for the model :py:data:`best_grid_model`"""

best_grid_result  = log0
""" ln(probability) result for the model :py:data:`best_grid_model`"""

best_age_range    = 0.0
""" Age or evolutionary parameter range on track with :py:data:`best_model_model`"""

# variables associated with the best model from the MCMC run:
best_MCMC_model   = None
""" best model from the MCMC run """

best_MCMC_params  = None
""" parameters for the model :py:data:`best_MCMC_model`"""

best_MCMC_result  = log0
""" ln(probability) result for the model :py:data:`best_MCMC_model`"""

# variables associated with the statistical parameters
statistical_model = None
""" model corresponding to statistical parameters """

statistical_params= None
""" parameters for the model :py:data:`statistical_model`"""

statistical_result= log0
""" ln(probability) result for the model :py:data:`statistical_model`"""

# other parameters
output_folder = None
""" folder in which to write the results """

pool          = None
""" pool from which to carry out parallel computations """

my_map        = None
""" pointer to the map function (either the parallel or sequential versions) """

accepted_parameters = []
""" list of parameters associated with accepted models """

rejected_parameters = []
""" list of parameters associated with rejected models """

nreject_classic = 0
""" Number of models rejected because of classic constraints """

nreject_seismic = 0
""" Number of models rejected because of seismic constraints """

nreject_prior   = 0
""" Number of models rejected based on priors """

class Distribution:

    """
    A class which represents a probability distribution, and can yield
    its value for a given input parameter, or provide a random realisation.

    .. note::
      Derived from a class originally written by G. Davies.
    """

    def __init__(self, _type, _values):
        """
        :param _type: type of probability function (current options include
                       "Gaussian", "Truncated_gaussian", "Uniform", "IMF1",
                       "IMF2", "Uninformative")
        :param _values: list of parameters relevant to the probability function

        :type _type: string
        :type _values: list of floats
        """

        self.type = _type
        """
        Type of probability function ("Uniform", "Gaussian", "Truncated_gaussian",
        "IMF1", "IMF2", or "Uninformative")
        """

        self.values = _values
        """List of parameters relevant to probability function"""

    def __call__(self, value):
        """
        Return ln of probability distribution function for a particular value.

        :param value: value at which to evaluate ln of probability distribution function.
        :type value: float

        :return: ln(distribution(value))
        :rtype: float

        .. note::
          The results are given to within an additive constant.
        """

        if self.type == "Gaussian":
            return -((value - self.values[0])/self.values[1])**2 / 2.0
        elif self.type == "Truncated_gaussian":
            if (abs(value - self.values[0]) < self.values[1]*self.values[2]):
                return -((value - self.values[0])/self.values[1])**2 / 2.0
            else:
                return log0
        elif self.type == "Uniform":
            if value < self.values[0] or value > self.values[1]:
                return log0
            else:
                return 0.0
        elif self.type == "IMF1":
            # values[0] < values[1] = mass limits
            # values[2] = exponent
            if (self.values[0] <= value <= self.values[1]):
                return -self.values[2]*math.log(value/self.values[0])
            else:
                return log0
        elif self.type == "IMF2":
            # values[0] < values[1] < values[2] = mass limits
            # values[3], values[4] = exponents
            if (self.values[0] <= value <= self.values[1]):
                return -self.values[3]*math.log(value/self.values[0])
            elif (self.values[1] < value <= self.values[2]):
                return -self.values[3]*math.log(self.values[1]/self.values[0]) \
                       -self.values[4]*math.log(value/self.values[1])
            else:
                return log0
        elif self.type == "Uninformative":
            return 0.0
        else:
            sys.exit("Unrecognised distribution: "+self.type)

    def realisation(self, size=None):
        """
        Return random values which statistically follow the probability distribution.

        :param size: shape of random variates
        :type size: int or tuple of ints

        :return: a set of random realisations
        :rtype: float
        """

        if self.type == "Gaussian":
            return np.random.normal(self.values[0],self.values[1],size=size)
        elif self.type == "Truncated gaussian":
            return truncnorm.rvs(-self.values[2],self.values[2],loc=self.values[0], \
                                 scale=self.values[1],size=size)
        elif self.type == "Uniform":
            return np.random.uniform(self.values[0],self.values[1],size=size)
        elif self.type == "IMF1":
            cnst2 = self.values[0]**(1.0-self.values[2])
            cnst1 = self.values[1]**(1.0-self.values[2]) - cnst2
            return (cnst1*np.random.uniform(0.0,1.0,size=size)+cnst2)**(1.0/(1.0-self.values[2]))
        elif self.type == "IMF2":
            rlim  = 1.0/(1.0 + ((1.0-self.values[3])/(1.0-self.values[4]))                       \
                  * self.values[1]**(self.values[4]-self.values[3])                              \
                  * (self.values[2]**(1.0-self.values[4])-self.values[1]**(1.0-self.values[4]))  \
                  / (self.values[1]**(1.0-self.values[3])-self.values[0]**(1.0-self.values[3])))
            cnst2 = self.values[0]**(1.0-self.values[3])
            cnst1 = (self.values[1]**(1.0-self.values[3]) - cnst2)/rlim
            cnst4 = self.values[1]**(1.0-self.values[4])
            cnst3 = (self.values[2]**(1.0-self.values[4]) - cnst4)/(1.0-rlim)
            def Finv(r):
                if (r <= rlim):
                    return (cnst1*r+cnst2)**(1.0/(1.0-self.values[3]))
                else:
                    return (cnst3*(r-rlim)+cnst4)**(1.0/(1.0-self.values[4]))
            vFinv = np.vectorize(Finv)
            return vFinv(np.random.uniform(0.0,1.0,size=size))
        elif self.type == "Uninformative":
            sys.exit("Unable to produce a realisation for an uninformative distribution")
        else:
            sys.exit("Unrecognised distribution: "+self.type)

    def re_centre(self,value):
        """
        Re-centre the probability distribution around the input value.

        :param value: new value around which to centre the distribution
        :type value: float
        """

        if self.type == "Gaussian":
            self.values[0] = value
        elif self.type == "Truncated_gaussian":
            self.values[0] = value
        elif self.type == "Uniform":
            dist = (self.values[1] - self.values[0])/2.0
            self.values[0] = value - dist
            self.values[1] = value + dist
        elif self.type == "IMF1":
            sys.exit("re_centre method not implemented for IMF1")
        elif self.type == "IMF2":
            sys.exit("re_centre method not implemented for IMF2")
        elif self.type == "Uninformative":
            # do nothing
            pass
        else:
            sys.exit("Unrecognised distribution: "+self.type)

    def re_normalise(self,value):
        """
        Re-normalise the probability distribution so that its characteristic
        width corresponds to the input value.

        :param value: new value around for the chacteristic width
        :type value: float
        """

        if self.type == "Gaussian":
            self.values[1] = value
        elif self.type == "Truncated_gaussian":
            self.values[1] = value
        elif self.type == "Uniform":
            centre = (self.values[0] + self.values[1])/2.0
            self.values[0] = centre - value
            self.values[1] = centre + value
        elif self.type == "IMF1":
            sys.exit("re_normalise method not implemented for IMF1")
        elif self.type == "IMF2":
            sys.exit("re_normalise method not implemented for IMF2")
        elif self.type == "Uninformative":
            sys.exit("Unable to renormalise an uninformative distribution")
        else:
            sys.exit("Unrecognised distribution: "+self.type)

    @property
    def mean(self):
        """
        Returns the mean value of the probability distribution.

        :return:  the mean value of the probability distribution
        :rtype: float
        """
        if self.type == "Gaussian":
            return self.values[0]
        elif self.type == "Truncated_gaussian":
            return self.values[0]
        elif self.type == "Uniform":
            return (self.values[0] + self.values[1])/2.0
        elif self.type == "IMF1":
            return ((1.0-self.values[2]) / (2.0-self.values[2]))  \
              * (self.values[1]**(2.0-self.values[2])-self.values[0]**(2.0-self.values[2])) \
              / (self.values[1]**(1.0-self.values[2])-self.values[0]**(1.0-self.values[2]))
        elif self.type == "IMF2":
            k = 1.0/((self.values[1]**(1.0-self.values[3])-self.values[0]**(1.0-self.values[3])) \
              / ((1.0-self.values[3])*self.values[0]**(-self.values[3]))                         \
              + (self.values[1]/self.values[0])**(-self.values[3])                               \
              * (self.values[2]**(1.0-self.values[4])-self.values[1]**(1.0-self.values[4]))      \
              / ((1.0-self.values[4])*self.values[1]**(-self.values[4])))
            return k*((self.values[1]**(2.0-self.values[3])-self.values[0]**(2.0-self.values[3])) \
              / ((2.0-self.values[3])*self.values[0]**(-self.values[3]))                         \
              + (self.values[1]/self.values[0])**(-self.values[3])                               \
              * (self.values[2]**(2.0-self.values[4])-self.values[1]**(2.0-self.values[4]))      \
              / ((2.0-self.values[4])*self.values[1]**(-self.values[4])))
        elif self.type == "Uninformative":
            return np.nan
        else:
            sys.exit("Unrecognised distribution: "+self.type)

    @property
    def error_bar(self):
        """
        Returns an error bar based on the distribution.  This
        does not necessarily correspond to the one-sigma value
        but rather to what is the most convenient value.

        :return: the error bar
        :rtype: float
        """
        if self.type == "Gaussian":
            return self.values[1]
        elif self.type == "Truncated_gaussian":
            return self.values[1]
        elif self.type == "Uniform":
            return (self.values[1] - self.values[0])/2.0
        elif self.type == "IMF1":
            k = (1.0-self.values[2]) \
              / (self.values[1]**(1.0-self.values[2])-self.values[0]**(1.0-self.values[2]))
            mean1 = (self.values[1]**(2.0-self.values[2])-self.values[0]**(2.0-self.values[2])) \
              * k / (2.0-self.values[2])
            mean2 = (self.values[1]**(3.0-self.values[2])-self.values[0]**(3.0-self.values[2])) \
              * k / (3.0-self.values[2])
            return math.sqrt(mean2 - mean1*mean1)
        elif self.type == "IMF2":
            k = 1.0/((self.values[1]**(1.0-self.values[3])-self.values[0]**(1.0-self.values[3])) \
              / ((1.0-self.values[3])*self.values[0]**(-self.values[3]))                         \
              + (self.values[1]/self.values[0])**(-self.values[3])                               \
              * (self.values[2]**(1.0-self.values[4])-self.values[1]**(1.0-self.values[4]))      \
              / ((1.0-self.values[4])*self.values[1]**(-self.values[4])))
            mean1 = k*((self.values[1]**(2.0-self.values[3])-self.values[0]**(2.0-self.values[3])) \
              / ((2.0-self.values[3])*self.values[0]**(-self.values[3]))                         \
              + (self.values[1]/self.values[0])**(-self.values[3])                               \
              * (self.values[2]**(2.0-self.values[4])-self.values[1]**(2.0-self.values[4]))      \
              / ((2.0-self.values[4])*self.values[1]**(-self.values[4])))
            mean2 = k*((self.values[1]**(3.0-self.values[3])-self.values[0]**(3.0-self.values[3])) \
              / ((3.0-self.values[3])*self.values[0]**(-self.values[3]))                         \
              + (self.values[1]/self.values[0])**(-self.values[3])                               \
              * (self.values[2]**(3.0-self.values[4])-self.values[1]**(3.0-self.values[4]))      \
              / ((3.0-self.values[4])*self.values[1]**(-self.values[4])))
            return math.sqrt(mean2 - mean1*mean1)
        elif self.type == "Uninformative":
            return np.nan
        else:
            sys.exit("Unrecognised distribution: "+self.type)

    @property
    def nparams(self):
        """
        Return the number of relevant parameters for a given distribution.

        :return: the number of relevant parameters
        :rtype: int
        """
        if self.type == "Gaussian":
            return 2
        elif self.type == "Truncated_gaussian":
            return 3
        elif self.type == "Uniform":
            return 2
        elif self.type == "IMF1":
            return 3
        elif self.type == "IMF2":
            return 5
        elif self.type == "Uninformative":
            return 0
        else:
            return 0

    def print_me(self):
        """Print type and parameters of probability distribution."""

        print(self.type + " " + str(self.values[0:self.nparams]))

    def to_string(self):
        """
        Produce nice string representation of the distribution.

        :return: nice string representation of the distribution
        :rtype: string
        """
        return self.type + " " + str(self.values[0:self.nparams])

class Prior_list:
    """
    A class which contains a list of priors as well as convenient methods for
    adding priors and for evaluating them.
    """

    def __init__(self):
        self.priors = []
        """A list of probability distributions which correspond to priors."""

    def add_prior(self, aPrior):
        """
        Add a prior to the list.

        :param aPrior: prior which is to be added to the list.
        :type aPrior: :py:class:`Distribution`
        """

        self.priors.append(aPrior)

    def realisation(self,size=None):
        """
        Return an array with realisations for each prior.  The last dimension
        will correspond to the different priors.

        :param size: shape of random variates (for each prior)
        :type size: int or tuple of ints

        :return: a set of realisations
        :rtype: numpy float array
        """

        if (size is None):
            return np.array([prior.realisation(size=size) for prior in self.priors])
        else:
            if (type(size) == int):
                output_size = (size,1)
            else:
                output_size = size + (1,)
            result = [prior.realisation(size=size).reshape(output_size) for prior in self.priors]
            return np.concatenate(result,axis=-1)

    def __call__(self, params):
        """
        Evaluate ln of priors for a list of parameters.

        :param params: list of parameters which intervene in the priors.
        :type params: array-like

        :return: ln(prior probability)
        :rtype: float

        .. warning::
          The list of parameters should contain the same number of elements
          as the number of priors.
        """
        print(params)
        assert (len(params) == len(self.priors)), "Incorrect number of parameters in call to Log_prior"
        lnprior = 0.0
        for aPrior, param in zip(self.priors, params):
            lnprior += aPrior(param)
        return lnprior

class Mode:
    """
    A class which describes an *observed* pulsation mode.
    """

    def __init__(self, _n, _l, _freq, _dfreq):
        """
        :param _n: radial order of observed mode
        :param _l: harmonic degree of observed mode.
        :param _freq: pulsation frequency (in :math:`\mathrm{\mu Hz}`).
        :param _dfreq: error bar on pulsation frequency (in :math:`\mathrm{\mu Hz}`).

        :type _n: int
        :type _l: int
        :type _freq: float
        :type _dfreq: float

        .. warning::
          Negative values are not accepted for ``_l``, ``_freq``, or ``_dfreq``.
        """

        # check inputs
        assert (_l >= 0),               "Modes with l < 0 do not exist"
        assert ((_n >= 0) or (_l > 0)), "Radial g-modes do not exist"
        assert (_freq >=0),             "I don't accept modes with negative frequencies"
        assert (_dfreq >=0),            "I don't accept modes with negative error bars"
        self.n = _n
        """Radial order of observed mode."""

        self.l = _l
        """Harmonic degree of observed mode."""

        self.freq = _freq
        """Pulsation frequency (in :math:`\mathrm{\mu Hz}`)."""

        self.dfreq = _dfreq
        """Error bar on pulsation frequency (in :math:`\mathrm{\mu Hz}`)."""

    def match(self,a_mode):
        """
        Check to see if input mode has the same (n,l) values as the current mode.

        :param a_mode: input mode which is being compared with current mode.
        :type a_mode: :py:class:`Mode`

        :return: ``True`` if the input mode has the same (n,l) values as the
          current mode.
        :rtype: boolean
        """

        if (self.l != a_mode.l): return False
        if (self.n != a_mode.n): return False
        return True

class Combination:
    """
    A class which contains indices and coefficients which intervene in:

    - linear combinations of frequencies
    - frequency ratios
    """

    def __init__(self):
        self.value = 0.0
        """Value of the frequency combination or ratio."""

        self.num = 0.0
        """Value of the frequency combination or numerator in a frequency ratio."""

        self.den = 0.0
        """Value of the denomenator in a frequency ratio."""

        self.num_index = []
        """Indices in a linear combination or numerator of a frequency ratio."""

        self.den_index = []
        """Indices in the denominator of a frequency ratio, otherwise empty."""

        self.num_coeff = []
        """Coefficients in a linear combination or numerator of a frequency ratio."""

        self.den_coeff = []
        """Coefficients in the denominator of a frequency ratio, otherwise empty."""

    def add_num(self,j,coeff):
        """
        Append the given index and coefficient to the list of numerator indices and
        coefficients.

        :param j: index of the mode
        :param coeff: coefficient used in the frequency combination

        :type j: int
        :type coeff: float
        """

        self.num_index.append(j)
        self.num_coeff.append(coeff)

    def add_den(self,j,coeff):
        """
        Append the given index and coefficient to the list of denominator indices and
        coefficients.

        :param j: index of the mode
        :param coeff: coefficient used in the frequency combination

        :type j: int
        :type coeff: float
        """

        self.den_index.append(j)
        self.den_coeff.append(coeff)

    def print_me(self):
        """ Print frequency combination. """

        print("value = "+str(self.value))
        print("*******************")
        print("num = "+str(self.num))
        for i, coeff in zip(self.num_index,self.num_coeff):
            print("  "+str(i)+"  "+str(coeff))
        print("*******************")
        print("den = "+str(self.den))
        for i, coeff in zip(self.den_index,self.den_coeff):
            print("  "+str(i)+"  "+str(coeff))

class Likelihood:
    """
    A class which described the likelihood function and allows users to evaluate it.
    """

    def __init__(self):
        self.modes = []
        """List of pulsation modes (of type :py:class:`Mode`)."""

        self.constraints = []
        """List of constraints which intervene in the likelihood function."""

        self.cov = None
        """Covariance matrix which intervenes when calculating frequency combinations."""

        self.invcov = None
        """Inverse of covariance matrix, :py:data:`Likelihood.cov`."""

        self.combinations = []
        """This contains indices and coefficients to frequency combinations and frequency ratios."""

        self.seismic_weight = 1.0
        """Absolute weight to be applied to seismic constraints"""

        self.classic_weight = 1.0
        """Absolute weight to be applied to classic constraints (incl. nu_max constraint)."""

        # the following variables are used for Fortran style optimisation
        self.nvalues = None
        """Array with the n values of the observed modes"""

        self.lvalues = None
        """Array with the l values of the observed modes"""

        self.fvalues = None
        """Array with the observed frequencies"""

        self.dfvalues = None
        """Array with the error bars on the observed frequencies"""

        self.values = None
        """1D float array with the value for each frequency combination"""

        self.ncoeff = None
        """
        2D int array with the number of terms for each frequency
        combination. The indices are:

          1. The type of term (0 = num, 1 = den)
          2. The index of the frequency combination
        """

        self.coeff = None
        """
        3D float array with the coefficients for each frequency
        combination. The indices are:

          1. The index of the term
          2. The type of term (0 = num, 1 = den)
          3. The index of the frequency combination
        """

        self.indices = None
        """
        3D int array with the mode indices for each frequency combination.
        The indices are:

          1. The index of the term
          2. The type of term (0 = num, 1 = den)
          3. The index of the frequency combination
        """

    def find_weights_new(self):
        """
        Find absolute weights for seismic and classic constraints based on options
        in ``AIMS_configure.py``.
        """
        nclassic = len(self.constraints)
        nseismic = len(self.combinations)
        config.find_weights(nclassic,nseismic)
        self.seismic_weight = config.seismic_weight
        self.classic_weight = config.classic_weight

    def find_weights(self):
        """
        Find absolute weights for seismic and classic constraints based on options
        in ``AIMS_configure.py``.
        """

        if (config.weight_option is None):
            self.seismic_weight = self.classic_weight = 1.0
            return

        if (config.weight_option == "Absolute"):
            self.seismic_weight = config.seismic_weight
            self.classic_weight = config.classic_weight
            return

        if (config.weight_option == "Relative"):
            nseismic = len(self.combinations)
            nclassic = len(self.constraints)
            self.seismic_weight = config.seismic_weight*float(nclassic)/float(nseismic)
            self.classic_weight = config.classic_weight
            return

        if (config.weight_option == "Reduced"):
            nseismic = len(self.combinations)
            nclassic = len(self.constraints)
            self.seismic_weight = 1.0/float(nclassic + nseismic - 1)
            self.classic_weight = 1.0/float(nclassic + nseismic - 1)
            return

        if (config.weight_option == "Reduced_bis"):
            nseismic = len(self.combinations)
            nclassic = len(self.constraints)
            self.seismic_weight = 1.0/float(nseismic - 1)
            self.classic_weight = 1.0/float(nclassic - 1)
            return

        sys.exit("Unknown weight_option: "+config.weight_option)


    def sort_modes(self):
        """
        Sort the modes.  The ordering will depend on the value of ``use_n`` from the
        ``AIMS_configure.py`` file.
        """

        if (config.use_n):
            self.modes.sort(key=attrgetter("l","n","freq"))
        else:
            self.modes.sort(key=attrgetter("l","freq","n"))

    def create_mode_arrays(self):
        """
        Create arrays with mode parameters (n, l, freq), which can
        be interfaced with fortran methods more easily.
        """

        size_obs = len(self.modes)
        self.lvalues = np.empty((size_obs,),dtype=model.ltype)
        self.nvalues = np.empty((size_obs,),dtype=model.ntype)
        self.fvalues = np.empty((size_obs,),dtype=model.ftype)
        self.dfvalues = np.empty((size_obs,),dtype=model.ftype)

        for i in range(size_obs):
            self.lvalues[i] = self.modes[i].l
            self.nvalues[i] = self.modes[i].n
            self.fvalues[i] = self.modes[i].freq
            self.dfvalues[i] = self.modes[i].dfreq

    def read_constraints(self,filename,factor=1.0):
        """
        Read a file with pulsation data and constraints.

        :param filename: name of file with pulsation data.
        :param factor: multiplicative factor for pulsation frequencies.  Can
          be used for conversions.

        :type filename: string
        :type factor: float
        """

        # read modes:
        self.modes = []

        obsfile = open(filename,"r")
        for line in obsfile:
            line = utilities.trim(line.strip())
            columns = line.split()
            if (len(columns) == 0): continue

            # distinguish between constraints and input frequencies by detecting
            # whether the first line is a number or a string:
            if (utilities.is_number(columns[0])):
                if (config.read_n):
                    if (len(columns) >= 4):
                        # l, n, f, df
                        self.modes.append(Mode(_n=int(columns[1]),_l=int(columns[0]), \
                             _freq=utilities.to_float(columns[2]),                    \
                             _dfreq=utilities.to_float(columns[3])*factor))
                else:
                    if (len(columns) >= 3):
                        # l, f, df
                        self.modes.append(Mode(_n=0, _l=int(columns[0]),    \
                             _freq=utilities.to_float(columns[1]),          \
                             _dfreq=utilities.to_float(columns[2])*factor))
            else:
                if (len(columns) < 3): continue
                if (utilities.is_number(columns[1])):
                    like.add_constraint((columns[0],           \
                        Distribution("Gaussian",               \
                        utilities.my_map(utilities.to_float,columns[1:]))))
                else:
                    like.add_constraint((columns[0],           \
                        Distribution(columns[1],               \
                        utilities.my_map(utilities.to_float,columns[2:]))))

        # don't forget to do this:
        obsfile.close()

        # print number of modes:
        print("I found %d modes in %s."%(len(self.modes),filename))

        # if the n values have not been provided, then they need to be guessed:
        if (not config.read_n): self.guess_n()

        # this needs to be done so that self.find_map() works correctly:
        self.sort_modes()
        self.create_mode_arrays()

    def add_constraint(self,constraint):
        """
        Add a supplementary constraint to the list of constraints.

        :param constraint: supplementary constraint
        :type constraint: (string, :py:class:`Distribution`)
        """

        (name,distribution) = constraint

        # look for single letter constraints and "translate" them to their full name:
        if (name == "T"): name = "Teff"
        if (name == "L"): name = "Luminosity"
        if (name == "R"): name = "Radius"
        if (name == "M"): name = "M_H"
        if (name == "G"): name = "log_g"
        self.constraints.append((name,distribution))

    def guess_dnu(self,with_n=False):
        """
        Guess the large frequency separation based on the radial modes.

        :param with_n: specifies whether to use the n values
          already stored with each mode, when calculating the
          large frequency separation.
        :type with_n: boolean

        :return: the large frequency separation
        :rtype: float
        """

        nmodes = len(self.modes)

        if (with_n):
            ind = np.arange(nmodes)
        else:

            # find index with mode order according to (l,freq):
            ind = np.lexsort(([mode.freq for mode in self.modes], \
                                 [mode.l for mode in self.modes]))  # sorts by l, then freq

            # find minimun frequency difference between radial modes as
            # a first estimate of the large frequency separation:
            min_diff = np.nan
            for i in range(1,nmodes):
                if (self.modes[ind[i]].l != 0): break
                diff = self.modes[ind[i]].freq - self.modes[ind[i-1]].freq
                # the following condition is "nan-resistant":
                if (not min_diff < diff): min_diff = diff

            # easy exit:
            if (min_diff != min_diff):
                print("WARNING: cannot find large frequency separation")
                return np.nan

        # refine large frequency separation using a least squares approach
        # on the radial modes, using temporary radial orders based on the
        # previous estimate of the large frequency separation.
        ntemp = 0
        one = n = nn = nu = nnu = 0.0
        for i in range(nmodes):
            if (self.modes[ind[i]].l != 0): break
            if ((i > 0) and (not with_n)):
                ntemp += max(1,int(round((self.modes[ind[i]].freq-self.modes[ind[i-1]].freq)/min_diff)))
            if (with_n): ntemp = self.modes[ind[i]].n
            cnst  = 1.0/self.modes[ind[i]].dfreq**2
            one   += cnst
            n     += cnst*ntemp
            nn    += cnst*ntemp*ntemp
            nu    += cnst*self.modes[ind[i]].freq
            nnu   += cnst*self.modes[ind[i]].freq*ntemp
        if (nn*one-n**2 == 0):
            return np.nan
        else:
            return (nnu*one-nu*n)/(nn*one-n**2)

    def guess_n(self):
        """
        Guess the radial order of the observed pulsations modes.

        This method uses the large frequency separation, as calculated with
        :py:meth:`guess_dnu`, to estimate the radial orders.  These orders are
        subsequently adjusted to avoid multiple modes with the same
        identification.  The resultant radial orders could be off by a constant
        offset, but this is not too problematic when computing frequency
        combinations or ratios.
        """

        # obtain large frequency separation
        # (and sort modes according to (l,freq))
        dnu = self.guess_dnu(with_n=False)

        # easy exit:
        if (np.isnan(dnu)): return

        # assign radial orders
        for mode in self.modes:
            mode.n = int(round(mode.freq/dnu - 1.0 - mode.l/2.0))

        # avoid having two modes with the same n value:
        ind = np.lexsort(([mode.freq for mode in self.modes], \
                             [mode.l for mode in self.modes]))  # sorts by l, then freq
        for i in range(len(self.modes)-2,0,-1):
            if ((self.modes[ind[i]].l == self.modes[ind[i+1]].l) and (self.modes[ind[i]].n == self.modes[ind[i+1]].n)):
                for j in range(i,0,-1):
                    if (self.modes[ind[j]].l != self.modes[ind[i]].l): break
                    self.modes[ind[j]].n -= 1

    def assign_n(self,my_model):
        """
        Assign the radial orders based on proximity to theoretical frequencies
        from an input model.

        :param my_model: input model
        :type my_model: :py:class:`model.Model`
        """

        mode_map, nmissing = self.find_map(my_model, False)

        # easy exit:
        if (nmissing > 0):
            print("WARNING: unable to assign radial orders due to unmatched modes")
            return

        # assign radial orders
        for i in range(len(self.modes)):
            self.modes[i].n = my_model.modes['n'][mode_map[i]]

    def clear_seismic_constraints(self):
        """
        This clears the seismic constraints.  Specifically, the list of seismic
        combinations, and associated covariance matrix and its inverse are
        reinitialised.
        """

        self.cov = None
        self.invcov = None
        self.combinations = []

    def add_seismic_constraint(self,string):
        """
        Add seismic contraints based on the keyword given in ``string``.

        :param string: keyword which specifies the type of constraint to be added.
          Current options include:

          - ``nu``: individual frequencies
          - ``nu0``: individual frequencies (radial modes only)
          - ``nu_min0``: radial mode with minimum frequency
          - ``r02``: :math:`r_{02}` frequency ratios
          - ``r01``: :math:`r_{01}` frequency ratios
          - ``r10``: :math:`r_{10}` frequency ratios
          - ``dnu``: individual large frequency separations (using all modes)
          - ``dnu0``: individual large frequency separations (using radial modes only)
          - ``d2nu``: second differences (using all modes)
          - ``avg_dnu``: average large frequency separation (using all modes)
          - ``avg_dnu0``: average large frequency separation (using radial modes only)

        :type string: string
        """

        if (string == "nu"):
            self.add_combinations(num_list=((0,0,1.0),))
            return

        if (string == "nu0"):
            self.add_combinations(num_list=((0,0,1.0),),target_ell=0)
            return

        if (string == "nu_min0"):
            self.add_nu_min_constraint(target_ell=0)
            return

        if (string == "r02"):
            self.add_combinations(num_list=((0,0,1.0),(-1,2,-1.0)),den_list=((0,1,1.0),(-1,1,-1.0)),target_ell=0)
            return

        if (string == "r01"):
            self.add_combinations(num_list=((-1,0,0.125),(-1,1,-0.5),(0,0,0.75),(0,1,-0.5),(1,0,0.125)), \
                                  den_list=((0,1,1.0),(-1,1,-1.0)),target_ell=0)
            return

        if (string == "r10"):
            self.add_combinations(num_list=((-1,0,-0.125),(0,-1,0.5),(0,0,-0.75),(1,-1,0.5),(1,0,-0.125)), \
                                  den_list=((1,-1,1.0),(0,-1,-1.0)),target_ell=1)
            return

        if (string == "dnu"):
            self.add_combinations(num_list=((0,0,1.0),(-1,0,-1.0)))
            return

        if (string == "dnu0"):
            self.add_combinations(num_list=((0,0,1.0),(-1,0,-1.0)),target_ell=0)
            return

        if (string == "d2nu"):
            self.add_combinations(num_list=((1,0,1.0),(0,0,-2.0),(-1,0,1.0)))
            return

        if (string == "avg_dnu"):
            self.add_dnu_constraint(l_targets=None)
            return

        if (string == "avg_dnu0"):
            self.add_dnu_constraint(l_targets=[0])
            return

        print("Unknown type of seismic constraint: "+string)
        print("Skipping ...")
        return

    def find_l_list(self, l_targets):
        """
        Find a list of l values with the following properties:

        - each l value only occurs once
        - each l value given in the parameter ``l_targets`` is
          in the result l list, except if there is 1 or less
          modes with this l value
        - if the parameter ``l_targets`` is ``None``, look for
          all l values with 2 or more modes associated with them

        :param l_targets: input list of l values
        :type l_targets: list of int

        :return: new list of l values with the above properties
        :rtype: list of int
        """

        # copy l_targets list to avoid modifying argument.  Make sure
        # each l value only occurs once in the copy.  Also, handle
        # the case l_targets is None:
        l_list = []
        if (l_targets is None):
            for mode in self.modes:
                if (mode.l not in l_list): l_list.append(mode.l)
        else:
            for l in l_targets:
                if (l not in l_list): l_list.append(l)

        # count the number of modes with different l values
        l_num = [0.0,]*len(l_list)
        for mode in self.modes:
            if (mode.l in l_list): l_num[l_list.index(mode.l)] += 1

        # remove l values from list if there are less than 2 modes with
        # that value (iterate backwards to avoid problems with shifting
        # indices).
        for i in range(len(l_list)-1,-1,-1):
            if (l_num[i] < 2): del l_num[i], l_list[i]
        return l_list

    def add_dnu_constraint_matrix(self, l_targets=[0]):
        """
        Add the large frequency separation as a contraint.  The coefficients
        are obtained via a least-squares approach.  The approach taken here
        has two advantages:

        1. Correlations between the large frequency separation and other
           seismic constraints will be taken into account.
        2. The same modes will be used in the same way, both for the
           observations and the models.

        :param l_targets: specifies for which l values the large frequency
          separation is to be calculated.  If ``None`` is supplied, all
          modes will be used.
        :type l_targets: list of int

        .. note::
          This uses a matrix approach and is therefore *not* the
          prefered method.
        """

        l_list = self.find_l_list(l_targets)

        # easy exit:
        if (len(l_list) == 0):
            print("WARNING: unable to find large frequency separation for given l targets.")
            print("Try including other l values in l_targets, or expanding your set of")
            print("observations.")
            return

        # set up linear system and invert it:
        n   = len(l_list) + 1
        mat = np.zeros((n,n),dtype=np.float64)
        v   = np.zeros(n,dtype=np.float64)
        for mode in self.modes:
            if (mode.l not in l_list): continue
            cnst = 1.0/(mode.dfreq*mode.dfreq)
            ndx  = l_list.index(mode.l) + 1
            mat[0,0]     += cnst*mode.n*mode.n
            mat[0,ndx]   += cnst*mode.n
            mat[ndx,0]   += cnst*mode.n
            mat[ndx,ndx] += cnst
            v[0]         += cnst*mode.freq*mode.n
            v[ndx]       += cnst*mode.freq
        invmat = np.linalg.inv(mat)

        # create associated mode combination
        a_combination = Combination()
        for j in range(len(self.modes)):
            if (self.modes[j].l not in l_list): continue
            ndx = l_list.index(self.modes[j].l) + 1
            a_combination.add_num(j,(invmat[0,0]*self.modes[j].n + invmat[0,ndx])/self.modes[j].dfreq**2)

        # find dnu:
        value  = np.dot(invmat[0,:],v)
        a_combination.value = value
        a_combination.num   = value
        a_combination.den   = 1.0

        print("Number of added seismic constraints:  1")
        self.combinations.append(a_combination)

    def add_dnu_constraint(self, l_targets=[0]):
        """
        Add the large frequency separation as a contraint.  The coefficients
        are obtained via a least-squares approach.  The approach taken here
        has two advantages:

        1. Correlations between the large frequency separation and other
           seismic constraints will be taken into account.
        2. The same modes will be used in the same way, both for the
           observations and the models.

        :param l_targets: specifies for which l values the large frequency
          separation is to be calculated.  If ``None`` is supplied, all
          modes will be used.
        :type l_targets: list of int

        .. note::
          This uses an analytical approach and is therefore the
          prefered method.
        """

        l_list = self.find_l_list(l_targets)

        # easy exit:
        if (len(l_list) == 0):
            print("WARNING: unable to find large frequency separation for given l targets.")
            print("Try including other l values in l_targets, or expanding your set of")
            print("observations.")
            return

        # find various summed quantities:
        n  = len(l_list)
        m  = np.zeros((n,2),dtype=np.float64)
        nn = 0.0
        for mode in self.modes:
            if (mode.l not in l_list): continue
            cnst = 1.0/(mode.dfreq*mode.dfreq)
            ndx  = l_list.index(mode.l)
            m[ndx,0] += cnst
            m[ndx,1] += cnst*mode.n
            nn += cnst*mode.n*mode.n

        det = nn
        for i in range(n): det -= m[i,1]*m[i,1]/m[i,0]

        # create associated mode combination
        a_combination = Combination()
        value = 0.0
        for j in range(len(self.modes)):
            if (self.modes[j].l not in l_list): continue
            ndx = l_list.index(self.modes[j].l)
            coeff = (self.modes[j].n - m[ndx,1]/m[ndx,0])/(det*self.modes[j].dfreq**2)
            a_combination.add_num(j,coeff)
            value += coeff * self.modes[j].freq

        # find dnu:
        a_combination.value = value
        a_combination.num   = value
        a_combination.den   = 1.0

        print("Number of added seismic constraints:  1")
        self.combinations.append(a_combination)

    def add_nu_min_constraint(self, target_ell=0, min_n=False):
        """
        Add the minimun frequency/mode of a specific ell value as a seismic
        constraint.  Typically, such constraints are used as an "anchor"
        when combined with constraints based on frequency ratios.

        :param target_ell: ell value of the minimum frequency/mode
        :param min_n: if ``False``, look for minimum observational
          frequency.  If ``True``, look for minimum radial order.

        :type target_ell: int
        :type min_n: boolean
        """

        ndx = -1
        min_val = np.nan
        for j in range(len(self.modes)):
            if (self.modes[j].l != target_ell): continue
            if (min_n):
                value = self.modes[j].n
            else:
                value = self.modes[j].freq

            # a nan-resistant condition:
            if (value >= min_val): continue
            min_val = value
            ndx = j

        if (ndx != -1):
            a_combination = Combination()
            a_combination.add_num(ndx,1.0)
            a_combination.value = self.modes[ndx].freq
            a_combination.num   = self.modes[ndx].freq
            a_combination.den   = 1.0
            self.combinations.append(a_combination)
            print("Number of added seismic constraints:  1")

    def add_combinations(self,num_list,den_list=[],target_ell=None):
        """
        This finds the indices of modes which intervene in a frequency combination or
        ratio, as specified by the mandatory and optional arguments.  These indices,
        the relevant coefficients, the numerator, the denominator, and the resultant
        value of the combination are stored in the :py:data:`combinations` variable.

        :param num_list: list of relative mode identifications and coefficients used
          to define a frequency combination or the numerator of a frequency ratio.
          This list contains tuples of the form (delta n, delta l, coeff).
        :param den_list: list of relative mode identifications and coefficients used
          to define the denominator of a frequency ratio.  If absent, then, it is
          assumed that a linear combination of frequencies is represented.  The form
          is the same as for ``num_list``.
        :param target_ell: this is used to impose a specific l value on the first
          selected mode.

        :type num_list: list of (int,int,float)
        :type den_list: list of (int,int,float)
        :type target_ell: int
        """

        # find indices of modes which follow the specified pattern
        n = 0
        for i in range(len(self.modes)):
            if ((target_ell is not None) and (self.modes[i].l != target_ell)): continue

            a_combination = Combination()
            skipMode = False

            num = 0.0
            for (n_diff, l_diff, coeff) in num_list:
                target_mode = Mode(self.modes[i].n+n_diff,self.modes[i].l+l_diff,0.0,0.0)
                for j in range(len(self.modes)):
                    if (target_mode.match(self.modes[j])):
                        a_combination.add_num(j,coeff)
                        num += coeff*self.modes[j].freq
                        break
                else:
                    skipMode = True
                    break

            if (skipMode): continue

            den = 0.0
            for (n_diff, l_diff, coeff) in den_list:
                target_mode = Mode(self.modes[i].n+n_diff,self.modes[i].l+l_diff,0.0,0.0)
                for j in range(len(self.modes)):
                    if (target_mode.match(self.modes[j])):
                        a_combination.add_den(j,coeff)
                        den += coeff*self.modes[j].freq
                        break
                else:
                    skipMode = True
                    break

            if (skipMode): continue

            if (len(den_list) > 0):
                a_combination.value = num/den
                a_combination.num   = num
                a_combination.den   = den
            else:
                a_combination.value = num
                a_combination.num   = num
                a_combination.den   = 1.0

            self.combinations.append(a_combination)
            n += 1
        print("Number of added seismic constraints: "+str(n))

    def find_covariance(self):
        """
        This prepares the covariance matrix and its inverse based on the frequency
        combinations in :py:data:`combinations`.

        .. warning::
          This method should be called *after* all of the methods which
          add to the list of frequency combinations.
        """

        n = len(self.combinations)
        self.cov = np.zeros((n,n),dtype=np.float64)

        for i in range(n):
            vec1 = self.find_vec(self.combinations[i])
            for j in range(n):
                vec2 = self.find_vec(self.combinations[j])
                for k in range(len(self.modes)):
                    self.cov[i][j] += vec1[k]*vec2[k]*self.modes[k].dfreq*self.modes[k].dfreq

        self.invcov = np.linalg.inv(self.cov)

    def find_vec(self,a_combination):
        """
        This finds a set of coefficients which intervene when constructing the
        coviance matrix for frequency combinations.

        :param a_combination: variable which specifies the frequency combination.
        :type a_combination: :py:class:`Combination`

        :return: the above set of coefficients
        :rtype: np.array
        """

        vec = np.zeros(len(self.modes),dtype=np.float64)
        if (len(a_combination.den_index) == 0):
            for j, coeff in zip(a_combination.num_index, a_combination.num_coeff):
                vec[j] += coeff
        else:

            num = 0.0
            for j, coeff in zip(a_combination.num_index, a_combination.num_coeff):
                vec[j] += coeff/a_combination.den
            for j, coeff in zip(a_combination.den_index, a_combination.den_coeff):
                vec[j] -= coeff*a_combination.value/a_combination.den

        return vec

    def create_combination_arrays(self):
        """
        Create array form of frequency combinations to be used with a
        fortran based routine for calculating the seismic chi^2 value.
        """

        # find appropriate array dimensions:
        ncomb = len(self.combinations)
        nmax  = 0
        for comb in self.combinations:
            n0 = len(comb.num_index)
            if (nmax < n0): nmax = n0
            n1 = len(comb.den_index)
            if (nmax < n1): nmax = n1

        self.values = np.empty((ncomb),dtype=model.ftype)
        self.ncoeff = np.empty((2,ncomb),dtype=np.int)
        self.coeff  = np.empty((nmax,2,ncomb),dtype=model.ftype)
        self.indices= np.empty((nmax,2,ncomb),dtype=np.int)

        for i in range(ncomb):
            self.values[i] = self.combinations[i].value

            self.ncoeff[0,i] = len(self.combinations[i].num_index)
            self.ncoeff[1,i] = len(self.combinations[i].den_index)

            n0 = len(self.combinations[i].num_index)
            n1 = len(self.combinations[i].den_index)
            self.coeff[0:n0,0,i] = np.array(self.combinations[i].num_coeff)
            self.coeff[0:n1,1,i] = np.array(self.combinations[i].den_coeff)

            self.indices[0:n0,0,i] = np.array(self.combinations[i].num_index)
            self.indices[0:n1,1,i] = np.array(self.combinations[i].den_index)

    def apply_constraints(self,my_model):
        """
        Calculate a :math:`\chi^2` value for the set of constraints (excluding
        seismic constraints based on mode frequencies).

        :param my_model: model for which the :math:`\chi^2` value is being calculated
        :type my_model: :py:class:`model.Model`

        :return: the :math:`\chi^2` value deduced from classic constraints
        :rtype: float
        """

        chi2 = 0.0
        for constraint in self.constraints:
            (string, distrib) = constraint
            chi2 += distrib(my_model.string_to_param(string))
        return chi2

    def find_map(self, my_model, use_n):
        """
        This finds a map which indicates the correspondance between observed
        modes and theoretical modes from ``my_model``.

        :param my_model: model for which the :math:`\chi^2` value is being calculated
        :param use_n: specify whether to use the radial order when finding the map
                     from observed modes to theoretical modes.  If ``False``, the map
                     is based on frequency proximity.

        :type m_model: :py:class:`model.Model`
        :type use_n: boolean

        :return: the correspondance between observed and theoretical modes from the
          above model, and the number of observed modes which weren't mapped onto
          theoretical modes
        :rtype: list of int, int

        .. note::
          - a value of -1 is used to indicate that no theoretical mode corresponds to
            a particular observed mode.
          - only zero or one observed mode is allowed to correspond to a theoretical
            mode
        """

        mode_map = np.empty((len(self.modes),),dtype=np.int)

        # sanity checks:
        if (len(self.modes) == 0): return mode_map, 0
        if (my_model.modes.shape[0] == 0): return mode_map, len(self.modes)

        if (use_n): # No, this is not a mistake - you wan't to use the local value of use_n, not the
                    # one from AIMS_configure.py.
            # NOTE: this assumes the observed modes are sorted according to (l,n)
            return aims_fortran.find_map_n(self.nvalues,self.lvalues,my_model.modes['n'],\
                                 my_model.modes['l'],mode_map)

        else:
            # NOTE: this assumes the observed modes are sorted according to (l,freq)

            # sorts by l, then freq
            ind = np.lexsort((np.copy(my_model.modes['freq']),
                                 np.copy(my_model.modes['l'])))

            # print(aims_fortran.find_map_freq(self.fvalues,self.lvalues,               \
            #                      my_model.modes['freq']*my_model.glb[model.ifreq_ref], \
            #                      my_model.modes['l'],ind,mode_map))
            return aims_fortran.find_map_freq(self.fvalues,self.lvalues,               \
                                 my_model.modes['freq']*my_model.glb[model.ifreq_ref], \
                                 my_model.modes['l'],ind,mode_map)

    def compare_frequency_combinations(self, my_model, mode_map, a=[]):
        """
        This finds a :math:`\chi^2` value based on a comparison of frequencies
        combinations, as defined in the :py:data:`combinations` variable.

        :param my_model: model for which the :math:`\chi^2` value is being calculated
        :param mode_map: a mapping which relates observed modes to theoretical ones
        :param a: parameters of surface correction terms

        :type my_model: :py:class:`model.Model`
        :type mode_map: list of int
        :type a: array-like

        :return: the :math:`\chi^2` value for the seismic constraints
        :rtype: float

        .. note::
          I'm assuming none of the modes are missing (i.e. that mode_map doesn't
          contain the value -1)
        """

        # sanity check:
        n = len(self.combinations)
        if (n == 0): return 0.0

        freq = my_model.get_freq(surface_option=config.surface_option,a=a) \
             * my_model.glb[model.ifreq_ref]

        dvalues = np.zeros((n,),dtype=model.ftype)
        dvalues = aims_fortran.compare_frequency_combinations(freq,mode_map,
                  self.values,self.ncoeff,self.coeff,self.indices,dvalues)

        return -np.dot(dvalues,np.dot(self.invcov,dvalues))/2.0

    def get_optimal_surface_amplitudes(self, my_model, mode_map):
        """
        Find optimal surface correction amplitude, for the surface correction
        specified by ``surface_option``.

        :param my_model: the model for which we're finding the surface correction amplitude
        :param mode_map: a mapping which relates observed modes to theoretical ones

        :type my_model: :py:class:`model.Model`
        :type mode_map: list of int

        :return: optimal surface correction amplitudes
        :rtype: np.array
        """

        n = len(self.modes)

        nsurf_free = nsurf
        if (config.surface_option == "Kjeldsen2008_2"): nsurf_free -= 1
        if (config.surface_option == "Sonoi2015_2"):    nsurf_free -= 1

        X = np.zeros((nsurf_free,n),dtype=np.float64)
        Y = np.zeros((n,),dtype=np.float64)

        for i in range(nsurf_free):
            a = np.zeros(nsurf,dtype=np.float64)
            a[i] = 1.0

            if (config.surface_option == "Kjeldsen2008_2"): a[-1] = my_model.b_Kjeldsen2008
            if (config.surface_option == "Sonoi2015_2"):    a[-1] = my_model.beta_Sonoi2015

            dsurf = my_model.get_surface_correction(config.surface_option,a)
            for j in range(n): X[i,j] = dsurf[mode_map[j]]/self.modes[j].dfreq

        for i in range(n):
            Y[i] = (self.modes[i].freq/my_model.glb[model.ifreq_ref] \
                 - my_model.modes['freq'][mode_map[i]]) / self.modes[i].dfreq

        mat = np.zeros((nsurf_free,nsurf_free),dtype=np.float64)
        b   = np.zeros((nsurf_free,),dtype=np.float64)
        for i in range(nsurf_free):
            b[i] = np.dot(X[i,:],Y)
            for j in range(i,nsurf_free):
                mat[i,j] = np.dot(X[i,:],X[j,:])
                mat[j,i] = mat[i,j]

        # if mat is singular, use np.linalg.lstsq(mat,b) instead
        result = np.linalg.solve(mat,b)

        # if need be append the fixed parameter:
        if (config.surface_option == "Kjeldsen2008_2"):
            result = np.append(result,my_model.b_Kjeldsen2008)
        if (config.surface_option == "Sonoi2015_2"):
            result = np.append(result,my_model.beta_Sonoi2015)

        return result


    def evaluate(self, my_model):
        """
        Calculate ln of likelihood function (i.e. a :math:`\chi^2` value) for a given
        model.

        :param my_model: model for which the :math:`\chi^2` value is being calculated
        :type my_model: :py:class:`model.Model`

        :return: the :math:`\chi^2` value, optionally the optimal surface amplitudes
           (depending on the value of :py:data:`AIMS_configure.surface_option`),
           integers to indicate whether the model was rejected due to classic
           or seismic contraints
        :rtype: float, np.array (optional), int, int

        .. note::
          This avoids model interpolation and can be used to gain time.
        """

        chi2 =  self.classic_weight*self.apply_constraints(my_model)
        reject_classic = 0
        if (chi2 < threshold): reject_classic = 1
        mode_map, nmissing = self.find_map(my_model, config.use_n and config.read_n)
        if (config.surface_option is None):
            if (nmissing > 0): return log0, reject_classic, 1
            chi2 += self.seismic_weight*self.compare_frequency_combinations(my_model,mode_map)
            return chi2, reject_classic, 0
        else:
            if (nmissing > 0):
                print(nmissing)
                return log0, [0.0]*nsurf, reject_classic, 1
            optimal_amplitudes = self.get_optimal_surface_amplitudes(my_model, mode_map)
            chi2 += self.seismic_weight*self.compare_frequency_combinations(my_model,mode_map,a=optimal_amplitudes)
            return chi2, optimal_amplitudes, reject_classic, 0

    def __call__(self, params):
        """
        Calculate ln of likelihood function (i.e. a :math:`\chi^2` value) for a given
        set of parameters.

        :param params: set of parmaeters for which the :math:`\chi^2` value is being
          calculated.
        :type params: array-like

        :return: the :math:`\chi^2` value
        :rtype: float

        .. note::
          If surface corrections are applied, then the last element(s)
          of params is/are assumed to be surface correction amplitude(s).
        """

        if (params is None): return log0
        if (config.interp_type == "Age"):
            my_model = model.interpolate_model(grid,params[0:ndims-nsurf],grid.tessellation,grid.ndx)
            if (my_model is None): return log0
            log_slope = 0.0
        elif (config.interp_type == "mHe"):
            my_model, slope = model.interpolate_model_mHe(grid,params[0:ndims-nsurf],grid.tessellation,grid.ndx)
            if (my_model is None): return log0
            log_slope = math.log(slope)
        else:
            sys.exit("An unexpected error occurred.  Please contact the authors of AIMS.")

        mode_map, nmissing = self.find_map(my_model, config.use_n)
        if (nmissing > 0): return log0
        chi2 = self.seismic_weight*self.compare_frequency_combinations(my_model,mode_map,a=params[ndims-nsurf:ndims])
        chi2 += self.classic_weight*self.apply_constraints(my_model)

        return chi2

    def is_outside(self, params):
        """
        Test to see if the given set of parameters lies outside the grid of
        models.  This is done by evaluate the probability and seeing if
        the result indicates this.

        :param params: input set of parameters
        :type params: array-like

        :return: ``True`` if the set of parameters corresponds to a point
           outside the grid.
        :rtype: boolean
        """

        # the following condition is "nan-resistant":
        return not (self(params) >= threshold)

class Probability:
    """
    A class which combines the priors and likelihood function, and allows the
    the user to evalute ln of the product of these.
    """

    def __init__(self,_priors,_likelihood):
        """
        :param _priors: input set of priors
        :param _likelihood: input likelihood function

        :type _priors: :py:class:`Prior_list`
        :type _likelihood: :py:class:`Likelihood`
        """

        self.priors = _priors
        """The set of priors."""

        self.likelihood = _likelihood
        """The likelihood function."""

    def evaluate(self, my_model):
        """
        Evalulate the ln of the product of the priors and likelihood function,
        i.e. the probability, for a given model, to within an additive
        constant.

        :param my_model: input model
        :type my_model: :py:class:`model.Model`

        :return: the ln of the probability, integers indicating if the model has bee
                 rejected based on classic constraints, seismic constraints, and/or priors
        :rtype: float, int, int, int

        .. note::
          This avoids model interpolation and can be used to gain time.
        """

        if (config.surface_option is None):
            result1, reject_classic, reject_seismic = self.likelihood.evaluate(my_model)
            result2 = self.priors(utilities.my_map(my_model.string_to_param,grid_params_MCMC))
        else:
            result1, surf_amplitudes, reject_classic, reject_seismic = self.likelihood.evaluate(my_model)
            result2 = self.priors(utilities.my_map(my_model.string_to_param,grid_params_MCMC)+list(surf_amplitudes))

        reject_prior = 0
        if (result2 < threshold): reject_prior = 1
        return result1 + result2, reject_classic, reject_seismic, reject_prior

    def __call__(self, params):
        """
        Evalulate the ln of the product of the priors and likelihood function,
        i.e. the probability, for a given model, to within an additive
        constant.

        :return: the ln of the probability
        :rtype: float

        :param params: input set of parameters
        :type params: array-like
        """

        result1 = self.likelihood(params)
        result2 = self.priors(params)
        return result1 + result2

def check_configuration():
    """
    Test the values of the variables in check_configuration to make
    sure they're acceptable.  If an unacceptable value is found, then
    this will stop AIMS and explain what variable has an erroneous
    value.
    """

    assert ((config.nwalkers%2) == 0), "nwalkers should be even.  Please modify AIMS_configure.py"
    return

def write_binary_data(infile,outfile):
    """
    Read an ascii file with a grid of models, and write corresponding binary file.

    :param infile: input ascii file name
    :param outfile: output binary file name

    :type infile: string
    :type outfile: string
    """

    grid = model.Model_grid()
    grid.read_model_list(infile)
    grid.tessellate()
    output = open(outfile,"wb")
    dill.dump(grid,output)
    output.close()
    grid.plot_tessellation()

def load_binary_data(filename):
    """
    Read a binary file with a grid of models. Option to find model in grid and
    output parameters (uncomment for use).

    :param filename: name of file with grid in binary format
    :type filename: string

    :return: the grid of models
    :rtype: :py:class:`model.Model_grid`
    """

    input_data = open(filename,"rb")
    grid = dill.load(input_data)
    input_data.close()
    retessellate = grid.remove_tracks(config.track_threshold)  # filter out incomplete tracks
    if (config.distort_grid):
        grid.distort_grid()
        retessellate = True
    if (retessellate or config.retessellate):
        grid.tessellate()

    # for track in grid.tracks:
    #     for model in track.models:
    #         # print model.glb
    #         # sys.exit()
    #         # print(model.glb[1]/constants.solar_mass)
    #         if (model.glb[1] == 3.1219100000000001E+033) & (model.glb[3] == 0.0140) & (model.glb[0] == 2101.7100000000000):
    #             print model.glb[1]/constants.solar_mass
    #             print model.glb
    #             print model.get_freq()
    #             model.write_file_simple('NGC_6819_test')
    #             print model.find_large_separation()
    #             print model.numax
    #             print model.FeH
    #             model.print_me()
    # sys.exit()

    if (grid.user_params != config.user_params):
        print("WARNING: mismatch between the user_params in the binary grid file")
        print("         and in AIMS_configure.py.  Will overwrite user_params.")
        print("  Binary grid file:  ",grid.user_params)
        print("  AIMS_configure.py: ",config.user_params)

        # manually correct relevant variables:
        config.user_params = grid.user_params
        model.nglb = 8 + len(config.user_params)
        model.nlin = 5 + len(config.user_params)
        model.ifreq_ref = 5 + len(config.user_params)
        model.iradius = 6 + len(config.user_params)
        model.iluminosity = 7 + len(config.user_params)
        model.init_user_param_dict()

    return grid

def write_list_file(filename):
    """
    Write list file from which to generate binary grid.  Various filters
    can be included to reduce the number of models.

    .. note::
      This code is intended for developpers not first time users.
    """

    # write results to file
    output = open(filename,"w")
    output.write(grid.prefix+" "+grid.postfix+"\n")
    for track in grid.tracks:
        #if (int(round(track.params[0]*100.0))%4 != 0): continue # impose step of 0.04 Msun
        if (int(round(track.params[1]*100.0))%10 != 0): continue # impose step of 0.1 on alpha_MLT
        for i in range(0,len(track.models)):
            amodel = track.models[i]  # skip every other model along evolutionary track
            output.write("%20s "%(amodel.name))                       # le nom du fichier
            output.write("%22.15e "%(amodel.glb[model.imass]))        # la masse (en g)
            output.write("%22.15e "%(amodel.glb[model.iradius]))      # le rayon (en cm)
            output.write("%22.15e "%(amodel.glb[model.iluminosity]))  # la luminosité (en erg/s)
            output.write("%17.15f "%(amodel.glb[model.iz0]))          # le Z0
            output.write("%17.15f "%(amodel.glb[model.ix0]))          # le X0
            output.write("%22.15e "%(amodel.glb[model.iage]))         # l'age (en Mans)
            output.write("%16.10f "%(amodel.glb[model.itemperature])) # la température effective (en K)
            output.write("%7.5f "%(amodel.string_to_param("alpha_MLT"))) # la longueur de mélange
            output.write("%17.15f "%(amodel.string_to_param("Zs")))   # le Zs
            output.write("%17.15f "%(amodel.string_to_param("Xs")))   # le Xs
            output.write("%17.15f "%(amodel.string_to_param("Zc")))   # le Zc
            output.write("%17.15f "%(amodel.string_to_param("Xc")))   # le Xc
            output.write("%22.15e\n"%(amodel.string_to_param("Tc")))  # la température centrale
    output.close()

def find_best_model():
    """
    Scan through grid of models to find "best" model for a given probability
    function (i.e. the product of priors and a likelihood function).
    """

    aux = range(len(grid.tracks))
    # find best models in each track, in a parallelised way:
    results = my_map(find_best_model_in_track,aux)

    global best_grid_model, best_grid_params, best_grid_result, best_age_range
    global accepted_parameters, rejected_parameters
    global nreject_classic, nreject_seismic, nreject_prior

    best_grid_result = log0
    best_grid_model  = None
    accepted_parameters = []
    rejected_parameters = []
    nreject_classic = 0
    nreject_seismic = 0
    nreject_prior   = 0
    for (result,model,accepted,rejected,age_range,rc,rs,rp) in results:
        accepted_parameters += accepted
        rejected_parameters += rejected
        nreject_classic += rc
        nreject_seismic += rs
        nreject_prior   += rp
        if result == float('inf'):
            result = log0
        if (result > best_grid_result):
            best_grid_result = result
            best_grid_model  = model
            best_age_range   = age_range

    print("Number of accepted models: %d"%(len(accepted_parameters)))
    print("Number of rejected models: %d"%(len(rejected_parameters)))
    print("  Rejections from classic constraints: %d"%(nreject_classic))
    print("  Rejections from seismic constraints: %d"%(nreject_seismic))
    print("  Rejections based on priors:          %d"%(nreject_prior))

    if (best_grid_model is None):
        print("ERROR:  Unable to find a model in the grid which matches your")
        print("        constraints.  Try extending the frequency spectra of")
        print("        your models before running AIMS.   Aborting.")
        sys.exit(1)

    best_grid_params  = utilities.my_map(best_grid_model.string_to_param,grid_params_MCMC)

    if (config.surface_option is not None):
        # print(best_grid_model)
        mode_map, nmissing = prob.likelihood.find_map(best_grid_model, config.use_n and config.read_n)
        if (nmissing > 0):
            # print(nmissing,best_grid_result)
            # sys.exit()
            print("An unexpected error occured.  Please contact the authors of AIMS.")
        optimal_amplitudes  = prob.likelihood.get_optimal_surface_amplitudes(best_grid_model, mode_map)
        best_grid_params = best_grid_params + list(optimal_amplitudes)

    best_grid_params = best_grid_params + utilities.my_map(best_grid_model.string_to_param,config.output_params)

    print("Best model:  "+str(best_grid_result)+" "+str(best_grid_params))

def find_best_model_in_track(ntrack):
    """
    Scan through an evolutionary track to find "best" model for :py:data:`prob`,
    the probability function (i.e. the product of priors and a likelihood function).

    :param ntrack: number of the evolutionary track
    :type ntrack: int

    :return: the ln(probability) value, the "best" model, the parameters of accepted models,
             the parameters of rejected models, the age range of the track, and the numbers
             of models rejected due to classic constraints, seismic contraints, and priors
    :rtype: (float, :py:class:`model.Model`, 2D float list, 2D float list, float, int, int, int)
    """

    best_result_local = log0
    best_model_local  = None
    rejected_parameters_local = []
    accepted_parameters_local = []
    reject_classic = 0
    reject_seismic = 0
    reject_prior   = 0
    nmodels = len(grid.tracks[ntrack].models)
    for i in range(nmodels):
    # for model in grid.tracks[ntrack].models:
        model = grid.tracks[ntrack].models[i]
        result, rc, rs, rp = prob.evaluate(model)
        reject_classic += rc
        reject_seismic += rs
        reject_prior   += rp
        if (config.interp_type == "mHe"):
            istart = max(i-1,0)
            istop  = min(i+1,nmodels-1)
            try: log_slope = math.log(abs((grid.tracks[ntrack].models[istop].string_to_param("Age")   \
                                    - grid.tracks[ntrack].models[istart].string_to_param("Age")) \
                                    /(grid.tracks[ntrack].models[istop].string_to_param("mHe")   \
                                    - grid.tracks[ntrack].models[istart].string_to_param("mHe"))))
            except: log_slope = 0.0
        else:
            log_slope = 0.0
        result += log_slope
        if (result > best_result_local):
            best_result_local = result
            best_model_local  = model
        if (result >= threshold):
            accepted_parameters_local.append(utilities.my_map(model.string_to_param,grid_params_MCMC))
        else:
            rejected_parameters_local.append(utilities.my_map(model.string_to_param,grid_params_MCMC))

    return (best_result_local, best_model_local, accepted_parameters_local, \
            rejected_parameters_local, grid.tracks[ntrack].age_range, reject_classic, \
            reject_seismic, reject_prior)

def init_walkers():
    """
    Initialise the walkers used in emcee.

    :return: array of starting parameters
    :rtype: np.array
    """

    # set up initial distributions according to the value of config.tight_ball:

    global tight_ball_distributions
    if (config.tight_ball):

        # find the best model:
        find_best_model()

        # create probability distributions for a tight ball:
        tight_ball_distributions = Prior_list() # just reuse the same class as for the priors

        #for param_name in grid_params_MCMC:
        for i in range(ndims):
            param_name = grid_params_MCMC_with_surf[i]
            # the star notation unpacks the tuple:
            if (param_name in config.tight_ball_range):
                new_distrib = Distribution(*config.tight_ball_range[param_name])
                new_distrib.re_centre(best_grid_params[i])
                tight_ball_distributions.add_prior(new_distrib)
            else:
                if ((param_name == "Age") or (param_name == "mHe")):
                    tight_ball_distributions.add_prior(Distribution("Gaussian",[best_grid_params[i],best_age_range/30.0]))
                else:
                    if (i >= ndims-nsurf):
                        # the parameter is an surface correction parameter
                        param_range = abs(best_grid_params[i])
                    else:
                        # the parameter is not an age surface correction parameter
                        param_range = grid.range(param_name)
                        param_range = (param_range[1]-param_range[0])/30.0
                    tight_ball_distributions.add_prior(Distribution("Gaussian",[best_grid_params[i],param_range]))

    else:
       # set the initial distributions to the priors
       tight_ball_distributions = prob.priors

    if (config.PT):
        p0 = np.zeros([config.ntemps, config.nwalkers, ndims])

        for k in range(config.ntemps):
            for j in range(config.nwalkers):
                params = None
                counter = 0
                while (prob.likelihood.is_outside(params)):
                    if (counter > config.max_iter):
                        sys.exit("ERROR: too many iterations to produce walker.  Aborting")
                    params = tight_ball_distributions.realisation()
                    counter+=1
                p0[k,j,:] = params

        # include the parameters of the best model as the first walker:
        if (config.tight_ball):
            for i in range(ndims): p0[0,0,i] = best_grid_params[i]

    else:
        p0 = np.zeros([config.nwalkers, ndims])

        for j in range(config.nwalkers):
            params = None
            counter = 0
            while (prob.likelihood.is_outside(params)):
                if (counter > config.max_iter):
                    sys.exit("ERROR: too many iterations to produce walker.  Aborting")
                params = tight_ball_distributions.realisation()
                counter+=1
            p0[j,:] = params

        # include the parameters of the best model as the first walker:
        if (config.tight_ball):
            for i in range(ndims): p0[0,i] = best_grid_params[i]

    return p0

def run_emcee():
    """
    Run the emcee program.

    :return: the ``emcee`` sampler for the MCMC run
    """

    print("Number of walkers:    "+str(config.nwalkers))
    print("Number of steps:      "+str(config.nsteps))

    if (config.PT):
        ''' v1 of emcee sampler '''
    #     print("Number of temp.:      "+str(config.ntemps))
    #     sampler = emcee.PTSampler(config.ntemps, config.nwalkers, ndims, prob.likelihood, prob.priors, pool=pool)
    #
    #     # initial burn-in:
    #     for p, lnprob, lnlike in tqdm(sampler.sample(p0, iterations = config.nsteps0),total=config.nsteps0): pass
    #
    #     # production run:
    #     sampler.reset()
    #     for p, lnprob, lnlike in tqdm(sampler.sample(p, lnprob0 = lnprob, lnlike0 = lnlike, iterations = config.nsteps),total=config.nsteps):
    #         pass
    #
    # else:
    #     sampler = emcee.EnsembleSampler(config.nwalkers, ndims, prob, pool=pool)
    #
    #     # initial burn-in:
    #     p, new_prob, state = sampler.run_mcmc(p0,config.nsteps0)
    #
    #     # production run:
    #     sampler.reset()
    #     p, new_prob, state = sampler.run_mcmc(p,config.nsteps)
    #
    # # Print acceptance fraction
    # print("Mean acceptance fraction: {0:.5f}".format(np.mean(sampler.acceptance_fraction)))
    # # Estimate the integrated autocorrelation time for the time series in each parameter.
    # try:
    #     autocorr_time =  sampler.get_autocorr_time(c=1.0)
    #     print("Autocorrelation time: "+str(autocorr_time))
    # except emcee.autocorr.AutocorrError:
    #     print("Autocorrelation time not available")
        ''' v2 of emcee sampler '''
        print "Number of temp.:      ", config.ntemps
        sampler = ptemcee.Sampler(config.nwalkers, ndims, prob.likelihood, prob.priors, ntemps=config.ntemps, threads=config.nprocesses)
        # initial burn-in:
        # print(prob.priors.priors.count,prob.priors.realisation.im_self,prob())
        # sys.exit()
        # for p, lnprob, lnlike in tqdm(sampler.sample(p0,adapt=True,iterations=config.nsteps0),total=config.nsteps0): pass
	    # Test for convergence after burn in. Taken from grd349/Hacks GitHub <--- at some point add in full version so process repeated until convergence
        # or breaks if insufficient convergence to solution after burn in.
        max_conv = 1
        rhat_accept = 1.05
        steps = config.nsteps0
        for i in range(max_conv):
            for p, lnprob, lnlike in tqdm(sampler.sample(p0,adapt=True,iterations=steps),total=steps): pass
            ''' Calculation of Rhat parameter (Gelman-Rubin diagnostic) to test walker convergence. '''
            med = np.median(sampler.chain[0,:,-steps:-1,:],axis=0)
            W = np.sum(np.var(sampler.chain[0,:,-steps:-1,:],axis=0))/len(med) # within chain variance
            mean_med = np.sum(med, axis=0)/len(med)
            B = (steps/(len(med)-1)) * np.sum((med - mean_med)**2)
            Var = (1 - 1.0/steps)*W + (B/steps)
            rhat = np.sqrt(Var/W)
            print "Rhat = ", rhat
            if np.all(rhat < rhat_accept):
                print "Sufficient convergence after burn-in (", steps,"): ", rhat," < ", rhat_accept
                break
            else:
                print "Insufficient convergence(", steps,"): ", rhat," > ", rhat_accept
                steps += config.add_steps

        # production run:
        sampler.reset()
        for p, lnprob, lnlike in tqdm(sampler.sample(p,adapt=True,storechain=True,iterations=config.nsteps),total=config.nsteps): pass

    else:
        sampler = emcee.EnsembleSampler(config.nwalkers, ndims, prob, pool=pool)
        max_conv = 5
        rhat_accept = 1.05
        steps = config.nsteps0
        for i in range(max_conv):
            # initial burn-in:
            p, new_prob, state = sampler.run_mcmc(p0,config.nsteps0)
            ''' Calculation of Rhat parameter (Gelman-Rubin diagnostic) to test walker convergence. '''
            med = np.median(sampler.chain[0,:,-steps:-1,:],axis=0)
            W = np.sum(np.var(sampler.chain[0,:,-steps:-1,:],axis=0))/len(med) # within chain variance
            mean_med = np.sum(med, axis=0)/len(med)
            B = (steps/(len(med)-1)) * np.sum((med - mean_med)**2)
            Var = (1 - 1.0/steps)*W + (B/steps)
            rhat = np.sqrt(Var/W)
            print "Rhat = ", rhat
            # conv = np.std(med, axis=0) / np.median(med, axis=0)
            # sys.exit()
            if np.all(rhat < rhat_accept):
                print "Sufficient convergence after burn-in (", steps,"): ", rhat," < ", rhat_accept
                break
            else:
                print "Insufficient convergence(", steps,"): ", rhat," > ", rhat_accept
                steps += config.add_steps

        # production run:
        sampler.reset()
        p, new_prob, state = sampler.run_mcmc(p,config.nsteps)

    # Print acceptance fraction
    print("Mean acceptance fraction: {0:.5f}".format(np.mean(sampler.acceptance_fraction)))
    # Estimate the integrated autocorrelation time for the time series in each parameter.
    print "Autocorrelation time: ", sampler.get_autocorr_time()

    return sampler

def find_blobs(samples):
    """
    Find blobs (i.e. supplementary output parameters) from a set of samples
    (i.e. for multiple models).

    :param samples: input set of samples
    :type samples: list/array of array-like

    :return: set of supplementary output parameters
    :rtype: np.array
    """

    blobs = my_map(find_a_blob,samples)
    return np.asarray(blobs)

def find_a_blob(params):
    """
    Find a blob (i.e. supplementary output parameters) for a given set of parameters
    (for one model).  The blob also includes the log(P) value as a first entry.

    :param params: input set of parameters
    :type params: array-like

    :return: list of supplementary output parameters
    :rtype: list of floats
    """

    if config.interp_type == "Age":
        my_model = model.interpolate_model(grid,params[0:ndims-nsurf],grid.tessellation,grid.ndx)
    elif config.interp_type == "mHe":
        my_model, slope = model.interpolate_model_mHe(grid,params[0:ndims-nsurf],grid.tessellation,grid.ndx)
    return utilities.my_map(my_model.string_to_param,config.output_params)

def write_samples(filename, labels, samples):
    """
    Write raw samples to a file.

    :param filename: name of file in which to write the samples
    :param labels:   names of relevant variables (used to write a header)
    :param samples:  samples for which statistical properties are calculated

    :type filename: string
    :type labels:   list of strings
    :type samples:  array-like
    """

    # write results to file
    (m,n) = np.shape(samples)
    output_file = open(filename,"w")
    output_file.write("#")
    for label in labels: output_file.write(' {0:22}'.format(label))
    output_file.write("\n ")
    for i in range(m):
        for j in range(n):
            output_file.write(' {0:22.15e}'.format(samples[i,j]))
        output_file.write("\n ")
    output_file.close()

def write_statistics(filename, labels, samples):
    """
    Write statistical properties based on a sequence of realisations to a file.
    The results include:

    - average values for each variable (statistical mean)
    - error bars for each variable (standard mean deviation)
    - correlation matrix between the different variables

    :param filename: name of file in which to write the statistical properties
    :param labels:   names of relevant variables
    :param samples:  samples for which statistical properties are calculated

    :type filename: string
    :type labels:   list of strings
    :type samples:  np.array
    """

    # initialisation
    (m,n) = np.shape(samples)
    average = samples.sum(axis=0, dtype=np.float64)/(1.0*m)
    covariance = np.zeros((n,n), dtype=np.float64)

    # find covariance matrix:
    for i in range(n):
        for j in range(i+1):
            covariance[i,j] = np.dot(samples[:,i]-average[i],samples[:,j]-average[j])/(1.0*m)
            covariance[j,i] = covariance[i,j]

    # write results to file
    output_file = open(filename,"w")
    output_file.write("Summary\n=======\n");
    for i in range(n):
        output_file.write('{0:25} {1:25.15e} {2:25.15e}\n'.format(labels[i], \
                          average[i],math.sqrt(covariance[i,i])))
    output_file.write("\nCorrelation matrix\n==================\n");

    output_file.write('{0:25} '.format(" "))
    for i in range(n):
        output_file.write('{0:25} '.format(labels[i]))
    output_file.write("\n")

    for i in range(n):
        output_file.write('{0:25} '.format(labels[i]))
        for j in range(n):
            output_file.write('{0:25.15e} '.format(covariance[i,j] \
                        /math.sqrt(covariance[i,i]*covariance[j,j])))
        output_file.write("\n")
    output_file.close()

def write_percentiles(filename, labels, samples):
    """
    Write percentiles based on a sequence of realisations to a file.
    The results include:

      +/- 2 sigma error bars (using the 2.5th and 97.5th percentiles)
      +/- 1 sigma error bars (using the 16th and 84th percentiles)
      the median value (i.e. the 50th percentile)

    :param filename: name of file in which to write the percentiles
    :param labels:   names of relevant variables
    :param samples:  samples for which statistical properties are calculated

    :type filename: string
    :type labels:   list of strings
    :type samples:  np.array
    """

    # n sigma values (from https://en.wikipedia.org/wiki/Normal_distribution)
    one_sigma = 0.682689492137
    two_sigma = 0.954499736104
    #three_sigma = 0.997300203937

    # initialisation
    (m,n) = np.shape(samples)
    percentiles = np.empty((n,7), dtype=np.float64)
    for i in range(n):
        percentiles[i,0] = np.percentile(samples[:,i],100.0*(0.5-two_sigma/2.0))
        percentiles[i,4] = np.percentile(samples[:,i],100.0*(0.5+two_sigma/2.0))
        percentiles[i,1] = np.percentile(samples[:,i],100.0*(0.5-one_sigma/2.0))
        percentiles[i,3] = np.percentile(samples[:,i],100.0*(0.5+one_sigma/2.0))
        percentiles[i,2] = np.percentile(samples[:,i],50.0)
        percentiles[i,5] = abs((percentiles[i,1]-percentiles[i,3])/2.)
        percentiles[i,6] = abs((percentiles[i,0]-percentiles[i,4])/2.)

    # write results to file
    output_file = open(filename,"w")
    output_file.write("Percentiles\n===========\n\n");
    output_file.write('{0:25} '.format("Quantity"))
    output_file.write('{0:25} '.format("-2sigma"))
    output_file.write('{0:25} '.format("-1sigma"))
    output_file.write('{0:25} '.format("Median"))
    output_file.write('{0:25} '.format("+1sigma"))
    output_file.write('{0:25}'.format("+2sigma"))
    output_file.write('{0:25} '.format("1sig"))
    output_file.write('{0:25}\n'.format("2sig"))

    for i in range(n):
        output_file.write('{0:25} '.format(labels[i]))
        for j in range(7):
            output_file.write('{0:25.15e} '.format(percentiles[i,j]))
        output_file.write('\n')
    output_file.close()

def write_LEGACY_summary(filename, KIC, labels, samples):
    """
    Write a one line summary of the statistical properties based on a
    sequence of realisations to a file.  The format matches that of the
    LEGACY project.

    The results include:

    - average values for each variable (statistical mean)
    - error bars for each variable (standard mean deviation)

    :param filename: name of file in which to write the statistical properties
    :param KIC: KIC number of the star
    :param labels:   names of relevant variables
    :param samples:  samples for which statistical properties are calculated

    :type filename: string
    :type KIC: string
    :type labels:   list of strings
    :type samples:  np.array
    """

    # initialisation
    (m,n) = np.shape(samples)
    average = samples.sum(axis=0, dtype=np.float64)/(1.0*m)
    uncertainties = np.zeros((n,), dtype=np.float64)

    # find covariance matrix:
    for i in range(n):
        uncertainties[i] = math.sqrt(np.dot(samples[:,i]-average[i],samples[:,i]-average[i])/(1.0*m))

    # write results to file
    output_file = open(filename,"w")
    output_file.write(KIC)

    print(n)
    print(labels)
    print(average)
    i = labels.index("Radius")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    i = labels.index("Mass")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    i = labels.index("log_g")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    i = labels.index("Rho")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    i = labels.index("Age")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    i = labels.index("Teff")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    i = labels.index("Fe_H")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    # MCore
    output_file.write(" -99.999 -99.999 -99.999")
    # RCore
    output_file.write(" -99.999 -99.999 -99.999")
    # Mbcz
    output_file.write(" -99.999 -99.999 -99.999")
    # Rbcz
    output_file.write(" -99.999 -99.999 -99.999")
    i = labels.index("Luminosity")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    i = labels.index("X")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    i = labels.index("Y")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    # Xsup
    output_file.write(" -99.999 -99.999 -99.999")
    # Ysup
    output_file.write(" -99.999 -99.999 -99.999")
    i = labels.index("Xc")
    output_file.write(" %f %f %f"%(average[i],uncertainties[i],-uncertainties[i]))
    # Ycen
    output_file.write(" -99.999 -99.999 -99.999")

    output_file.write("\n")
    output_file.close()

def write_readme(filename, elapsed_time):
    """
    Write parameters relevant to this MCMC run.

    :param filename: name of file in which to write the statistical properties
    :type filename: string
    """

    # various format related strings:
    boolean2str = ("False","True")
    str_decimal = "{0:40}{1:d}\n"
    str_string  = "{0:40}{1:40}\n"
    str_float   = "{0:40}{1:22.15e}\n"

    output_file = open(filename,"w")

    output_file.write(string_to_title("Observational constraints"))
    output_file.write(str_decimal.format("Number of modes",len(prob.likelihood.modes)))
    output_file.write("# NOTE: the radial orders may have been recalculated (see Radial orders section)\n")
    output_file.write("{0:4} {1:4} {2:17} {3:17}\n".format("l","n","freq","dfreq"))
    for mode in prob.likelihood.modes:
        output_file.write("{0:4d} {1:4d} {2:17.10e} {3:17.10e}\n".format( \
                           mode.l, mode.n, mode.freq, mode.dfreq))
    for (param,distrib) in prob.likelihood.constraints:
        output_file.write(str_string.format(param,distrib.to_string()))
    output_file.write(str_string.format("Seismic constraints",str(config.seismic_constraints)))
    output_file.write(str_string.format("Surface option",config.surface_option))
    if (config.surface_option == "Kjeldsen2008"):
        output_file.write(str_float.format("Surface exponent, b",config.b_Kjeldsen2008))
    if (config.surface_option == "Sonoi2015"):
        output_file.write(str_float.format("Surface parameter, beta",config.beta_Sonoi2015))


    output_file.write(string_to_title("Priors"))
    for name,distribution in zip(grid_params_MCMC_with_surf, prob.priors.priors):
        output_file.write(str_string.format("Prior on "+name,distribution.to_string()))

    nseismic = len(prob.likelihood.combinations)
    nclassic = len(prob.likelihood.constraints)
    output_file.write(string_to_title("Weighting"))
    output_file.write(str_string.format("Weight option",config.weight_option))
    output_file.write(str_float.format("Absolute seismic weight",prob.likelihood.seismic_weight))
    output_file.write(str_float.format("Absolute classic weight",prob.likelihood.classic_weight))
    if (nclassic != 0):
        output_file.write(str_float.format("Relative seismic weight",prob.likelihood.seismic_weight*float(nseismic)/float(nclassic)))
    else:
        output_file.write(str_float.format("Relative seismic weight",np.nan))
    output_file.write(str_float.format("Relative classic weight",prob.likelihood.classic_weight))

    output_file.write(string_to_title("Radial orders"))
    output_file.write(str_string.format("Radial orders from input file",boolean2str[config.read_n]))
    output_file.write(str_string.format("Radial orders from best model",boolean2str[config.assign_n]))
    output_file.write(str_string.format("Use radial orders in mode map",boolean2str[config.use_n]))

    output_file.write(string_to_title("The grid and interpolation"))
    output_file.write(str_string.format("Binary grid",config.binary_grid))
    output_file.write(str_string.format("User defined parameters",str([x[0] for x in config.user_params])))
    output_file.write(str_string.format("Scaled age interpolation",boolean2str[config.scale_age]))
    output_file.write(str_string.format("Distort grid prior to tessellation",boolean2str[config.distort_grid]))
    if (config.distort_grid):
        output_file.write("Distortion matrix:")
        output_file.write(str(grid.distort_mat))

    output_file.write(string_to_title("EMCEE parameters"))
    output_file.write(str_string.format("With parallel tempering",boolean2str[config.PT]))
    if (config.PT):
        output_file.write(str_decimal.format("Number of temperatures",config.ntemps))
    output_file.write(str_decimal.format("Number of walkers",config.nwalkers))
    output_file.write(str_decimal.format("Number of burn-in steps",config.nsteps0))
    output_file.write(str_decimal.format("Number of production steps",config.nsteps))
    output_file.write(str_decimal.format("Thinning parameter",config.thin))

    output_file.write(string_to_title("Initialisation"))
    output_file.write(str_string.format("Initialise with a tight ball",boolean2str[config.tight_ball]))
    if (config.tight_ball):
        for i in range(ndims):
            name = grid_params_MCMC_with_surf[i]
            output_file.write(str_string.format("Tight ball distrib. on "+name, \
                tight_ball_distributions.priors[i].to_string()))

    output_file.write(string_to_title("Output"))
    output_file.write(str_string.format("List of output parameters",str(config.output_params)))
    output_file.write(str_string.format("Write OSM files",boolean2str[config.with_osm]))
    output_file.write(str_string.format("Plot walkers",boolean2str[config.with_walkers]))
    output_file.write(str_string.format("Plot echelle diagrams",boolean2str[config.with_echelle]))
    output_file.write(str_string.format("Plot histograms",boolean2str[config.with_histograms]))
    output_file.write(str_string.format("Plot triangle plots",boolean2str[config.with_triangles]))
    output_file.write(str_string.format("List of plot extensions",str(config.plot_extensions)))
    output_file.write(str_string.format("List of triangle plot extensions",str(config.tri_extensions)))

    output_file.write(string_to_title("Miscellaneous"))
    output_file.write(str_string.format("With parallelisation",boolean2str[config.parallel]))
    if (config.parallel):
        output_file.write(str_decimal.format("Number of processes",config.nprocesses))
    output_file.write(str_float.format("Execution time (s)",elapsed_time))
    output_file.write(str_decimal.format("Number of accepted models",len(accepted_parameters)))
    output_file.write(str_decimal.format("Number of rejected models",len(rejected_parameters)))
    output_file.write(str_decimal.format("  Classic rejections",nreject_classic))
    output_file.write(str_decimal.format("  Seismic rejections",nreject_seismic))
    output_file.write(str_decimal.format("  Prior-based rejections",nreject_prior))

    output_file.close()

def write_combinations(filename,samples):
    """
    Produce a list of linear combinations of grid models (based on interpolation) corresponding
    to the provided model parameters.

    :param filename: name of the file to which to write the model combinations
    :param samples: set of model parameters for which we would like to obtain the grid models
                   and interpolation coefficients

    :type filename: string
    :type samples: np.array
    """

    output_file = open(filename,"w")

    output_file.write("{0:s} {1:e}\n".format(grid.prefix,constants.G))
    for params in samples:
        results = model.find_combination(grid,params[0:ndims-nsurf])
        if config.interp_type == "Age":
            my_model = model.interpolate_model(grid,params[0:ndims-nsurf],grid.tessellation,grid.ndx)
        if config.interp_type == "mHe":
            my_model, slope = model.interpolate_model_mHe(grid,params[0:ndims-nsurf],grid.tessellation,grid.ndx)

        if (results is None): continue  # filter out combinations outside the grid
        output_file.write("{0:d} {1:.15e} {2:.15e} {3:.15e} {4:.5f} {5:.5f} {6:.15e} {7:.5f}\n".format(
                          len(results), my_model.glb[model.imass], my_model.glb[model.iradius],
                          my_model.glb[model.iluminosity],my_model.glb[model.iz0],
                          my_model.glb[model.ix0],my_model.glb[model.iage],
                          my_model.glb[model.itemperature]))
        for (coef,model_name) in results:
           output_file.write("{0:.15f} {1:s}\n".format(coef, model_name))
        output_file.write("\n")

    output_file.close()

def write_model(my_model,my_params,my_result,model_name,extended=False):
    """
    Write text file with caracteristics of input model.

    :param my_model: model for which we're writing a text file
    :param my_params: parameters of the model
    :param my_result: ln(P) value obtained for the model
    :param model_name: name used to describe this model.  This is also used
      when naming the text file.
    :param extended: if set to True, all of the theoretical modes are
                     saved in the text file, including those not matched
                     to observations

    :type my_model: :py:class:`model.Model`
    :type my_params: array-like
    :type my_result: float
    :type model_name: string
    """

    # sanity check
    if (my_model is None): return

    # initialisations:
    str_float   = "{0:40}{1:22.15e}\n"
    param_names = grid_params_MCMC_with_surf + config.output_params
    mode_map, nmissing = prob.likelihood.find_map(my_model, config.use_n)
    mode_map = list(mode_map)

    output_file = open(os.path.join(output_folder,model_name+"_model.txt"),"w")
    output_file.write(string_to_title("Model: "+model_name))
    output_file.write(str_float.format("ln(P)",my_result))
    for name,value in zip(param_names,my_params):
        output_file.write(str_float.format(name,value))

    output_file.write(string_to_title("Frequencies"))
    if (config.surface_option is None):
        output_file.write("{0:4} {1:4} {2:17} {3:17} {4:17}\n".format("l","n","freq_theo", \
                          "freq_obs", "dfreq_obs"))
        for i in range(my_model.modes.size):
            if (i in mode_map):
                j = mode_map.index(i)
                output_file.write("{0:4d} {1:4d} {2:17.10e} {3:17.10e} {4:17.10e}\n".format( \
                           my_model.modes['l'][i], my_model.modes['n'][i], \
                           my_model.modes['freq'][i]*my_model.glb[model.ifreq_ref], \
                           prob.likelihood.modes[j].freq, prob.likelihood.modes[j].dfreq))
            else:
                if (extended):
                    output_file.write("{0:4d} {1:4d} {2:17.10e}\n".format( \
                               my_model.modes['l'][i], my_model.modes['n'][i], \
                               my_model.modes['freq'][i]*my_model.glb[model.ifreq_ref]))

    else:
        freq = my_model.get_freq(surface_option=config.surface_option, a=my_params[ndims-nsurf:ndims])
        output_file.write("{0:4} {1:4} {2:17} {3:17} {4:17} {5:17}\n".format("l","n", \
                          "freq_theo", "freq+surf. corr.", "freq_obs", "dfreq_obs"))
        for i in range(my_model.modes.size):
            if (i in mode_map):
                j = mode_map.index(i)
                output_file.write("{0:4d} {1:4d} {2:17.10e} {3:17.10e} {4:17.10e} {5:17.10e}\n".format( \
                           my_model.modes['l'][i], my_model.modes['n'][i], \
                           my_model.modes['freq'][i]*my_model.glb[model.ifreq_ref], \
                           freq[i]*my_model.glb[model.ifreq_ref], prob.likelihood.modes[j].freq, \
                           prob.likelihood.modes[j].dfreq))
            else:
                if (extended):
                    output_file.write("{0:4d} {1:4d} {2:17.10e} {3:17.10e}\n".format( \
                               my_model.modes['l'][i], my_model.modes['n'][i], \
                               my_model.modes['freq'][i]*my_model.glb[model.ifreq_ref], \
                               freq[i]*my_model.glb[model.ifreq_ref]))

    output_file.close()

def string_to_title(string):
    """
    Create fancy title from string.

    :param string: string from which the title is created.
    :type string: string

    :return: the fancy string title
    :rtype: string
    """

    n = len(string)
    ntitle = 80
    nhalf  = (ntitle - n - 6)//2
    result = "\n"+"#"*ntitle + "\n"
    result += "#"*(ntitle - nhalf - n - 6)
    result += "   "+string+"   "+"#"*nhalf+"\n"
    result += "#"*ntitle + "\n"
    return result

def write_osm_frequencies(filename, my_model):
    """
    Write file with frequencies for Optimal Stellar Model (OSM), written
    by R. Samadi.

    :param filename: name of file which will contain the frequencies
    :type filename: string

    :param my_model: model from which are derived the radial orders
    :type my_model: :py:class:`model.Model`

    .. note::
      Written by B. Herbert.
    """

    # initialisations:
    mode_map, nmissing = prob.likelihood.find_map(my_model, config.use_n)

    output_file = open(os.path.join(config.output_osm,filename+"_freq.dat"),"w")
    output_file.write("# {0:4} {1:4} {2:17} {3:17} \n".format("n", "l",\
                          "nu[muHz]", "sigma [muHz]"))
    for i in range(len(mode_map)):
        j = mode_map[i]
        output_file.write("{0:4d} {1:4d} {2:17.10e} {3:17.10e}\n".format( \
                           my_model.modes['n'][j], my_model.modes['l'][j], \
                           prob.likelihood.modes[i].freq, \
                           prob.likelihood.modes[i].dfreq))

    output_file.close()


def write_osm_don(filename, my_model):
    """
    Write file with choice of physical ingredients to be used by CESAM or CESTAM
    and OSM.

    :param filename: name of file which will contain the physical ingredients
    :type filename: string

    :param my_model: model from which which is derived various physical constraints/settings
    :type my_model: :py:class:`model.Model`

    .. note::
      Written by B. Herbert.
    """

    # test mixing length parameter:
    if "alpha_MLT" in model.user_params_index:
        alpha = my_model.string_to_param("alpha_MLT")
    else:
        print("WARNING: no 'alpha_MLT' parameter in grid.  Using")
        print("         solar-calibrated value in OSM files.")
        alpha = 1.6 # approximate solar calibrated value

    # initialisation:
    mode_map, nmissing = prob.likelihood.find_map(my_model, config.use_n)
    output_file = open(os.path.join(config.output_osm,filename+"_seismic.don"),"w")
    output_file.write("&NL_CESAM \n")
    output_file.write(" NOM_CHEMIN = '/home/hg_central/cestam/SUN_STAR_DATA/eos_opal2005/, \n")
    output_file.write(" NOM_CTES = 'ctes_94_ags09', \n")
    output_file.write(" NOM_DES = 'no_des', \n")
    output_file.write(" NOM_OUTPUT = 'osc_adia', \n")
    output_file.write(" N_MAX = 6000, \n")
    output_file.write(" PRECISION = 'sa', \n")
    output_file.write(" / \n")
    output_file.write("&NL_MASS \n")
    output_file.write(" MTOT = %17e ,\n" %(my_model.string_to_param("Mass")))
    output_file.write(" NOM_PERTM = 'pertm_ext', \n")
    output_file.write(" MDOT = 0.0, \n")
    output_file.write(" / \n")
    output_file.write("&NL_EVOL \n")
    output_file.write(" AGEMAX = %17e , \n" %(my_model.string_to_param("Age")))
    output_file.write(" ARRET = 'else', \n")
    output_file.write(" DTLIST = 1000000000000.0, \n")
    output_file.write(" LOG_TEFF = %17e , \n" %(my_model.string_to_param("log_Teff")))
    output_file.write(" NB_MAX_MODELES = 10000, \n")
    output_file.write(" HE_CORE = -1.0, \n")
    output_file.write(" T_STOP = 30000000000.0, \n")
    output_file.write(" X_STOP = -1.0 ,\n")
    output_file.write(" / \n")
    output_file.write("&NL_CHIM \n")
    output_file.write(" GRILLE_FIXE = F, \n")
    output_file.write(" NOM_ABON = 'solaire_ags09', \n")
    output_file.write(" MODIF_CHIM = F, \n")
    output_file.write(" GARDE_XISH = F, \n")
    output_file.write(" X0 = %.15f , \n" %(my_model.string_to_param("X")))
    output_file.write(" Y0 = %.15f , \n" %(my_model.string_to_param("Y")))
    output_file.write(" ZSX0 = %.15f , \n" %(my_model.string_to_param("zsx_0")))
    output_file.write(" / \n")
    output_file.write("&NL_CONV \n")
    output_file.write(" NOM_CONV = 'conv_jmj', \n")
    output_file.write(" ALPHA = %17.15f , \n" %(alpha))
    output_file.write(" OVSHTS = 0.0, \n")
    output_file.write(" OVSHTI = 0.0, \n")
    output_file.write(" JPZ = T, \n")
    output_file.write(" CPTURB = 0.0, \n")
    output_file.write(" LEDOUX = F ,\n")
    output_file.write(" / \n")
    output_file.write("&NL_DIFF \n")
    output_file.write(" DIFFUSION = F, \n")
    output_file.write(" NOM_DIFFM = 'diffm_mp', \n")
    output_file.write(" NOM_DIFFT = 'difft_nu', \n")
    output_file.write(" D_TURB = 0.0, \n")
    output_file.write(" RE_NU = 0.0, \n")
    output_file.write(" NOM_FRAD = 'no_frad' ,\n")
    output_file.write(" / \n")
    output_file.write("&NL_ROT \n")
    output_file.write(" W_ROT = 0.0, \n")
    output_file.write(" UNIT = 'kms/s', \n")
    output_file.write(" NOM_DIFFW = 'diffw_0', \n")
    output_file.write(" NOM_DIFFMG = 'diffmg_0', \n")
    output_file.write(" NOM_THW = 'rot_0', \n")
    output_file.write(" NOM_PERTW = 'pertw_0', \n")
    output_file.write(" P_PERTW = 6.5e47, \n")
    output_file.write(" LIM_JPZ = T, \n")
    output_file.write(" NOM_DES_ROT = 'no_des', \n")
    output_file.write(" TAU_DISK = 0.0, \n")
    output_file.write(" P_DISK = 0.0, \n")
    output_file.write(" W_SAT = 8.0, \n")
    output_file.write(" / \n")
    output_file.write("&NL_ETAT \n")
    output_file.write(" NOM_ETAT = 'etat_opal5Z', \n")
    output_file.write(" F_EOS = 'eos_opal5_cur_0.013067.bin', ' ', ' ', ' ', ' ', ' ', ' ', ' ', \n")
    output_file.write(" / \n")
    output_file.write("&NL_OPA \n")
    output_file.write(" NOM_OPA = 'opa_yveline', \n")
    output_file.write(" F_OPA = 'opa_yveline_ags09.bin', ' ', ' ', ' ', ' ', ' ', ' ', ' ', \n")
    output_file.write(" / \n")
    output_file.write("&NL_NUC \n")
    output_file.write(" NOM_NUC = 'ppcno3a9', \n")
    output_file.write(" NOM_NUC_CPL = 'NACRE_LUNA', \n")
    output_file.write(" MITLER = F , \n")
    output_file.write(" / \n")
    output_file.write("&NL_ATM \n")
    output_file.write(" NOM_ATM = 'lim_atm', \n")
    output_file.write(" NOM_TDETAU = 'edding', \n")
    output_file.write(" TAU_MAX = 20.0, \n")
    output_file.write(" LIM_RO = F , \n")
    output_file.write(" / \n")
    output_file.close()

def write_osm_xml(filename,my_params, my_model):
    """
    Write file with classic constraints for OSM

    :param filename: name of file with classic constraints
    :type filename: string

    :param my_model: model used in deriving some of the constraints
    :type my_model: :py:class:`model.Model`

    .. note::
      Originally written by B. Herbert.  Includes some modifications.
    """

    # test mixing length parameter:
    if "alpha_MLT" in model.user_params_index:
        alpha = my_model.string_to_param("alpha_MLT")
    else:
        alpha = 1.6 # approximate solar calibrated value

    # initialisations:
    mode_map, nmissing = prob.likelihood.find_map(my_model, config.use_n)

    biblio_osm = {"T":"teff","Teff":"teff","L":"l","Luminosity":"l","R":"r","Radius":"r","G":"logg","log_g":"logg"}

    config_osm = etree.Element("config")

    # write free parameters:
    append_osm_parameter(config_osm, "agemax", my_model.string_to_param("Age"), \
                         20.0, 5.0, (20.0, 13819.0))
    append_osm_parameter(config_osm, "mtot", my_model.string_to_param("Mass"), \
                         0.03, 3.0, (0.5, 5.0))
    append_osm_parameter(config_osm, "zsx0", my_model.string_to_param("zsx_0"), \
                         0.0001, 5.0, (-0.3, 0.3))
    append_osm_parameter(config_osm, "alpha", alpha, 0.03, 5.0, (1.0, 2.5))

    # treat surface option:
    if config.surface_option == "Kjeldsen2008":
        append_osm_parameter(config_osm, "se_a", my_params[ndims-nsurf] \
            *(my_model.numax/my_model.glb[model.ifreq_ref])**(config.b_Kjeldsen2008 - 1.0), \
            0.0001, 5.0, (-1.0, 0.0))

    if config.surface_option == "Kjeldsen2008_scaling":
        append_osm_parameter(config_osm, "se_a", my_params[ndims-nsurf] \
            *(my_model.numax/my_model.glb[model.ifreq_ref])**(my_model.b_Kjeldsen2008 - 1.0), \
            0.0001, 5.0, (-1.0, 0.0))

    if config.surface_option == "Kjeldsen2008_2":
        append_osm_parameter(config_osm, "se_a", my_params[ndims-nsurf] \
            *(my_model.numax/my_model.glb[model.ifreq_ref])**(my_params[ndims-nsurf+1] - 1.0), \
            0.0001, 5.0, (-1.0, 0.0))
        append_osm_parameter(config_osm, "se_b", my_params[ndims-nsurf+1], \
            0.1, 5.0, (1.0, 6.0))

    if "Sonoi2015" in config.surface_option:
        append_osm_parameter(config_osm, "se_a", my_params[ndims-nsurf] \
            *(my_model.glb[model.ifreq_ref]/my_model.numax), \
            0.0001, 5.0, (-1.0, 0.0))

    if config.surface_option == "Sonoi2015_2":
        # Parameter a has already been introduced above.
        # Hence, only introduce parameter b:
        append_osm_parameter(config_osm, "se_b", my_params[ndims-nsurf+1], \
            0.1, 5.0, (1.0, 8.0))

    for (param,distrib) in prob.likelihood.constraints:
        if (param == "M_H"):
            target_osm = etree.Element("target")
            target_osm.set("name", "log_zsx_s")
            value_osm = etree.SubElement(target_osm, "value")
            value_osm.text = "%22.15e"%(distrib.mean + np.log10(constants.solar_z/constants.solar_x))
            sigma_osm = etree.SubElement(target_osm, "sigma")
            sigma_osm.text = "%22.15e"%(distrib.error_bar)
            config_osm.append(target_osm)
        else:
            if param in biblio_osm:
                target_osm = etree.Element("target")
                target_osm.set("name", biblio_osm[param])
                value_osm = etree.SubElement(target_osm, "value")
                value_osm.text = "%22.15e"%(distrib.mean)
                sigma_osm = etree.SubElement(target_osm, "sigma")
                sigma_osm.text = "%22.15e"%(distrib.error_bar)
                config_osm.append(target_osm)
            else:
                print("WARNING: unused constraint in OSM file: %s"%(param))

    seismic_constraints_osm = etree.Element("seismic_constraints")
    fich_osm = etree.SubElement(seismic_constraints_osm, "file")
    fich_osm.text = filename+"_freq.dat"
    types_osm = etree.SubElement(seismic_constraints_osm, "types")
    types_osm.text = "nu"
    matching_osm = etree.SubElement(seismic_constraints_osm, "matching")
    matching_osm.text = "order"
    config_osm.append(seismic_constraints_osm)

    setting_osm = etree.Element("settings")
    levmar_osm = etree.SubElement(setting_osm, "levmar")
    maxiter_osm = etree.SubElement(levmar_osm, "maxiter")
    maxiter_osm.text = "%d"%(50)
    chi2min_osm = etree.SubElement(levmar_osm, "chi2min")
    chi2min_osm.text = "%15.8e"%(0.1)
    ftol_osm = etree.SubElement(levmar_osm, "ftol")
    ftol_osm.text = "%15.8e"%(0.001)
    autostep_osm = etree.SubElement(levmar_osm, "autostep")
    autostep_osm.text = "%d"%(0)
    cov_cdtnb_thr_osm = etree.SubElement(levmar_osm, "cov_cdtnb_thr")
    cov_cdtnb_thr_osm.text = "%15.8e"%(1e13)
    hess_cdtnb_thr_osm = etree.SubElement(levmar_osm, "hess_cdtnb_thr")
    hess_cdtnb_thr_osm.text = "%15.8e"%(1e13)

    modes_osm = etree.SubElement(setting_osm, "modes")

    if config.surface_option == "Kjeldsen2008":
        append_osm_surface_effects(modes_osm, "kb2008", my_model.numax, (my_params[ndims-nsurf] \
            *(my_model.numax/my_model.glb[model.ifreq_ref])**(config.b_Kjeldsen2008 - 1.0), \
            config.b_Kjeldsen2008))

    if config.surface_option == "Kjeldsen2008_scaling":
        append_osm_surface_effects(modes_osm, "kb2008", my_model.numax, (my_params[ndims-nsurf] \
            *(my_model.numax/my_model.glb[model.ifreq_ref])**(my_model.b_Kjeldsen2008 - 1.0), \
            my_model.b_Kjeldsen2008))

    if config.surface_option == "Kjeldsen2008_2":
        append_osm_surface_effects(modes_osm, "kb2008", my_model.numax, (my_params[ndims-nsurf] \
            *(my_model.numax/my_model.glb[model.ifreq_ref])**(my_params[ndims-nsurf+1] - 1.0), \
            my_params[ndims-nsurf+1]))

    if config.surface_option == "Sonoi2015":
        append_osm_surface_effects(modes_osm, "lorentz2", my_model.numax, (my_params[ndims-nsurf] \
            *(my_model.glb[model.ifreq_ref]/my_model.numax) , config.beta_Sonoi2015))

    if config.surface_option == "Sonoi2015_scaling":
        append_osm_surface_effects(modes_osm, "lorentz2", my_model.numax, (my_params[ndims-nsurf] \
            *(my_model.glb[model.ifreq_ref]/my_model.numax) , my_model.beta_Sonoi2015))

    if config.surface_option == "Sonoi2015_2":
        append_osm_surface_effects(modes_osm, "lorentz2", my_model.numax, (my_params[ndims-nsurf] \
            *(my_model.glb[model.ifreq_ref]/my_model.numax) , my_params[ndims-nsurf+1]))

    oscprog_osm = etree.SubElement(modes_osm, "oscprog")
    oscprog_osm.text = "adipls"

    models_osm = etree.SubElement(setting_osm, "models")
    dy_dz_osm = etree.SubElement(models_osm, "dy_dz")
    dy_dz_osm.text = "%15.8e"%(0.44)
    yp_osm = etree.SubElement(models_osm, "yp")
    yp_osm.text = "%15.8e"%(0.2475)
    zp_osm = etree.SubElement(models_osm, "zp")
    zp_osm.text = "%15.8e"%(0.0)
    start_osm = etree.SubElement(models_osm, "start")
    start_osm.text = "pms"
    cf_osm = etree.SubElement(models_osm, "cf")
    cf_osm.text = "%15.8e"%(8e-5)
    config_osm.append(setting_osm)

    output_file = open(os.path.join(config.output_osm,filename+"_seismic.xml"),"w")
    output_file.write(etree.tostring(config_osm, pretty_print=True))
    output_file.close()

def append_osm_parameter(config_osm, name, value, step, rate, bounds):
    """
    Add a parameter in xlm format in the file with the classic
    constraints for OSM.

    :param config_osm: XLM etree element to which to add the parameter
    :type config_osm: lxml.etree._Element

    :param name: name of the parameter
    :type name: string

    :param value: value of the parameter
    :type value: float

    :param step: parameter step (this intervenes when numerically calculating
                 derivatives with respect to this parameter)
    :type step: float

    :param rate: parameter rate (this corresponds to a tolerance on this
                 parameter)
    :type rate: float

    :param bounds: bounds on the parameter
    :type bounds: float tuple
    """

    parameter_osm = etree.Element("parameter")
    parameter_osm.set("name", name)
    value_osm = etree.SubElement(parameter_osm, "value")
    value_osm.text = "%15.8e"%(value)
    step_osm = etree.SubElement(parameter_osm, "step")
    step_osm.text = "%15.8e"%(step)
    rate_osm = etree.SubElement(parameter_osm, "rate")
    rate_osm.text = "%15.8e"%(rate)
    bounds_osm = etree.SubElement(parameter_osm, "bounds")
    bounds_osm.text = "%15.8e , %15.8e"%(bounds[0], bounds[1])
    config_osm.append(parameter_osm)

def append_osm_surface_effects(modes_osm, name, numax, values):
    """
    Add a method with which to calculate surface effects to the
    OSM contraint file.

    :param modes_osm: XML element to which to add the surface effects method
    :type modes_osm: lxml.etree._Element

    :param name: name of the method
    :type name: string

    :param numax: value of numax
    :type numax: float

    :param values: values which intervene in the method
    :type values: float tuple
    """

    surf_effects_osm = etree.SubElement(modes_osm, "surface_effects")
    formula_osm = etree.SubElement(surf_effects_osm, "formula")
    formula_osm.text = name
    params_surf_osm = etree.SubElement(surf_effects_osm, "parameters")
    params_surf_osm.text = "%15.8e , %15.8e"%(values[0],values[1])
    numax_osm = etree.SubElement(surf_effects_osm, "numax")
    numax_osm.text = "%15.8e"%(numax)

def echelle_diagram(my_model,my_params,model_name):
    """
    Write text file with caracteristics of input model.

    :param my_model: model for which we're writing a text file
    :param my_params: parameters of the model
    :param model_name: name used to describe this model.  This is also used
      when naming the text file.

    :type my_model: :py:class:`model.Model`
    :type my_params: array-like
    :type model_name: string
    """

    # sanity checks
    if (my_model is None): return
    if (len(prob.likelihood.modes) == 0): return

    # initialisations:
    lmax = int(max([mode.l for mode in prob.likelihood.modes]))
    msize = 6
    mode_map, nmissing = prob.likelihood.find_map(my_model, config.use_n)
    dnu = prob.likelihood.guess_dnu(with_n=True)
    if (dnu < 1.0):
        dnu_str = "{0:7.2e}".format(dnu)
    elif (dnu < 10.0):
        dnu_str = "{0:4.2f}".format(dnu)
    elif (dnu < 100.0):
        dnu_str = "{0:4.1f}".format(dnu)
    elif (dnu < 1000.0):
        dnu_str = "{0:5.1f}".format(dnu)
    else:
        dnu_str = "{0:7.2e}".format(dnu)

    nu_obs = [None]*(lmax+1)
    error_obs = [None]*(lmax+1)
    nu_theo = [None]*(lmax+1)
    style = ["d","s","^","o"]  # this should be modified if lmax is modified

    freq = my_model.get_freq(surface_option=None)*my_model.glb[model.ifreq_ref]
    if (config.surface_option is not None):
        nu_theo_surf = [None]*(lmax+1)
        freq_surf = my_model.get_freq(surface_option=config.surface_option, \
                    a=my_params[ndims-nsurf:ndims]) \
                  * my_model.glb[model.ifreq_ref]

    # obtain frequency sets (and limits)
    nu_min = 1e300
    nu_max =-1e300
    for l in range(lmax+1):
        nu_obs[l] = np.array([mode.freq for mode in prob.likelihood.modes if mode.l == l])
        error_obs[l] = np.array([mode.dfreq for mode in prob.likelihood.modes if mode.l == l])
        nu_theo[l] = np.array([freq[mode_map[j]] for j in range(len(mode_map)) \
                      if (mode_map[j] != -1) and (my_model.modes['l'][mode_map[j]] == l)])
        if (config.surface_option is not None):
            nu_theo_surf[l] = np.array([freq_surf[mode_map[j]] for j in range(len(mode_map)) \
                               if (mode_map[j] != -1) and (my_model.modes['l'][mode_map[j]] == l)])

        if (nu_obs[l].shape[0] > 0):
            nu_min = min((nu_min,np.nanmin(nu_obs[l])))
            nu_max = max((nu_max,np.nanmax(nu_obs[l])))

        if (nu_theo[l].shape[0] > 0):
            nu_min = min((nu_min,np.nanmin(nu_theo[l])))
            nu_max = max((nu_max,np.nanmax(nu_theo[l])))

    nu_diff = nu_max - nu_min
    nu_min -= 0.05*nu_diff
    nu_max += 0.05*nu_diff

    # make the echelle diagram:
    # NOTE: (b)lue, (g)reen, (r)ed, (c)yan, (m)agenta, (y)ellow, (k) black, (w)hite

    plt.figure()
    plt.plot((dnu,dnu),(nu_min,nu_max),"k:")
    for l in range(lmax+1):
        plt.errorbar(nu_obs[l]%dnu, nu_obs[l], xerr=error_obs[l], fmt="r"+style[l], markersize=msize)

    if (config.surface_option is not None):
        for l in range(lmax+1):
            plt.plot(nu_theo_surf[l]%dnu, nu_theo_surf[l], "c"+style[l], markersize=msize)

    for l in range(lmax+1):
        plt.plot(nu_theo[l]%dnu, nu_theo[l], "b"+style[l], markersize=msize)

    plt.title("Model: "+model_name)
    plt.xlabel(r"Reduced frequency, $\nu$ mod "+dnu_str+r" $\mu$Hz")
    plt.ylabel(r"Frequency, $\nu$ (in $\mu$Hz)")
    plt.xlim((0.0,dnu*1.4))
    plt.ylim((nu_min,nu_max))

    # create a legend
    handles = [mpatches.Patch(color='red')]
    labels  = [r'$\nu_{\mathrm{obs}}$']
    handles.append(mpatches.Patch(color='blue'))
    labels.append(r'$\nu_{\mathrm{theo}}$')
    if (config.surface_option is not None):
        handles.append(mpatches.Patch(color='cyan'))
        labels.append(r'$\nu_{\mathrm{theo}}^{\mathrm{surf}}$')
    for l in range(lmax+1):
        handles.append(mlines.Line2D([],[], color='white', marker=style[l], markersize=msize, \
                                     linestyle='None',drawstyle='steps-mid'))
        labels.append(r'$\ell='+"{0:d}".format(l)+r'$')
    plt.legend(handles, labels)

    # save plot
    for ext in config.plot_extensions:
        plt.savefig(os.path.join(output_folder,"echelle_"+model_name.replace(" ", "_")+"."+ext))

    plt.clf()

def plot_walkers(samples, labels, filename, nw=3):
    """
    Plot individual walkers.

    :param samples: samples from the emcee run
    :param labels: labels for the different dimensions in parameters space
    :param filename: specify name of file in which to save plots of walkers.
    :param nw: number of walkers to be plotted

    :type samples: np.array
    :type labels: list of strings
    :type filename: string
    :type nw: int

    .. warning::
      This method must be applied before the samples are reshaped,
      and information on individual walkers lost.
    """

    plt.figure()
    itr = np.arange(config.nsteps)
    for i in range(nw):
        for j in range(ndims):
            #plt.subplot(ndims*100+nw*10+1+i+j*nw)
            plt.subplot(ndims,nw,1+i+j*nw)
            plt.title(labels[j])
            plt.plot(itr,samples[i,:,j],'-')
    for ext in config.plot_extensions:
        plt.savefig(filename+ext)

    plt.clf()

def plot_distrib_iter(samples, labels, folder):
    """
    Plot individual distribution of walkers as a function of iterations.

    :param samples: samples from the emcee run
    :param labels: labels for the different dimensions in parameters space
    :param folder: specify name of file in which to save plots of walkers.

    :type samples: np.array
    :type labels: list of strings
    :type folder: string

    .. warning::
      This method must be applied before the samples are reshaped,
      and information on individual walkers lost.
    """


    mid_values  = np.empty((config.nsteps,),dtype=np.float64)
    yfill       = np.empty((2*config.nsteps,),dtype=np.float64)
    xfill       = np.array(list(range(config.nsteps))+list(range(config.nsteps-1,-1,-1)))
    for i in range(ndims):
        plt.figure()
        for j in range(config.nsteps):
            # print(samples[:,j,i])
            mid_values[j] = np.percentile(samples[:,j,i],50.0)
            yfill[j] = np.percentile(samples[:,j,i],25.0)
            yfill[-j-1] = np.percentile(samples[:,j,i],75.0)
        plt.fill(xfill,yfill,"c")
        plt.plot(xfill[:config.nsteps],mid_values,"b")
        plt.title(labels[i])
        plt.xlabel(r"Iteration, $n$")
        plt.ylabel(r"Walker distribution")
        for ext in config.plot_extensions:
            plt.savefig(os.path.join(output_folder,"distrib_iter_"+grid_params_MCMC_with_surf[i]+"."+ext))
        plt.clf()

def plot_histograms(samples, names, fancy_names, truths=None):
    """
    Plot a histogram based on a set of samples.

    :param samples: samples form the emcee run
    :param names: names of the quantities represented by the samples.  This will
      be used when naming the file with the histogram
    :param fancy_names: name of the quantities represented by the samples. This
      will be used as the x-axis label in the histogram.
    :param truths: reference values (typically the true values or some other
      important values) to be added to the histograms as a vertical line

    :type samples: np.array
    :type names: list of strings
    :type fancy_names: list of strings
    :type truths: list of floats
    """

    for i in range(len(names)):
        plt.figure()
        # print(i)
        # print(samples[:,i])
        n, bins, patches = plt.hist(samples[:,i],50,density=True,histtype='bar')
        if (truths is not None):
            ylim = plt.ylim()
            plt.plot([truths[i],truths[i]],ylim,'g-')
            plt.ylim(ylim) # restore limits, just in case
        plt.xlabel(fancy_names[i])
        for ext in config.plot_extensions:
            plt.savefig(os.path.join(output_folder,"histogram_"+names[i]+"."+ext))

    plt.clf()

def interpolation_tests(filename):
    """
    Carry out various interpolation tests and write results in
    binary format to file.

    :param filename: name of file in which to write test results.
    :type filename: string

    .. note::
      The contents of this file may be plotted using methods from
      ``plot_interpolation_test.py``.
    """

    grid = load_binary_data(config.binary_grid)
    titles = utilities.my_map(model.string_to_latex,grid.grid_params+("Age",))
    results_age1 = [track.test_interpolation(1) for track in grid.tracks]
    results_age2 = [track.test_interpolation(2) for track in grid.tracks]
    results_track, ndx1, ndx2, tessellation = grid.test_interpolation()
    output = open(filename,"wb")
    dill.dump([grid.ndim+1, model.nglb, titles, grid.grid, ndx1, ndx2, tessellation, \
               results_age1, results_age2, results_track],output)
    output.close()

def plot_frequencies(grid):
    """
    Plot frequencies along an evolutionary track
    """

    # arbitrarily choose one of the tracks:
    ntrack = 0
    track = grid.tracks[ntrack]

    nmin, nmax, lmin, lmax = track.find_mode_range()
    # nmin = 15
    # nmax = 25

    title = r"$M = %.2f M_{\odot}$, $\log(Z) = %.3f$"%(track.params[0],track.params[1])

    nages = 1000
    ages, freqs = track.find_modes(nmax,0)
    ages_ext = np.linspace(min(ages),max(ages),nages)
    freqs_ext = np.empty((nages,nmax-nmin+1,2),dtype=float)
    freqs_ext_dim = np.empty((nages,nmax-nmin+1,2),dtype=float)
    for i in range(nages):
        aModel = track.interpolate_model(ages_ext[i])
        for n in range(nmin,nmax+1):
            for l in range(2):
                freqs_ext[i,n-nmin,l] = aModel.find_mode(n,l)
                freqs_ext_dim[i,n-nmin,l] = aModel.find_mode(n,l)*aModel.glb[model.ifreq_ref]

    for n in range(nmin,nmax+1):
        if n == nmin:
            label0 = "radial"
            label1 = r"$\ell = 1$"
        else:
            label0 = label1 = None
        ages, freqs = track.find_modes(n,0)
        plt.plot(ages,freqs,"b.")
        plt.plot(ages_ext,freqs_ext[:,n-nmin,0],"b-",label=label0)
        ages, freqs = track.find_modes(n,1)
        plt.plot(ages,freqs,"r.")
        plt.plot(ages_ext,freqs_ext[:,n-nmin,1],"r-",label=label1)
    plt.legend(fontsize=12,loc=2)
    plt.xlabel(r"Stellar age (in Myrs)",fontsize=12)
    plt.ylabel(r"Frequency, $\omega/\sqrt{GM/R^3}$",fontsize=12)
    plt.title(title,fontsize=15)
    # plt.ylim(21.5,40)
    plt.savefig("freq_non_dim.pdf")
    plt.clf()
    plt.close()

    for n in range(nmin,nmax+1):
        if n == nmin:
            label0 = "radial"
            label1 = r"$\ell = 1$"
        else:
            label0 = label1 = None
        ages, freqs = track.find_modes_dim(n,0)
        plt.plot(ages,freqs,"b.")
        plt.plot(ages_ext,freqs_ext_dim[:,n-nmin,0],"b-",label=label0)
        ages, freqs = track.find_modes_dim(n,1)
        plt.plot(ages,freqs,"r.")
        plt.plot(ages_ext,freqs_ext_dim[:,n-nmin,1],"r-",label=label1)
    plt.legend(fontsize=10)
    plt.xlabel(r"Stellar age (in Myrs)",fontsize=15)
    plt.ylabel(r"Frequency, $\nu$ (in $\mu$Hz)",fontsize=15)
    plt.title(title,fontsize=15)
    plt.savefig("freq_dim.pdf")

    sys.exit(0)


if __name__ == "__main__":
    """
    AIMS = Asteroseismic Inference on a Massive Scale
    """

    # initialise the user-defined parameter dictionaries
    model.init_user_param_dict()

    # this if for writing binary data
    if (config.write_data):
        write_binary_data(config.list_grid,config.binary_grid)
        sys.exit(0)

    # this if for testing the interpolation
    if (config.test_interpolation):
        interpolation_tests(config.interpolation_file)
        sys.exit(0)

    # check number of arguments
    assert (len(sys.argv) > 1), "Usage: AIMS.py observations_file"

    # check parameters in AIMS_configure.py
    check_configuration()

    print("AIMS = Asteroseismic Inferences on a Massive Scale")
    t0 = time.time()

    # create output folder:
    output_folder = os.path.join(config.output_dir,os.path.basename(sys.argv[1]))
    if (os.path.exists(output_folder)):
        if (os.path.isfile(output_folder)):
            sys.exit('Unable to overwrite file "%s" with folder'%(output_folder))
    else:
        os.makedirs(output_folder)

    # create output folder for OSM:
    if (config.with_osm):
        if (os.path.exists(config.output_osm)):
            if (os.path.isfile(config.output_osm)):
                sys.exit('Unable to overwrite file "%s" with folder'%( \
                         config.output_osm))
        else:
            os.makedirs(config.output_osm)

    # seed random number generator (NOTE: this is not thread-safe):
    np.random.seed()

    # define likelihood function:
    like = Likelihood()
    like.read_constraints(sys.argv[1],factor=1.0)
    print("Estimated large frequency separation (in microHz): "+str(like.guess_dnu(with_n=True)))
    utilities.my_map(like.add_seismic_constraint,config.seismic_constraints)
    like.find_covariance()
    like.create_combination_arrays()
    like.find_weights()

    # load grid and associated quantities
    grid = load_binary_data(config.binary_grid)
    # plot_frequencies(grid) # this will stop the program
    grid_params_MCMC = grid.grid_params + (config.interp_type,)
    grid_params_MCMC_with_surf = grid_params_MCMC \
                               + model.get_surface_parameter_names(config.surface_option)
    nsurf            = len(model.get_surface_parameter_names(config.surface_option))
    ndims            = len(grid_params_MCMC_with_surf)

    # define priors:
    priors = Prior_list()
    for param_name in grid_params_MCMC_with_surf:
        # the star notation unpacks the tuple:
        if (param_name in config.priors):
            priors.add_prior(Distribution(*config.priors[param_name]))
        else:
            if (config.tight_ball):
                priors.add_prior(Distribution("Uninformative",[]))
            else:
                if (param_name in model.get_surface_parameter_names(config.surface_option)):
                    sys.exit("Please define prior for surface effects in configuration file.")
                else:
                    priors.add_prior(Distribution("Uniform",grid.range(param_name)))

    # combine the above:
    prob = Probability(priors,like)

    # titles, labels, etc.
    labels     = ["ln(P)"]+utilities.my_map(model.string_to_latex,grid_params_MCMC_with_surf)
    labels_big = labels + utilities.my_map(model.string_to_latex,config.output_params)
    names_big  = ("lnP",)+grid_params_MCMC_with_surf + config.output_params

    # start pool of parallel processes
    # NOTE: multiprocessing.pool.TheadPool doesn't duplicate memory
    #       like Pool, but can actually slow down execution (even
    #       compared to non-parallel execution).
    if (config.parallel):
        pool = Pool(processes = config.nprocesses)
        my_map = pool.map
    else:
        pool = None
        my_map = utilities.my_map

    # initialised walkers:
    p0 = init_walkers()

    if (config.assign_n):
        if (best_grid_model is None): find_best_model()
        like.clear_seismic_constraints()
        like.assign_n(best_grid_model)
        like.sort_modes()
        like.create_mode_arrays()
        utilities.my_map(like.add_seismic_constraint,config.seismic_constraints)
        like.find_covariance()
        like.create_combination_arrays()
        like.find_weights()

    # run emcee:
    sampler = run_emcee()

    # Collect results:
    if (config.PT):
        samples = sampler.chain[0,...]
        # lnprob  = sampler.lnprobability[0,...]
        lnprob  = sampler.logprobability[0,...]
    else:
        samples = sampler.chain
        # lnprob  = sampler.lnprobability
        lnprob  = sampler.logprobability

    # Write file with parameters used in this run
    elapsed_time = time.time() - t0
    write_readme(os.path.join(output_folder,"README"),elapsed_time)

    # Various diagnostic which must be done before reshaping the samples:
    if (config.with_walkers):
        plot_walkers(samples, labels[1:], os.path.join(output_folder,"walkers."), nw = 3)

    if (config.with_distrib_iter):
        plot_distrib_iter(samples, labels[1:], output_folder)

    # Reshape the samples and obtain auxiliary quantities:
    # NOTE: choosing order='C' leads to much better sub-sampling since the
    #       the thin_samples array is much more likely to draw from all of the
    #       walkers rather than a subset. Otherwise, if order='F' and if
    #       nwalkers is a multiple of thin, than only 1 out thin walkers
    #       is used (and this leads to poor sampling, repititions, etc.).
    samples = samples.reshape((config.nwalkers*config.nsteps,ndims),order='C')
    lnprob  = lnprob.reshape((config.nwalkers*config.nsteps,1),order='C')
    samples = np.concatenate((lnprob,samples),axis=1)

    thin_samples = samples[0::config.thin,:]
    blobs = find_blobs(thin_samples[:,1:])
    samples_big = np.concatenate((thin_samples,blobs),axis=1)

    # end pool of parallel processes
    if (config.parallel):
        pool.close()
        pool.join()

    # write various text files:
    write_samples   (os.path.join(output_folder,"samples.txt"),     labels,        samples)
    write_samples   (os.path.join(output_folder,"samples_big.txt"), labels_big,    samples_big)
    write_statistics(os.path.join(output_folder,"results.txt"),     names_big[1:], samples[:,1:])
    write_statistics(os.path.join(output_folder,"results_big.txt"), names_big[1:], samples_big[:,1:])
    write_percentiles(os.path.join(output_folder,"percentiles.txt"),     \
                      names_big[1:], samples[:,1:])
    write_percentiles(os.path.join(output_folder,"percentiles_big.txt"), \
                      names_big[1:], samples_big[:,1:])
    if (config.with_combinations):
        write_combinations(os.path.join(output_folder,"combinations.txt"),samples[0::config.thin_comb,1:])

    # write best grid model
    if (best_grid_model is not None):
        write_model(best_grid_model,best_grid_params,best_grid_result, \
                   "best_grid",extended=config.extended_model)
        write_combinations(os.path.join(output_folder,"combinations_best_grid.txt"),[best_grid_params])
        if (config.with_echelle):
            echelle_diagram(best_grid_model,best_grid_params,"Best grid")

    # find best MCMC model
    ndx_max = lnprob.argmax()
    best_MCMC_result = lnprob[ndx_max,0]
    best_MCMC_params = samples[ndx_max,1:]
    best_MCMC_params.reshape(ndims)
    if config.interp_type == "Age":
        best_MCMC_model  = model.interpolate_model(grid,best_MCMC_params[0:ndims-nsurf],grid.tessellation,grid.ndx)
    elif config.interp_type == "mHe":
        best_MCMC_model, slope  = model.interpolate_model_mHe(grid,best_MCMC_params[0:ndims-nsurf],grid.tessellation,grid.ndx)
    best_MCMC_params = list(best_MCMC_params) + utilities.my_map(best_MCMC_model.string_to_param,config.output_params)
    write_model(best_MCMC_model,best_MCMC_params,best_MCMC_result, \
               "best_MCMC",extended=config.extended_model)
    write_combinations(os.path.join(output_folder,"combinations_best_MCMC.txt"),[best_MCMC_params])
    if (config.with_echelle):
        echelle_diagram(best_MCMC_model,best_MCMC_params,"Best MCMC")

    # find statistical model
    statistical_params = samples[:,1:].sum(axis=0, dtype=np.float64)/(1.0*config.nwalkers*config.nsteps)
    statistical_params.reshape(ndims)
    statistical_result = prob(statistical_params)
    print(statistical_result)
    if config.interp_type == "Age":
        statistical_model  = model.interpolate_model(grid,statistical_params[0:ndims-nsurf],grid.tessellation,grid.ndx)
    elif config.interp_type == "mHe":
        statistical_model, slope  = model.interpolate_model_mHe(grid,statistical_params[0:ndims-nsurf],grid.tessellation,grid.ndx)
    if (statistical_model is not None):
        statistical_params = list(statistical_params) + utilities.my_map(statistical_model.string_to_param,config.output_params)
        write_model(statistical_model,statistical_params,statistical_result, \
                   "statistical",extended=config.extended_model)
        write_combinations(os.path.join(output_folder,"combinations_statistical.txt"),[statistical_params])
        if (config.with_echelle):
            echelle_diagram(statistical_model,statistical_params,"statistical")

    # write OSM files:
    if (config.with_osm):
        write_osm_frequencies(os.path.basename(sys.argv[1]), best_MCMC_model)
        write_osm_don(os.path.basename(sys.argv[1]), best_MCMC_model)
        write_osm_xml(os.path.basename(sys.argv[1]), best_MCMC_params, best_MCMC_model)

    # write one line summary for the LEGACY project
    #write_LEGACY_summary(os.path.join(output_folder,"LEGACY_summary.txt"), \
    #                     sys.argv[1][3:], names_big[1:], samples_big[:,1:])

    # make various plots:
    if (best_grid_model is None):
        plot_histograms(samples[:,0:1],["lnP"],["ln(P)"])
    else:
        plot_histograms(samples[:,0:1],["lnP"],["ln(P)"],truths=[best_grid_result])

    if (config.with_histograms):
        plot_histograms(samples_big[:,1:],names_big[1:],labels_big[1:],truths=best_grid_params)

    if (config.with_rejected):
        if (len(rejected_parameters) >= ndims-nsurf):
            fig = corner.corner(np.array(rejected_parameters), labels=labels[1:ndims-nsurf+1])
            for ext in config.tri_extensions:
                fig.savefig(os.path.join(output_folder,"rejected."+ext))
                plt.close('all')
        if (len(accepted_parameters) >= ndims-nsurf):
            fig = corner.corner(np.array(accepted_parameters), labels=labels[1:ndims-nsurf+1])
            for ext in config.tri_extensions:
                fig.savefig(os.path.join(output_folder,"accepted."+ext))
                plt.close('all')

    if (config.with_triangles):
        fig = corner.corner(samples[:,1:], truths=best_grid_params, labels=labels[1:])
        for ext in config.tri_extensions:
            fig.savefig(os.path.join(output_folder,"triangle."+ext))
            plt.close('all')
        fig = corner.corner(samples_big[:,1:], truths=best_grid_params, labels=labels_big[1:])
        for ext in config.tri_extensions:
            fig.savefig(os.path.join(output_folder,"triangle_big."+ext))
            plt.close('all')
