#!/usr/bin/env python
# $Id: constants.py
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
A module which contains the following physical constants:

+------------------------------+--------------------------------------+-------------------------------------+
| **Name of variable**         | **Quantity it describes**            | **Units**                           |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_radius`     | the solar radius                     | :math:`\mathrm{cm}`                 |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_mass`       | the solar mass                       | :math:`\mathrm{g}`                  |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_luminosity` | the solar luminosity                 | :math:`\mathrm{g.cm^2.s^{-3}}`      |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_temperature`| the solar effective temperature      | :math:`\mathrm{K}`                  |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_dnu`        | the solar large frequency separation | :math:`\mathrm{\mu Hz}`             |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_numax`      | the solar frequency at maximum power | :math:`\mathrm{\mu Hz}`             |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_cutoff`     | the solar cutoff frequency           | :math:`\mathrm{\mu Hz}`             |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`G`                | the gravitational constant           | :math:`\mathrm{cm^3.g^{-1}.s^{-2}}` |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_x`          | the solar hydrogen content           | dimensionless                       |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`solar_z`          | the solar metallicity content        | dimensionless                       |
+------------------------------+--------------------------------------+-------------------------------------+
| :py:const:`A_FeH`            | multiplicative constant in           | dimensionless                       |
|                              | :math:`\mathrm{[M/H]=A_{FeH}[Fe/H]}` |                                     |
+------------------------------+--------------------------------------+-------------------------------------+

.. note::
  These values can be edited according to the latest discoveries.  As
  good practise, it is helpful to include the relevant reference.
"""

__docformat__ = 'restructuredtext'

solar_radius     = 6.9599e10    # Allen (1973)
""" the solar radius in :math:`\mathrm{cm}` """

solar_mass       = 1.98919e33   # Same as in Model S
""" the solar mass in :math:`\mathrm{g}` """

solar_luminosity = 3.844e33
""" the solar luminosity in :math:`\mathrm{g.cm^2.s^{-3}}` """

solar_temperature= 5777.0
""" the solar temperature in :math:`\mathrm{K}` """

solar_dnu        = 135.1        # solar delta nu value (in \mu Hz), Huber et al. (2011)
""" the solar large frequency separation in :math:`\mathrm{\mu Hz}` """

solar_numax      = 3090.0       # solar nu_max value (in \mu Hz), Huber et al. (2011)
""" the solar frequency at maximum power in :math:`\mathrm{\mu Hz}` """

solar_cutoff     = 5300.0       # Jimenez et al. (2011) (see Balmforth & Gough 1990, Fossat et al. 1992)
""" the solar cut-off frequency separation in :math:`\mathrm{\mu Hz}` """

G                = 6.67428e-8   # CODATA 2006
#G                = 6.6716823e-8 # the gravitational constant in cm^3.g^-1.s^-2 (CoRoT/ESTA value)
""" the gravitational constant in :math:`\mathrm{cm^3.g^{-1}.s^{-2}}` """

A_FeH            = 1.0
"""
multiplicative constant which intervenes in the
conversion from metal content to iron content
"""

solar_x          = 0.7336  # Grevesse & Noels (1993)
#solar_x          = 0.7345  # Grevesse & Sauval (1998)
#solar_x          = 0.7392  # Asplund et al. (2005)
#solar_x          = 0.7381  # Asplund et al. (2009)
""" the solar hydrogen content """

solar_z          = 0.0179  # Grevesse & Noels (1993)
#solar_z          = 0.0169  # Grevesse & Sauval (1998)
#solar_z          = 0.0122  # Asplund et al. (2005)
#solar_z          = 0.0134  # Asplund et al. (2009)
""" the solar metallicity content """
