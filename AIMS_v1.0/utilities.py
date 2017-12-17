#!/usr/bin/env python
# $Id: utilities.py
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
A module which contains various utility methods for handling strings and floats.
"""

__docformat__ = 'restructuredtext'


def to_float(s):
    """
    Convert a string to a float.

    :param s: string which will be converted to a float
    :type s: string

    :return: the corresponding float
    :rtype: float

    .. note:: 
      This method allows "d" and "D" as an exponent (i.e. for
      Fortran style numbers).
    """

    return float(s.replace("d","e").replace("D","e"))

def is_number(s):
    """
    Test a string to see if it is a number.

    :param s: string which is being tested
    :type s: string

    :return: ``True`` if ``s`` is a number, and ``False`` otherwise
    :rtype: boolean

    .. note::
      This method allows "d" and "D" as an exponent (i.e. for
      Fortran style numbers).
    """

    try:
        float(s.replace("d","e").replace("D","e"))
        return True
    except ValueError:
        return False

def trim(s):
    """
    Return a string with comments (starting with "#") removed.

    :param s: the string for which we would like to remove comments.
    :type s: string

    :return: the string without comments
    :rtype: string
    """
    ndx = s.find("#")
    if (ndx == -1): ndx = len(s)
    return s[:ndx]

