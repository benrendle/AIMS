#!/usr/bin/env python
# coding: utf-8
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

import numpy as np

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

def sparse_print(filename,mat):
    """
    Print a sparse matrix (for debug purposes only):

    :param filename: name of the file in which to print the matrix
    :type filename: string

    :param mat: the matrix to be printed
    :type mat: numpy array
    """

    output_file = open(filename,"w")
    it = np.nditer(mat, flags=['multi_index'])
    while not it.finished:
       if (it[0] != 0.0): output_file.write("%22.15e %s\n"%(it[0], it.multi_index))
       it.iternext()
    output_file.close()

def my_map(fct,lst):
    """
    Systematically applies a function to a list of items.  This deals with
    the python3 behaviour of map which returns a map object rather than a
    list.

    :param fct: the function to be applied to each element of a list
    :type fct: function

    :param lst: the list to which is applied the function
    :type lst: list
    """
    return list(map(fct,lst))
