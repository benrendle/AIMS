#!/usr/bin/env python

"""
A simple utility for correcting binary AIMS grids
"""

import sys
import dill
import numpy as np
import matplotlib.pyplot as plt

def write_binary_data(grid,outfile):
    """
    Read an ascii file with a grid of models, and write corresponding binary file.
    
    :param grid: grid object to be written to a file
    :param outfile: output binary file name

    :type grid: AIMS grid
    :type outfile: string
    """

    output = open(outfile,"wb")
    dill.dump(grid,output)
    output.close()

def load_binary_data(filename):
    """
    Read a binary file with a grid of models.
    
    :param filename: name of file with grid in binary format
    :type filename: string

    :return: the grid of models
    :rtype: :py:class:`model.Model_grid`
    """

    input_data = open(filename,"rb")
    grid = dill.load(input_data)
    input_data.close()
    return grid

if __name__ == "__main__":

    # check number of arguments
    assert (len(sys.argv) > 2), "Usage: test_grid.py input_grid output_grid"

    grid = load_binary_data(sys.argv[1])
    grid.user_params = (("alpha_MLT", r'Mixing length parameter, $%s\alpha_{\mathrm{MLT}}%s$'), \
               ("Zs", r'Surface metallicity, $%sZ_c%s$'), \
               ("Xs", r'Surface hydrogen, $%sX_c%s$'), \
               ("Zc", r'Central metallicity, $%sZ_c%s$'), \
               ("Xc", r'Central hydrogen, $%sX_c%s$'), \
               ("r_BCZ",r'Radius at BCZ, $%sr_{\mathrm{BCZ}}%s$'), \
               ("tau_BCZ",r'Acoustic depth of BCZ, $%s\tau_{\mathrm{BCZ}}%s$'), \
               ("tau",r'Acoustic radius, $%s\tau%s$'), \
               ("alpha_OV", r'Overshoot parameter, $%s\alpha_{\mathrm{OV}}%s$'))
    write_binary_data(grid,sys.argv[2])

