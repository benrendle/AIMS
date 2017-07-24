import dill
import sys
import copy
import random
import numpy as np
import matplotlib.pyplot as plt

import plt_interp_MS as pli


def generate_grids(params, param_dict, ttls, grd, res, dim):
    ax_title = copy.deepcopy(ttls)
    for key in param_dict:
        del ax_title[key]
        const_param = ttls[key][:ttls[key].find(',')]

        for value in param_dict[key]:
            new_grid = list(copy.deepcopy(grd))
            for i in xrange(len(new_grid)):
                new_grid[i] = list(new_grid[i])

            new_results = list(copy.deepcopy(res))

            for location, item in reversed(list(enumerate(grd))):
                if item[key] != value:
                    del new_grid[location]
                    del new_results[location]

            for item in new_grid:
                del item[key]

            new_grid = np.array(new_grid)
            plot_title = "%s = %.4f Avg radial error (nincr = 1)" % (const_param, value)
            grid_title = "%s = %.4f" % (const_param, value)
            pli.surface2D(1, new_results, 0, dim, "avg", plot_title, 1, ax_title, params)
            pli.plot_grid(new_grid, dim - 1, ax_title, grid_title)


if __name__ == "__main__":
    filename = sys.argv[1]
    data = open(filename, "r")
    [ndim, nglb, titles, grid, ndx1, ndx2, tessellation, results_age1, \
     results_age2, results_track] = dill.load(data)
    data.close()

    print "\n\n\nThe available parameters are:"
    for i in xrange(len(titles) - 1):
        title = titles[i]
        comma = title.find(',')
        print "%d: %s" % (i, title[:comma])
    print "\nSelect 2 parameters to plot (the others will be kept constant)"
    print "Hint: Type the number of the parameter from the list above"

    input1 = raw_input("Parameter 1: ")
    input2 = raw_input("Parameter 2: ")

    param1 = int(input1)
    param2 = int(input2)

    assert param1 != param2, "Parameters must not be identical"
    assert param1 < len(titles) -1, "Please choose a number from the list"

    dimensions = len(titles) - 3

    other_dim_vals = {}
    for track in grid:
        for index in xrange(len(titles)-1):
            if index not in (param1, param2):
                if index not in other_dim_vals:
                    other_dim_vals[index] = []
                value = track[index]
                if value not in other_dim_vals[index]:
                    other_dim_vals[index].append(value)

    if len(other_dim_vals) > 1:
        print "\nCannot support more than 3 free parameters at this time"
        sys.exit(1)
    else:
        for key in other_dim_vals:
            length = len(other_dim_vals[key])
            if length > 10:
                print "\nThe expected number of graphs plotted is %d." % length
                print "Do you wish to select a subset of these graphs to plot or exit program?"
                cont = raw_input("Enter 1 to continue and 0 to exit: ")
                cont = int(cont)
                cont = bool(cont)

                if cont:
                    random.seed()
                    raw_num = raw_input("Enter the number of combinations to plot (maximum: %d): " % length)
                    num_to_plot = int(raw_num)
                    while len(other_dim_vals[key]) > num_to_plot:
                        index = random.randint(0, len(other_dim_vals[key]) - 1)
                        del other_dim_vals[key][index]

            else:
                cont = True

    if cont:
        generate_grids((param1, param2), other_dim_vals, titles, grid, results_age1, ndim)
        plt.show()
    else:
        sys.exit(0)
