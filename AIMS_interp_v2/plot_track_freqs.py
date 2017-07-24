import dill
import sys
import matplotlib.pyplot as plt
import itertools
import numpy as np

import model as mdl


def find_high_low_deviation(model_grid, num=1):
    """
    return lowest and highest variance tracks and their location in the grid
    """
    highest_variance = np.zeros(num, dtype=float)
    lowest_variance = np.array([1000 for i in xrange(num)], dtype=float)
    high_track_id = np.empty(num, dtype=int)
    low_track_id = np.empty(num, dtype=int)

    for index, track in enumerate(model_grid.tracks):
        mode_dict = {}  # dict containing the radial order as key and a list of the frequencies with that order
        total_variance = 0.

        for model in track.models:  # extract frequencies from each model in track and assign them to a dict of lists
            modes, freqs = find_frequencies(model)
            for mode, freq in itertools.izip(modes, freqs):
                if mode not in mode_dict:
                    mode_dict[mode] = []
                mode_dict[mode].append(freq)

        for key in mode_dict:  # get total variance of this track
            variance = np.var(mode_dict[key])
            total_variance += variance

        for i in xrange(num - 1, -1, -1):  # find and order the highest and lowest variance tracks
            if total_variance > highest_variance[i]:
                if i == num - 1:
                    highest_variance[i] = total_variance
                    high_track_id[i] = index
                else:
                    highest_variance[i], highest_variance[i+1] = total_variance, highest_variance[i]
                    high_track_id[i], high_track_id[i+1] = index, high_track_id[i]

            if total_variance < lowest_variance[i]:
                if i == num - 1:
                    lowest_variance[i] = total_variance
                    low_track_id[i] = index
                else:
                    lowest_variance[i], lowest_variance[i+1] = total_variance, lowest_variance[i]
                    low_track_id[i], low_track_id[i+1] = index, low_track_id[i]

    return high_track_id, highest_variance, low_track_id, lowest_variance


def plot_tracks(model_grid, track_list, offset=0):
    for x, track_id in enumerate(track_list):
        plot_modes(model_grid.tracks[track_id], x + offset)
        print "figure plotting COMPLETE"


def plot_modes(track, fig_num):
    max_mode = get_max_mode(track)
    fig = plt.figure(fig_num)
    print "plotting figure %d" % fig_num

    params = track.grid_params
    param_vals = track.params
    title_array = []
    for par_index in xrange(len(params)):
        title_array.append(params[par_index]+" = "+str(param_vals[par_index]))

    string_title = ", ".join(title_array)
    plt.title(string_title)
    for i in xrange(max_mode):
        age = []
        frequency = []
        for model in track.models:
            mode_n, mode_freq = find_frequencies(model)
            if len(mode_n) > i:
                age.append(model.glb[0])
                frequency.append(mode_freq[i])
                plt.plot(age, frequency)
                plt.xlabel("Age, Myr")
                plt.ylabel("Frequency")


def get_max_mode(track):
    max_mode = len(track.models[0].modes)
    for model in track.models:
        if len(model.modes) > max_mode:
            max_mode = len(model.modes)
    return max_mode


def find_frequencies(model):
    mode_n = []
    mode_freq = []

    for mode in model.modes:
        mode_n.append(mode['n'])
        mode_freq.append(mode['freq'])

    return mode_n, mode_freq


if __name__ == "__main__":
    filename = sys.argv[1]
    data = open(filename, "r")
    grid = dill.load(data)
    data.close()
    if len(sys.argv) > 2:
        plot_number = int(sys.argv[2])
    else:
        plot_number = 1

    choose_tracks = find_high_low_deviation(grid, plot_number)
    worst_id = choose_tracks[0]
    worst_variance = choose_tracks[1]
    best_id = choose_tracks[2]
    best_variance = choose_tracks[3]
    print "Models with highest variance: %s Variance: %s" % (worst_id, worst_variance)
    print "Models with lowest variance: %s Variance: %s" % (best_id, best_variance)
    plot_tracks(grid, np.append(worst_id, best_id))
    plt.show()
