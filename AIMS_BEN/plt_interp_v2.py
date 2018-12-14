import sys
import dill
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as interp
import colormaps
import constants
import pandas as pd
import sys

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
plt.rcParams["font.family"] = "serif"

'''
error_ndx:
[0] = maximum error on the radial modes
[1] = RMS error on the radial modes
[2] = RMS error on the radial modes near :math:`\\nu_{\mathrm{max}}`
[3] = maximum error on the non radial modes
[4] = RMS error on the non radial modes
[5] = RMS error on the non radial modes near :math:`\\nu_{\mathrm{max}}`
[6] = Age
[7] = Mass
[8] = Teff
[9] = Z0
[10] = X0
[11] = -
[12] = Xc
[13] = Period Spacing
[14] = Reference Frequency
[15] = Radius
[16] = Luminosity
'''

def interp_scatter(p,results,error_ndx,truncate=0,tpe="max"):#,a,tpe="max",title=None,truncate=0):
    df = pd.DataFrame()
    n = len(results)
    x = []
    y = []
    z1 = []
    z2 = []
    mu = []
    for i in xrange(n):
        start = truncate
        stop  = results[i].shape[0] - truncate
        if (tpe == "max"):
            value = np.nanmax(results[i][start:stop,ndim+error_ndx])
        elif (tpe == "avg"):
            num = den = 0.0
            for j in xrange(start,stop):
                if (not math.isnan(results[i][j,ndim+error_ndx])):
                    num += results[i][j,ndim+error_ndx]**2
                    den += 1.0
            if (den > 0.0):
                value = math.sqrt(num/den)
            else:
                value = np.nan
        else:
            print("ERROR: unrecognised type: ",tpe)
            sys.exit(1)
        if (value > 0.0):
	    # print value
            z1.append(value)	# math.log10
            z2.append(math.log10(value))#/constants.solar_radius))
            x.append(results[i][0,0])
            y.append(results[i][0,1])
            # mu.append(results[i][0,24])

    x = np.array(x,dtype = np.float64)
    y = np.array(y,dtype = np.float64)
    z1 = np.array(z1,dtype = np.float64)
    z2 = np.array(z2,dtype = np.float64)
    df['x'] = x
    df['y'] = y
    df['z1'] = z1
    df['z2'] = z2
    # df['mu'] = mu

    return df


if __name__ == "__main__":
    """
    Plot interpolation tests.
    """


    # assert (len(sys.argv) > 1), "Usage: plot_interpolation_test.py data_file"

    # filename = sys.argv[1]
    # filename_ms = 'interp_MS_v3.9'
    filename_ms = '/home/bmr135/bison/Sim2/AIMS_Gael/interp_MS_v3.9.1'
    filename_rgb = 'interp_RGB_v3.10'

    input_data = open(filename_ms,"r")
    [ndim, nglb, titles, grid, ndx1, ndx2, tessellation, results_age1_ms, \
        results_age2_ms, results_track_ms] = dill.load(input_data)
    results_age_ms = [results_age1_ms, results_age2_ms]
    input_data.close()

    input_data = open(filename_rgb,"r")
    [ndim, nglb, titles, grid, ndx1, ndx2, tessellation, results_age1_rgb, \
        results_age2_rgb, results_track_rgb] = dill.load(input_data)
    results_age_rgb = [results_age1_rgb, results_age2_rgb]
    input_data.close()

    a = interp_scatter(1,results_age1_ms,2,truncate=1,tpe="max")
    b = interp_scatter(1,results_age2_ms,2,truncate=1,tpe="max")
    c = interp_scatter(1,results_age1_rgb,2,truncate=1,tpe="max")
    d = interp_scatter(1,results_age2_rgb,2,truncate=1,tpe="max")

    a = a.sort_values(['z2'])
    b = b.sort_values(['z2'])
    c = c.sort_values(['z2'])
    d = d.sort_values(['z2'])

    # print(np.median(a['z2']/b['z2']))
    # print(np.median(a['mu']/b['mu']))
    # # print(np.median(c['z2']/d['z2']))
    # # print(np.sqrt(2))
    # sys.exit()

    max = np.max([np.max(a['z2']),np.max(b['z2']),np.max(c['z2']),np.max(d['z2'])])
    min = np.min([np.min(a['z2']),np.min(b['z2']),np.min(c['z2']),np.min(d['z2'])])
    dd = min-max

    ca = np.where(c['z2'] > -1.745)
    cb = np.where(c['z2'] > -1.854)
    da = np.where(d['z2'] > -1.745)
    db = np.where(d['z2'] > -1.854)
    # sys.exit()
    ''' Determine RGB threshold contours '''
    xi, yi = np.linspace(d['x'].min(),d['x'].max(),200), np.linspace(d['y'].min()+0.04,d['y'].max()+0.04,200)
    xi, yi = np.meshgrid(xi,yi)
    rbf = interp.Rbf(d['x'],d['y'],d['z2'],function='linear')
    zi = rbf(xi,yi)
    # plt.figure()
    # plt.close()
    '''  '''

    fig, ((ax,ax1,ax4),(ax2,ax3,ax5)) = plt.subplots(2,3,sharex='col',sharey='row',gridspec_kw={"width_ratios" : [5,5,0.2]},figsize=(9,4))
    col = ax.scatter(a['x'], a['y'], c=a['z2'], cmap=colormaps.parula,vmin=min,vmax=max,s=75)
    ax.set_title(r'Main Sequence',fontsize=12)
    ax.set_xlim(0.75,2.25)
    ax.set_ylim(-2.54,-1.46)
    ax2.scatter(b['x'], b['y'], c=b['z2'], cmap=colormaps.parula,vmin=min,vmax=max,s=75)
    ax2.set_ylim(-2.54,-1.46)

    ax1.scatter(c['x'], c['y'], c=c['z2'], cmap=colormaps.parula,vmin=min,vmax=max,s=75)
    ax1.set_xlim(0.75,2.25)
    # ax1.scatter(c['x'].iloc[cb], c['y'].iloc[cb], facecolors='none', edgecolors='k')
    # ax1.scatter(c['x'].iloc[ca], c['y'].iloc[ca], facecolors='none', edgecolors='m')
    ax1.set_title(r'Red Giant Branch',fontsize=12)
    cont = ax1.contourf(xi,yi,zi,100,alpha=0.)
    cont.collections[-18].set_color('m')
    cont.collections[-18].set_alpha(1.)
    cont.collections[-22].set_color('k')
    cont.collections[-22].set_alpha(1.)


    ax3.scatter(d['x'], d['y'], c=d['z2'], cmap=colormaps.parula,vmin=min,vmax=max,s=75)
    cont = ax3.contourf(xi,yi,zi,100,alpha=0.)
    cont.collections[-25].set_color('m')
    cont.collections[-25].set_alpha(1.)
    cont.collections[-29].set_color('k')
    cont.collections[-29].set_alpha(1.)

    ax4.axis('off')
    ax5.axis('off')

    plt.subplots_adjust(hspace=0.07, wspace=0.05,top=0.91,bottom=0.15)
    cbar_ax = fig.add_axes([.88, 0.13, 0.02, .79])
    cbar_ax.yaxis.tick_right()
    cbar_ax.axes.get_xaxis().set_ticks([])
    cbar = fig.colorbar(col,cax=cbar_ax,ticks=[-1.5,-1.75,-2.0,-2.25,-2.50,-2.75,-3.0,-3.25,-3.5,-3.75])
    cbar_ax.hlines(1+((max--1.745)/dd),0,1,colors='m',linewidth=2)
    cbar_ax.hlines(1+((max--1.854)/dd),0,1,colors='k',linewidth=2)
    fig.text(0.5, 0.03, r'Mass [M$_{\odot}$]', ha='center', va='center',fontsize=12)
    fig.text(0.025, 0.5, r'Metallicity [log$_{10}$($Z_{0}$)]', ha='center', va='center', rotation='vertical',fontsize=12)
    fig.text(0.0525, 0.725, r'Single Increment', ha='center', va='center', rotation='vertical',fontsize=8)
    fig.text(0.0525, 0.325, r'Double Increment', ha='center', va='center', rotation='vertical',fontsize=8)
    fig.text(0.975, 0.5, r'log$_{10}$(max. RMS error)', ha='center', va='center', rotation=270,fontsize=12)
    ax.xaxis.set_ticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2])
    ax1.xaxis.set_ticks([0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2])
    ax.yaxis.set_ticks([-2.5,-2.25,-2.0,-1.75,-1.5])
    ax2.yaxis.set_ticks([-2.5,-2.25,-2.0,-1.75,-1.5])

    plt.show()
    fig.savefig('model_interp.pdf', bbox_inches='tight')
