# Track interpolation testing plotting program

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
import matplotlib.lines as lines
import matplotlib.lines as mlines
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib.patches import Ellipse
from matplotlib.ticker import NullFormatter
import seaborn as sns
import sys
import dill
import math
import scipy.interpolate as interpol
import colormaps
import utilities

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams.update({'font.size': 12})
plt.rcParams["font.family"] = "serif"

''' Track Parameters '''
m = 1.47
x = 0.740
z = 0.0057
''' Echelle comes from data point 4 in original file '''
'''
MS - 1.47, 0.740, 0.0057, 272
RGB - 1.19, 0.731, 0.0100, 612
'''
''' Delaunay tessellation models formatting to save files (formatting
is time consuming therefore save out and read in later on). '''

dy = pd.read_csv('/media/bmr135/SAMSUNG/AIMS-interp-testing2/Delaunay_MS_Models_mHe.txt',names=['model'],delimiter=r'\s+')
df = pd.read_csv('/media/bmr135/SAMSUNG/AIMS-interp-testing2/Delaunay_MS_Mod_vals_mHe.txt',\
                    names=['n','mass','rad','lumo','z','x','age','teff'],delimiter=r'\s+')


dy['mass'], dy['x'], dy['z'], dy['n'] = 0, 0, 0, 0
dy['mass']=dy['model'].map(lambda x: x.split('M')[1][:4])
dy['x']=dy['model'].map(lambda x: x.split('X')[1][:5])
dy['z']=dy['model'].map(lambda x: x.split('Z')[1][:6])
dy['n']=dy['model'].map(lambda x: x.split('-')[1][1:4])

dy0 = dy['n'].iloc[::3]
dy0 = dy0.reset_index(drop=True)
dy1 = dy['n'].iloc[1::3]
dy1 = dy1.reset_index(drop=True)
dy2 = dy['n'].iloc[2::3]
dy2 = dy2.reset_index(drop=True)

df['mass'] = pd.to_numeric(df['mass'],errors='coerce')
df = df.dropna()
df = df.reset_index(drop=True)
orig = df.iloc[::4]
orig = orig.reset_index(drop=True)
orig['n'] = dy0
mod1 = df.iloc[1::4]
mod1 = mod1.reset_index(drop=True)
mod1['n'] = dy.iloc[1::2]
mod1['n'] = dy1
mod2 = df.iloc[2::4]
mod2 = mod2.reset_index(drop=True)
mod2['n'] = dy.iloc[2::3]
mod2['n'] = dy2
interp = df.iloc[3::4]
interp = interp.reset_index(drop=True)

df1, df2, df3, df4 = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
df1 = orig[(orig['mass'] == m) & (orig['x'] == x) & (orig['z'] == z)]
a = orig[(orig['mass'] == m) & (orig['x'] == x) & (orig['z'] == z)].index.get_values()


''' Plotting Function with GridSpec '''
fig = plt.figure(figsize=(12,8))
gs = gridspec.GridSpec(3,1,wspace=0.1,top=0.95,bottom=0.075,hspace=0.275)

# gs.update(wspace=0.02, hspace=0.03)


''' HRD with Residuals '''

if len(df1) > 0:
    df2 = mod1.iloc[a]
    df3 = mod2.iloc[a]
    df4 = interp.iloc[a]
    df1 = df1.reset_index(drop=True) # Original Track
    df2 = df2.reset_index(drop=True) # Interpolation Track 1
    df3 = df3.reset_index(drop=True) # Interpolation Track 2
    df4 = df4.reset_index(drop=True) # Interpolated Track
    if len(df1) > len(df2):
        diff = len(df1) - len(df2)
        df1 = df1[:-diff]

    dT, dL = pd.DataFrame(), pd.DataFrame()
    dT['res_T'] = np.zeros(len(df1))
    dT['res_T'] = (df1['teff'] - df4['teff'])/df1['teff']
    dL['res_L'] = np.zeros(len(df1))
    dL['res_L'] = (df1['lumo'] - df4['lumo'])/df1['lumo']

    gs00 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1], height_ratios=[1.0,0.2], width_ratios=[0.9,0.1], hspace=0.1, wspace=0.05)

    ax = plt.subplot(gs00[0,0])
    axT = plt.subplot(gs00[1,0])
    axL = plt.subplot(gs00[0,1])

    tmax = max(list(map(max,[df1['teff'], df2['teff'], df3['teff'], df4['teff']])))
    tmin = min(list(map(min,[df1['teff'], df2['teff'], df3['teff'], df4['teff']])))
    lmax = max(list(map(max,[df1['lumo'], df2['lumo'], df3['lumo'], df4['lumo']])))
    lmin = min(list(map(min,[df1['lumo'], df2['lumo'], df3['lumo'], df4['lumo']])))

    for k in range(len(df4)):
        ax.add_line(lines.Line2D([df4['teff'][k],df2['teff'][k]], [df4['lumo'][k],df2['lumo'][k]-.5], \
                            linewidth=1, color='k', linestyle='--', axes=ax, alpha=0.3))
    for k in range(len(df4)):
        ax.add_line(lines.Line2D([df4['teff'][k],df3['teff'][k]], [df4['lumo'][k],df3['lumo'][k]+.5], \
                            linewidth=1, color='k', linestyle='--', axes=ax, alpha=0.3))

    ax.plot(df1['teff'],df1['lumo'],label=r'Orig. Track',alpha=0.5,linewidth=2)
    ax.scatter(df2['teff'],df2['lumo']-.5,label=r'M:%s, Z:%s, X:%s'%(df2['mass'][0],df2['z'][0],df2['x'][0]),marker='.',s=75)
    # ax.plot(df2['teff'],df2['lumo']-.5,alpha=0.3,label='_nolegend_',linewidth=2)
    ax.scatter(df3['teff'],df3['lumo']+.5,label=r'M:%s, Z:%s, X:%s'%(df3['mass'][0],df3['z'][0],df3['x'][0]),marker='.',s=75)
    # ax.plot(df3['teff'],df3['lumo']+.5,alpha=0.3,label='_nolegend_',linewidth=2)
    ax.scatter(df4['teff'],df4['lumo'],label=r'Interp. Models',color='r',marker='3',linewidths=2,s=75)
    ax.scatter(df1['teff'].iloc[4], df1['lumo'].iloc[4], facecolors='none', edgecolors='r', linewidths=3, s=150)
    # ax.set_xlabel(r'T$_{\rm{eff}}$ [K]')
    ax.set_ylabel(r'L [L$_{\odot}$]')
    ax.set_xlim(tmax+50,tmin-50)
    ax.set_ylim(lmin-1,lmax+1)
    # ax.set_title(r'M = %s M$_{\odot}$, Z = %s, X = %s' %(df1['mass'][0],df1['z'][0],df1['x'][0]),fontsize=20)
    ax.legend(prop={'size':9},loc=4,ncol=2)

    axT.scatter(df1['teff'],dL['res_L'],marker='+',s=75)
    axL.scatter(dT['res_T']*100,df1['lumo'],marker='+',s=75)

    # no labels
    nullfmt = NullFormatter()
    ax.xaxis.set_major_formatter(nullfmt)
    # axL.yaxis.set_major_formatter(nullfmt)
    axT.set_xlim(ax.get_xlim())
    axL.set_ylim(ax.get_ylim())
    axL.get_yaxis().set_ticklabels([])
    axT.set_xlabel(r'T$_{\rm{eff}}$ [K]')#,fontsize=20)
    axT.set_ylabel(r'$\Delta$L/L')
    axL.set_xlabel(r'$\Delta$T/T')
axL.xaxis.set_major_locator(plt.MaxNLocator(3))
axL.set_xlim(-0.5,0.5)
axT.yaxis.set_major_locator(plt.MaxNLocator(3))
axT.set_ylim(-0.01,0.01)
ax.text(0.025, 0.9, '(B)', horizontalalignment='center',\
      verticalalignment='center', transform=ax.transAxes)


''' Interpolation Test Result for Grid '''

gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios=[2.0,0.1], wspace=0.05)
axSc = plt.subplot(gs01[0,0])
axCol = plt.subplot(gs01[0,1])


filename = '/home/bmr135/AIMS/AIMS_BEN/interp_MS_v3.9'
input_data = open(filename,"r")
[ndim, nglb, titles, grid, ndx1, ndx2, tessellation, results_age1, \
    results_age2, results_track] = dill.load(input_data)
results_age = [results_age1, results_age2]
input_data.close()

def plot_partition_tessellation(grid, ndx1, ndx2, tessellation):
    """
    Make an interactive tessellation plot based on the supplied partition
    on the grid. Clicking on the blue dots will produce a 2D slice showing
    track interpolation errors for the associated evolutionary track.

    :param grid: array containing basic grid parameters (excluding age)
    :param ndx1: list with the indices of the first part of the partition.
    :param ndx2: list with the indices of the second part of the partition.
    :param tessellation: grid tessellation associated with ``ndx2``

    :type grid: np.array
    :type ndx1: list of int
    :type ndx2: list of int

    .. warning::
      This only works for two-dimensional tessellations.
    """

    # remember ndim includes the age dimension which is not included
    # in the tessellation:
    if (ndim != 3):
        print("Only able to plot the tessellation in two dimensions.")
        return

    # find bounds:
    xmin = np.nanmin(grid[:,0])
    xmax = np.nanmax(grid[:,0])
    ymin = np.nanmin(grid[:,1])
    ymax = np.nanmax(grid[:,1])
    dx = xmax-xmin
    dy = ymax/ymin
    xmin -= dx*0.03
    xmax += dx*0.03
    ymin /= dy**0.05
    ymax *= dy**0.05

    # axSc.triplot(grid[ndx1,0],grid[ndx1,1],tessellation.simplices.copy(),alpha=0.2,color='grey')
    axSc.triplot(grid[ndx2,0],grid[ndx2,1],tessellation.simplices.copy(),alpha=0.2,color='grey')
    axSc.plot(grid[ndx2,0],grid[ndx2,1],'o',alpha=0.2,color='grey')
    # axSc.plot(0.75,-2.244,'rd', alpha=0.2)
    axSc.plot(0.75,-2.494,'rd', alpha=0.2)
    # axSc.plot(0.77,-2.494,'rd', alpha=0.2)
    axSc.plot(2.25,-2.494,'rd', alpha=0.2)
    axSc.plot(2.25,-2.244,'rd', alpha=0.2)
    axSc.plot(2.25,-1.523,'rd', alpha=0.2)
    axSc.plot(2.23,-1.523,'rd', alpha=0.2)
    axSc.plot(2.21,-1.523,'rd', alpha=0.2)
    axSc.plot(2.19,-1.523,'rd', alpha=0.2)
    axSc.plot(2.17,-1.523,'rd', alpha=0.2)
    axSc.plot(2.15,-1.523,'rd', alpha=0.2)

# print results_track

p = 1
results = results_track
error_ndx = 0
tpe = 'max'
title = 'Max. radial error between tracks'
truncate = 1
n = len(results)
x1 = []
y1 = []
z1 = []
z2 = []
for i in range(len(results)):
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
        z1.append(value)	# math.log10
        z2.append(math.log10(value))
        x1.append(results[i][0,0])
        y1.append(results[i][0,1])
        # print results[i][0,0],results[i][0,1], value
x1 = np.array(x1,dtype = np.float64)
y1 = np.array(y1,dtype = np.float64)
z1 = np.array(z1,dtype = np.float64)
z2 = np.array(z2,dtype = np.float64)
min = np.min(z2)
max = np.max(z2)
# xi, yi = np.linspace(x1.min(),x1.max(),200), np.linspace(y1.min(),y1.max(),200)
# xi, yi = np.meshgrid(xi,yi)
# rbf = interpol.Rbf(x1,y1,z2,function='linear')
# zi = rbf(xi,yi)
# track = axSc.contourf(xi,yi,zi,100,cmap=colormaps.parula)

# elle = Ellipse(xy=[m,np.log10(z)], width=0.015, height=0.03, fill=None, linewidth=2, color='r')
# circle1 = plt.Circle((m, np.log10(z)), .0075, color='k', fill=False,linewidth=3)
# axSc.add_artist(elle)
plot_partition_tessellation(grid, ndx1, ndx2, tessellation)
col = axSc.scatter(x1, y1, c=z2, cmap=colormaps.parula,vmin=min,vmax=max,alpha=1.1)
axSc.scatter(m, np.log10(z), facecolors='none', edgecolors='r', linewidths=3, s=160)
axSc.set_xlabel(titles[0])#,fontsize=20)
axSc.set_ylabel(titles[1])#,fontsize=20)
cbar = plt.colorbar(col,cax=axCol)

# axCol = plt.cb()
# axSc.set_xticks(fontsize=15)
# axSc.set_yticks(fontsize=15)

# cax = cb.ax
# cax.text(3.5,0.7,r"$\log_{10}$(Max. error)",rotation=270,fontsize=20)
# cax.tick_params(labelsize=15)
# cbar = Colorbar(ax = axCol, mappable = track, orientation = 'vertical')#, ticklocation = 'right')
cbar.set_label(r"$\log_{10}$(max. RMS error)", labelpad=20, rotation=270)
# if (title is not None): axSc.set_title(title,fontsize=15)
axSc.text(0.025, 0.9, '(A)', horizontalalignment='center',\
      verticalalignment='center', transform=axSc.transAxes)


''' Echelle Diagram '''

# theo = pd.read_csv('') # Theoretical Frequencies
filename = '/media/bmr135/SAMSUNG/AIMS-interp-testing2/M147_MS_272' # Original Frequencies
filename1 = '/media/bmr135/SAMSUNG/AIMS-interp-testing2/Interp_Freqs_MS/M1.47.X0.740.Z0.0057-atm-0272' # Interpolated Frequencies
freqfile = open(filename)
freqfile.readline() # skip head
mode_temp = []
l, nn, fre = [], [], []
for line in freqfile:
    line = line.strip()
    columns = line.split()
    n = int(columns[1])
    freq = utilities.to_float(columns[2])
    mode_temp.append((n,int(columns[0]),freq,utilities.to_float(columns[4]))) # CLES
    l.append(int(columns[0]))
    nn.append(n)
    fre.append(freq)

freqfile1 = open(filename1)
freqfile1.readline() # skip head
dnu = float(freqfile1.readline())
mode_temp1 = []
for line in freqfile1:
    line = line.strip()
    columns = line.split()
    n = int(columns[0])
    freq1 = utilities.to_float(columns[2])
    mode_temp1.append((n,int(columns[1]),freq1)) # ,utilities.to_float(columns[4]))) # CLES



axE = plt.subplot(gs[2])

l0, l1, l2 = [], [], []
for i in mode_temp:
    if i[1] == 0: l0.append(i)
    if i[1] == 1: l1.append(i)
    if i[1] == 2: l2.append(i)
l01, l11, l21 = [], [], []
for i in mode_temp1:
    if i[1] == 0: l01.append(i)
    if i[1] == 1: l11.append(i)
    if i[1] == 2: l21.append(i)
n = 8
m = 20
shift = 5
# axE.scatter((mode_temp[n][2]+shift)%dnu,mode_temp[n][2],label='Original: $l=0$',color='blue')
# axE.scatter((mode_temp1[n][2]+shift)%dnu,mode_temp1[n][2],color='orange',label='Interp.: $l=0$',facecolors='none',s=140,linewidths=1.5)

axE.scatter((l0[n][2]+shift)%dnu,l0[n][2],label='Original: $l=0$',color='blue')
axE.scatter((l01[n][2]+shift)%dnu,l01[n][2],color='orange',label='Interp.: $l=0$',facecolors='none',s=140,linewidths=2.)
axE.scatter((l1[n][2]+shift)%dnu,l1[n][2],label='Original: $l=1$',marker='D',color='blue')
axE.scatter((l11[n][2]+shift)%dnu,l11[n][2],color='orange',label='Interp.: $l=1$',facecolors='none',s=140,linewidths=2.,marker='D')
axE.scatter((l2[n][2]+shift)%dnu,l2[n][2],label='Original: $l=2$',marker='^',color='blue')
axE.scatter((l21[n][2]+shift)%dnu,l21[n][2],color='orange',label='Interp.: $l=2$',facecolors='none',s=140,linewidths=2.,marker='^')

# k=0
# for i in mode_temp:
#     if (k>n) & (k<m): axE.scatter((i[2]+shift)%dnu,i[2],color='blue')
#     k+=1
k=0
for i in l0:
    if (k>n)  & (k<m): axE.scatter((i[2]+shift)%dnu,i[2],color='blue')
    k+=1
k=0
for i in l1:
    if (k>n-1)  & (k<m): axE.scatter((i[2]+shift)%dnu,i[2],marker='D',color='blue')
    k+=1
k=0
for i in l2:
    if (k>n-1)  & (k<m): axE.scatter((i[2]+shift)%dnu,i[2],marker='^',color='blue')
    k+=1
# k=0
# for i in mode_temp1:
#     if (k>n) & (k<m): axE.scatter((i[2]+shift)%dnu,i[2],color='orange',facecolors='none',s=140,linewidths=2.)
#     k+=1
k=0
for i in l01:
    if (k>n)  & (k<m): axE.scatter((i[2]+shift)%dnu,i[2],color='orange',facecolors='none',s=140,linewidths=1.5)
    k+=1
k=0
for i in l11:
    if (k>n-1)  & (k<m): axE.scatter((i[2]+shift)%dnu,i[2],color='orange',facecolors='none',s=140,linewidths=1.5,marker='D')
    k+=1
k=0
for i in l21:
    if (k>n-1)  & (k<m): axE.scatter((i[2]+shift)%dnu,i[2],color='orange',facecolors='none',s=140,linewidths=1.5,marker='^')
    k+=1

# axE.legend()
axE.text(0.025, 0.9, '(C)', horizontalalignment='center',\
      verticalalignment='center', transform=axE.transAxes)


axE.set_xlabel(r"Reduced frequency, $\nu$ mod "+str(dnu)+r" $\mu$Hz")
axE.set_ylabel(r"Frequency, $\nu$ (in $\mu$Hz)")
axE.set_xlim(0.0,dnu)
axE.legend(prop={'size':9})
print(dnu)

fig.savefig('track_interp_MS_alt.jpeg', bbox_inches='tight', dpi=52)
plt.show()
