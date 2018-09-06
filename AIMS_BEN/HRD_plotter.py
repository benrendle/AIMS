''' HRD Plotter from AIMS input list file '''

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import colormaps
import constants

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

df = pd.read_csv('/home/bmr135/bison/Sim2/AIMS_Gael/CLES_MS_v3.2',names=['model','mass','radius','lumo','met','X','age','teff','mHe','Xc','pi'],skiprows=1,delimiter=r'\s+')
df1 = pd.read_csv('/home/bmr135/bison/Sim2/AIMS_Gael/CLES_RGB_v3.2',names=['model','mass','radius','lumo','met','X','age','teff','mHe','Xc','pi'],skiprows=1,delimiter=r'\s+')
# print(df['model'].iloc[0][17:19])
df = df[df['X'] == 0.731]
df = df.reset_index(drop=True)
df1 = df1[df1['X'] == 0.731]
df1 = df1.reset_index(drop=True)
# print(df1)

cm = colormaps.parula
plt.gca().set_color_cycle([cm(1.*i/75) for i in range(75)])
plt.xlabel(r'T$_{\mathrm{eff}}$ [K]',fontsize=20)
plt.ylabel(r'Luminosity [$\log_{10}\left(L/L_{\odot}\right)$]',fontsize=20)

# print(df['model'].iloc[-1][:5])
f=0
d=0
x = []
for i in range(len(df)-1):
    x = np.append(x,i)
    if (df['mass'].iloc[i] != df['mass'].iloc[i+1]):
        if f==0:
            plt.plot(df['teff'][x],np.log10(df['lumo'][x]/constants.solar_luminosity),label = df['model'].iloc[i][:5] if (f/15.0).is_integer() else '_nolegend_',lw=2,zorder=-1)
            x = []
        if f>0:
            plt.plot(df['teff'][x],np.log10(df['lumo'][x]/constants.solar_luminosity),label = df['model'].iloc[i+1][:5] if ((f+1)/15.0).is_integer() else '_nolegend_',lw=2,zorder=-1)
            x = []
        # print(f/15.0)
        f+=1
# print(f)
x = []
for i in range(len(df1)-1):
    x = np.append(x,i)
    if (df1['mass'].iloc[i] != df1['mass'].iloc[i+1]):
        plt.plot(df1['teff'][x],np.log10(df1['lumo'][x]/constants.solar_luminosity),label='_nolegend_')
        x = []
        d+=1
# print(d)




ax = plt.gca()
ax.legend(loc=3,prop={'size':15})#,frameon=False)
ax.set_xlim(ax.get_xlim()[::-1])
ax.tick_params(labelsize=20)
plt.tight_layout()
plt.show()
