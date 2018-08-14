''' Plotting of interpolated HRD from AIMS '''

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.ticker import NullFormatter
import seaborn as sns
import sys

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams.update({'font.size': 20})
# matplotlib.rcParams.update({'axes.labelsize': 5})

def HRD(df1,df2,m,x,z,ext,save):
    ''' Plot HRD '''
    tmax = max(list(map(max,[df1['teff'], df2['teff']])))
    tmin = min(list(map(min,[df1['teff'], df2['teff']])))
    lmax = max(list(map(max,[df1['lumo'], df2['lumo']])))
    lmin = min(list(map(min,[df1['lumo'], df2['lumo']])))
    dT, dL = pd.DataFrame(), pd.DataFrame()
    dT['res_T'] = np.zeros(len(df1))
    dT['res_T'] = (df1['teff'] - df2['teff'])/df1['teff']
    dL['res_L'] = np.zeros(len(df1))
    dL['res_L'] = (df1['lumo'] - df2['lumo'])/df1['lumo']

    ''' Original Tracks '''
    names = ['num','age','logTe','logL','Xc','Yc','Lg/L','rhoc','Tc','rad','MnC','logg','ZXs', \
            'Dif','alpha','over','Rconv','Z','X','Mass','DPl1','Mover','mHe']
    t1 = pd.read_csv('/home/ben/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.98.X0.679.Z0.0332-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')
    t2 = pd.read_csv('/home/ben/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.98.X0.713.Z0.0175-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')
    t3 = pd.read_csv('/home/ben/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.98.X0.732.Z0.0090-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')
    # t2 = pd.read_csv('/home/bmr135/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.97.X0.679.Z0.0332-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')
    # t3 = pd.read_csv('/home/bmr135/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.99.X0.679.Z0.0332-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.35, 0.55
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_T = [left, bottom-0.23, width, 0.2]
    rect_L = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(figsize=(12, 6))

    axHRD = plt.axes(rect_scatter)
    axT = plt.axes(rect_T)
    axL = plt.axes(rect_L)#,sharey=axHRD)

    # for i in range(len(df1)):
    #     axHRD.add_line(lines.Line2D([df1['teff'][i],df2['teff'][i]], [df1['lumo'][i],df2['lumo'][i]+0.5], \
    #                         linewidth=1, color='k', linestyle='--', axes=axHRD, alpha=0.3))

    # the scatter plot:
    axHRD.scatter(df1['teff'],df1['lumo'],label=r'Original',marker='+')
    axHRD.scatter(df2['teff'],df2['lumo']+0.5,label=r'Interpolated',marker='3')
    axHRD.plot(10**t1['logTe'],10**t1['logL'],color='m',alpha=0.3,label='Track')
    # axHRD.plot(10**t1['logTe'],10**t1['logL']+0.5,color='m',alpha=0.3,label='Track')
    axHRD.plot(10**t2['logTe'],10**t2['logL'],color='m',linestyle='--',alpha=0.3,label='Track')
    # axHRD.plot(10**t2['logTe'],10**t2['logL']+0.5,color='m',linestyle='--',alpha=0.3,label='Track')
    axHRD.plot(10**t3['logTe'],10**t3['logL'],color='m',linestyle='-.',alpha=0.3,label='Track')
    # axHRD.plot(10**t3['logTe'],10**t3['logL']+0.5,color='m',linestyle='-.',alpha=0.3,label='Track')
    axHRD.set_xlim(tmax+50,tmin-50)
    axHRD.set_ylim(lmin-1,lmax+1)

    axT.scatter(df1['teff'],dL['res_L'],marker='+')
    axL.scatter(dT['res_T'],df1['lumo'],marker='+')

    # no labels
    nullfmt = NullFormatter()
    axHRD.xaxis.set_major_formatter(nullfmt)
    # axL.yaxis.set_major_formatter(nullfmt)

    axT.set_xlim(axHRD.get_xlim())
    axL.set_ylim(axHRD.get_ylim())
    # axHRD.set_yscale('log')
    # axL.set_yscale('log')

    axL.get_yaxis().set_ticklabels([])


    axT.set_xlabel(r'T$_{\rm{eff}}$ [K]',fontsize=20)
    axT.set_ylabel(r'$\Delta$L')
    axL.set_xlabel(r'$\Delta$T$_{\rm{eff}}$')
    axHRD.set_ylabel(r'L [L$_{\odot}$]',fontsize=20)
    axHRD.set_title(r'M = %s M$_{\odot}$, Z = %s, X = %s' %(df1['mass'][0],df1['z'][0],df1['x'][0]),fontsize=20)
    # axHRD.legend(prop={'size':15},framealpha=0.5)
    # plt.savefig(ext+save+str(m)+'_X'+str(x)+'_Z'+str(z)+'.png')
    plt.show()

def lumo_age(df1,df2,m,x,z,ext,save):
    ''' Plot luminosity as a function of age '''
    ''' Original Tracks '''
    names = ['num','age','logTe','logL','Xc','Yc','Lg/L','rhoc','Tc','rad','MnC','logg','ZXs', \
            'Dif','alpha','over','Rconv','Z','X','Mass','DPl1','Mover','mHe']
    t1 = pd.read_csv('/home/ben/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.98.X0.679.Z0.0332-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')
    # t2 = pd.read_csv('/home/bmr135/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.98.X0.713.Z0.0175-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')
    # t3 = pd.read_csv('/home/bmr135/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.98.X0.732.Z0.0090-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')
    t2 = pd.read_csv('/home/ben/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.97.X0.679.Z0.0332-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')
    t3 = pd.read_csv('/home/ben/Dropbox/AIMS_Work_Files/AIMS_Interp_Tests/M0.99.X0.679.Z0.0332-sumN.txt',names=names,skiprows=3,delimiter=r'\s+')


    fig, ax = plt.subplots(1)
    for i in range(len(df1)):
        ax.add_line(lines.Line2D([df1['age'][i],df2['age'][i]], [df1['lumo'][i],df2['lumo'][i]+0.5], \
                            linewidth=1, color='k', linestyle='--', axes=ax, alpha=0.3))
    ax.scatter(df1['age'],df1['lumo'],label=r'Original',marker='+')
    ax.scatter(df2['age'],df2['lumo']+.5,label=r'Interpolated',marker='3')
    ax.plot(t1['age']*10**-6,10**t1['logL'],color='m',alpha=0.3,label='Track')
    # ax.plot(t1['age']*10**-6,10**t1['logL']+0.5,color='m',alpha=0.3,label='Track')
    ax.plot(t2['age']*10**-6,10**t2['logL'],color='m',linestyle='--',alpha=0.3,label='Track')
    # ax.plot(t2['age']*10**-6,10**t2['logL']+0.5,color='m',linestyle='--',alpha=0.3,label='Track')
    ax.plot(t3['age']*10**-6,10**t3['logL'],color='m',linestyle='-.',alpha=0.3,label='Track')
    # ax.plot(t3['age']*10**-6,10**t3['logL']+0.5,color='m',linestyle='-.',alpha=0.3,label='Track')
    ax.set_xlabel(r'Age [Myr]')
    ax.set_ylabel(r'L [L$_{\odot}$]')
    plt.legend(prop={'size':15})
    plt.tight_layout()
    # plt.savefig(ext+save+str(m)+'_X'+str(x)+'_Z'+str(z)+'.png')
    plt.show()

if __name__ == '__main__':

    if sys.argv[1] == '1':
        ''' MS linear interp models '''
        df = pd.read_csv('/home/bmr135/Dropbox/combinations_Delaunay_MS.txt',names=['n','mass','rad','lumo','z','x','age','teff'],delimiter=r'\s+')
        s_hrd = 'MS_HRD/M'
        s_la = 'MS_Lumo_Age/M'
        ext = '/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/results/Interp_HRDs/'
        # m = np.linspace(0.9,1.5,61)
        # x = [0.679,0.713,0.732,0.741,0.746]
        # z = [0.0332,0.0175,0.0090,0.0046,0.0023]

        m = [0.98]
        x = [0.679]
        z = [0.0332]

        df['mass'] = pd.to_numeric(df['mass'],errors='coerce')
        df = df.dropna()
        df = df.reset_index(drop=True)
        for i in range(len(m)):
            df1 = pd.DataFrame()
            for j in range(len(x)):
                df1 = df[(df['mass'] == m[i]) & (df['x'] == x[j]) & (df['z'] == z[j])]
                if len(df1) > 0:
                    df1 = df1.reset_index(drop=True)
                    orig = df1.iloc[::2]
                    orig = orig.reset_index(drop=True)
                    interp = df1.iloc[1::2]
                    interp = interp.reset_index(drop=True)
                    if len(orig) > len(interp):
                        diff = len(orig) - len(interp)
                        orig = orig[:-diff]
                    # HRD(orig,interp,m[i],x[j],z[j],ext,s_hrd)
                    lumo_age(orig, interp,m[i],x[j],z[j],ext,s_la)
        #interp.to_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/interp_MS_mods',index=False)


    if sys.argv[1] == '2':
        ''' RGB linear interp models '''
        df = pd.read_csv('/home/bmr135/Dropbox/combinations_Delaunay_RGB.txt',names=['n','mass','rad','lumo','z','x','age','teff'],delimiter=r'\s+')
        s_hrd = 'RGB_HRD/M'
        s_la = 'RGB_Lumo_Age/M'
        ext = '/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/results/Interp_HRDs/'
        # m = np.linspace(0.9,1.5,61)
        # x = [0.679,0.713,0.732,0.741,0.746]
        # z = [0.0332,0.0175,0.0090,0.0046,0.0023]

        m = [0.98]
        x = [0.679]
        z = [0.0332]

        df['mass'] = pd.to_numeric(df['mass'],errors='coerce')
        df = df.dropna()
        df = df.reset_index(drop=True)
        for i in range(len(m)):
            df1 = pd.DataFrame()
            for j in range(len(x)):
                df1 = df[(df['mass'] == m[i]) & (df['x'] == x[j]) & (df['z'] == z[j])]
                if len(df1) > 0:
                    df1 = df1.reset_index(drop=True)
                    orig = df1.iloc[::2]
                    orig = orig.reset_index(drop=True)
                    interp = df1.iloc[1::2]
                    interp = interp.reset_index(drop=True)
                    if len(orig) > len(interp):
                        diff = len(orig) - len(interp)
                        orig = orig[:-diff]
                    HRD(orig,interp,m[i],x[j],z[j],ext,s_hrd)
                    lumo_age(orig, interp,m[i],x[j],z[j],ext,s_la)
        #interp.to_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/results/interp_RGB_mods',index=False)


    if sys.argv[1] == '3':
        ''' Delaunay tessellation models formatting to save files (formatting
        is time consuming therefore save out and read in later on). '''
        ext = '/media/bmr135/SAMSUNG/AIMS-interp-testing2/Interp_HRDs_RGB/'
        # ext = '/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/results/Interp_HRDs/'
        # s_hrd = 'MS_HRD/M'
        s_hrd = 'RGB_mHe_HRD/M'

        # dy = pd.read_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/Delaunay_MS_Models.txt',names=['model'],delimiter=r'\s+')
        # dy = pd.read_csv('/home/bmr135/git_AIMS/AIMS/AIMS_BEN/Delaunay_MS_Models_mHe.txt',names=['model'],delimiter=r'\s+')
        dy = pd.read_csv('/media/bmr135/SAMSUNG/AIMS-interp-testing2/Delaunay_RGB_Models_mHe.txt',names=['model'],delimiter=r'\s+')
        # df = pd.read_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/Delaunay_MS_Mod_vals.txt',\
        #                     names=['n','mass','rad','lumo','z','x','age','teff'],delimiter=r'\s+')
        df = pd.read_csv('/media/bmr135/SAMSUNG/AIMS-interp-testing2/Delaunay_RGB_Mod_vals_mHe.txt',\
                           names=['n','mass','rad','lumo','z','x','age','teff'],delimiter=r'\s+')
        # df = pd.read_csv('/media/bmr135/SAMSUNG/AIMS-interp-testing2/Delaunay_MS_Mod_vals_mHe.txt',\
        #                     names=['n','mass','rad','lumo','z','x','age','teff'],delimiter=r'\s+')

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

        ''' Save out split values if desired '''
        # orig.to_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/Del_RGB_orig_mods', index=False)
        # mod1.to_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/Del_RGB_A_mods', index=False)
        # mod2.to_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/Del_RGB_B_mods', index=False)
        # orig.to_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/Del_MS_orig_mods', index=False)
        # orig.to_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/Del_MS_A_mods', index=False)
        # orig.to_csv('/home/bmr135/git_AIMS/George_Smetana/AIMS-Interp/Del_MS_B_mods', index=False)


        m = np.linspace(0.75,2.25,76)
        x = [0.691,0.716,0.731,0.740,0.745]
        z = [0.0300,0.0175,0.0100,0.0057,0.0032]
        # m = [1.17]
        # x = [0.746]
        # z = [0.0023]
        # print(df)

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

        for i in range(len(m)):
            df1, df2, df3, df4 = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
            for j in range(len(x)):
                df1 = orig[(orig['mass'] == m[i]) & (orig['x'] == x[j]) & (orig['z'] == z[j])]
                a = orig[(orig['mass'] == m[i]) & (orig['x'] == x[j]) & (orig['z'] == z[j])].index.get_values()
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

                    # fig, ax = plt.subplots(1)
                    left, width = 0.1, 0.65
                    bottom, height = 0.35, 0.55
                    bottom_h = left_h = left + width + 0.02
                    rect_scatter = [left, bottom, width, height]
                    rect_T = [left, bottom-0.23, width, 0.2]
                    rect_L = [left_h, bottom, 0.2, height]

                    plt.figure(figsize=(12, 6))
                    ax = plt.axes(rect_scatter)
                    axT = plt.axes(rect_T)
                    axL = plt.axes(rect_L)#,sharey=axHRD)

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

                    ax.plot(df1['teff'],df1['lumo'],label=r'Orig. Track',alpha=0.5)
                    ax.scatter(df2['teff'],df2['lumo']-.5,label=r'M:%s, Z:%s, X:%s'%(df2['mass'][0],df2['z'][0],df2['x'][0]),marker='.')
                    ax.plot(df2['teff'],df2['lumo']-.5,alpha=0.3,label='_nolegend_')
                    ax.scatter(df3['teff'],df3['lumo']+.5,label=r'M:%s, Z:%s, X:%s'%(df3['mass'][0],df3['z'][0],df3['x'][0]),marker='.')
                    ax.plot(df3['teff'],df3['lumo']+.5,alpha=0.3,label='_nolegend_')
                    ax.scatter(df4['teff'],df4['lumo'],label=r'Interp. Models',color='r',marker='3')
                    # ax.set_xlabel(r'T$_{\rm{eff}}$ [K]')
                    ax.set_ylabel(r'L [L$_{\odot}$]')
                    ax.set_xlim(tmax+50,tmin-50)
                    ax.set_ylim(lmin-1,lmax+1)
                    ax.set_title(r'M = %s M$_{\odot}$, Z = %s, X = %s' %(df1['mass'][0],df1['z'][0],df1['x'][0]),fontsize=20)
                    ax.legend(prop={'size':10})

                    axT.scatter(df1['teff'],dL['res_L'],marker='+')
                    axL.scatter(dT['res_T'],df1['lumo'],marker='+')

                    # no labels
                    nullfmt = NullFormatter()
                    ax.xaxis.set_major_formatter(nullfmt)
                    # axL.yaxis.set_major_formatter(nullfmt)
                    axT.set_xlim(ax.get_xlim())
                    axL.set_ylim(ax.get_ylim())
                    # axHRD.set_yscale('log')
                    # axL.set_yscale('log')
                    axL.get_yaxis().set_ticklabels([])
                    axT.set_xlabel(r'T$_{\rm{eff}}$ [K]',fontsize=20)
                    axT.set_ylabel(r'$\Delta$L')
                    axL.set_xlabel(r'$\Delta$T$_{\rm{eff}}$')

                    # plt.tight_layout()
                    plt.savefig(ext+'M'+str(m[i])+'_X'+str(x[j])+'_Z'+str(z[j])+'.png')
                    # plt.show()
                    plt.close()
		            # sys.exit()



        '''
        Requirements:
        - Plot orginal track plus tracks interpolated from
        - Plot HRD scatter of interpolated points
        - Plot scatter of models for each interpolation
        - Join model points to respective interpolated models
        '''
