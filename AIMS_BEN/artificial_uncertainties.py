''' Frequency uncertainties based on 16Cyg A distribution '''

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import argrelextrema

n = [13,14,15,16,17,18,19,20,21,22,23,24,25,26]
l = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
nu = [1495.00,1598.69,1700.91,1802.31,1904.61,2007.58,2110.91,2214.22,2317.32, \
      2420.90,2525.07,2629.20,2733.61,2838.40]
sig = [0.07,0.07,0.08,0.07,0.06,0.05,0.04,0.05,0.05,0.08,0.16,0.18,0.46,0.78]
numax = 2200

x = np.linspace(1400,3000,1601)
x3 = (x-1400) * (x-numax)**2 * 1e-9

# plt.figure()
# plt.scatter(nu,sig)
# plt.plot(x,x3)
# plt.show()

df = pd.read_csv('/home/bmr135/git_AIMS/AIMS/AIMS_BEN/M127_MS_386t',delimiter=r'\s+')
# print(df)
nmx = 1332.62326228

x = np.linspace(100,2500,2401)
x3 = ((df['nu_theo']-nmx) * (2*df['nu_theo']-nmx)**2 * 1e-9)/30 + 0.14

# plt.figure()
# plt.scatter(df['nu_theo'],x3)
# # plt.plot(x,x3)
# plt.show()

df['nu_err'] = x3

df0 = df[df['l']==0]
df1 = df[df['l']==1]
df1 = df1.reset_index(drop=True)
df2 = df[df['l']==2]
df2 = df2.reset_index(drop=True)

a=[]
for i in range(len(df0)-1):
    if df0['nu_err'][i] > df0['nu_err'][i+1]:
        a.append(i)
b=a[0]
for i in range(0,b,1):
    df0['nu_err'][i] = df0['nu_err'][i] + (df0['nu_err'][b] - df0['nu_err'][i])*1.25

a=[]
for i in range(len(df1)-1):
    if df1['nu_err'][i] > df1['nu_err'][i+1]:
        a.append(i)
b=a[0]
for i in range(0,b,1):
    df1['nu_err'][i] = df1['nu_err'][i] + (df1['nu_err'][b] - df1['nu_err'][i])*1.25

a=[]
for i in range(len(df2)-1):
    if df2['nu_err'][i] > df2['nu_err'][i+1]:
        a.append(i)

for i in range(0,a[0],1):
    df2['nu_err'][i] = df2['nu_err'][i] + (df2['nu_err'][a[0]] - df2['nu_err'][i])*1.25

df3 = pd.concat([df0,df1,df2],ignore_index=True)
# print(df3)

plt.figure()
plt.scatter(df3['nu_theo'],x3)
plt.scatter(df3['nu_theo'],df3['nu_err'])
plt.show()

df3.to_csv('/home/bmr135/git_AIMS/AIMS/AIMS_BEN/MS_test_v386',index=False,sep=' ')
