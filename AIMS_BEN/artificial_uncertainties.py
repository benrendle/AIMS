''' Frequency uncertainties based on 16Cyg A distribution '''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

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

df = pd.read_csv('/home/bmr135/git_AIMS/AIMS/AIMS_BEN/MS_test',delimiter=r'\s+')
print(df)
nmx = 1811.511
