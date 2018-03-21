''' Script to remove models that have no He-core masses '''

import numpy as np
import pandas as pd

df = pd.read_csv('/home/bmr135/git_AIMS/AIMS/AIMS_BEN/CLES_RGB_v3',delimiter=r'\s+',names=['ext','m','r','l','Z','X','a','t','me','xc','p'],skiprows=1)
df = df[df['me'] > 0.0]
df.to_csv('/home/bmr135/git_AIMS/AIMS/AIMS_BEN/CLES_RGB_v3.1', index=False, sep = ' ', header=False)

with open('/home/bmr135/git_AIMS/AIMS/AIMS_BEN/CLES_RGB_v3.1',"r+") as f:
	content = f.read()
	line = "/home/bmr135/GridGiantsClesV0.3/models_grad_rad_under/	.freq"
	f.seek(0,0)
	f.write(line.rstrip('\r\n') + '\n' + content)
