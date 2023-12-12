
# from typing import List, Tuple
from __future__ import print_function
import RobClust
import matlab

#import sys
#sys.path.append("/Applications/MATLAB_R2022a.app/extern/engines/python/dist/matlab/engine/")
#import matlab.engine
# needed for files ops
import numpy as np
import pandas as pd



# import matplotlib.pyplot as plt
# 
# import seaborn as sns; sns.set()  
# 
# from sklearn.linear_model import Ridge
# from sklearn.base import RegressorMixin, BaseEstimator, clone
# 
# from sklearn.metrics import r2_score
# import statsmodels.api as sm
# from statsmodels.api import add_constant
# 
# from collections import defaultdict
# import re
# from os import path


my_Rob = RobClust.initialize()


# Change the inputs
p1=1
p2=0
kmax = 2
clust_fun = 'RLGA'

# import csv
# Import dataset
name2 = pd.read_csv('PY_name_orig.csv')

XX = name2['AuctionPrice']
yy = name2['Aggregated-Positive-InfraMarginalRent-EUR']

# Scale the values
X = XX/np.max(np.absolute(XX))
y = yy/np.max(np.absolute(yy))




# make arrays available to MATLAB
AA = matlab.double(X)
bb = matlab.double(y)

# 
nn = AA.size[1]
pp = AA.size[0]

# Apply the robust clustering methodology with RobClust()
robclustout = my_Rob.RobClust(bb, AA, nn, pp, kmax, clust_fun, p1, p2)

# struct are returned to Python as dictionaries!
print(type(robclustout))

# extract label field from dict.
labelspy=robclustout.get('labels')

# print type and size (debug)
print(type(labelspy))
print(labelspy.size)

# convert to array with numpy
labelnp=np.asarray(labelspy)

# extract label field from dict.
linecoefpy=robclustout.get('linecoef')

# convert to array with numpy
linecoefnp=np.asarray(linecoefpy)

# scalars are treated as floats
alphaoptpy=robclustout.get('alphaopt')
coptpy=robclustout.get('copt')
koptpy=robclustout.get('kopt')
okclusterpy=robclustout.get('okcluster')
slopetestpy=robclustout.get('slopetest')


print(robclustout, sep='\n')
input("Press Enter to continue...")


my_Rob.terminate()