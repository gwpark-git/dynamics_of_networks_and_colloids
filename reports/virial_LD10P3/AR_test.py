
from numpy import *
import matplotlib.pyplot as plt

import statsmodels.api as sm
# from statsmodels.sandbox.regression.predstd import wls_prediction_std


index_st = 2000
RF_NP0400 = loadtxt('RF_NP0400_NC25.dat')[index_st:,1]




# signal = np.ones(20) + 1e-6*np.random.randn(20)
# ar_mod = AR(signal)
# ar_res = ar_mod.fit(4)

# ar_res.predict(4, 60)
