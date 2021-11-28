# Frontiers fig 1 and 2
import numpy as np
import pandas as pd

import azapy as az


sdate = pd.to_datetime("2012-01-01")
edate = pd.to_datetime('today')
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdir = "../MkTdata"

mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                     dir=mktdir, force=False) 

#=============================================================================
coef = np.ones(3)
coef = coef / coef.sum()

cr1 = az.LSSDAnalyzer(coef, mktdata)

print("\nFrontiers evaluations\n")

opt = {'title': "Portfolio frontiers", 'tangent': True}

saveto = "./frontiers_1.png"
rft = cr1.viewFrontiers(musharpe=0., randomport=100, options=opt, saveto=saveto)


saveto = "./frontiers_2.png"
_ = cr1.viewFrontiers(data=rft, fig_type='Sharpe_RR', saveto=saveto)