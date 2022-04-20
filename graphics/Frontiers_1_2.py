# Frontiers fig 1 and 2
import numpy as np
import azapy as az
 
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

sdate = "2012-01-01"
edate = "2021-07-27"

mktdir = "../MkTdata"

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

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