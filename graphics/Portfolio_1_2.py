# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 09:47:11 2021

@author: mirce
"""
# Examples
import pandas as pd

import azapy as az

#=============================================================================
# Collect some market data
sdate = pd.to_datetime("2012-01-01")
edate = pd.to_datetime('today')
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "./MkTdata"

# force=True read from alphavantage server
# force=False read from local directory if data exists
mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                     dir=mktdir, force=False) 


# Setup CVaR parameters
alpha = [0.91, 0.90, 0.85]
coef = [0.1, 0.1, 0.8]

#=============================================================================
# Compute C-Sharpe optimal portfolio
p4 = az.Port_CVaR(mktdata, pname='CVaRPort')
 
import time
tic = time.perf_counter()

port4 = p4.set_model(mu=0., alpha=alpha, coef=coef, method='ecos')   

toc = time.perf_counter()
print(f"time get_port: {toc-tic}")

ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)

p5 = az.Port_InvVol(mktdata, pname='Portfolio')    

import time
tic = time.perf_counter()

port5 = p5.set_model()   

toc = time.perf_counter()
print(f"time get_port: {toc-tic}")

ww = p5.get_weights()
ffile = 'graphics/Portfolio_1.png'
df = p5.port_view(saveto=ffile)
#fig1 = df.plot().get_figure()

hfile = 'graphics/Portfolio_2.png'
df = p5.port_view(fancy=True, saveto=hfile)


