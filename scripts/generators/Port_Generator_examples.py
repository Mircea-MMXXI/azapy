# Examples
import numpy as np

import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'

symb = ['GLD', 'TLT', 'IHI', 'VGT', 'OIH',
        'XAR', 'XBI', 'XHE', 'XHS', 'XLB',
        'XLE', 'XLF', 'XLI', 'XLK', 'XLU', 
        'XLV', 'XLY', 'XRT', 'SPY', 'ONEQ', 
        'QQQ', 'DIA', 'ILF', 'XSW', 'PGF', 
        'IDV', 'JNK', 'HYG', 'SDIV', 'VIG', 
        'SLV', 'AAPL', 'MSFT', 'AMZN', 'GOOG', 
        'IYT', 'VGI', 'IWM', 'BRK-B', 'ITA']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir, 
                     verbose=False)

#==============================================================================
# build CorrClusterSelector
# correlation threshold 
# i.e., assets with correlation larget than this value are clustered together
corr_threshold = 0.98
# tenor of relevant rate of return 
# its is set to quarterly across all components - since this is a quarterly 
# rebalanced portfolio 
freq = 'Q'
ccs = az.CorrClusterSelector(corr_threshold=corr_threshold, freq=freq)

# build a DualMomentumSelctor
# maximum number of selected symbol
nw = 5
# minimum number of symbols with positive momentum 
#   for a full capital allocation - in our case
#   roughly 70% of the initial number of symbols (moderat aggresive)
ths = np.floor(len(symb) * 0.7)
dms = az.DualMomentumSelector(nw=nw, threshold=ths)

# buid a CVaR optimizer
# equaly weithed mixture of CVaR measures 
alpha = [0.95, 0.9]
# calibration period
hlength = 1.
# optimizatin type - maximization of CVaR-Sharpe ratio
rtype = 'Sharpe'
mu0 = 0.01
cvar = az.CVaRAnalyzer(alpha=alpha, freq=freq, hlength=hlength,
                       rtype=rtype, mu0=mu0)

# build the ModelPilpeline
model = az.ModelPipeline([ccs, dms, cvar])

# backtesting
# set the histoffset withthe same value of hlength 
# (historical data needed for the first period calibration)
histoffset = hlength
fixoffset = -1
pp = az.Port_Generator(mktdata, 
                       freq=freq, 
                       fixoffset=fixoffset, 
                       histoffset=histoffset)
# build the historical portfolio
port = pp.set_model(model)
# other info
_ = pp.port_view()
print(f"portfolio drawdown\n{pp.port_drawdown(fancy=True)}")
print("portfolio perfarmance including all components\n"
      f"{pp.port_perf(fancy=True)}")
print(f"portfolio anual returns\n{pp.port_annual_returns()}")
print(f"portfolio monthly returns\n{pp.port_monthly_returns()}")
print("portfolio returns per invesment period\n"
      f"{pp.port_period_returns(fancy=True)}")
print("portfolio performance for each investment period\n" 
      f"{pp.port_period_perf(fancy=True)}")
print(f"portfolio weights for each investment period\n{pp.ww}")

print(f"{pp.port_perf(fancy=True).iloc[0]} * 100")


