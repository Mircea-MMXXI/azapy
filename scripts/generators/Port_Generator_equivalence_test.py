# Examples
import pandas as pd
import time

import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# Collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = '2023-03-30'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT', 'OIH']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

tt = {}
portx = {}
port = {}

freq = 'M'
hlength = 1.25
fixoffset = 0

#==============================================================================
pnamex = 'InvVolx'
aa = az.InvVolEngine(freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'InvVol'
tic = time.perf_counter()
pp = az.Port_InvVol(mktdata, freq=freq, fixoffset=fixoffset, 
                    histoffset=hlength,pname=pname)
port[pname] = pp.set_model(hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
pnamex = 'InvVarx'
aa = az.InvVarEngine(freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'InvVar'
tic = time.perf_counter()
pp = az.Port_InvVar(mktdata, freq=freq, fixoffset=fixoffset, 
                    histoffset=hlength, pname=pname)
port[pname] = pp.set_model(hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
pnamex = 'InvDDx'
aa = az.InvDDEngine(freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'InvDD'
tic = time.perf_counter()
pp = az.Port_InvDD(mktdata, freq=freq, fixoffset=fixoffset, 
                   histoffset=hlength, pname=pname)
port[pname] = pp.set_model(hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
pnamex = 'Kellyx'
aa = az.KellyEngine(freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'Kelly'
tic = time.perf_counter()
pp = az.Port_Kelly(mktdata, freq=freq, fixoffset=fixoffset, 
                   histoffset=hlength, pname=pname)
port[pname] = pp.set_model(hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
ww = pd.Series(1, index=symb)

pnamex = 'ConstWx'
aa = az.ConstWEngine(ww)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'ConstW'
tic = time.perf_counter()
pp = az.Port_ConstW(mktdata, freq=freq, fixoffset=fixoffset, 
                    histoffset=hlength,pname=pname)
port[pname] = pp.set_model(ww)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
alpha = [0.95, 0.9]
coef = [1] * len(alpha)
rtype = 'Sharpe'
mu0 = 0

pnamex = 'CVaRx'
aa = az.CVaRAnalyzer(alpha=alpha, coef=coef, rtype=rtype, mu0=mu0, 
                     freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'CVaR'
tic = time.perf_counter()
pp = az.Port_CVaR(mktdata, freq=freq, fixoffset=fixoffset, 
                  histoffset=hlength, pname=pname)
port[pname] = pp.set_model(alpha=alpha, coef=coef, rtype=rtype, mu0=mu0, 
                           hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
alpha = [0.90, 0.80]
coef = [1] * len(alpha)
rtype = 'Sharpe'
mu0 = 0

pnamex = 'SMCRx'
aa = az.SMCRAnalyzer(alpha=alpha, coef=coef, rtype=rtype, mu0=mu0, 
                     freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'SMCR'
tic = time.perf_counter()
pp = az.Port_SMCR(mktdata, freq=freq, fixoffset=fixoffset, 
                  histoffset=hlength, pname=pname)
port[pname] = pp.set_model(alpha=alpha, coef=coef, rtype=rtype, mu0=mu0, 
                           hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
alpha = [0.65, 0.6]
coef = [1] * len(alpha)
rtype = 'Sharpe'
mu0 = 0

pnamex = 'EVaRx'
aa = az.EVaRAnalyzer(alpha=alpha, coef=coef, rtype=rtype, mu0=mu0, 
                     freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'EVaR'
tic = time.perf_counter()
pp = az.Port_EVaR(mktdata, freq=freq, fixoffset=fixoffset, 
                  histoffset=hlength, pname=pname)
port[pname] = pp.set_model(alpha=alpha, coef=coef, rtype=rtype, mu0=mu0, 
                           hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
alpha = [-0.01, 0.0, +0.01]
coef = [1] * len(alpha)
detrended = True
rtype = 'Sharpe'
mu0 = 0

pnamex = 'BTADx'
aa = az.BTADAnalyzer(alpha=alpha, coef=coef, detrended=detrended, 
                     rtype=rtype, mu0=mu0, freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'BTAD'
tic = time.perf_counter()
pp = az.Port_BTAD(mktdata, freq=freq, fixoffset=fixoffset, 
                  histoffset=hlength, pname=pname)
port[pname] = pp.set_model(alpha=alpha, coef=coef, detrended=detrended, 
                           rtype=rtype, mu0=mu0, hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
alpha = [-0.01, 0.0, +0.01]
coef = [1] * len(alpha)
detrended = True
rtype = 'Sharpe'
mu0 = 0

pnamex = 'BTSDx'
aa = az.BTSDAnalyzer(alpha=alpha, coef=coef, detrended=detrended, 
                     rtype=rtype, mu0=mu0, freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'BTSD'
tic = time.perf_counter()
pp = az.Port_BTSD(mktdata, freq=freq, fixoffset=fixoffset, 
                  histoffset=hlength, pname=pname)
port[pname] = pp.set_model(alpha=alpha, coef=coef, detrended=detrended, 
                           rtype=rtype, mu0=mu0, hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
maxorder = 5
coef = [1/maxorder] * maxorder
rtype = 'Sharpe'
mu0 = 0

pnamex = 'MADx'
aa = az.MADAnalyzer(coef=coef, rtype=rtype, mu0=mu0, 
                    freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'MAD'
tic = time.perf_counter()
pp = az.Port_MAD(mktdata, freq=freq, fixoffset=fixoffset, 
                 histoffset=hlength, pname=pname)
port[pname] = pp.set_model(coef=coef, rtype=rtype, mu0=mu0, hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
maxorder = 3
coef = [1/maxorder] * maxorder
rtype = 'Sharpe'
mu0 = 0

pnamex = 'LSDx'
aa = az.LSDAnalyzer(coef=coef, rtype=rtype, mu0=mu0, 
                    freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'LSD'
tic = time.perf_counter()
pp = az.Port_LSD(mktdata, freq=freq, fixoffset=fixoffset, 
                 histoffset=hlength, pname=pname)
port[pname] = pp.set_model(coef=coef, rtype=rtype, mu0=mu0, hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
rtype = 'Sharpe'
mu0 = 0

pnamex = 'MVx'
aa = az.MVAnalyzer(rtype=rtype, mu0=mu0, freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'MV'
tic = time.perf_counter()
pp = az.Port_MV(mktdata, freq=freq, fixoffset=fixoffset, 
                histoffset=hlength, pname=pname)
port[pname] = pp.set_model(rtype=rtype, mu0=mu0, hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
rtype = 'Sharpe'
mu0 = 0

pnamex = 'SDx'
aa = az.SDAnalyzer(rtype=rtype, mu0=mu0, freq=freq, hlength=hlength)
mod = az.ModelPipeline([aa])

ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
                        histoffset=hlength, pname=pnamex)
tic = time.perf_counter()
portx[pnamex] = ppx.set_model(mod)
toc = time.perf_counter()
ttx = toc - tic
print(f"time {pnamex}': {ttx}")

pname = 'SD'
tic = time.perf_counter()
pp = az.Port_SD(mktdata, freq=freq, fixoffset=fixoffset, 
                histoffset=hlength, pname=pname)
port[pname] = pp.set_model(rtype=rtype, mu0=mu0, hlength=hlength)
toc = time.perf_counter()
tt0 = toc - tic
print(f"time {pname}': {tt0}")

tt[pname] = [ttx, tt0]
ppc = az.Port_Simple([portx[pnamex], port[pname]])
ppc.set_model()
_ = ppc.port_view_all(componly=True)

#==============================================================================
# uncomment the lines below if you want to include Gini risk measure
# the computations may take some time

# rtype = 'Sharpe'
# mu0 = 0

# pnamex = 'GINIx'
# aa = az.GINIAnalyzer(rtype=rtype, mu0=mu0, freq=freq, hlength=hlength)
# mod = az.ModelPipeline([aa])

# ppx = az.Port_Generator(mktdata, freq=freq, fixoffset=fixoffset, 
#                         histoffset=hlength, pname=pnamex)
# tic = time.perf_counter()
# portx[pnamex] = ppx.set_model(mod)
# toc = time.perf_counter()
# ttx = toc - tic
# print(f"time {pnamex}': {ttx}")

# pname = 'GINI'
# tic = time.perf_counter()
# pp = az.Port_GINI(mktdata, freq=freq, fixoffset=fixoffset, 
#                   histoffset=hlength, pname=pname)
# port[pname] = pp.set_model(rtype=rtype, mu0=mu0, hlength=hlength)
# toc = time.perf_counter()
# tt0 = toc - tic
# print(f"time {pname}': {tt0}")

# tt[pname] = [ttx, tt0]
# ppc = az.Port_Simple([portx[pnamex], port[pname]])
# ppc.set_model()
# _ = ppc.port_view_all(componly=True)

#==============================================================================
# collect all


pts = pd.concat(list(port.values()), axis=1)
#pts.to_csv("NewPortAll.csv")
ppall = az.Port_Simple(list(port.values()))
ppall.port_view_all(componly=True)
ppall.port_perf(componly=True)

ttd = pd.DataFrame.from_dict(tt, orient='index')
ttd.columns = ['Model', 'Compact']
ttd['Diff'] = ttd['Model'] - ttd['Compact']
print(f"comparison of time executions\n{ttd}")

