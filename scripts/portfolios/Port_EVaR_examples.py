# Examples
import pandas as pd 
import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# Collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT', 'OIH']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#==============================================================================
# Define mEVaR measure parameters alpha and coef
alpha = [0.75, 0.65]
coef = [1.] * len(alpha)
hlength = 1.25
portname = 'mEVaR'

# set Port_EVaR class
p4 = az.Port_EVaR(mktdata, pname=portname)
 
#==============================================================================
# Beyond this point any section can be run independently 
#==============================================================================
# Sharpe optimal portfolio for 0 risk free rate
rtype = 'Sharpe'
mu0 = 0.
port4 = p4.set_model(alpha=alpha, coef=coef, rtype=rtype, mu0=mu0, 
                     hlength=hlength)   

# plots
_ = p4.port_view(title=portname + "-Sharpe", ylabel="price ($)")
_ = p4.port_view_all(title=portname + " Portfolio", ylabel="relative move")

# performance monitoring
performance = p4.port_perf()
drawdowns = p4.port_perf()
aret = p4.port_annual_returns()
qret = p4.port_quarterly_returns()
mret = p4.port_monthly_returns()
pret = p4.port_period_returns()
with pd.option_context('display.max_columns', None):
    print(f"Performance\n{performance.round(4)}")
    print(f"Portfolio Historical Drawdowns\n{drawdowns.round(4)}")
    print(f"Portfolio Annual Returns\n{aret.round(4)}")
    print(f"Portfolio Quarterly Returns\n{qret.round(4)}")
    print(f"Portfolio Monthly Returns\n{mret.round(4)}")
    print(f"Portfolio Period Returns\n{pret.round(2)}")

# accounting information
ww = p4.get_weights()
nshares = p4.get_nshares()
accinfo = p4.get_account()
with pd.option_context('display.max_columns', None):
    print(f"Portfolio Historical Weights\n{ww.round(4)}")
    print(f"Portfolio Numbers of Shares\n{nshares}")
    print(f"Portfolio Rolling Accounting Information\n{accinfo.round(0)}")

#==============================================================================
# compare several standard strategies with equal weighted portfolio
# MaxDiverse - maximum diversified portfolio
# MinRisk - minimum risk portfolio
# Sharpe - maximum generalized Sharpe portfolio
# InvNrisk - optimal-risk portfolio with same risk as EWP
# InvNdrr - optimal-diversified portfolio with 
#           same expected rate of return as EWP
# InvNdiverse - optimal-diversified portfolio with 
#               same diversification factor as EWP
rtypes = ['MaxDiverse', 'MinRisk', 'Sharpe', 
          'InvNrisk', 'InvNdrr', 'InvNdiverse']

port = []
for rtype in rtypes:
    print(f"rtype {rtype}")
    port4 = p4.set_model(alpha=alpha, coef=coef, rtype=rtype, hlength=hlength, method='excp')
    port4.columns = [rtype]
    port.append(port4)
    
# add EWP (1/N  portfolio)
p5 = az.Port_ConstW(mktdata, pname="1/N")
port5 = p5.set_model()
port.append(port5)
    
# compare
pp = az.Port_Simple(port)
_ = pp.set_model()
_ = pp.port_view_all(componly=True, 
                     title=portname + " Portfolios - Relative Performance", 
                     xlabel="year")
# compare performances
perfs = pp.port_perf(componly=True)
print(f"Portfolio Performances\n{perfs.round(4)}")
arets = pp.port_annual_returns(withcomp=True, componly=True)
print(f"Annual Returns\n{arets.round(4) * 100}")

#==============================================================================
# Other examples
# Optimal-risk portfolio for fixed aversion factor
rtype = 'RiskAverse'
aversion = 0.4
port4 = p4.set_model(alpha=alpha, coef=coef, rtype=rtype, aversion=aversion, 
                     hlength=hlength)   

# plots
_ = p4.port_view(title=portname + " aversion = " + str(aversion), 
                 ylabel="price ($)")
_ = p4.port_view_all(title=portname + " Portfolio", ylabel="relative move")

# performance monitoring
performance = p4.port_perf()
drawdowns = p4.port_perf()
aret = p4.port_annual_returns()
qret = p4.port_quarterly_returns()
mret = p4.port_monthly_returns()
pret = p4.port_period_returns()
with pd.option_context('display.max_columns', None):
    print(f"Performance\n{performance.round(4)}")
    print(f"Portfolio Historical Drawdowns\n{drawdowns.round(4)}")
    print(f"Portfolio Annual Returns\n{aret.round(4)}")
    print(f"Portfolio Quarterly Returns\n{qret.round(4)}")
    print(f"Portfolio Monthly Returns\n{mret.round(4)}")
    print(f"Portfolio Period Returns\n{pret.round(2)}")

# accounting information
ww = p4.get_weights()
nshares = p4.get_nshares()
accinfo = p4.get_account()
with pd.option_context('display.max_columns', None):
    print(f"Portfolio Historical Weights\n{ww.round(4)}")
    print(f"Portfolio Numbers of Shares\n{nshares}")
    print(f"Portfolio Rolling Accounting Information\n{accinfo.round(0)}")
    
#==============================================================================
# Optimal-risk portfolio for targeted expected rate of return
rtype = 'Risk'
mu = 0.06
port4 = p4.set_model(alpha=alpha, coef=coef, rtype=rtype, mu=mu, 
                     hlength=hlength)   

# plots
_ = p4.port_view(title=portname + " - Optimal Risk for mu = " + str(mu), 
                 ylabel="price ($)")
_ = p4.port_view_all(title=portname + " Portfolio", ylabel="relative move")

# performance monitoring
performance = p4.port_perf()
drawdowns = p4.port_perf()
aret = p4.port_annual_returns()
qret = p4.port_quarterly_returns()
mret = p4.port_monthly_returns()
pret = p4.port_period_returns()
with pd.option_context('display.max_columns', None):
    print(f"Performance\n{performance.round(4)}")
    print(f"Portfolio Historical Drawdowns\n{drawdowns.round(4)}")
    print(f"Portfolio Annual Returns\n{aret.round(4)}")
    print(f"Portfolio Quarterly Returns\n{qret.round(4)}")
    print(f"Portfolio Monthly Returns\n{mret.round(4)}")
    print(f"Portfolio Period Returns\n{pret.round(2)}")

# accounting information
ww = p4.get_weights()
nshares = p4.get_nshares()
accinfo = p4.get_account()
with pd.option_context('display.max_columns', None):
    print(f"Portfolio Historical Weights\n{ww.round(4)}")
    print(f"Portfolio Numbers of Shares\n{nshares}")
    print(f"Portfolio Rolling Accounting Information\n{accinfo.round(0)}")
    
#==============================================================================
# Optimal-diversified portfolio for targeted expected rate of return
rtype = 'Diverse'
mu = 0.06
port4 = p4.set_model(alpha=alpha, coef=coef, rtype=rtype, mu=mu, 
                     hlength=hlength)   

# plots
_ = p4.port_view(title=portname + " - Optimal Diverse. for mu = " + str(mu), 
                 ylabel="price ($)")
_ = p4.port_view_all(title=portname + " Portfolio", ylabel="relative move")

# performance monitoring
performance = p4.port_perf()
drawdowns = p4.port_perf()
aret = p4.port_annual_returns()
qret = p4.port_quarterly_returns()
mret = p4.port_monthly_returns()
pret = p4.port_period_returns()
with pd.option_context('display.max_columns', None):
    print(f"Performance\n{performance.round(4)}")
    print(f"Portfolio Historical Drawdowns\n{drawdowns.round(4)}")
    print(f"Portfolio Annual Returns\n{aret.round(4)}")
    print(f"Portfolio Quarterly Returns\n{qret.round(4)}")
    print(f"Portfolio Monthly Returns\n{mret.round(4)}")
    print(f"Portfolio Period Returns\n{pret.round(2)}")

# accounting information
ww = p4.get_weights()
nshares = p4.get_nshares()
accinfo = p4.get_account()
with pd.option_context('display.max_columns', None):
    print(f"Portfolio Historical Weights\n{ww.round(4)}")
    print(f"Portfolio Numbers of Shares\n{nshares}")
    print(f"Portfolio Rolling Accounting Information\n{accinfo.round(0)}")
    
#==============================================================================