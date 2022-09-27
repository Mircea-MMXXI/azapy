import numpy as np
import azapy as az
print(f"azapy version {az.getVersion()} >= 1.0.2")

#=============================================================================
# Collect market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2022-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ',]

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Define mSMCR measure parameters alpha and coef
alpha = np.array([0.90, 0.85])
coef = np.full(len(alpha), 1/len(alpha))

# set Port class
p4 = az.Port_SMCR(mktdata, pname='mSMCRPort')
 
#=============================================================================
# Compute MaxDiverse optimal portfolio (example)
rtype = "MaxDiverse"
port4 = p4.set_model(alpha=alpha, coef=coef, rtype=rtype)   

ww = p4.get_weights()
p4.port_view(title="mSMCR Max Diversified Portfolio", ylabel="price ($)")
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()   
p4.get_account(fancy=True)

#=============================================================================
# compare several standard strategies with 1/N (EWP - equal weighted portfolio)
# MaxDiverse - maximum diversified portfolio
# MinRisk - minimum risk portfolio
# Sharpe - maximum generalized Sharpe portfolio
# InvNrisk - efficient portfolio with same risk as EWP
# InvNrr - efficient-diversified portfolio with 
#          same expected rate of return as EWP
# InvNdiverse - efficient-diversified portfolio with 
#               same diversification factor as EWP
rtype = ["MaxDiverse", "MinRisk", "Sharpe", 
         "InvNrisk", "InvNrr", "InvNdiverse"]

port = []
for rty in rtype:
    port4 = p4.set_model(alpha=alpha, coef=coef, rtype=rty)  
    port4.columns = [rty]
    port.append(port4)
    
# add EWP (1/N  portfolio)
p5 = az.Port_ConstW(mktdata, pname="1/N")
port5 = p5.set_model()
port.append(port5)
    
# compare
pp = az.Port_Simple(port)
pp.set_model()
pp.port_view_all(componly=True, 
                 title="mSMCR Portfolios - Relative Performance", 
                 xlabel="year")
# compare performances
print(pp.port_perf(componly=True))
