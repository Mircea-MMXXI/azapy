# Compare BTAD for L=1, alpha1=0 and detrended RR with first order mMAD 
# they should be identical up to machine precision 
import pandas as pd
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# set BTAD for alpha1=0 (default) detrended
alpha = [0.]
coef = [1.]
cr1 = az.BTADAnalyzer(alpha, coef, mktdata, detrended=True)
# set MAD first order (default)
coef = [1]
cr2 = az.MADAnalyzer(coef, mktdata)

# collect reference data 
_ = cr1.getWeights(mu=0.)
RR_ = cr1.RR
sharpe_ = cr1.sharpe

# dict rtype: mu 
mus = {'Risk': RR_, 'MinRisk': 0., 'Sharpe': 0., 'Sharpe2': 0., 'InvNrisk': 0.,
       'RiskAverse': sharpe_}

# loop over all rtype's
for rtype in mus.keys():
    ww1 = cr1.getWeights(mu=mus[rtype], rtype=rtype)
    RR1 = cr1.RR
    risk1 = cr1.risk
    
    ww2 = cr2.getWeights(mu=mus[rtype], rtype=rtype)
    RR2 = cr2.RR
    risk2 = cr2.risk
    
    print("\nComparisons: BTAD (alpha0=0) vs. MAD (first order) "
          + f"rtype = {rtype} \n")
    print(f"risk: BTAD {risk1} MAD {risk2} Diff {risk1-risk2}")
    print(f"RR: BTAD {RR1} MAD {RR2} Diff {RR1 - RR2}")
    ww = pd.DataFrame({'BTAD': ww1, 'MAD': ww2, 'Diff': ww1 - ww2})
    print(f"weights:\n{ww}")