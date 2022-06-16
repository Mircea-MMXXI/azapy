# Compare BTSD for L=1, alpha1=0 and detrended RR with first order mLSSD
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
# set BTSD for alpha=0 (default) detrended
alpha = [0.]
coef = [1.]
cr1 = az.BTSDAnalyzer(alpha, coef, mktdata, detrended=True)
# set LSSD first order (default)
coef = [1]
cr2 = az.LSSDAnalyzer(coef, mktdata)

# collect reference data 
_ = cr1.getWeights(mu=0.)
RR_ = cr1.RR
sharpe_ = cr1.sharpe

# rtype collection
rtype_collection = ['Risk', 'MinRisk', 'InvNrisk', 'RiskAverse', 'Sharpe', 
                    'Sharpe2' ]

# loop over all rtype's
for rtype in rtype_collection:
    ww1 = cr1.getWeights(rtype=rtype, mu=RR_, aversion=sharpe_)
    RR1 = cr1.RR
    risk1 = cr1.risk
    
    ww2 = cr2.getWeights(rtype=rtype, mu=RR_, aversion=sharpe_)
    RR2 = cr2.RR
    risk2 = cr2.risk
    
    print("\nComparisons: BTSD (alpha0=0) vs. LSSD (first order) "
          + f"rtype = {rtype} \n")
    print(f"risk: BTSD {risk1} LSSD {risk2} Diff {risk1-risk2}")
    print(f"RR: BTSD {RR1} LSSD {RR2} Diff {RR1 - RR2}")
    ww = pd.DataFrame({'BTSD': ww1, 'LSSD': ww2, 'Diff': ww1 - ww2})
    print(f"weights:\n{ww}")
