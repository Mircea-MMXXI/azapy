# Compare BTSD for L=1, alpha1=0 and detrended RR with first order mLSD
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
cr1 = az.BTSDAnalyzer(mktdata=mktdata, detrended=True)
# set LSD first level (default)
cr2 = az.LSDAnalyzer(mktdata=mktdata)

# collect reference data for Sortino optimal portfolio
_ = cr1.getWeights()
RR_ = cr1.RR
sharpe_ = cr1.sharpe

# rtype collection
rtype_collection = ['Risk', 'MinRisk', 'InvNrisk', 'RiskAverse', 'Sharpe', 
                    'Sharpe2' ]

# loop over all rtype's and compare
for rtype in rtype_collection:
    ww1 = cr1.getWeights(rtype=rtype, mu=RR_, aversion=sharpe_)
    RR1 = cr1.RR
    risk1 = cr1.risk
    
    ww2 = cr2.getWeights(rtype=rtype, mu=RR_, aversion=sharpe_)
    RR2 = cr2.RR
    risk2 = cr2.risk
    
    print("\nComparisons: BTSD (alpha0=0) vs. LSD (first level) "
          + f"rtype = {rtype} \n")
    print(f"risk: BTSD {risk1} LSD {risk2} Diff {risk1-risk2}")
    print(f"RR: BTSD {RR1} LSD {RR2} Diff {RR1 - RR2}")
    ww = pd.DataFrame({'BTSD': ww1, 'LSD': ww2, 'Diff': ww1 - ww2})
    print(f"weights:\n{ww}")
