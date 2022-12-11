# Compare BTAD for L=1, alpha1=0 and detrended RR with first order mMAD 
# they should be identical up to machine precision 
import pandas as pd
import azapy as az
print(f"azapy version {az.version()} >= 1.1.0", flush=True)

#=============================================================================
# Collect some market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = '2021-07-27'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# set BTAD for alpha1=0 (default) detrended
cr1 = az.BTADAnalyzer(mktdata=mktdata, detrended=True)
# set MAD first level (default)
cr2 = az.MADAnalyzer(mktdata=mktdata)

# collect reference data for Omega optimal portfolio
_ = cr1.getWeights()
RR_ = cr1.RR
sharpe_ = cr1.sharpe

# rtype collection
rtype_collection = ['Risk', 'MinRisk', 'InvNrisk', 'RiskAverse', 'Sharpe', 
                    'Sharpe2', 'Diverse', 'Diverse2', 'InvNdiverse', 'InvNdrr']

# loop over all rtype's and compare
for rtype in rtype_collection:
    ww1 = cr1.getWeights(rtype=rtype, mu=RR_, aversion=sharpe_)
    RR1 = cr1.RR
    risk1 = cr1.risk
    
    ww2 = cr2.getWeights(rtype=rtype, mu=RR_, aversion=sharpe_)
    RR2 = cr2.RR
    risk2 = cr2.risk
    
    print("\nComparisons: BTAD (alpha1=0) vs. MAD (first level) "
          + f"rtype = {rtype} \n")
    print(f"risk: BTAD {risk1:f} MAD {risk2:f} Diff {risk1-risk2:f}")
    print(f"RR: BTAD {RR1:f} MAD {RR2:f} Diff {RR1 - RR2:f}")
    ww = pd.DataFrame({'BTAD': ww1, 'MAD': ww2, 'Diff': ww1 - ww2})
    print(f"weights:\n{ww.round(4)}")
