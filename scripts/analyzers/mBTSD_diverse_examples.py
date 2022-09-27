import numpy as np
import pandas as pd
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
# Define mBTSD measure parameters alpha and coef
alpha = np.array([-0.01, 0.0, 0.01])
# equal weighted risk mixture
coef = np.full(len(alpha), 1/len(alpha))

# build the analyzer object
cr1 = az.BTSDAnalyzer(alpha, coef, mktdata)

#=============================================================================
# plot diverse-efficient frontier
# a choice for targeted expected rate of return
mu = cr1.muk.mean()
cr1.getWeights("Diverse", mu=mu)

rr = np.linspace(np.min(cr1.muk), np.max(cr1.muk), num=10)
div = np.zeros(len(rr))
rate = np.zeros(len(rr))
for i in range(len(rr)):
    cr1.getWeights("Diverse", mu=rr[i])
    div[i] = cr1.diverse
    rate[i] = cr1.RR
    
res = pd.DataFrame({"RR": rr, "Div": div, "Rate": rate})
_ = res.plot("Rate", "Div")

#=============================================================================
# maximum diversified portfolio
ww1 = cr1.getWeights("MaxDiverse")

# manual computation (must have identical results)
ww2 = cr1.getWeights("Diverse", mu=cr1.muk.min())

ww_comp = pd.DataFrame({"Direct": ww1, 
                        "Manual": ww2, 
                        "Diff": ww1 - ww2})
print(f"compare diverse computations\n {ww_comp}")

#=============================================================================
# diverse-efficient portfolio with the same diversification factor as EWP

ww1 = cr1.getWeights("InvNdiverse")
div1 = cr1.diverse

ww2 = np.full(len(symb), 1/len(symb))
div2 = cr1.getDiversification(ww2)

print(f"diversification factor comp\n {div1} {div2}, {div2-div1}")

#=============================================================================
# diverse-efficient portfolio with the same diversification factor as 
# a benchmark portfolio

# benchmark portfolio weights (arbitrary in this example)
ww = range(1, len(symb) + 1) 
ww = ww / np.sum(ww)
ww_diverse = cr1.getDiversification(ww)

ww1 = cr1.getWeights("InvNdiverse", ww0=ww)
ww1_diverse = cr1.diverse

print(f"diverse comp\n benchmark: {ww_diverse} efficient: {ww1_diverse}\
      diff: {ww1_diverse -  ww1_diverse}")

ww_comp = pd.DataFrame({"Benchmark": ww, "Diverse": ww1})
print(f"compare diverse computations\n {ww_comp}")

#=============================================================================
# diverse-efficient portfolio with same expected rate of return as EWP

ww1 = cr1.getWeights("InvNrr")
RR1 = cr1.RR
ww1 = cr1.ww

# manual computation (must have identical results)
RR2 = cr1.muk.mean()
wwm = cr1.getWeights("MaxDiverse")
RRm = cr1.RR
dd = 1 if RR2 >= RRm else -1
ww2 = cr1.getWeights('Diverse', mu=RR1, d=dd)

print(f"RR comp\n EWP {RR2} MaxDivers {RRm} InvNrr {RR1}")
 
comp = pd.DataFrame({"InvNrr": ww1, "manual": ww2, "Diff": ww2 - ww1})
print(f"InvNrr weights comparison\n{comp}")

#=============================================================================
# diverse-efficient portfolio with same expected rate of return as 
# a benchmark portfolio

# benchmark portfolio weights (arbitrary in this example)
ww = range(1, len(symb) + 1) 
ww = ww / np.sum(ww)
ww_diverse = cr1.getDiversification(ww)
ww_RR = cr1.RR

ww1 = cr1.getWeights("InvNrr", ww0=ww)
ww1_diverse = cr1.diverse
ww1_RR = cr1.RR

print(f"diverse comp\n benchmark: {ww_diverse} efficient: {ww1_diverse}\
      diff: {ww1_diverse -  ww1_diverse}")
      
print(f"RR comp\n benchmark: {ww_RR} efficient: {ww1_RR}\
      diff: {ww1_RR -  ww1_RR}")

ww_comp = pd.DataFrame({"Benchmark": ww, "Diverse": ww1})
print(f"compare diverse computations\n {ww_comp}")

#=============================================================================
# examples of viewFrontiers

# 10 (random in this example) additional portfolios to be added to the plot
rng = np.random.RandomState(42)
addp = {}
for i in range(10):
    addp['p' + str(i+1)] = rng.dirichlet([0.5] * len(symb))
addport = pd.DataFrame().from_dict(addp, 'index', columns=symb)

# default views
# expected rate of retun vs risk 
fdat = cr1.viewFrontiers()
# Sharpe vs expected rate of return
_ = cr1.viewFrontiers(fig_type='Sharpe_RR', data=fdat)
# diversification factor vs expected rate of return 
# note that diverse-efficient and diverse-inefficient frontiers are missing
# the dotted lines are the efficient and inefficient frontiers
_ = cr1.viewFrontiers(fig_type='Diverse_RR', data=fdat)


# custom views
opt = {'tangent': True, 'title': 'Custom view 1', 'minrisk_label': 'mRx', 
       'sharpe_label': 'sharpe', 'addport_label': True, 'xlabel': "R-R"}
fd1 = cr1.viewFrontiers(minrisk=True, efficient=30, inefficient=30, 
                        maxdiverse=True, 
                        diverse_efficient=30, diverse_inefficient=30,
                        invNdiverse=True, invNrr=True,
                        randomport=10,
                        options=opt, addport=addport)
# note that the options values are not propagated to next plot 
_ = cr1.viewFrontiers(fig_type='RR_risk', title="Custom view 2", 
                      invNdiverse_label=None, invNrr_label=None, data=fd1)
_ = cr1.viewFrontiers(fig_type='Sharpe_RR', xlabel="RR", 
                      invNdiverse_label=None, data=fd1)
_ = cr1.viewFrontiers(fig_type='Diverse_RR', invNrisk_label=None, data=fd1)