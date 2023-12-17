# Examples
import numpy as np
import pandas as pd
import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = '2021-07-27'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#==============================================================================
# Define mSMCR measure parameters alpha and coef
alpha = np.array([0.85, 0.75])
# equal weighted risk mixture
coef = np.full(len(alpha), 1/len(alpha))
# set now the title of the frontiers plots
title_plot = 'mSMCR frontiers'
hlength = 3.25
method = 'ecos' # default choice

# build the analyzer object
cr1 = az.SMCRAnalyzer(alpha, coef, mktdata, hlength=hlength, method=method)

#==============================================================================
# Beyond this point any section can be run independently 
#==============================================================================
print("\n******************************************************************\n")
print("\n*** Risk of a given portfolio ***")
print("---we choose a random portfolio---")
ww = np.random.dirichlet([0.5] * len(symb))

risk = cr1.getRisk(ww)
status = cr1.status
RR = cr1.RR
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
comp_time = cr1.time_level1
print(f"Risk comp time {comp_time:f}\n "
      f"Portfolio parameters:\nweights {ww.round(4)}\n"
      f"expected rate of return {RR:f}\n"
      f"risk {risk:f}\n"
      f"primary risk comp   {primary_risk.round(6)}\n"
      f"secondary risk comp {secondary_risk.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("\n*** Diversification + Risk of a given portfolio ***")
print("---we choose a random portfolio---")
ww = np.random.dirichlet([0.5] * len(symb))

diverse = cr1.getDiversification(ww)
status = cr1.status
risk = cr1.risk
RR = cr1.RR
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
comp_time = cr1.time_level1
print(f"Diverse + Risk comp time {comp_time:f}\n "
      f"Portfolio parameters:\nweights {ww.round(4)}\n"
      f"expected rate of return {RR:f}\n"
      f"diversification factor {diverse:f}\n"
      f"risk {risk:f}\n"
      f"primary risk comp   {primary_risk.round(6)}\n"
      f"secondary risk comp {secondary_risk.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Optimal risk portfolio for targeted expected rate of return ***")
rtype = 'Risk'
mu = 0.04
ww = cr1.getWeights(rtype, mu)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
comp_time = cr1.time_level1

print(f"rtype {rtype} for mu {mu} "
      f"computation status {status} time {comp_time:f}\n"
      f"optimal weights:\n{ww.round(4)}\n"
      f"expected rate of return {RR:f}\n"
      f"risk {risk:f}\n"
      f"primary risk comp   {primary_risk.round(6)}\n"
      f"secondary risk comp {secondary_risk.round(6)}\n")

print("=== test - compute risk for portfolio with optimal weights ===")
print(f"optimal weights\n{ww}")
risk_test = cr1.getRisk(ww)
status_test = cr1.status
RR_test = cr1.RR
risk_test = cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp

prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"rtype {rtype} for mu {mu} computation status {status}\n"
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Minimum risk porfolio ***")
rtype = 'MinRisk'
ww = cr1.getWeights(rtype)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
comp_time = cr1.time_level1
print(f"rtype {rtype} computation status {status} comp time {comp_time}\n")

print("=== test - compute optimal risk portfolio for "
      "mu = min component expected rate of return ===")
# results should be identical
rtype_test = 'Risk'
mu = max(cr1.muk.min(), 0)
ww_test = cr1.getWeights(rtype_test, mu)
status_test = cr1.status
RR_test = cr1.RR
risk_test= cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
comp_time_test = cr1.time_level1
print(f"test rtype {rtype_test} computation status {status_test} "
      f"time {comp_time_test:f}\n")

weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"optimal weigts\n{weights.round(4)}\n" + 
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test}\n"
      f"primary risk comp\n{prc}\n"
      f"secondary risk comp\n{src}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Sharpe optimal portfolio - max Sharpe ratio ***")
rtype = 'Sharpe' 
mu0 = 0. # 0. risk free rate (default value)
ww = cr1.getWeights(rtype, mu0=mu0)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
sharpe = cr1.sharpe
comp_time = cr1.time_level1
print(f"rtype {rtype} computation status {status} comp time {comp_time}\n")

print("=== test1 - compute risk for portfolio with Sharpe weights ===")
risk_test = cr1.getRisk(ww)
status_test = cr1.status
ww_test = cr1.ww
RR_test = cr1.RR
risk_test = cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
comp_time_test = cr1.time_level1
sharpe_test = (RR_test - mu0) / risk_test

prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"test computation status {status_test} comp time {comp_time_test:f}\n\n"
      f"optimal weigts\n{weights.round(4)}\n" 
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"sharpe {sharpe:f} test {sharpe_test:f} diff {sharpe - sharpe_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

print("=== test2 - compute optimal risk portfolio for "
      "mu equal to Sharpe portfolio expected rate of return ===")
ww_test = cr1.getWeights('Risk', mu=RR)
status_test = cr1.status
RR_test = cr1.RR
risk_test = cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
comp_time_test = cr1.time_level1
sharpe_test = (RR_test - mu0) / risk_test
print(f"test computation status {status_test} time {comp_time_test:f}\n")

weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"optimal weigts\n{weights.round(4)}\n" 
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"sharpe {sharpe:f} test {sharpe_test:f} diff {sharpe - sharpe_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Sharpe optimal portfolio - min inverse Sharpe ratio ***")
rtype = 'Sharpe2' 
mu0 = 0. # 0. risk free rate (default value)
ww = cr1.getWeights(rtype, mu0=mu0)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
sharpe = cr1.sharpe
comp_time = cr1.time_level1
print(f"rtype {rtype} computation status {status} comp time {comp_time:f}\n")

print("=== test - compare Sharpe with Sharpe2 ===")
rtype_test = 'Sharpe' 
mu0_test = mu0
ww_test = cr1.getWeights(rtype_test, mu0=mu0_test)
status_test = cr1.status
RR_test = cr1.RR
risk_test = cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
sharpe_test = cr1.sharpe
comp_time_test = cr1.time_level1

weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"rtype {rtype_test} computation status {status_test} "
      f"comp time {comp_time_test:f}\n\n"
      f"optimal weigts\n{weights.round(4)}\n"
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"sharpe {sharpe:f} test {sharpe_test:f} diff {sharpe - sharpe_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Optimal risk portfolio for fixed risk-aversion factor ***")
# set the aversion factor equal to Sharpe ratio for mu0=0.
# compute Sharpe portfolio for mu0=0. (default)
rtype_test = 'Sharpe'
ww_test = cr1.getWeights(rtype_test)
status_test = cr1.status
RR_test = cr1.RR
risk_test = cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
sharpe_test = cr1.sharpe
comp_time_test = cr1.time_level1

# actual computation
rtype = 'RiskAverse'
aversion = np.abs(sharpe_test)
print(f"aversion = Sharpe ratio = {aversion:f}")
ww = cr1.getWeights(rtype, aversion=aversion)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
comp_time = cr1.time_level1
print(f"rtype {rtype} computation status {status} comp time {comp_time:f}\n")

print("=== test - compare optimal risk portfolio for aversion factor equal to "
      "Sharpe ratio ===\n=== (must return the Sharpe portfolio) ===")

weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})

print(f"test rtype {rtype_test} computation status {status_test} " 
      f"comp time {comp_time_test:f}\n\n" 
      f"optimal weigts\n{weights.round(4)}\n" 
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n" 
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n" 
      f"primary risk comp\n{prc.round(6)}\n" 
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Optimal risk portfolio "
      "with same risk as a benchmark portfolio ***")
print("\t--------------------------------------------------------------------"
      "\n"
      "\tNote: If the benchmark portfolio risk is greater than the risk\n"
      "\tof single asset portfolio with the highest expected rate of return,\n"
      "\tthen the InvNrisk portfolio defaults to this single asset portfolio."
      "\n"
      "\t--------------------------------------------------------------------"
      "\n")
ww0 = np.random.dirichlet([0.5] * len(symb))
# for equal weighted portfolio uncomment the line below
# ww0 = np.full(len(symb), 1/len(symb))
print(f"benchmark portfolio weights {ww0.round(4)}")

rtype = 'InvNrisk'
ww = cr1.getWeights(rtype, ww0=ww0)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
comp_time = cr1.time_level1
print(f"rtype {rtype} computation status {status} comp time {comp_time:f}\n")

print("=== test1 - compare the risk with the benchmark ===")
symb_max = cr1.muk.idxmax()
ww_s = pd.Series(0., index=symb)
ww_s[symb_max] = 1.
risk_s = cr1.getRisk(ww_s)
risk_test = cr1.getRisk(ww0)
if risk_s < risk_test:
    print(f"benchmark port risk {risk_test:f} smaller than {risk_s:f}\n"
          f"default to single asset portfolio {symb_max}")
    risk_test = risk_s
print(f"risk {risk:f} benchmark risk {risk_test:f} diff {risk - risk_test:f}")


print("\n=== test2 - compare with the optimal risk portfolio for "
      "mu = InvNrisk port expected rate of return ===\n"
      "=== must be the same (up to precision) ===")
rtype_test = 'Risk'
mu_test = RR
ww_test = cr1.getWeights(rtype_test, mu=mu_test)
status_test = cr1.status
RR_test = cr1.RR
risk_test = cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
comp_time_test = cr1.time_level1

weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"rtype {rtype_test:} computation status {status_test} "
      f"comp time {comp_time_test:f}\n\n"
      f"optimal weights\n{weights.round(4)}\n"
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Optimal diversified portfolio for targeted "
      "expected rate of return ***")
rtype = 'Diverse'
mu = 0.04
ww = cr1.getWeights(rtype, mu)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
diverse = cr1.diverse
comp_time = cr1.time_level1
print(f"rtype {rtype} computation status {status} comp time {comp_time:f}\n")

print("=== test - compute risk/diversification for a portfolio with weights "
      "equal to the optimal weights ===")
print(f"optimal weights\n{ww.round(4)}")
diverse_test = cr1.getDiversification(ww)
status_test = cr1.status
RR_test = cr1.RR
risk_test = cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp

prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"diversification factor {diverse:f} test {diverse_test:f} "
      f"diff {diverse - diverse_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Optimal diversified portfolio for targeted "
      "expected rate of return (alternative) ***")
rtype = 'Diverse2'
mu = 0.04
ww = cr1.getWeights(rtype, mu)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
diverse = cr1.diverse
comp_time = cr1.time_level1
print(f"rtype {rtype} computation status {status} comp time {comp_time:f}\n")

print("=== test - compare rtype 'Diverse2' with 'Diverse' ===")
rtype_test = 'Diverse'
ww_test = cr1.getWeights(rtype_test, mu)
status_test = cr1.status
RR_test = cr1.RR
risk_test = cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
diverse_test = cr1.diverse
comp_time_test = cr1.time_level1
print(f"rtype {rtype_test} computation status {status_test} "
      f"comp time {comp_time_test:f}\n")

weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"optimal weights\n{weights.round(4)}\n"
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"diversification factor {diverse:f} test {diverse_test:f} "
      f"diff {diverse - diverse_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Maximum diversified portfolio ***")
rtype = 'MaxDiverse'
ww = cr1.getWeights(rtype)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
diverse = cr1.diverse
comp_time = cr1.time_level1
print(f"rtype {rtype} computation status {status} comp time {comp_time:f}\n")

print("=== test - compute optimal diversified portfolio for "
      "mu = min component expected rate of return ===")
# results should be identical
rtype_test = 'Diverse'
mu = max(cr1.muk.min(), 0)
ww_test = cr1.getWeights(rtype_test, mu)
status_test = cr1.status
RR_test = cr1.RR
risk_test= cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
diverse_test = cr1.diverse
comp_time_test = cr1.time_level1
print(f"test rtype {rtype_test} computation status {status_test} "
      f"time {comp_time_test:f}\n")

weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"optimal weigts\n{weights.round(4)}\n" + 
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"diversification factor {diverse:f} test {diverse_test:f} "
      f"diff {diverse - diverse_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("******************************************************************\n")
print("\n*** Optimal diversified portfolio with same diversification "
      "factor as a benchmark portfolio ***")
ww0 = np.random.dirichlet([0.5] * len(symb))
# for equal weighted portfolio uncomment the line below
# ww0 = np.full(len(symb), 1/len(symb))
print(f"benchmark portfolio weights {ww0.round(4)}")

rtype = 'InvNdiverse'
ww = cr1.getWeights(rtype, ww0=ww0)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
comp_time = cr1.time_level1
diverse = cr1.diverse
print(f"rtype {rtype} computation status {status} comp time {comp_time:f}\n")

print("=== test1 - compare the diversification factors ===")
diverse_test = cr1.getDiversification(ww0)
print(f"diversfication {diverse:f} benchmark port {diverse_test:f} "
      f"diff {diverse - diverse_test:f}\n")

print("=== test2 - compare with optimal diversified portfolio for "
      "mu = InvNdiverse portfolio expected rate of return ===")
rtype_test = 'Diverse'
mu_test = RR
ww_test = cr1.getWeights(rtype_test, mu_test)
status_test = cr1.status
RR_test = cr1.RR
risk_test= cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
diverse_test = cr1.diverse
comp_time_test = cr1.time_level1 
print(f"test rtype {rtype_test} computation status {status_test} "
      f"time {comp_time_test:f}\n")

weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"optimal weights\n{weights.round(4)}\n" + 
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"diversification factor {diverse:f} test {diverse_test:f} "
      f"diff {diverse - diverse_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("******************************************************************\n")
print("\n*** Optimal diversified portfolio with same expected rate of return "
      "as a benchmark portfolio ***")
ww0 = np.random.dirichlet([0.5] * len(symb))
# for equal weighted portfolio uncomment the line below
# ww0 = np.full(len(symb), 1/len(symb))
print(f"benchmark portfolio weights {ww0.round(4)}")

rtype = 'InvNdrr'
ww = cr1.getWeights(rtype, ww0=ww0)
status = cr1.status
RR = cr1.RR
risk = cr1.risk
primary_risk = cr1.primary_risk_comp
secondary_risk = cr1.secondary_risk_comp
comp_time = cr1.time_level1
diverse = cr1.diverse
print(f"rtype {rtype} computation status {status} comp time {comp_time:f}\n")

print("=== test1 - compare the portfolios expected rate of return ===")
_ = cr1.getRisk(ww0)
RR_test = cr1.RR
print(f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n")

print("=== test2 - compare with optimal diversified portfolio for "
      "mu = benchmark portfolio expected rate of return ===")
mu = np.dot(ww0, cr1.muk)
# first compute MaxDiverese portfolio expected rate of return
rtype_test = 'MaxDiverse'
ww_test = cr1.getWeights(rtype_test)
status_test = cr1.status
comp_time_test = cr1.time_level1
print(f"rtype {rtype_test} comp time {comp_time_test:f}")
if status_test == 0:
    mu_maxd = cr1.RR
    # close the branch of the frontiers
    rtype_test = 'Diverse'
    if mu > mu_maxd:
        ww_test = cr1.getWeights(rtype_test, mu, d=1)
        comp_time_test += cr1.time_level1
        print(f"rtype {rtype_test} d=1 comp time {cr1.time_level1:f}\n"
              f"test comp time {comp_time_test:f}")
    elif mu < mu_maxd:
        ww_test = cr1.getWeights(rtype_test, mu, d=-1)
        comp_time_test += cr1.time_level1
        print(f"rtype {rtype_test} d=-1 comp time {cr1.time_level1:f}\n"
              f"test comp time {comp_time_test:f}")
status_test = cr1.status
RR_test = cr1.RR
risk_test= cr1.risk
primary_risk_test = cr1.primary_risk_comp
secondary_risk_test = cr1.secondary_risk_comp
diverse_test = cr1.diverse
print(f"test rtype {rtype_test} computation status {status_test} "
      f"time {comp_time_test:f}\n")

weights = pd.DataFrame({'ww': ww, 'test': ww_test, 'diff': ww - ww_test})
prc = pd.DataFrame({'primary': primary_risk, 'test': primary_risk_test,
                   'diff': primary_risk - primary_risk_test})
src = pd.DataFrame({'secondary': secondary_risk, 'test': secondary_risk_test,
                   'diff': secondary_risk - secondary_risk_test})
print(f"optimal weigts\n{weights.round(4)}\n" + 
      f"expected rate of return {RR:f} test {RR_test:f} diff {RR - RR_test:f}\n"
      f"diversification factor {diverse:f} test {diverse_test:f} "
      f"diff {diverse - diverse_test:f}\n"
      f"risk {risk:f} test {risk_test:f} diff {risk - risk_test:f}\n"
      f"primary risk comp\n{prc.round(6)}\n"
      f"secondary risk comp\n{src.round(6)}\n")

#==============================================================================
print("\n******************************************************************\n")
print("*** Frontiers evaluations - standard view***")
opt = {'title': title_plot, 'tangent': True}
print("\n expected rate of return vs risk representation")
rft = cr1.viewFrontiers(options=opt)
print("\n Sharpe vs expected rate of return representation")
_ = cr1.viewFrontiers(data=rft, fig_type='Sharpe_RR', options=opt)
print("\n diverification factor vs expected rate of return")
_ = cr1.viewFrontiers(data=rft, fig_type='Diverse_RR', options=opt)

#==============================================================================
print("\n******************************************************************\n")
print("*** Frontiers evaluations - custom view***")
# 10 (random in this example) additional portfolios to be added to the plot
rng = np.random.RandomState(42)
addp = {}
for i in range(10):
    addp['p' + str(i+1)] = rng.dirichlet([0.5] * len(symb))
addport = pd.DataFrame().from_dict(addp, 'index', columns=symb)

opt = {'tangent': True, 'title': title_plot, 'minrisk_label': 'mRx', 
       'sharpe_label': 'sharpe', 'addport_label': True, 'xlabel': "RofR"}
print("\n expected rate of return vs risk representation")
fd1 = cr1.viewFrontiers(minrisk=True, efficient=20, inefficient=20, 
                        maxdiverse=True, 
                        diverse_efficient=20, diverse_inefficient=20,
                        invNdiverse=True, invNdrr=True,
                        randomport=10,
                        options=opt, addport=addport)
print("\n Sharpe vs expected rate of return representation")
_ = cr1.viewFrontiers(fig_type='Sharpe_RR', 
                      invNdiverse_label=None, data=fd1, options=opt)
print("\n diverification factor vs expected rate of return")
_ = cr1.viewFrontiers(fig_type='Diverse_RR', 
                      invNrisk_label=None, data=fd1, options=opt)

#==============================================================================
print("\n******************************************************************\n")
print("*** Example of rebalancing positions for a Sharpe strategy ***")
# set Sharpe strategy
rtype = 'Sharpe' 
mu0 = 0. # 0. risk free rate (default value)
ww = cr1.getWeights(rtype, mu0=mu0, verbose=True)

# assumed existing positions and cash
ns = pd.Series(100, index=symb)
cash = 0.

# new positions and rolling info
# optimization strategie
rtype = 'Sharpe'
mu0 = 0. # risk free rate

pos = cr1.getPositions(nshares=ns, cash=cash)
print(f" New position report\n {pos}")

#==============================================================================
print("\n******************************************************************\n")
print("*** Speed comparisons for different methods ***")
# may take some time to complete
# to run please uncomment the lines below
# methods = cr1.methods
# # remove 'interior_point' if exists - it is painfully slow
# if 'interior-point' in methods:
#     methods.remove('interior-point')
# rtypes = cr1.rtypes
# mu = 0.04
# mu0 = 0.
# aversion = 0.6
# ewp = np.full(len(symb), 1/len(symb))

# res_time = pd.DataFrame(0., index=rtypes, columns=methods)
# res_RR = pd.DataFrame(0., index=rtypes, columns=methods)

# for method_ in methods:
#     for rtype_ in rtypes:
#         cr1.set_method(method_)
#         _ = cr1.getWeights(rtype_, mu=mu, mu0=mu0, aversion=aversion, ww0=ewp)
#         print(f"method {method_} rtype {rtype_} status {cr1.status}")
#         res_time.loc[rtype_, method_] = \
#             cr1.time_level1 if cr1.status == 0 else np.nan
#         res_RR.loc[rtype_, method_] = cr1.RR if cr1.status == 0 else np.nan

# print(f"\nComputation time (s) per method per rtype\n{res_time.round(6)}\n")
# print(f"expected rate of return\n{res_RR.round(4)}\n")
# #restore the initial method
# cr1.set_method(method)

#==============================================================================

