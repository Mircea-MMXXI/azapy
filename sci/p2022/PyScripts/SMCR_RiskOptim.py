import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mixture SMCR measure
# confidence levels
alpha = [0.99, 0.975, 0.95]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)

# create the SMCRAnalyzer class calculator
calc = az.SMCRAnalyzer(alpha, coef, mktdata)

# compute optimal portfolio weights
# rtype='Risk' -> minimization of risk
# mu=0.035 -> portfolio expected rate of return set to 3.5%
ww = calc.getWeights(rtype='Risk', mu=0.035)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mixture SMCR
mSMCR = calc.risk
# SMCR components
SMCR = calc.primary_risk_comp.copy()
# associated SMVaR values
SMVaR = calc.secondary_risk_comp.copy()