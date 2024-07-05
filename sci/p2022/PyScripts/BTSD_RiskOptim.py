import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# define the parameters of mixture BTSD
# threshold levels
alpha = [-0.01, 0, 0.01]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)

# create the BTSDAnalyzer class calculator
calc = az.BTSDAnalyzer(alpha, coef, mktdata)

# compute optimal portfolio weights
# rtype='Risk' -> minimization of risk mBTSD
# mu=0.038 -> portfolio expected rate of return set to 3.8%
ww = calc.getWeights(rtype='Risk', mu=0.038)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio mBTSD
mBTSD = calc.risk
# BTSD components
BTSD = calc.primary_risk_comp.copy()