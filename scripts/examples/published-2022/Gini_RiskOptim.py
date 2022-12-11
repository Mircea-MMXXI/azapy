import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
mktdata = az.readMkT(symb)

# create the GINIAnalyzer class calculator
calc = az.GINIAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='Risk' -> minimization of risk
# mu=0.02 -> portfolio expected rate of return set to 2%
ww = calc.getWeights(rtype='Risk', mu=0.02)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# Gini dispersion
Gini = calc.risk