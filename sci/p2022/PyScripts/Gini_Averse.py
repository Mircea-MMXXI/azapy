import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# create the GINIAnalyzer class calculator
calc = az.GINIAnalyzer(mktdata)

# compute optimal portfolio weights
# aversion=3 -> risk aversion coefficient set to 3
# rtype='RiskAverse' -> optimal for fixed risk averse
ww = calc.getWeights(rtype='RiskAverse', aversion=3)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# Gini dispersion
Gini = calc.risk