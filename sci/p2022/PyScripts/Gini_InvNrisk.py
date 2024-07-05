import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# create the GINIAnalyzer class calculator
calc = az.GINIAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='InvNrisk' -> efficient portfolio with
#       the same Gini as equal weighted portfolio
ww = calc.getWeights(rtype='InvNrisk')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# Gini dispersion
Gini = calc.risk