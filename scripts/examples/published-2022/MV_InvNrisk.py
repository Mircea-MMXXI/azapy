import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
mktdata = az.readMkT(symb)

# create the MVAnalyzer class calculator
calc = az.MVAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='InvNrisk' -> optimal portfolio
#   with same variance as equal weighted portfolio
ww = calc.getWeights(rtype='InvNrisk')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio VAR
VAR = calc.risk