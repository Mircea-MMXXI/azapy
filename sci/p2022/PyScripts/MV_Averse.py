import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# create the MVAnalyzer class calculator
calc = az.MVAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='RiskAverse' -> optimal with fixed risk-aversion
# aversion=3 -> risk-aversion factor set to 3
ww = calc.getWeights(rtype='RiskAverse', aversion=3)

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio VAR
VAR = calc.risk