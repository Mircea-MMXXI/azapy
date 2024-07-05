import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
mktdata = az.readMkT(symb)

# create the MVAnalyzer class calculator
calc = az.MVAnalyzer(mktdata)

# compute optimal portfolio weights
# rtype='Risk' -> minimization of risk
# mu=0.035 -> portfolio expected rate of return set to 3.5%
ww = calc.getWeights(mu=0.035, rtype='Risk')

# other quantities of interest:
# portfolio expected rate of return
rate = calc.RR
# portfolio VAR
VAR = calc.risk