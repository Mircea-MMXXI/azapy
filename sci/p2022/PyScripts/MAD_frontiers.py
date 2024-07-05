import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
# set the end date of the time series
# default start date is sdate='2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, edate=edate)

# define the parameters of mMAD measure
# coefficients - the normalization is done internally
# implicitly the mMAD level is len(coef)
coef = [1] * 3

# create the MADAnalyzer class calculator
calc = az.MADAnalyzer(coef, mktdata)

# plot expected rate of return vs mMAD
# set some plot features
opt = {'title': "mMAD Port",
       'xlabel': "mMAD : [L=3, equal weighted]",
       'tangent': True}
fig1 = '../figs/mMADFrontier1.png'
rft = calc.viewFrontiers(musharpe=0, randomport=100,
                         options=opt, saveto=fig1)

# plot mMAD-Sharpe vs expected rate of return
# set some plot features
opt2 = {'title': "mMAD Port",
        'ylabel': "mMAD-Sharpe"}
fig2 = '../figs/mMADFrontier2.png'
# reuse the previous plot computations (data=rft)
_ = calc.viewFrontiers(data=rft, fig_type='Sharpe_RR',
                       options=opt2, saveto=fig2)
