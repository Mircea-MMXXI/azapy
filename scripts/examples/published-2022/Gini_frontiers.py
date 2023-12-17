import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
# set the end date of the time series
# default start date is sdate='2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, edate=edate)

# create the GINIAnalyzer class calculator
calc = az.GINIAnalyzer(mktdata)

# plot expected rate of return vs Gini
# set some plot features
opt = {'title': "Gini Port",
       'xlabel': "Gini",
       'tangent': True}
fig1 = '../figs/GINIFrontier1.png'
rft = calc.viewFrontiers(musharpe=0, randomport=100,
                         options=opt, saveto=fig1)

# plot Gini-Sharpe expected rate of return
# set some plot features
opt2 = {'title': "Gini Port",
        'ylabel': "Gini-Sharpe"}
fig2 = '../figs/GINIFrontier2.png'
# reuse the previous plot computations (data=rft)
_ = calc.viewFrontiers(data=rft, fig_type='Sharpe_RR',
                       options=opt2, saveto=fig2)
