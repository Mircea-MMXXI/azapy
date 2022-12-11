import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
# set the end date of the time series
# default start date is sdate='2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, edate=edate)

# create the MVAnalyzer class calculator
calc = az.MVAnalyzer(mktdata)

# plot expected rate of return vs MV
# set some plot features
opt = {'title': "MV Port",
       'xlabel': "VAR",
       'tangent': True}
fig1 = '../figs/MVFrontier1.png'
rft = calc.viewFrontiers(musharpe=0, randomport=100,
                         options=opt, saveto=fig1)

# plot MV-Sharpe vs expected rate of return
# set some plot features
opt2 = {'title': "MV Port",
        'ylabel': "MV-Sharpe"}
fig2 = '../figs/MVFrontier2.png'
# reuse the previous plot computations (data=rft)
_ = calc.viewFrontiers(data=rft, fig_type='Sharpe_RR',
                       options=opt2, saveto=fig2)
