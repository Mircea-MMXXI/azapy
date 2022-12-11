import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# collect market data
# set the end date of the time series
# default start date is sdate='2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, edate=edate)

# define the parameters of the mixture BTSD
# threshold levels
alpha = [-0.01, 0, 0.01]
# coefficients - the normalization is done internally
coef = [1] * len(alpha)


# create the BTSDAnalyzer class calculator
calc = az.BTSDAnalyzer(alpha, coef, mktdata)

# plot expected rate of return vs mBTSD
# set some plot features
opt = {'title': "mBTSD Port",
       'xlabel': "mBTSD : " + str(alpha) + ", equal weigthed",
       'tangent': True}
fig1 = '../figs/mBTSDFrontier1.png'
rft = calc.viewFrontiers(musharpe=0, randomport=100,
                         options=opt, saveto=fig1)
                         
# plot Sortino vs expected rate of return
# set some plot features
opt2 = {'title': "mBTSD Port",
        'ylabel': "Sortino"}
fig2 = '../figs/mBTSDFrontier2.png'
# reuse the previous plot computations (data=rft)
_ = calc.viewFrontiers(data=rft, fig_type='Sharpe_RR',
                       options=opt2, saveto=fig2)
