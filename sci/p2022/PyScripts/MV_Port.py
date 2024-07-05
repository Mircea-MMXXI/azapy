import azapy as az

# define portfolio components
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# collect market data
# set the start and end dates of the time series
sdate = '2012-01-01'
edate = '2021-07-27'
mktdata = az.readMkT(symb, sdate=sdate, edate=edate)

# create Port_MV calculator (quarterly rebalanced)
# pname="MV Port" - portfolio name
calc = az.Port_MV(mktdata, pname="MV Port")

# set the model parameters:
#    MV-Sharpe optimal portfolio (default rtype='Sharpe')
#    risk free rate set to 0 (default mu0=0)
# returns the portfolio time series (backtesting)
port = calc.set_model()

# plot the portfolio time series + EMA30 and EMA200
fig = '../figs/Port_MV_Sharpe.png'
_ = calc.port_view(saveto=fig)
