# Portfolio fig 
import pandas as pd

import azapy as az

sdate = pd.to_datetime("2012-01-01")
edate = pd.to_datetime('today')
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "../MkTdata"

mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                     dir=mktdir, force=False) 


p5 = az.Port_InvVol(mktdata, pname='Portfolio')    
port5 = p5.set_model()   

saveto = './Portfolio_1.png'
_ = p5.port_view(saveto=saveto)

# saveto = './Portfolio_2.png'
# _ = p5.port_view(fancy=True, saveto=saveto)
