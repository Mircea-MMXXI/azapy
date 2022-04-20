# Portfolio fig 
import azapy as az

symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

sdate = "2012-01-01"
edate = "2021-07-27"

mktdir = "../MkTdata"

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)


p5 = az.Port_InvVol(mktdata, pname='Portfolio')    
port5 = p5.set_model()   

saveto = './Portfolio_1.png'
_ = p5.port_view(saveto=saveto)

# saveto = './Portfolio_2.png'
# _ = p5.port_view(fancy=True, saveto=saveto)
