import azapy as az

mktdir = "../work/MkTdata"
#mktdir = '../work/empty'
        
sos = az.update_all_mktdata(mktdir, verbose=True)
    