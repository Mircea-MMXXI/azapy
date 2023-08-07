# Examples
import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# Collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = '2021-07-27'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#==============================================================================
# buid a fixing schedule
fixing_schedule = az.schedule_simple(sdate, edate, 
                                     freq='M', noffset=-5, fixoffset=0)

#==============================================================================
# compute weights 
# using uniformly distributed random vectors in n-simplex
puniv = az.UniversalEngine(mktdata, schedule=fixing_schedule)
ww = puniv.getWeights(mc_paths=100, nr_batches=10, mc_seed=42, verbose=True)
print(f"weights\n{ww}")

# using Flat Dirichlet distribution (equivalent to above uniform distribution)
dirichlet_alpha = [1] * len(symb)
puniv = az.UniversalEngine(mktdata, schedule=fixing_schedule, 
                           dirichlet_alpha=dirichlet_alpha)
ww = puniv.getWeights(mc_paths=100, nr_batches=10, mc_seed=42, verbose=True)
print(f"weights\n{ww}")

# using Dirichlet distribution with all alpha equal to 1/2
dirichlet_alpha = [0.5] * len(symb)
puniv = az.UniversalEngine(mktdata, schedule=fixing_schedule, 
                           dirichlet_alpha=dirichlet_alpha)
ww = puniv.getWeights(mc_paths=100, nr_batches=10, mc_seed=42, verbose=True)
print(f"weights\n{ww}")

# using Dirichlet distribution with all alpha equal to the inverse number of 
# portfolio symbols
dirichlet_alpha = [1 / len(symb)] * len(symb)
puniv = az.UniversalEngine(mktdata, schedule=fixing_schedule, 
                           dirichlet_alpha=dirichlet_alpha)
ww = puniv.getWeights(mc_paths=100, nr_batches=10, mc_seed=42, verbose=True)
print(f"weights\n{ww}")