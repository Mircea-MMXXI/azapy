# Example of how to call schedule functions
import azapy as az

sdate = '2015-01-01'
edate = '2020-12-31'

sims = az.schedule_simple(sdate=sdate, edate=edate)
print(f"simple schedule\n{sims}")

srol = az.schedule_roll(sdate=sdate, edate=edate)
print(f"rolling history schedule\n{srol}")

soff = az.schedule_offset(sdate=sdate, edate=edate)
print(f"simple schedule with offset start\n{soff}")
