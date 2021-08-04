# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 15:45:43 2021

@author: mirce
"""
import pandas as pd
import azapy as az

sdate = pd.to_datetime('2015-01-01')
edate = pd.to_datetime('2020-12-31')

sims = az.schedule_simple(sdate=sdate, edate=edate)
print(f"simple schedule\n {sims}")

srol = az.schedule_roll(sdate=sdate, edate=edate)
print(f"roll schedule\n {srol}")
