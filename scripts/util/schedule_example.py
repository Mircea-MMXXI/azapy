# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 15:45:43 2021

@author: mirce
"""
import pandas as pd
import azapy as az

sdate = pd.to_datetime('2015-01-01')
edate = pd.to_datetime('2020-12-31')

az.simple_schedule(sdate=sdate, edate=edate)
az.schedule_roll(sdate=sdate, edate=edate)
