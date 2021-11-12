# Example of how to call gambilngKelly function

import azapy as az

# 3 simultaneous independent games 
p = [0.55, 0.6, 0.65]

ww = az.gamblingKelly(p)

# bet sizes for each game in percentage
print(f"bet sizes as fraction of capital (in percent)\n{ww}")

# percentage of the total capital invested in each round
print(f"total fraction of capital invested in all games (in precent): {ww.sum()}")
