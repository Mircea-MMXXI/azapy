import azapy as az

mktdir = "../../MkTdata"

ercode = az.update_MkTdata(mktdir)

# error code 200 is OK (see docs)
print(f"error code: {ercode}")

