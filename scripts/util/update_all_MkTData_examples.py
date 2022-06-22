import azapy as az

mktdir = "../../MkTdata_test"

ercode = az.update_all_MkTData(mktdir)

# error code 200 is OK (see docs)
print(f"error code: {ercode}")

