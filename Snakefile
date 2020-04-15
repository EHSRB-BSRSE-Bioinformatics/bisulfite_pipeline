import re
import os
import pandas as pd
from glob import glob

configfile: "config.json"

SAMPLES_WT = pd.read_csv("metadata.csv").query('desc == "WT"')["SRR"].values
SAMPLES_TKO = pd.read_csv("metadata.csv").query('desc == "TKO"')["SRR"].values

print("WT samples: ")
print(SAMPLES_WT)
print("TKO samples: ")
print(SAMPLES_TKO)

