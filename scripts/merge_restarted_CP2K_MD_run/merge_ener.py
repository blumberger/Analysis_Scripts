import pandas as pd
import numpy as np
import sys
import subprocess

_, merge_fold, curr_fold = sys.argv

def read_ener_csv(fp):
    df = pd.read_csv(fp, names=("Step", "Time", "Kin", "Temp", "Pot", "Etot", "CPU"), skiprows=[0], delim_whitespace=True)
    return df

df = read_ener_csv(f"{merge_fold}/run-1.ener")
last_line=subprocess.check_output(f"tail -1 {curr_fold}/run-1.ener", shell=True).decode("utf-8")
last_step, last_time  = last_line.split()[:2]
last_step = int(last_step)
last_time = float(last_time)

df['Step'] += last_step
df['Time'] += last_time
np_df = df.to_numpy()
np.savetxt("new_file.tmp", np_df[1:], fmt=('%10i', '%19.6f', '%19.9f', '%19.9f',
                                       '%19.9f', '%19.9f', '%19.9f'))
