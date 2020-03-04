import matplotlib.pyplot as plt
import os


if not os.path.isdir("snapshots"):
   os.makedirs("snapshots")

f, a = plt.subplots()

temp_data = log_data.data.collected_csv_data[1]

timestep = L_dump.data.metadata['timestep']

temp_data = temp_data[temp_data['Step'] <= timestep]

a.plot(temp_data['Step']*1e-3, temp_data['Temp'], 'r--')
a.set_ylabel("Temperature [K]", fontsize=27)
a.set_xlabel("Time [ps]", fontsize=27)


f2, a2 = plt.subplots(nrows=2)
ang_dist.data.plot(axes=a2, label=f"Time: {timestep*1e-3} ps")

count = str(count).zfill(5)

f.savefig(f"snapshots/temp_frame_{count}.png")
f2.savefig(f"snapshots/ang_dist_frame_{count}.png")

plt.close("all")

