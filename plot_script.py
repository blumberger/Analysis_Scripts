import matplotlib.pyplot as plt
import os


if not os.path.isdir("snapshots"):
   os.makedirs("snapshots")


# Plot the temperature
f, a = plt.subplots(figsize=(12, 6))

temp_data = log_data.data.collected_csv_data[1]

timestep = L_dump.data.metadata['timestep']

temp_data = temp_data[temp_data['Step'] <= timestep]

a.plot(temp_data['Step']*1e-3, temp_data['Temp'], 'r--')
a.set_ylabel("Temperature [K]", fontsize=24)
a.set_xlabel("Time [ps]", fontsize=24)
plt.tight_layout()


# Plot the angular distributions
f2, a2 = plt.subplots(nrows=2, figsize=(6,6))
ang_dist.data.plot(axes=a2, label=f"Time: {timestep*1e-3} ps")
plt.tight_layout()


count = str(count).zfill(5)


# Plot the RDF
f3, a3 = plt.subplots(figsize=(6, 6))
a3.plot(rdf_data.data.radii, rdf_data.data.RDF)
a3.set_ylabel("RDF", fontsize=24)
a3.set_xlabel(r"R [$\AA$]", fontsize=24)
plt.tight_layout()


f.savefig(f"snapshots/temp_frame_{count}.png")
f2.savefig(f"snapshots/ang_dist_frame_{count}.png")
f3.savefig(f"snapshots/RDF_frame_{count}.png")
plt.close("all")

