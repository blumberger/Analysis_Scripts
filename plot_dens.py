import matplotlib.pyplot as plt

dens.data.append_csvs()
plt.plot(dens.data.collected_csv_data[0]['Step']*1e-3,
         dens.data.collected_csv_data[0]['Density'])

plt.xlabel("Time [ps]", fontsize=27)
plt.ylabel(r"Density [g cm$^{-3}$]", fontsize=27)
plt.savefig("Density.png")
plt.close()
