import matplotlib.pyplot as plt

plt.plot(rdf_data.data.radii, rdf_data.data.RDF, 'k--')

plt.xlabel(r"R [$\AA$]", fontsize=27)
plt.ylabel(r"RDF", fontsize=27)
plt.show()
