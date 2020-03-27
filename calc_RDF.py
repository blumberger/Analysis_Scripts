folder = "~/Data/Pentacene_Crystal_Melt/xyz_md/"
filepath = "$folder/snapshot0.dat"


read $filepath lammps_dump as data
set system data to pentacene

calc RDF from data as rdf_data
