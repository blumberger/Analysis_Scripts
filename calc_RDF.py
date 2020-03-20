folder = "~/Data/Pentacene_Crystal_Melt/xyz_md/"
filepath = "$folder/snapshot0.dat"
#filepath = "/home/oem/Downloads/RDF/xyz_md_300K/snapshot0.dat"


read $filepath lammps_dump as data
set system data to pentacene
#data['coordinate_wrapping'] = "wrapped"

calc RDF from data as rdf_data
