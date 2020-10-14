#!/usr/bin/env python3from src.calc import general_calc as gen_type

# -*- coding: utf-8 -*-
"""
Holds the dictionaries that point to the functions/classes that carry out loading, writing and calculating functions/classes.
"""

# Calculator functions/classes
from src.calc import pvecs as pvec_lib
from src.calc import NN
from src.calc import angular_distribution as ang_dist
from src.calc import density as dens
from src.calc import RDF as rdf
from src.calc import create_psf as psf_calc
from src.calc import crystallinity
from src.calc import couplings as coupl
from src.calc import layers
from src.calc import rotation as rotate
from src.calc import tmp_calc as tmp

# File handling functions/classes
from src.io_utils import general_io as gen_io
from src.io_utils import CP2K_inp_files as CP2K_inp
from src.io_utils import xyz_files as xyz
from src.io_utils import lammps
from src.io_utils import json_files as json
from src.io_utils import csv_files
from src.io_utils import massif_files as M_files
from src.io_utils import param_files
from src.io_utils import pseudo_hamiltonian as psu_ham
from src.io_utils import numpy_files
from src.io_utils import esp

# Plot handling functions/classes
from src.plot import coupling_networks as coup_net


load_fncs = {
             'cp2k_inp': CP2K_inp.Read_INP, 'xyz': xyz.XYZ_File,
             'json': json.read_json, 'lammps_log': lammps.Lammps_Log_File,
             'txt': gen_io.DataFileStorage, 'lammps_data': lammps.Lammps_Data_File,
             'lammps_dump': lammps.Lammps_Dump, 'massif_file': M_files.Massif_File,
             'params': param_files.Params, "pseudo_ham": psu_ham.Pseudo_Ham,
             "esp": esp.ESP_File, "lammps_input": lammps.Lammps_Input_File,
}
write_fncs = {
              'cp2k_inp': "self.write", 'xyz': xyz.Write_XYZ_File,
              'json': json.write_json, 'csv': csv_files.Write_CSV,
              "psf": gen_io.Write_File, "mol_layers": json.write_json,
              "numpy": numpy_files.Numpy, "lammps_input": lammps.Write_Lammps_Input,
              "txt": gen_io.Write_File,  #"lammps_dump": lammps.Write_Lammps_Dump,
}
calc_fncs = {
             'pvecs': pvec_lib.PVecs, 'NN': NN.NN, 'density': dens.Density,
             'angular_dist': ang_dist.Angular_Dist, 'RDF': rdf.RDF,
             'psf_file': psf_calc.Create_PSF, 'crystallinity': crystallinity.Crystallinity,
             'couplings': coupl.Couplings, "coupling_connections": coupl.Coupling_Connections,
             "long_ax_rotation": rotate.Long_Ax_Rot, "mol_layers": layers.Molecular_Layers,
             "tmp": tmp.TMP, "all_aom_couplings": coupl.Calc_All_AOM_Couplings,
             "et_rates": coupl.Calc_ET_Rates,
}
plot_fncs = {
             'coupling_network': coup_net.Coupling_Networks,
}
