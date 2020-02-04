import RDF
import angle_dist

import os
import matplotlib.pyplot as plt
import re



do_RDF = False
do_angles = True




folders = ["/scratch/mellis/AMORPHOUS_PENTACENE/1ns_quenched_crashedAfter209ps_NVT_equilibriation",
           "/scratch/mellis/AMORPHOUS_PENTACENE/",
           ]

filenames = ["1ns_quenching_crashed_after_209ps.xyz", "crystal2.xyz"]

def find_timelabel(filename):
   """
   Will find the quenching time from the filename
   """
   time_label = re.findall("[0-9]+[a-z]s", filename)
   if len(time_label) >= 1:
      time_label = time_label[0]
   else:
      time_label = ""
      
   if '10ns' in time_label and 'new' in filename:
       time_label = '10ns_new'

   return time_label


if do_RDF:
   f, a = plt.subplots()
   for xyz_folder, xyz_filename in zip(folders, filenames):
      xyz_filepath = os.path.join(xyz_folder, xyz_filename)
      if not os.path.isfile(xyz_filepath):
          raise SystemExit("Can't find the file: '%s'!" % xyz_filepath)
      
      
      time_label = find_timelabel(xyz_filename)
      if not time_label: time_label = "Crystal"
      else:              time_label += " Amorphous"

      _, _, _ = RDF.plot_RDF_from_file(xyz_filepath, a=a, num_ats_per_mol=36,
                                       nbins=60, legend_lab=time_label)
         
   f.legend()
      
      #f.savefig("/homes/mellis/Documents/Graphs/Amorphous_Pentacene/RDFs/%s.jpg" % title)
   
   a.set_xlim([0, 50])
   plt.show() 

if do_angles:
   f, a = plt.subplots(3)
   for xyz_folder, xyz_filename in zip(folders, filenames):
      xyz_filepath = os.path.join(xyz_folder, xyz_filename)
      if not os.path.isfile(xyz_filepath):
          raise SystemExit("Can't find the file: '%s'!" % xyz_filepath)
   
      time_label = find_timelabel(xyz_filename)
      if not time_label: time_label = "Crystal"
      else:              time_label += " Amorphous"

      angle_dist.plot_angle_dist_from_file(xyz_filepath, a, nbins=70, legend_label=time_label)

   for ax in a:
      ax.set_ylim([0, 12])
   f.legend()
   plt.show()
