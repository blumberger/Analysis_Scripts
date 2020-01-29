import RDF
import angle_dist



do_RDF = False
do_angles = True




folders = ["/scratch/mellis/AMORPHOUS_PENTACENE/10ns_new_quenched_NVT_equilibration",
           "/scratch/mellis/AMORPHOUS_PENTACENE/100ps_quenched_NVT_equilibriation",
           "/scratch/mellis/AMORPHOUS_PENTACENE/10ns_quenched_NVT_equilibriation",
           "/scratch/mellis/AMORPHOUS_PENTACENE/1ns_quenched_crashedAfter209ps_NVT_equilibriation",
           "/scratch/mellis/AMORPHOUS_PENTACENE/1ns_quenched_crashedAfter209ps_NVT_equilibriation",
           ]

filenames = ["10ns_quenched_(new).xyz", "100ps_quench.xyz", "10nsNPT_quench.xyz",
             "1ns_quenching_crashed_after_209ps.xyz", "crystal.xyz"]


if do_RDF:
   for xyz_folder, xyz_filename in zip(folders, filenames):
      xyz_filepath = os.path.join(xyz_folder, xyz_filename)
      if not os.path.isfile(xyz_filepath):
          raise SystemExit("Can't find the file: '%s'!" % xyz_filepath)
      
      f, a, title = plot_RDF_from_file(xyz_filepath)
      f.savefig("/homes/mellis/Documents/Graphs/Amorphous_Pentacene/RDFs/%s.jpg" % title)
   

if do_angles:
   for xyz_folder, xyz_filename in zip(folders, filenames):
      xyz_filepath = os.path.join(xyz_folder, xyz_filename)
      if not os.path.isfile(xyz_filepath):
          raise SystemExit("Can't find the file: '%s'!" % xyz_filepath)
   
      f, a, title = angle_dist.plot_angle_dist_from_file(xyz_filepath)
      plt.savefig("/homes/mellis/Documents/Graphs/Amorphous_Pentacene/Angle_Distributions/%s.jpg" % title)

