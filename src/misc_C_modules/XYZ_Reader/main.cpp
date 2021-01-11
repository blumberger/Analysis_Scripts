#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <unordered_map>

#include "XYZ_File.hpp"
#include "Veloc_File.hpp"
#include "input_file.hpp"
#include "CP2K_files.hpp"


namespace fs = std::filesystem;

// Will write the positions, xyz files, for the clustered atoms
void write_cluster_pos(std::string fp, xyz::XYZ_File &FullXYZ, std::vector<std::set<int>> const &cluster_inds, int apm);
// Will write the velocities, VEL.init files, for the clustered atoms
void write_cluster_vel(std::string fp, vel::Veloc_File &FullVel, std::vector<std::set<int>> const &cluster_inds, int apm);
// Will write a fixed atoms file -letting CP2K know which atoms should be fixed in place
void write_fixed_at_file(std::string fp, std::vector<int> fixed_ats, const int fixed_at_list_lim);


int main() {

    /* 
        Init
    */
    // Parse any input
    inp::INP_Parser INP;
    INP.read("instruct.txt");

    auto folder_out = INP.get_string("out_folder", "");
    auto folder_in = INP.get_string("in_folder");

    // Make directory
	std::cerr << "Folder out: " << folder_out << std::endl;
	if (folder_out == "") {
		std::cerr << "No folder to output files to..." << std::endl;
		exit(1);
	}

    if (not fs::exists(folder_out)) {
        if (fs::create_directories(folder_out.c_str()) == false) {
            std::cerr << "Couldn't make directory '" << folder_out << "'." << std::endl;
            exit(1);
        }
    } else {
        fs::remove_all(folder_out);
        fs::create_directories(folder_out);
    }

    fs::copy("instruct.txt", folder_out + "/instruct.txt");




    xyz::XYZ_File Pos;
	xyz::XYZ_File OrigPos;
    vel::Veloc_File Vel;
    vel::Veloc_File OrigVel;

	std::vector<int> files_to_read;
	if (INP.get_string("pos_filename", "") != "" or 
		INP.data_deps.find(inp::POS_XYZ) != INP.data_deps.end())
	{
		files_to_read.push_back(inp::POS_XYZ);
	}
	if (INP.get_string("vel_filename", "") != "" or 
		INP.data_deps.find(inp::VEL_INP) != INP.data_deps.end())
	{
		files_to_read.push_back(inp::VEL_INP);
	}
	if (INP.get_string("AOM_single_mol", "") != "" or
		INP.data_deps.find(inp::AOM_INP) != INP.data_deps.end()) 
	{
		files_to_read.push_back(inp::AOM_INP);
	}

	/*
	 * Read any data files we need to for the methods' dependencies
	*/
	std::cerr << "\nReading any required data...";
	for (auto dep: files_to_read) {
		switch (dep) {
			// Read pos
			case inp::POS_XYZ: {
				auto fni_pos = INP.get_string("pos_filename");
				std::cerr << "\n\t* Postions: '"<<folder_in+"/"+fni_pos<<"'";
				Pos.read(folder_in + "/" + fni_pos);
				break;
			}

			// Read Vel
			case inp::VEL_INP: {
				auto fni_vel = INP.get_string("vel_filename");
				std::cerr << "\n\t* Velocities: '"<<folder_in+"/"+fni_vel<<"'";
				Vel.read(folder_in + "/" + fni_vel);
				break;
			}

			// Read AOM
			case inp::AOM_INP: {
			    auto fni_AOM = INP.get_string("AOM_single_mol", "");
				break;
			}

			default:
				std::cerr << "\n\nDon't understand a data dependency in the input file. Please check it." << std::endl;
				exit(1);
		}
	} std::cerr << "\nFinished reading input data\n" << std::endl;


	/*
	 * Define some global variables
	 */
	// Mol properties
	int apm = INP.get_int("atoms_per_molecule");
	int nmol = xyz::get_nmol_from_atoms_atoms_per_mol(Pos.natoms, apm);
	// Cluster array
	std::vector<std::set<int>> clusters;
	// Mol nums to track which atoms we're working with
	std::vector<int> step1mol_nums(nmol, 0);
	for (int i=0; i<nmol; i++) step1mol_nums[i] = i;
	std::vector<std::vector<int>> current_mol_nums(Pos.nsteps, step1mol_nums);
	// Original file data.
	OrigPos.set_data(Pos.xyz, Pos.cols);
	OrigVel = Vel.copy();


	/*
	 * Carry out the methods
	*/ 
    // Loop over the methods requested in the input file and complete them in sequence
	for (auto method: INP.methods) {
		auto method_name = inp::method_codes_rev[method];
		std::cerr << "\n\nCarrying out task: '" << method_name << "'" << std::endl;


		// Choose which method to act on the data
		switch (method) {

			// This will slice the position data by the center of mass (COM).
			//  molecules will not be split up -if they aren't already.
			case (inp::SLICE_POS):
			{
				struct inp::slice_instruct slices;

				// Init slices dict
				slices.zmin = (double) INP.get_num("zmin", -99999);
				slices.zmax = (double) INP.get_num("zmax",  99999);
	    		slices.ymin = (double) INP.get_num("ymin", -99999);
				slices.ymax = (double) INP.get_num("ymax",  99999);
				slices.xmin = (double) INP.get_num("xmin", -99999);
				slices.xmax = (double) INP.get_num("xmax",  99999);

				// Find indices of slices
				std::vector<std::vector<int>> slice_inds;
			    auto act_ats = Pos.slice_by_COM(apm, slices, slice_inds, 
					     					   INP.get_string("mol_nums_slice_file", ""),
						     				   INP.get_string("COM_file", ""));

				// Edit the current_mol_num array to contain current mol nums -from slice
				std::vector<int> buffer;	
				for (int i=0; i<current_mol_nums.size(); i++) {
					buffer.resize(current_mol_nums[i].size());
					buffer = current_mol_nums[i];
					current_mol_nums[i].resize(slice_inds[i].size());
					for (int j=0; j<slice_inds[i].size(); j++) {
						current_mol_nums[i][j] = buffer[slice_inds[i][j]]; 
					}
				}
				
				// Actually to the indexing
				Pos = Pos.index_atoms(act_ats);
				Vel = Vel.index_arr(act_ats[0]);

				// Write files
				std::string at_fp = folder_out+"/"+INP.get_string("slice_filepath", "");
				if (at_fp != folder_out+"/") Pos.write(at_fp);
				at_fp = folder_out+"/"+INP.get_string("slice_vel_filepath", "");
				if (at_fp != folder_out+"/") Vel.write(at_fp);
				
				break;
			}

			// This will calculate the center of mass of the molecules.
			//  Will output the centers of mass in the specified file.
			case (inp::CALC_COM):
			{
				Pos = Pos.get_COMs(apm);
				auto out_path = folder_out+"/"+INP.get_string("COM_filepath", "");
				if (out_path != folder_out+"/") { std::cerr << "\tWriting COM positions\n";  Pos.write(out_path);}
				break;
			}

			// This will attempt to cluster positions by using a density based clustering method.
			//  Outputs will be of the original atomic positions with each cluster in a separate xyz file.
			case (inp::CLUSTER_POS):
			{
				// A reduced version of DBSCAN (without boundary values) is the default clustering algorithm
				auto str_clust_type = INP.get_string("cluster_type", "naive");
				int cluster_min_points = INP.get_int("cluster_min_points", 10);
				double cluster_cutoff = INP.get_num("cluster_cutoff");
				
				clusters = Pos.get_clusters(cluster_cutoff, str_clust_type, cluster_min_points);	

				// Save the data
				auto fp = folder_out + "/" + INP.get_string("cluster_filepath", "");
				if (fp != folder_out + "/") { std::cerr << "\tWriting All Cluster Positions\n"; write_cluster_pos(fp, OrigPos, clusters, apm);}
				fp = folder_out + "/" + INP.get_string("cluster_vel_filepath", "");
				if (fp != folder_out + "/")	{ std::cerr << "\tWriting All Cluster Velocities\n"; write_cluster_vel(fp, OrigVel, clusters, apm);}
				
				break;
			}

			// Will set the current positions to the cluster selected in instruct.txt
			case (inp::SELECT_CLUSTER): 
			{
				// Init params and error checking
				auto cluster_ind = INP.get_int("cluster_selection");
				if (clusters.size() == 0) {
					std::cerr << "\n\nPlease first use the 'CALC_CLUSTERS' method before selecting a cluster.\n\n" << std::endl;
					exit(1);
				}

				// Get the cluster mol indices and conver to atomic indices
				std::vector<int> cluster_inds (clusters[cluster_ind].begin(), clusters[cluster_ind].end());
				std::vector<std::vector<int>> at_inds = {xyz::mol_inds_to_at_inds(cluster_inds, apm)};	
				// Update the current mol nums
				for (auto i=0; i<Pos.nsteps; i++) {
					current_mol_nums[i].resize(cluster_inds.size());
					current_mol_nums[i] = cluster_inds;
				}
				// Do the indexing
				Pos = OrigPos.index_atoms(at_inds);
				Vel = OrigVel.index_arr(at_inds[0]);
				// Write files
				auto fp = folder_out+"/"+INP.get_string("selected_cluster_filepath", "");
				if (fp != folder_out+"/") { std::cerr << "\tWriting Selected Cluster Positions\n"; Pos.write(fp); }
				fp = folder_out+"/"+INP.get_string("selected_cluster_vel_filepath", "");
				if (fp != folder_out+"/") { std::cerr << "\tWriting Selected Cluster Velocities\n"; Vel.write(fp); }

				break;
			}

			// Will wrap get the nearest neighbour mols by center of mass
			case (inp::CALC_NN): 
			{
				// Get nearest neighbours
				double NN_cutoff = INP.get_num("NN_cutoff");
				xyz::XYZ_File COM = Pos.get_COMs(apm);
		        auto NN_ats = OrigPos.get_nearest_COM_to_pos(COM, apm,
															 current_mol_nums[0],
															 NN_cutoff);
				/* Write files */
				// Write the nearest neighbours
				auto fp = folder_out+"/"+INP.get_string("NN_wrapped_filepath", "");
				if (fp != folder_out+"/") OrigPos.write(fp, NN_ats);
		        
				// Write the FIXED_ATOMS file
				fp = folder_out+"/"+INP.get_string("fixed_at_filepath", "");
				if (fp != folder_out+"/") {
					std::vector<int> fixed_ats (NN_ats[0].size(), 0);
					for (int i=0; i<NN_ats[0].size(); i++) fixed_ats[i] = i;
					write_fixed_at_file(fp, fixed_ats, 35);
				}

				// Write the full system -pos from the slice after the pos of fixed ats
				fp = folder_out+"/"+INP.get_string("fixed_and_slice_filepath", "");
				auto fp_v = folder_out+"/"+INP.get_string("fixed_and_slice_vel_filepath", "");
				bool write = (fp != folder_out+"/") or (fp_v != folder_out+"/");
				if (write == true) {
					// MD_ats holds the index of nearest neigh atoms and then slice ats
					std::vector<int> MD_ats(NN_ats[0].size()+Pos.natoms, 0);
					for (int i=0; i<NN_ats[0].size(); i++) MD_ats[i] = NN_ats[0][i];
					auto slice_inds = xyz::mol_inds_to_at_inds(current_mol_nums[0], apm);
					for (int i=0; i<slice_inds.size(); i++) MD_ats[i+NN_ats[0].size()] = slice_inds[i];
					std::vector<std::vector<int>> write_inds = {MD_ats};

					// Write the files
					if (fp != folder_out+"/") OrigPos.write(fp, write_inds);
					if (fp_v != folder_out+"/") OrigVel.write(fp_v, write_inds[0]);
				}
				break;
			}

			// Let the user know the method hasn't been coded up yet and exit
			default:
			{
				std::cerr << "\n\n! No Code Written For '"<<method_name<<"' !\n\n"<<std::endl;
				exit(1);
			}
		}
	}

	auto fp = folder_out + "/" + "initial_pos.xyz";
	OrigPos.write(fp);
	std::cerr << "\n\nAll Finished\n\n";

    return 0;
}



/*
 * Will write each cluster of atoms in an individual xyz file.
 *
 * Inputs:
 *		fp <std::string> => The filepath template (including folder) pointing to where we should save the data.
 *		FullXYZ <xyz::XYZ_File> => The xyy data for the full structure.
 *		cluster_inds <std::vector<std::set<int>> => The indices of each cluster.
 *		apm <int> => The number of atoms per molecule
*/
void write_cluster_pos(std::string fp, xyz::XYZ_File &FullXYZ, std::vector<std::set<int>> const &cluster_inds, int apm) {
	if (FullXYZ.empty) return;

	auto newfp = fp;

	int count = 0;
	for (auto cluster: cluster_inds) {
		size_t pos = fp.find("*");

		if (pos != -1) {
			newfp.replace(pos, 1, std::to_string(count));

			std::vector<int> inds(cluster.begin(), cluster.end());
			auto tmp_at_inds = xyz::mol_inds_to_at_inds(inds, apm);
			std::vector<std::vector<int>> at_inds = {tmp_at_inds};
			FullXYZ.write(newfp, at_inds);
			
			newfp = fp;
		} 
		else { 
			std::cerr << "\n\nPlease put a '*' in the `cluster_filepath` where you would like the cluster number to go"; exit(1); 
		}
		count++;
	}
}



/*
 * Will write the velocities in the VEL.init format.
 *
 * Inputs:
 *		fp <std::string> => The filepath template (including folder) pointing to where we should save the data.
 *		FullVel <vel::Veloc_File> => The velocity data for the full structure.
 *		cluster_inds <std::vector<std::set<int>> => The indices of each cluster.
 *		apm <int> => The number of atoms per molecule
*/
void write_cluster_vel(std::string fp, vel::Veloc_File &FullVel, std::vector<std::set<int>> const &cluster_inds, int apm) {
	if (FullVel.empty) return;
	
	int count = 0;
	for (auto cluster: cluster_inds) {
		size_t pos = fp.find("*");
		auto newfp = fp;

		if (pos != -1) {
			newfp.replace(pos, 1, std::to_string(count));

			std::vector<int> inds(cluster.begin(), cluster.end());
			auto at_inds = xyz::mol_inds_to_at_inds(inds, apm);
			
			FullVel.write(newfp, at_inds);
			
			newfp = fp;
		} 
		else { 
			std::cerr << "\n\nPlease put a '*' in the `cluster_filepath` where you would like the cluster number to go"; exit(1); 
		}
		count++;
	}
}


/*
 * Will write the fixed atoms file. This tells CP2K which atoms to feeze in place (probably NN wrapped ones)
 *
 * Inputs:
 *		fp <std::string> => The filepath template (including folder) pointing to where we should save the data.
 *		fixed_ats <std::vector<int>> => The atoms to fix in place
 *		fixed_at_list_lim <int> => The num of items in each list in the fixed at file (a sensible value is atoms per mol)
 */
void write_fixed_at_file(std::string fp, std::vector<int> fixed_ats, const int fixed_at_list_lim) {
    std::fstream file;
    file.open(fp, std::ios::out);
    file << "&FIXED_ATOMS";

    for (auto i=0; i<fixed_ats.size(); i++) {
		if (i % fixed_at_list_lim == 0) file << "\n\tLIST ";
        file << fixed_ats[i] << " ";
    }

    file << "\n&END FIXED_ATOMS";            
    file.close();
}
//
//    /* 
//        Read and Parse Files
//    */
//    // Parse the XYZ file
//    std::cerr << "Start" << std::endl;
//    xyz::XYZ_File Pos;
//    if (fno["all_at"] != "")
//        Pos.write(fno["all_at"]);
//
//    vel::Veloc_File Vel;
//    Vel.read(folder_in + "/" + fni_vel);
//
//
//    /* 
//        Perform Operations on Data
//    */
//    // Slice the data
//    std::vector<std::vector<int>> mol_nums;
//    act_ats = Pos.slice_by_COM(apm, slices, mol_nums, fno["mol_nums"], fno["COM"]);
//    int num_act_ats = act_ats[0].size();
//    int num_tot_ats = Pos.natoms;
//    int num_inact_ats = num_tot_ats - num_act_ats;
//
//    // Get the inactive atom indices
//    // Total Atoms - active ones
//    inact_ats.resize(Pos.nsteps);
//    for (auto istep=0; istep<Pos.nsteps; istep++) {
//        int count = 0; int j = 0;
//        inact_ats[istep].resize(num_inact_ats);
//
//        for (auto i=0; i<Pos.natoms; i++) {
//            if (i == act_ats[istep][count]) {
//                count++;
//            } else {
//                inact_ats[istep][j] = i;
//                j++;
//            }
//        }
//
//    }
//
//    // Get the molecules with COM within a defined distance of the active molecules.
//    if (do_NN) {
//        auto Active_Pos_File = Pos.index_atoms(act_ats);
//        auto COM_File = Active_Pos_File.get_COMs(apm);
//        NN_ats = Pos.get_nearest_COM_to_pos(COM_File, mol_nums[0], apm, NN_cutoff);
//        Pos.write(fno["NN"], NN_ats);
//
//        if (fno["MD"] != "") {
//            // Populate an MD_ats array with all the indices
//            std::vector<int> tmp (NN_ats[0].size() + act_ats[0].size(), 0);
//            fixed_ats.resize(1);
//            fixed_ats[0].resize(NN_ats[0].size());
//            for (auto i=0; i<NN_ats[0].size(); i++) {
//                tmp[i] = NN_ats[0][i];
//                fixed_ats[0][i] = i;
//            }
//            auto count=NN_ats[0].size();
//            for (auto i=0; i<act_ats[0].size(); i++) {
//                tmp[count] = act_ats[0][i];
//                count ++;
//            }
//            MD_ats.resize(1);
//            MD_ats[0] = tmp;
//            Pos.write(fno["MD"], MD_ats);
//
//            std::fstream mol_num_file;
//            mol_num_file.open(fno["mol_nums"], std::ios::out);
//            for (auto istep=0; istep<NN_ats.size(); istep++) {
//                mol_num_file << istep << ":\n";
//                for (auto i=(int) NN_ats[istep].size()/apm; i<MD_ats[0].size()/apm; i++) {
//                    mol_num_file << i << " ";
//                }
//                mol_num_file << "\n";
//            }   
//            mol_num_file.close();
//        }
//    } 
//
//    // Otherwise the MD ats are all atoms and fixed ats are inactive atoms
//    else {
//        MD_ats.resize(Pos.nsteps);
//        for (auto istep=0; istep<Pos.nsteps; istep++) {
//            MD_ats[istep].resize(Pos.natoms);
//            for (auto iat=0; iat<Pos.natoms; iat++) {
//                MD_ats[istep][iat] = iat;
//            }
//        }
//        fixed_ats = inact_ats;
//    }
//
//
//
//    /*
//        Write Any Data Required.
//    */
//    // Write the files
//    if (fno["all_at"] != "")
//        Pos.write(fno["act_at"], act_ats);
//    if (fno["veloc_MD"] != "")
//        Vel.write(fno["veloc_MD"], MD_ats[0]);
//    if (fno["veloc_act"] != "")
//        Vel.write(fno["veloc_act"], act_ats[0]);
//
//    // Set the active atoms in an array
//    if (fno["fix_at"] != "" or fno["inact_at"] != "") {
//        // Create mask of bools: shape (natoms, ); true = active; false = inactive
//        std::vector<bool> act_at_mask (Pos.natoms, false);
//        for (auto iat : act_ats[0])
//            act_at_mask[iat] = true;
//
//
//        // Write the fixed atoms file.
//        if (fno["fix_at"] != "") {
//        } 
//
//        // If required write the inactive atoms
//        std::cerr << "\n\n";
//        if (fno["inact_at"] != "") {
//            std::vector<std::vector<int>> inactive_atoms_arr;
//            inactive_atoms_arr.resize(1);
//            inactive_atoms_arr[0].resize(num_inact_ats);
//            int count = 0;
//            for (auto iat=0; iat<Pos.natoms; iat++){
//                if (not act_at_mask[iat]) {
//                    inactive_atoms_arr[0][count] = iat;
//                    count++;
//                }
//            }
//
//            // Finally write the inactive atoms.
//            Pos.write(fno["inact_at"], inactive_atoms_arr);
//
//        }
//
//        std::cerr << "Num inactive atoms: " << num_inact_ats << std::endl;
//        std::cerr << "Num Active atoms:   " << num_act_ats << std::endl;
//        std::cerr << "Num Total atoms:    " << Pos.natoms << std::endl;
//        std::cerr << "Num Active Mols:    " << num_act_ats / apm << std::endl;
//        if (do_NN)
//            std::cerr << "Num NN atoms:       " << NN_ats[0].size() << std::endl;
//            std::cerr << "Num System Mols:    " << (num_act_ats + NN_ats[0].size()) / apm << std::endl;
//        if (num_inact_ats + num_act_ats != Pos.natoms) {
//            std::cerr << "\n\nERROR: # Inactive Ats != # Active Ats\n\n" << std::endl;
//        }
//    }
//
//    // Write the DECOMP file if requested.
//    if (fno["decomp"] != "") {
//        if (do_NN) {
//            std::fstream decomp_file;
//            decomp_file.open(fno["decomp"], std::ios::out);
//            decomp_file << "&ENERGY_DECOMP\n\tINDEX_MOL_DECOMP ";
//            int first_mol = (int) NN_ats[0].size() / apm;
//            int num_mols = (int) act_ats[0].size() / apm;
//            for (auto i=first_mol+1; i<1+num_mols+first_mol; i++)
//                decomp_file << i << " ";
//            decomp_file << "\n\tNUM_ACTIVE_ATOMS " << act_ats[0].size();
//            decomp_file << "\n&END ENERGY_DECOMP";
//            decomp_file.close();
//        } else {
//            std::cerr << "\n\nNot Yet Implemented: DECOMP file without the nearest neighbour slicing.\n\n";
//        }
//    }
//
//    // Write the AOM file if required
//    if (fno["AOM"] != "") {
//        if (fni_AOM == "") {
//            std::cerr << std::endl << std::endl;
//            std::cerr << "Couldn't write AOM file as I don't know what coeffs to use." << std::endl;
//            std::cerr << "Please point me to a file that has the coeffs for a single ";
//            std::cerr << "molecule. This should be in the correct CP2K AOM format and";
//            std::cerr << "have the correct order of atoms." << std::endl << std::endl;
//        } else {
//            if (do_NN) {
//                std::fstream AOM_in_file, AOM_out_file;
//                AOM_in_file.open(fni_AOM, std::ios::in);
//                AOM_out_file.open(fno["AOM"], std::ios::out);
//
//                std::string inert_at = "XX   1    0   0.0   0.0";
//                for (auto i=0; i<NN_ats[0].size()-1; i++)
//                    AOM_out_file << inert_at << "\n";
//                AOM_out_file << inert_at;
//
//                // Read the aom file
//                std::string single_mol, line;
//                while( getline(AOM_in_file, line))
//                    single_mol = single_mol + "\n" + line;
//
//                for (auto i=0; i<(int) act_ats[0].size()/36; i++)
//                    AOM_out_file << single_mol;
//
//                AOM_out_file.close();
//                AOM_in_file.close();
//            } else {
//                std::cerr << "\n\nNot Yet Implemented: AOM_COEFF file without the nearest neighbour slicing.\n\n";
//            }
//        }
//    }



















//    // Set some variables
//	// fno is a dict holding output filenames
//    std::unordered_map<std::string, std::string> fno = {};
//    fno["inact_at"] = INP.get_string("inactive_ats_filename", "");
//    fno["MD"] = INP.get_string("MD_system", "");
//    fno["act_at"] = INP.get_string("active_ats_filename", "");
//    fno["all_at"] = INP.get_string("all_ats_filename", "");
//    fno["veloc_MD"] = INP.get_string("new_veloc_MD_filename", "");
//    fno["veloc_act"] = INP.get_string("new_veloc_active_filename", "");
//    fno["fix_at"] = INP.get_string("fixed_atoms_filename", "");
//    fno["decomp"] = INP.get_string("decomp_filename", "");
//    fno["AOM"] = INP.get_string("AOM_filename", "");
//    fno["mol_nums"] = INP.get_string("mol_nums_filename", "");
//    fno["COM"] = INP.get_string("COM_filename", "");
//    fno["NN"] = INP.get_string("NN_filename", "");
//
//    // Add folderpath
//    for (const auto &element : fno) {
//        if (element.second != "") {
//            std::string last_elm = folder_out.substr(folder_out.size()-1, 1);
//            if (last_elm != "/")
//                fno[element.first] = folder_out + "/" + element.second;
//            else
//                fno[element.first] = folder_out + element.second;
//        }
//	  }
