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


namespace fs = std::filesystem;

bool do_NN = false;
double NN_cutoff = -1;


int main() {

    /* 
        Init
    */
    // Parse any input
    inp::INP_Parser INP;
    INP.read("instruct.txt");

    // Set some variables
    std::unordered_map<std::string, std::string> fno = {};
    fno["inact_at"] = INP.get_string("inactive_ats_filename", "");
    fno["MD"] = INP.get_string("MD_system", "");
    fno["act_at"] = INP.get_string("active_ats_filename", "");
    fno["all_at"] = INP.get_string("all_ats_filename", "");
    fno["veloc"] = INP.get_string("new_veloc_filename", "");
    fno["fix_at"] = INP.get_string("fixed_atoms_filename", "");
    fno["decomp"] = INP.get_string("decomp_filename", "");
    fno["AOM"] = INP.get_string("AOM_filename", "");
    fno["mol_nums"] = INP.get_string("mol_nums_filename", "");
    fno["COM"] = INP.get_string("COM_filename", "");
    fno["NN"] = INP.get_string("NN_filename", "");
    auto folder_out = INP.get_string("out_folder", "");

    // Add folderpath
    for (const auto &element : fno) {
        if (element.second != "") {
            std::string last_elm = folder_out.substr(folder_out.size()-1, 1);
            if (last_elm != "/")
                fno[element.first] = folder_out + "/" + element.second;
            else
                fno[element.first] = folder_out + element.second;
        }
    }

    auto folder_in = INP.get_string("in_folder");
    auto fni_pos = INP.get_string("pos_filename");
    auto fni_vel = INP.get_string("vel_filename");
    auto fni_AOM = INP.get_string("AOM_single_mol", "");

    // Slice by COM
    struct inp::slice_instruct slices;
    slices.zmin = (double) INP.get_num("zmin", -9999999.);
    slices.zmax = (double) INP.get_num("zmax", 9999999.);
    slices.ymin = (double) INP.get_num("ymin", -9999999.);
    slices.ymax = (double) INP.get_num("ymax", 9999999.);
    slices.xmin = (double) INP.get_num("xmin", -9999999.);
    slices.xmax = (double) INP.get_num("xmax", 9999999.);

    // System Propeties
    int apm = INP.get_int("atoms_per_molecule", 0);
    if (apm == 0){
        std::cerr << "\n\nAtoms per molecule not set in the instruct.txt file." << std::endl;
        std::cerr << "Please set it using the syntax: `atoms_per_molecule = <int>`\n\n" << std::endl;
        throw -1;
    }

    if (INP.get_num("NN_cutoff", -1) != -1 and fno["NN"] != "") {
        do_NN = true;
        NN_cutoff = (double) INP.get_num("NN_cutoff");
    }

    // Make directory
    if (not fs::exists(folder_out)) {
        if (fs::create_directories(folder_out.c_str()) == false) {
            std::cerr << "Couldn't make directory '" << folder_out << "'." << std::endl;
            throw -1;
        }
    } else {
        fs::remove_all(folder_out);
        fs::create_directories(folder_out);
    }
    std::vector<std::vector<int>> NN_ats;


    /* 
        Read and Parse Files
    */
    // Parse the XYZ file
    std::cout << "Start" << std::endl;
    xyz::XYZ_File Pos;
    Pos.read(folder_in + "/" + fni_pos);
    if (fno["all_at"] != "")
        Pos.write(fno["all_at"]);

    vel::Veloc_File Vel;
    Vel.read(folder_in + "/" + fni_vel);


    /* 
        Perform Operations on Data
    */
    // Slice the data
    std::vector<std::vector<int>> mol_nums;
    auto act_atoms = Pos.slice_by_COM(apm, slices, mol_nums, fno["mol_nums"], fno["COM"]);
    // Get the molecules with COM within a defined distance of the active molecules.
    if (do_NN) {
        auto Active_Pos_File = Pos.index_atoms(act_atoms);
        auto COM_File = Active_Pos_File.get_COMs(apm);
        NN_ats = Pos.get_nearest_COM_to_pos(COM_File, mol_nums[0], apm, NN_cutoff);
        Pos.write(fno["NN"], NN_ats);

        if (fno["MD"] != "") {
            // Populate an array with all the 
            std::vector<int> tmp (NN_ats[0].size() + act_atoms[0].size(), 0);

            for (auto i=0; i<NN_ats[0].size(); i++) {
                tmp[i] = NN_ats[0][i];
            }
            auto count=NN_ats[0].size();
            for (auto i=0; i<act_atoms[0].size(); i++) {
                tmp[count] = act_atoms[0][i];
                count ++;
            }
            std::vector<std::vector<int>> MD_ats(1, tmp);
            Pos.write(fno["MD"], MD_ats);
        }
    }

    /*
        Write Any Data Required.
    */
    // Write the files
    if (fno["all_at"] != "")
        Pos.write(fno["act_at"], act_atoms);
    if (fno["veloc"] != "")
        Vel.write(fno["veloc"], act_atoms[0]);

    // Set the active atoms in an array
    if (fno["fix_at"] != "" or fno["inact_at"] != "") {
        // Create mask of bools: shape (natoms, ); true = active; false = inactive
        std::vector<bool> act_at_mask (Pos.natoms, false);
        for (auto iat : act_atoms[0])
            act_at_mask[iat] = true;

        // Total Atoms - active ones
        int inact_atoms = Pos.natoms - act_atoms[0].size();
        // Write the fixed atoms file.
        if (fno["fix_at"] != "") {
            std::fstream fix_at_file;
            fix_at_file.open(fno["fix_at"], std::ios::out);
            fix_at_file << "&FIXED_ATOMS\n\tLIST ";

            // Write the FIXED_ATOMS file
            if (do_NN) 
                for (auto i=0; i<NN_ats[0].size(); i++)
                    fix_at_file << i+1 << " ";

            else {
                for (auto iat=0; iat<Pos.natoms; iat++) {
                    if (not act_at_mask[iat]) {
                        fix_at_file << iat+1 << " ";
                    }
                }
            }
            fix_at_file << "\n&END FIXED_ATOMS";            
            fix_at_file.close();
        } 

        // If required write the inactive atoms
        std::cout << "\n\n";
        if (fno["inact_at"] != "") {
            std::vector<std::vector<int>> inactive_atoms_arr;
            inactive_atoms_arr.resize(1);
            inactive_atoms_arr[0].resize(inact_atoms);
            int count = 0;
            for (auto iat=0; iat<Pos.natoms; iat++){
                if (not act_at_mask[iat]) {
                    inactive_atoms_arr[0][count] = iat;
                    count++;
                }
            }

            // Finally write the inactive atoms.
            Pos.write(fno["inact_at"], inactive_atoms_arr);

        }

        std::cout << "Num inactive atoms: " << inact_atoms << std::endl;
        std::cout << "Num Active atoms:   " << act_atoms[0].size() << std::endl;
        std::cout << "Num Total atoms:    " << Pos.natoms << std::endl;
    }

    // Write the DECOMP file if requested.
    if (fno["decomp"] != "") {
        if (do_NN) {
            std::fstream decomp_file;
            decomp_file.open(fno["decomp"], std::ios::out);
            decomp_file << "&ENERGY_DECOMP\n\tINDEX_MOL_DECOMP ";
            int first_mol = (int) NN_ats[0].size() / apm;
            int num_mols = (int) act_atoms[0].size() / apm;
            for (auto i=first_mol+1; i<1+num_mols+first_mol; i++)
                decomp_file << i << " ";
            decomp_file << "\n\tNUM_ACTIVE_ATOMS " << act_atoms[0].size();
            decomp_file << "\n&END ENERGY_DECOMP";
            decomp_file.close();
        } else {
            std::cerr << "\n\nNot Yet Implemented: DECOMP file without the nearest neighbour slicing.\n\n";
        }
    }

    // Write the AOM file if required
    if (fno["AOM"] != "") {
        if (fni_AOM == "") {
            std::cerr << std::endl << std::endl;
            std::cerr << "Couldn't write AOM file as I don't know what coeffs to use." << std::endl;
            std::cerr << "Please point me to a file that has the coeffs for a single ";
            std::cerr << "molecule. This should be in the correct CP2K AOM format and";
            std::cerr << "have the correct order of atoms." << std::endl << std::endl;
        } else {
            if (do_NN) {
                std::fstream AOM_in_file, AOM_out_file;
                AOM_in_file.open(fni_AOM, std::ios::in);
                AOM_out_file.open(fno["AOM"], std::ios::out);

                std::string inert_at = "XX   1    0   0.0   0.0";
                for (auto i=0; i<NN_ats[0].size()-1; i++)
                    AOM_out_file << inert_at << "\n";
                AOM_out_file << inert_at;

                // Read the aom file
                std::string single_mol, line;
                while( getline(AOM_in_file, line))
                    single_mol = single_mol + "\n" + line;

                for (auto i=0; i<(int)act_atoms[0].size()/36; i++)
                    AOM_out_file << single_mol;

                AOM_out_file.close();
                AOM_in_file.close();
            } else {
                std::cerr << "\n\nNot Yet Implemented: AOM_COEFF file without the nearest neighbour slicing.\n\n";
            }
        }
    }

    return 0;
}
