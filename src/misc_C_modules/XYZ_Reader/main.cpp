#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include "Periodic_table.hpp"
#include "XYZ_File.hpp"



std::string folder = "/home/matt/Helford/Documents/Code_bits_and_bobs/MD_Analysis_Scripts";
std::string filename = "unwrapped_slab.xyz";



int main() {

    // Parse the XYZ file
    xyz::XYZ_File Pos;
    Pos.read(folder + "/" + filename);

    // Slice the COM by some dims
    xyz::slice_instruct slice_instruct;
    slice_instruct.zmin = 15;
    slice_instruct.zmax = 29;
    slice_instruct.xmin = -100;
    slice_instruct.xmax = 132;
    slice_instruct.ymin = -100;
    slice_instruct.ymax = 132;
    auto at_inds = Pos.slice_by_COM(36, slice_instruct);
    Pos.write("Slice.xyz", at_inds);

    return 0;
}
