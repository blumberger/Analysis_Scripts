#ifndef VELOC_FILE_HEADER
#define VELOC_FILE_HEADER

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

namespace vel {
	class Veloc_File {
		private:
			std::string line;
			std::vector<double> vec3;

	        // Will open a file
	        void open_read() {
	            // Open file to read
	            file.open(fp, std::ios::in);

	            if (! file.is_open()) {
	            	std::cerr << "Can't Find VELOC File!" << std::endl;
	            	std::cerr << "Filepath: " << fp << "\n" << std::endl;
	                throw "No File Found!";
	            }

	        }

	        // Will get the number of lines in the file
	        void get_natom() {
	        	natoms = 0;
	        	while (getline(file, line)) {
	        		natoms++;
	        	}
	        	file.close();

	        	open_read();
	        }


		protected:
			std::string word;

			void parse_lines() {
				int iat=0;
				while (getline(file, line)) {
					std::stringstream ss;
					ss << line;

					for (int i=0; i<3; i++) {
						ss >> word;
						vel[iat][i] = std::stod(word);
					}
					
					iat++;
				}
			}

			// Will get the number of molecules from the num atoms and atoms per molecule.
	        // This should just be a simple divide but this checks if the 2 integers divide perfectly.
	        int get_nmol_from_atoms_atoms_per_mol(int natom, int natom_per_molecule) {
	            // Get number of molecules
	            float nmolF = natoms / natom_per_molecule;
	            int nmol = (int) nmolF;
	            if (nmol != nmolF) {
	                std::cerr << "Incorrect number of atoms per molecule" << std::endl << std::endl;
	                std::cerr << "\t* Num Atoms: " << natoms << std::endl;
	                std::cerr << "\t* Num Atoms Per Mol: " << natom_per_molecule << std::endl;
	                std::cerr << "\t* Num Mol: " << nmolF << std::endl;
	                throw "Incorrect number of atoms per molecule";
	            }

	            return nmol;
	        }

		public:
			std::fstream file;
			std::string fp="";
			int natoms = -1;
			std::vector<std::vector<double>> vel;

			// Will read a VELOC.init file
			void read(std::string filepath) {
				fp = filepath;
				open_read();
				get_natom();

				vel.resize(natoms);
				for (auto iat=0; iat<natoms; iat++) {
					vel[iat].resize(3);
				}

				parse_lines();

				file.close();
			}

 			// Will return a 4D array of shape (nmol, nat_per_mol, ndim)
	        // Will reshape the xyz array to group by molecule
	        std::vector<std::vector<std::vector<double>>> atoms_to_mols(int natom_per_molecule) {
	            // Create the array
	            auto nmol = get_nmol_from_atoms_atoms_per_mol(natoms, natom_per_molecule);
	            auto imol=0;
	            auto iapm=0;
	            auto iat=0;

	            std::vector<std::vector<std::vector<double>>> mol_data;
                mol_data.resize(nmol);

                for (imol=0; imol<nmol; imol++) {
                    mol_data[imol].resize(natom_per_molecule);

                    for (iapm=0; iapm<natom_per_molecule; iapm++) {
                        mol_data[imol][iapm].resize(3);
                        
                        iat = (imol * natom_per_molecule) + iapm;
                        mol_data[imol][iapm] = vel[iat];
                    }
                }

	            return mol_data;
	        }

	        // Will return the array indexed by an integer array of shape (natom)
	        Veloc_File index_arr(std::vector<int> at_inds) {
	        	std::vector<double> vec3 (3, 0);
	        	std::vector<std::vector<double>> new_vel_arr (at_inds.size(), vec3);

	        	int count = 0;
	        	for (auto i : at_inds)	{
	        		new_vel_arr[count] = vel[i];
	        	}

	        	Veloc_File new_vel;
	        	new_vel.natoms = at_inds.size();
	        	new_vel.vel = new_vel_arr;
	        	return new_vel;
	        }

	        // Will write a VELOC init file
	        void write(std::string fp, std::vector<int> at_inds={}) {
	        	std::fstream out_file;
	        	out_file.open(fp, std::ios::out);

	        	// If we don't pass any at_inds -standard behaviour
	        	if (at_inds.size() == 0) {
		        	for (auto each_vel: vel) {
		        		for (auto dim: each_vel) 
		        			out_file << "     " << dim;
		        		out_file << "\n";
		        	}
		        }

		        // If we pass some at inds only print those
		        else {
		        	for (auto iat=0; iat<at_inds.size(); iat++) {
		        		out_file << "  " << std::to_string(vel[iat][0]) 
		        			 	 << "  " << std::to_string(vel[iat][1]) 
		        			 	 << "  " << std::to_string(vel[iat][2]) << "\n";
		        	}
		        }

	        	out_file.close();
	        }
	};
}

#endif