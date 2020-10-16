#ifndef XYZ_File_header
#define XYZ_File_header

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

namespace xyz {

	const double unset = -921983741.719823712;
	struct slice_instruct {
		double xmin = unset;
		double xmax = unset;
		double ymin = unset;
		double ymax = unset;
		double zmin = unset;
		double zmax = unset;
	};

	class XYZ_File {
	    private:

	        std::string line;
	        std::vector<std::string> line_col;
	        std::vector<double> line_xyz;

	        // Will open a file
	        void open_read() {
	            // Open file to read
	            file.open(fp, std::ios::in);

	            if (! file.is_open()) {
	                throw "No File Found!";
	            }
	        }

	        // Will get how many
	        void get_num_ats() {
	            // Number of atoms is on the first line of the file:
	            std::string nat_str;
	            getline(file, nat_str);
	            natoms = (int) std::stoi(nat_str);
	        }

	        // Will get how many steps are in the xyz file
	        void get_nstep() {
	            std::string tmp;
	            while (getline(file, tmp)) {
	                nlines++;
	            }
	            file.close();

	            float nstepsF = nlines / ((float) natoms + 2.);
	            int nstepsI = (int) nstepsF;
	            if (nstepsF != nstepsI) {
	                throw "Corrupted XYZ File!";
	            }
	            nsteps = nstepsI;

	            // Re-open the file again
	            open_read();
	        }

	        // Initialise the vectors
	        void init_arrays() {
	            xyz.resize(nsteps);
	            cols.resize(nsteps);
	            for (auto i=0; i<nsteps; i++) {
	                xyz[i].resize(natoms);
	                cols[i].resize(natoms);

	                for (auto j=0; j<natoms; j++) {
	                    xyz[i][j].resize(ndata_cols);
	                    cols[i][j].resize(nmeta_cols);
	                }
	            }

	            line_col.resize(nmeta_cols);
	            line_xyz.resize(ndata_cols);
	        }

	        // Will determine how many data columns there are and how many 
	        //  metadata cols (element name etc...)
	        void get_num_data_cols() {
	            bool is_float[100];
	            int count = 0;
	            double tmp;
	            std::stringstream ss;
	            std::string word;

	            for (auto i=0; i<100; i++) { is_float[i] = false; }

	            // Convert to a stringstream
	            ss << line;

	            // Get the doubles in string
	            count = 0;
	            while (ss >> word) {
	                count++;
	                try {
	                    tmp = std::stod(word);
	                    is_float[count] = true;
	                } catch (std::invalid_argument& e) {
	                    is_float[count] = false;
	                }
	            }
	            nall_cols = count;

	            // Get how many doubles from the end there are
	            ndata_cols=0;
	            for (auto i=nall_cols; i>=0; i--) {
	                if (! is_float[i]) {
	                    break;
	                }
	                ndata_cols++;
	            }
	            nmeta_cols = nall_cols - ndata_cols;
	        }

	        // Will parse a line in the data section
	        void parse_line(int istep, int iatom) {
	            std::stringstream ss;
	            std::string word;
	            auto i=0;

	            ss << line;
	            word = "";
	            for (i=0; i<nmeta_cols; i++) {
	                ss >> cols[istep][iatom][i];
	            }

	            for (i=0; i<ndata_cols; i++) {
	                ss >> xyz[istep][iatom][i];
	            }

	        }

	        // Fill the column and xyz data array
	        void fill_arrays() {
	            auto istep = 0;
	            auto iat = 0;

	            // Skip the first 2 lines
	            getline(file, line);
	            getline(file, line);

	            // Get some metadata used to parse the file
	            getline(file, line);
	            get_num_data_cols();

	            // Allocate the memory for the xyz data and columns
	            init_arrays();
	            parse_line(0, 0);

	            for (iat=1; iat<natoms; iat++) {
	                if (getline(file, line)) {
	                    parse_line(0, iat);
	                } else {
	                    std::cerr << "Declaration of natom at the top of the file is incorrect!" << std::endl;
	                    throw "Corrupted XYZ File";
	                }
	            }

	            for (istep=1; istep<nsteps; istep++) {
	                // Ignore the title lines
	                getline(file, line);
	                getline(file, line);

	                // Loop over atoms
	                for (iat=0; iat<natoms; iat++) {
	                    if (getline(file, line)) {
	                        parse_line(0, iat);
	                    } else {
	                        std::cerr << "Declaration of natom at the top of the file is incorrect!" << std::endl;
	                        throw "Corrupted XYZ File";
	                    }
	                }
	            }

	        }

	    protected:

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

	        // Set a mask to tell which atoms have a particular dimension within a range
	        void set_at_mask(std::vector<std::vector<bool>> &mask, int dim,
	        				 double max=-9999999999, double min=9999999999999) {

	        	double val;
	        	for (auto istep=0; istep<nsteps; istep++) {
	        		for (auto iat=0; iat<natoms; iat++) {
	        			val = xyz[istep][iat][dim];

	        			if ((val > min) and (val < max)) {
	        				mask[istep][iat] = true;
	        			}
	        		}
	        	}
	        }

	    public:
	        std::fstream file;
	        std::string fp="";
	        int natoms=-1;
	        int nlines=-1;
	        int nsteps=-1;
	        int ndata_cols=-1;
	        int nmeta_cols=-1;
	        int nall_cols=-1;
	        std::vector<std::vector<std::vector<double>>> xyz;
	        std::vector<std::vector<std::vector<std::string>>> cols;

	        void parse(std::string line) {

	            std::cout << line << std::endl;
	        }

	        // Will read the xyz file.
	        void read(std::string filepath) {
	            fp = filepath;

	            open_read();

	            // Get the number of atoms and steps
	            get_num_ats();
	            nlines = 1;
	            get_nstep();


	            fill_arrays();

	            file.close();
	        }

	        // Will write the data to a xyz file.
	        void write(std::string filepath, std::vector<std::vector<int>> at_inds={{}}) {
	        	file.open(filepath, std::ios::out);
	        	if (nsteps == -1) {
	        		std::cerr << "Steps not properly set" << std::endl;
	        		throw "Steps not properly set";
	        	}
	        	if (natoms == -1) {
	        		std::cerr << "Atoms not properly set" << std::endl;
	        		throw "Atoms not properly set";
	        	}
	        	if (nmeta_cols == -1) {
	        		std::cerr << "Cols not properly set" << std::endl;
	        		throw "Cols not properly set";
	        	}
	        	if (ndata_cols == -1) {
	        		std::cerr << "xyz data not properly set" << std::endl;
	        		throw "xyz data not properly set";
	        	}

	        	// Write all atoms -standard behaviour
	        	if (at_inds[0].size() == 0)
	        	{
		        	for (auto istep=0; istep<nsteps; istep++) {
		        		file << std::to_string(natoms) << "\n\n";

		        		for (auto iat=0; iat<natoms; iat++) {
		        			for (auto col : cols[istep][iat])
		        				file << col << "\t";

		        			for (auto data : xyz[istep][iat])
		        				file << data << "\t";
		        			file << "\n";
		        		}
		        	}
		        } 

		        // Only write the atoms we want in the at_inds array
		        else {
		        	for (auto istep=0; istep<at_inds.size(); istep++) {
		        		int ats = at_inds[istep].size();
		        		file << std::to_string(ats) << "\n\n";

		        		for (auto iat: at_inds[istep]) {
		        			for (auto col : cols[istep][iat])
		        				file << col << "\t";

		        			for (auto data : xyz[istep][iat])
		        				file << data << "\t";
		        			file << "\n";
		        		}
		        	}
		        }
	        	file.close();
	        }

	        void set_xyz(std::vector<std::vector<std::vector<double>>> &xyz_in)
	        {
	        	xyz = xyz_in;

	        	nsteps = xyz_in.size();
	        	natoms = xyz_in[0].size();
	        	ndata_cols = xyz_in[0][0].size();
	        }

	        void set_cols(std::vector<std::vector<std::vector<std::string>>> &cols_in)
	        {
	        	cols = cols_in;

	        	nsteps = cols_in.size();
	        	natoms = cols_in[0].size();
	        	nmeta_cols = cols_in[0][0].size();
	        }

	        // Will set the xyz and cols data.
	        void set_data(std::vector<std::vector<std::vector<double>>> &xyz_in,
	        			  std::vector<std::vector<std::vector<std::string>>> &cols_in)
	        {
	        	set_xyz(xyz_in);
	        	set_cols(cols_in);
	        }




	        // Will return a 3D array of shape (nstep, nmol, nat_per_mol, ncol)
	        // Will reshape the cols array to group by molecule
	        std::vector<std::vector<std::vector<std::vector<std::string>>>> cols_to_mols(int natom_per_molecule) {
	            // Create the array
	            auto nmol = get_nmol_from_atoms_atoms_per_mol(natoms, natom_per_molecule);
	            auto istep=0;
	            auto imol=0;
	            auto iapm=0;
	            auto iat=0;
	            auto idim=0;

	            std::vector<std::vector<std::vector<std::vector<std::string>>>> mol_data;
	            mol_data.resize(nsteps);
	            for (istep=0; istep<nsteps; istep++) {
	                mol_data[istep].resize(nmol);

	                for (imol=0; imol<nmol; imol++) {
	                    mol_data[istep][imol].resize(natom_per_molecule);

	                    for (iapm=0; iapm<natom_per_molecule; iapm++) {
	                        mol_data[istep][imol][iapm].resize(nmeta_cols);
	                        
	                        iat = (imol * natom_per_molecule) + iapm;
	                        for (idim=0; idim<nmeta_cols; idim++) {
	                            mol_data[istep][imol][iapm][idim] = cols[istep][iat][idim];
	                        }
	                    }
	                }
	            }

	            return mol_data;
	        }

	        // Will return a 4D array of shape (nstep, nmol, nat_per_mol, ndim)
	        // Will reshape the xyz array to group by molecule
	        std::vector<std::vector<std::vector<std::vector<double>>>> atoms_to_mols(int natom_per_molecule) {
	            // Create the array
	            auto nmol = get_nmol_from_atoms_atoms_per_mol(natoms, natom_per_molecule);
	            auto istep=0;
	            auto imol=0;
	            auto iapm=0;
	            auto iat=0;
	            auto idim=0;

	            std::vector<std::vector<std::vector<std::vector<double>>>> mol_data;
	            mol_data.resize(nsteps);
	            for (istep=0; istep<nsteps; istep++) {
	                mol_data[istep].resize(nmol);

	                for (imol=0; imol<nmol; imol++) {
	                    mol_data[istep][imol].resize(natom_per_molecule);

	                    for (iapm=0; iapm<natom_per_molecule; iapm++) {
	                        mol_data[istep][imol][iapm].resize(ndata_cols);
	                        
	                        iat = (imol * natom_per_molecule) + iapm;
	                        for (idim=0; idim<ndata_cols; idim++) {
	                            mol_data[istep][imol][iapm][idim] = xyz[istep][iat][idim];
	                        }
	                    }
	                }
	            }

	            return mol_data;
	        }

	        // Will get the columns for just 1 molecule of shape (nat_per_mol, ndim)
	        std::vector<std::vector<std::string>> get_1mol_col(int atoms_per_molecule, int imol=0, int istep=0) {
	            auto nmol = get_nmol_from_atoms_atoms_per_mol(natoms, atoms_per_molecule);
	            std::vector<std::vector<std::string>> mol_col;
	            auto iapm=0;
	            auto idim=0;
	            auto iat=0;

	            // Grab the first mol
	            mol_col.resize(atoms_per_molecule);
	            for (iapm=0; iapm<atoms_per_molecule; iapm++) {
	                mol_col[iapm].resize(nmeta_cols);
	                iat = (atoms_per_molecule * imol) + iapm;

	                for (idim=0; idim<nmeta_cols; idim++) {
	                    mol_col[iapm][idim] = cols[istep][iat][idim];
	                }
	            }

	            return mol_col;
	        }

	        // Will get which indices are within constrainst defined by slice_instruct
	        std::vector<std::vector<int>> get_inds_of_slice(struct slice_instruct slice_instruct) {
	        	std::vector<std::vector<int>> at_inds;
	        	double xmin, xmax, ymin, ymax, zmin, zmax;
	        	xmin = -9999999999; xmax = 9999999999;
	        	ymin = -9999999999; ymax = 9999999999;
	        	zmin = -9999999999; zmax = 9999999999;

	        	if (slice_instruct.xmin != unset) xmin = slice_instruct.xmin;
	        	if (slice_instruct.xmax != unset) xmax = slice_instruct.xmax;
	        	if (slice_instruct.ymin != unset) ymin = slice_instruct.ymin;
	        	if (slice_instruct.ymax != unset) ymax = slice_instruct.ymax;
	        	if (slice_instruct.zmin != unset) zmin = slice_instruct.zmin;
	        	if (slice_instruct.zmax != unset) zmax = slice_instruct.zmax;

	        	std::vector<bool> false_vec (natoms, false);
	        	std::vector<std::vector<bool>> mask (nsteps, false_vec);
				double valx, valy, valz;
	        	for (auto istep=0; istep<nsteps; istep++) {
	        		for (auto iat=0; iat<natoms; iat++) {
	        			valx = xyz[istep][iat][0];
	        			valy = xyz[istep][iat][1];
	        			valz = xyz[istep][iat][2];

	        			if (valx < xmax and valx > xmin 
	        		     	and valy < ymax and valy > ymin
	        		     	and valz < zmax and valz > zmin) {
	        				mask[istep][iat] = true;
	        			}
	        		}
	        	}

	        	at_inds.resize(nsteps);
	        	for (auto istep=0; istep<nsteps; istep++) {
	        		int nelm = 0;
	        		for (auto i : mask[istep]) 
	        			nelm = nelm + (int) i;

	        		at_inds[istep].resize(nelm);

	        		int count = 0;
	        		for (auto iat=0; iat<natoms; iat++) {
	        			if (mask[istep][iat]) {
	        				at_inds[istep][count] = iat;
	        				count++;
	        			}
	        		}
	        	}

	        	return at_inds;
	        }
	    
	        // Get the center of masses with shape (nstep, nmol, ndim)
	        //  This will return an XYZ_File.
	        XYZ_File get_COMs(int atoms_per_molecule) {
	            auto istep=0;
	            auto imol=0;
	            auto idim=0;
	            auto iapm=0;
	            auto nmol = get_nmol_from_atoms_atoms_per_mol(natoms, atoms_per_molecule);
	            auto mol_crds = atoms_to_mols(atoms_per_molecule);
	            auto mol_cols = get_1mol_col(atoms_per_molecule);

				// Get the atomic masses in each mol
				auto tot_mass = 0.0;
				auto mass = 0.0;
	   			std::vector<double> at_masses(atoms_per_molecule); 
				auto count = 0;
				for (auto col: mol_cols) {
					if (atomic_masses.find(col[0]) != atomic_masses.end()) {
						mass = atomic_masses[col[0]];
						at_masses[count] = mass;
						count++;	
						tot_mass = tot_mass + mass;
					} else {
						std::cerr << "I don't know the atomic mass of " << col[0] << "!" << std::endl;
						throw "Unkown Atom Abbreviation";
					}
				}

	            std::vector<std::vector<std::vector<double>>> COM;
	            std::vector<std::vector<std::vector<std::string>>> COM_cols;
				double mean[3] = {0.0, 0.0, 0.0};

				COM.resize(nsteps);
				COM_cols.resize(nsteps);
				for (istep=0; istep<nsteps; istep++) {
					COM[istep].resize(nmol);
					COM_cols[istep].resize(nmol);

					for (imol=0; imol<nmol; imol++) {
						COM[istep][imol].resize(ndata_cols);
						COM_cols[istep][imol].resize(1);

						// Set the mean to 0
						for (idim=0; idim<ndata_cols; idim++)
							mean[idim] = 0.0;
						
						// Loop over atoms and find mean
						for (iapm=0; iapm<atoms_per_molecule; iapm++) {
							for (idim=0; idim<ndata_cols; idim++) {
								mean[idim] = mean[idim] + (mol_crds[istep][imol][iapm][idim] * at_masses[iapm]);
							}
						}

						// Divide by mass and set the COM
						for (idim=0; idim<ndata_cols; idim++) {
							COM[istep][imol][idim] = mean[idim] / tot_mass;
						}
						COM_cols[istep][imol][0] = "XX";

					}
				}



				xyz::XYZ_File COM_File;
				COM_File.set_data(COM, COM_cols);

	            return COM_File;
	        }

	        // Will slice an xyz file data by the molecule and return the atom indices of the slice
	        std::vector<std::vector<int>> slice_by_COM(int atoms_per_molecule,
	         										   struct slice_instruct slice_instruct)
	        {
			    // Calculate COMs
			    XYZ_File COM_File = get_COMs(36);

			    // Get the molecular indices within the slice
			    auto mol_inds = COM_File.get_inds_of_slice(slice_instruct);

			    // Convert mol index to atom index
			    std::vector<std::vector<int>> at_inds;
			    at_inds.resize(mol_inds.size());

			    for (int istep=0; istep<mol_inds.size(); istep++) {
			    	int new_size = mol_inds[istep].size()*atoms_per_molecule;
			    	int count = 0;
			    	at_inds[istep].resize(new_size);


			    	for (auto imol: mol_inds[istep]) {
			    		for (int iapm=0; iapm<atoms_per_molecule; iapm++) {
			    			int iat = (atoms_per_molecule * imol) + iapm;

			    			at_inds[istep][count] = iat;
			    			count++;
			    		}
			    	}
			    }
			    return at_inds;
	        }
	};
}



#endif