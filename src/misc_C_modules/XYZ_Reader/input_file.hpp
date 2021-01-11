#ifndef INP_PARSER_HEADER
#define INP_PARSER_HEADER

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <set>
#include <vector>

std::unordered_map <int, std::string> reverse_str2int (std::unordered_map <std::string, int> inDict) {
	std::unordered_map <int, std::string> outDict;
	for (auto key_val: inDict) {
		outDict[key_val.second] = key_val.first;
	}
	return outDict;
}


namespace inp {

	enum var_type {typeINT, typeFLOAT, typeSTR, typeARR};
	enum methods {SLICE_POS, CALC_COM, CLUSTER_POS, SELECT_CLUSTER, CALC_NN};
	enum data_files {POS_XYZ, VEL_INP, AOM_INP};


	std::unordered_map <int, std::vector<int>> data_dep_map = {
		{CALC_NN, {POS_XYZ}},
		{SELECT_CLUSTER, {}},
		{CLUSTER_POS, {POS_XYZ}},
		{CALC_COM, {POS_XYZ}},
		{SLICE_POS, {POS_XYZ}}
	};
	std::unordered_map <std::string, int> method_codes = {
		{"WRAP_NEAREST_NEIGHBOURS", CALC_NN},
		{"SELECT_CLUSTER", SELECT_CLUSTER},
		{"CALC_CLUSTERS", CLUSTER_POS},
		{"CALC_COM", CALC_COM},
		{"SLICE_POS", SLICE_POS}
	};
	std::unordered_map <int, std::string> method_codes_rev = reverse_str2int(method_codes);

	const double unset = -921983741.71982323423712;
	struct slice_instruct {
		double xmin = unset;
		double xmax = unset;
		double ymin = unset;
		double ymax = unset;
		double zmin = unset;
		double zmax = unset;
	};

	class INP_Parser {

		private:

		protected:

			void parse() {
				std::string line;
				std::cerr << "|------------------------" << std::endl;
				std::cerr << "|Input File Type Rundown:";
				while (getline(file, line)) {
					if (line.empty())  continue;

					parse_line(line);
				}
				std::cerr << "\n|------------------------\n" << std::endl;
				empty = false;
			}

			std::vector<std::string> split_string(std::string const &str, std::string const &delim) {
				std::stringstream ss(str);

				// Get num of occurances of delimeter
				int n = 0, len_delim = delim.size();
				for (auto i=0; i<=str.size()-len_delim; i++)
					if (str.substr(i, len_delim) == delim) n++;
				n++;

				// Populate the vector with the tokens
				std::vector<std::string> splitter(n, "") ;
				const char delim_cstr = *delim.c_str();
				for (auto i=0; i<n; i++)
					std::getline(ss, splitter[i], delim_cstr);

				return splitter;
			}

			inline void parse_line(std::string &line) {
				std::string first_elm;

				// Ignore comment lines
				if (line.find("#") != std::string::npos) return;

				// Tokenise into vector split by '='
				auto splitter = split_string(line, "=");

				// Ignore lines with no equals sign
				if (splitter.size() < 2) return;

				// Set the data and var_name variables
				std::string var_name = splitter[0], space = " ";
				std::string last_elm = var_name.substr(var_name.size()-1, 1);
				if (last_elm == space) var_name = var_name.substr(0, var_name.size()-1);

				std::string data = splitter[1];
				for (auto i=2; i<splitter.size(); i++)
					data = data + "=" + splitter[i];

				auto type = get_type(data);
				if (type == typeFLOAT) {
					doubleVals[var_name] = std::stod(data);	
					std::cerr << "\n|\t* " << doubleVals[var_name] << "\t\t -> float";
				}
				else if (type == typeINT) {
					intVals[var_name] = std::stoi(data); 
					std::cerr << "\n|\t* " << intVals[var_name] << "\t\t -> int";
				}
				else if (type == typeSTR) {
					first_elm = data.substr(0, 1);
					if (first_elm == space) data = data.substr(1, data.size()-1);
					stringVals[var_name] = data; 
					std::cerr << "\n|\t* " << stringVals[var_name] << " -> string";
				}
				else if (type == typeARR and var_name == "methods") {
					std::cerr << "\n\t* "<< data <<" methods array";
					parse_methods(data);				}
			}

			/*
			 * Will parse an array from the input file and add values to the dicts
			 *
			 * This will only parse non-nested arrays.
			*/
			void parse_methods(std::string str) {
				int startInd=0; int endInd=0;
				for (auto i=0; i<str.size(); i++) {
					auto c = str.substr(i, 1);
					if (c == "(") startInd = i;
					if (c == ")") endInd = i;
				}
				
				std::vector<std::string> str_methods = split_string(str.substr(startInd+1,
																			   endInd-startInd-1),
																    ",");

				for (auto &method_name: str_methods) {
					// Clean up the string
					method_name.erase(std::remove(method_name.begin(), method_name.end(), ' '),
												  method_name.end());
					std::transform(method_name.begin(), method_name.end(), method_name.begin(),
								   ::toupper);

					// Add the method and any data dependencies
					if (method_codes.find(method_name) != method_codes.end()) {
						auto method = method_codes[method_name];
						for (auto i: data_dep_map[method]) {
							data_deps.insert(i);
						}
						methods.push_back(method);

	                } else {
						std::cerr << "\n\n\nUnknown method: '"<<method_name<<"'"<<std::endl;
						std::cerr << "\n\n\n\n" << std::endl;
						exit(1);
					}
				}
			}

			inline int get_type(std::string str) {
				try {
					double tmpD = std::stod(str);
					int tmpI = (int) tmpD;

					// Check if is int
					if (tmpI == tmpD and str.find(".") == std::string::npos) {
						auto test_str = std::to_string(tmpI);
						str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
						if (str == test_str) return typeINT;
						else				  return typeSTR;
					}

					// Check if is float
					else {
						auto test_str = std::to_string(tmpD);
						str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
						str.erase(std::remove(str.begin(), str.end(), '0'), str.end());
						test_str.erase(std::remove(test_str.begin(), test_str.end(), '0'), test_str.end());

						if (test_str != str) return typeSTR;					
						return typeFLOAT;
					}
				} 

				catch (const std::invalid_argument& e) {
					auto Ostr = str;
					str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
					if (str.find("(") == 0) {
						if (str.find(")") == str.length()-1) {
							return typeARR;
						}
					}
					return typeSTR;
				}
			}


		public:
			std::string fp;
			std::fstream file;

			std::unordered_map<std::string, std::string> stringVals {};
			std::unordered_map<std::string, double> doubleVals {};
			std::unordered_map<std::string, int> intVals {};

			std::vector<int> methods={};
			std::set<int> data_deps;

			bool empty = true;


			void read(std::string filepath) {
				fp = filepath;
				file.open(fp, std::ios::in);

				std::cerr << "Opened " << fp << std::endl;

				parse();

				file.close();
			}

			// Will get a value if it exists if not throw and error
			std::string get_string(std::string var_name, std::string default_=")(THROW_ERROR)(*&$S}}~\"") {
				if (stringVals.find(var_name) != stringVals.end()) {
					return stringVals[var_name];
				}

				else {
					if (default_ == ")(THROW_ERROR)(*&$S}}~\"") {
						std::cerr << "\n\n\n\nNo string named '" << var_name << "' in the input file\n\n" << std::endl;
						exit(1);
					} else {
						return default_;
					}
				}
			}

			// Will get a value if it exists if not throw and error
			double get_double(std::string var_name, double default_= unset) {
				std::cerr << var_name << " \n"; // << doubleVals.find(var_name) << " " << doubleVals.end() << std::endl;
				if (doubleVals.find(var_name) != doubleVals.end()) {
					return doubleVals[var_name];
				}

				else {
					if (default_ == unset) {
						std::cerr << "\n\n\n\nNo double named '" << var_name << "' in  the input file\n\n" << std::endl;
						exit(1);
					} else {
						return default_;
					}
				}
			}

			// Will get a value if it exists if not throw and error
			double get_int(std::string var_name, int default_= (int) unset) {
				if (intVals.find(var_name) != intVals.end()) {
					return intVals[var_name];
				}

				else {
					if (default_ == (int) unset) {
						std::cerr << "\n\n\n\nNo int named '" << var_name << "'\n\n" << std::endl;
						exit(1);
					} else {
						return default_;
					}
				}
			}

			// Will get a value if it exists if not throw and error
			double get_num(std::string var_name, double default_= unset) {
				if (doubleVals.find(var_name) != doubleVals.end())
					return doubleVals[var_name];

				else if (intVals.find(var_name) != intVals.end())
					return (double) intVals[var_name];

				else {
					if (default_ == unset) {
						std::cerr << "\n\n\n\nNo number named '" << var_name << "'\n\n" << std::endl;
						exit(1);
					} else {
						return default_;
					}
				}
			}

			// // A work in progress to parse a general type from input file.
			// // Just for me to learn about Template types...
			// template <typename T>
			// inline void parse_line_util(std::stringstream &ss, T &variable) {
			// 	bool comment_break = true;
			// 	bool last_loop = false;
			// 	bool tripped_up = false;
			// 	std::string word;

			// 	// Iterate over the words in the line and decide what to do
			// 	while (ss >> word) {
			// 		// Ignore the equals that sets the value
			// 		if (word == "=") continue;

			// 		// If there is a comment break here.
			// 		auto comm_ind = word.find("#");
			// 		if (comment_break and comm_ind != std::string::npos) {
			// 			if (comm_ind == 0) break;

			// 			word = word.substr(0, comm_ind);
			// 			last_loop = true;
			// 		}

			// 		// Do the parsing
			// 		std::string temp_type = typeid(variable).name();
			// 		std::cerr << word << ",  " << temp_type << std::endl;

			// 		// Parse double
			// 		if (temp_type == "d") {
			// 			try {
			// 				variable = std::stod(word); break;
			// 			} catch (const std::invalid_argument& e) {
			// 				tripped_up = true; break; }
			// 		} 

			// 		// Parse int
			// 		else if (temp_type == "i") {
			// 			try {
			// 				variable = std::stoi(word); break;
			// 			} catch (const std::invalid_argument& e) {
			// 				tripped_up = true; break;  }
			// 		}

			// 		// Parse string
			// 		else if (temp_type == "NSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE") {
			// 			// variable = word;
			// 		}

			// 		// Unknown type
			// 		else {
			// 			std::cerr << "\n\n\n\nI don't recognise type '" << temp_type << "'.";
			// 			std::cerr << "\n\n Change the `parse_line_util` function in the `input_function.hpp` to change this.";
			// 			std::cerr << s\n\ntd::endl;
			// 			exit(1);

			// 		}

			// 		if (tripped_up) {
			// 			std::cerr << "\n\n\n\nI don't understand the value '" << word << "'.";
			// 			std::cerr << "\n\n\n Please change this to the type '" << temp_type << "'" << std::endl;
			// 		}

			// 		// Break if we reached a comment string
			// 		if (last_loop) break;
			// 	}
			// }

	};

}


#endif

//
//
//	/*
//	 * A curiosity -not actually used anywhere. It was my attempt at defining a general array type that could store multiple types.
//	*/
//	class Arr{
//		private:
//			std::vector<int> valTypes;
//			std::vector<double> Fvals;
//			std::vector<int> Ivals;
//			std::vector<std::string> Svals;
//
//			std::vector<int> inds;
//
//			int Icount=0;
//			int Scount=0;
//			int Fcount=0;
//			int totItems=0;
//
//		public:
//			void append (std::string val) {
//				Svals.push_back(val);
//				valTypes.push_back(typeSTR);
//				inds.push_back(Scount);
//				Scount++;
//				totItems++;
//			}
//			void append (double val) {
//				Fvals.push_back(val);
//				valTypes.push_back(typeFLOAT);
//				inds.push_back(Fcount);
//				Fcount++;
//				totItems++;
//			}
//			void append (float val) {
//				Fvals.push_back( (double) val);
//				valTypes.push_back(typeFLOAT);
//				inds.push_back(Fcount);
//				Fcount++;
//				totItems++;
//			}
//			void append (int val) {
//				Ivals.push_back(val);
//				valTypes.push_back(typeINT);
//				inds.push_back(Icount);
//				Icount++;
//				totItems++;
//			}
//
//
//			/*
//			 * Will print all values in the array
//			*/
//			void show() {
//				int i;
//				std::cerr << "[";
//				for (i=0; i<totItems-1; i++) {
//					switch (valTypes[i]) {
//						case typeSTR:
//							std::cerr << Svals[inds[i]] << ", ";
//							break;
//						case typeINT:
//							std::cerr << Ivals[inds[i]] << ", ";
//							break;
//						case typeFLOAT:
//							std::cerr << Fvals[inds[i]] << ", ";
//							break;
//						default:
//							std::cerr << "Broken!";
//							exit(1);
//					}
//				}
//				i = totItems-1;
//				switch (valTypes[i]) {
//					case typeSTR:
//						std::cerr << Svals[inds[i]];
//						break;
//					case typeINT:
//						std::cerr << Ivals[inds[i]];
//						break;
//					case typeFLOAT:
//						std::cerr << Fvals[inds[i]];
//						break;
//					default:
//						std::cerr << "Broken!";
//						exit(1);
//				}
//				std::cerr << "]";
//			}
//	};


