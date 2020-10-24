#ifndef INP_PARSER_HEADER
#define INP_PARSER_HEADER

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <vector>



namespace inp {

	enum var_type {typeINT, typeFLOAT, typeSTR};

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
				while (getline(file, line)) {
					if (line.empty())  continue;

					parse_line(line);
				}
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
				if (type == typeFLOAT)
					doubleVals[var_name] = std::stod(data);
				else if (type == typeINT)
					intVals[var_name] = std::stoi(data);
				else if (type == typeSTR)
					first_elm = data.substr(0, 1);
					if (first_elm == space) data = data.substr(1, data.size()-1);
					stringVals[var_name] = data;
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
					return typeSTR;
				}
			}


		public:
			std::string fp;
			std::fstream file;

			std::unordered_map<std::string, std::string> stringVals {};
			std::unordered_map<std::string, double> doubleVals {};
			std::unordered_map<std::string, int> intVals {};


			void read(std::string filepath) {
				fp = filepath;
				file.open(fp, std::ios::in);

				std::cout << "Opened " << fp << std::endl;

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
						std::cerr << "No string named '" << var_name << "'" << std::endl;
						throw 0;
					} else {
						return default_;
					}
				}
			}

			// Will get a value if it exists if not throw and error
			double get_double(std::string var_name, double default_= unset) {
				std::cout << var_name << " \n"; // << doubleVals.find(var_name) << " " << doubleVals.end() << std::endl;
				if (doubleVals.find(var_name) != doubleVals.end()) {
					return doubleVals[var_name];
				}

				else {
					if (default_ == unset) {
						std::cerr << "No string named '" << var_name << "'" << std::endl;
						throw 0;
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
						std::cerr << "No string named '" << var_name << "'" << std::endl;
						throw 0;
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
						std::cerr << "No string named '" << var_name << "'" << std::endl;
						throw 0;
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
			// 		std::cout << word << ",  " << temp_type << std::endl;

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
			// 			std::cerr << "I don't recognise type '" << temp_type << "'.";
			// 			std::cerr << " Change the `parse_line_util` function in the `input_function.hpp` to change this.";
			// 			std::cerr << std::endl;
			// 			throw "Exit";

			// 		}

			// 		if (tripped_up) {
			// 			std::cerr << "I don't understand the value '" << word << "'.";
			// 			std::cerr << "\n Please change this to the type '" << temp_type << "'" << std::endl;
			// 		}

			// 		// Break if we reached a comment string
			// 		if (last_loop) break;
			// 	}
			// }

	};

}


#endif