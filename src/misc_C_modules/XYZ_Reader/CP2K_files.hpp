#ifndef CP2K_FILE_HEADER
#define CP2K_FILE_HEADER

#include <string>
#include <iostream>
#include <vector>

std::string prep_string(const std::string &str) {
	std::string out = str;

	// Uppercase string
	transform(out.begin(), out.end(), out.begin(), ::toupper);

	// Remove Spaces either side.
	auto start_ind = 0;
	const char space = *" ";
	const char *wrd = str.c_str();
	for (auto i=0; i<str.length(); i++) {
		if (wrd[i] == space) continue;
		start_ind = i;
		break;
	}

	auto end_ind = str.length();
	for (auto i=str.length()-1; i>=0; i--) {
		if (wrd[i] == space) continue;
		end_ind = i+1;
		break;
	}

	return out.substr(start_ind, end_ind-start_ind);;
}

namespace cp2k {

	class Tag {
		private:
		protected:
		public:
			std::string name;
			std::string data;

			Tag(const std::string &in_name, const std::string &in_data) {
				name = in_name;
				data = in_data;
			}

			std::string print_tag(int nspaces=0) {
				std::string spaces = "";
				for (auto i=0; i<nspaces; i++) spaces += " ";
				std::string tag;
				tag = name + "   " + data;
				return tag;
			}
	};

	class Section {
		private:
		protected:
			Section create_section(const std::string &sect_name) {
				std::string name = prep_string(sect_name);
				Section Sect = Section(name);
				return Sect;
			}

			Tag create_tag(const std::string &tag_name,
						   const std::string &data) {
				std::string name = prep_string(tag_name);
				Tag tag = Tag(name, data);
				return tag;
			}

		public:
			std::string name;
			std::vector<Section> sections;
			std::vector<Tag> tags;
			unsigned int nSect=0;
			unsigned int nTag=0;

			Section(const std::string &sect_name) {
				name = sect_name;
			}

			void add_section(const std::string &sect_name) {
				Section Sect = create_section(sect_name);
				sections.push_back(Sect);
				nSect++;
			}

			void add_tag(const std::string &tag_name,
						 const std::string &data) {
				Tag tag = create_tag(tag_name, data);
				tags.push_back(tag);
				nTag++;
			}

			std::string construct_file(int nspaces=0, int indent_level=3) {
				std::string file_txt = "";
				std::string spaces="", indent="";
				for (auto i=0; i<nspaces; i++) spaces += " ";
				for (auto i=0; i<nspaces+indent_level; i++) indent += " ";

				file_txt = spaces + "&" + name + "\n";

				for (auto tag: tags)
					file_txt = file_txt + indent + tag.print_tag(nspaces) + "\n";

				nspaces += indent_level;
				for (auto child: sections)
					file_txt = file_txt + "\n" + child.construct_file(nspaces);

				file_txt = file_txt + spaces + "&END " + name + "\n";
				return file_txt;
			}

			bool find_section(const std::string &sect_name, Section &Sect) {
				auto search_str = prep_string(sect_name);

				for (auto i=0; i<nSect; i++) {
					if (sections[i].name == search_str) {
						Sect = sections[i];
						return true;
					}

					else {
						bool is_there = sections[i].find_section(search_str, Sect);
						if (is_there == true) return true;
					}
				}

				return false;
			}

	};

	class INP {
		private:

		protected:

			Section create_section(const std::string &sect_name) {
				auto name = prep_string(sect_name);
				auto Sect = Section(name);
				return Sect;
			}

			Tag create_tag(const std::string &tag_name,
							   const std::string &data) {
				auto name = prep_string(tag_name);
				Tag tag = Tag(name, data);
				return tag;
			}

		public:
			std::vector<Section> sections;
			std::vector<Tag> tags;
			unsigned int nSect=0;
			unsigned int nTag=0;

			// void add_tag(const std::string &tag_name,
			// 			 const std::string &data,
			// 		     const std::string &sect_to_change="") {

			// 	// If we add to root
			// 	if (sect_to_change == "") {
			// 		Tag tag = create_tag(tag_name, data);
			// 		tags.push_back(tag);
			// 		nTag++;
			// 	}

			// 	// Add to another section.
			// 	else {
			// 		Section Sect("");
			// 		bool is_there = find_section(sect_to_change, Sect);
			// 		if (not is_there) {
			// 			std::cerr << "No section named: '" << sect_to_change << "'.\n\n";
			// 		} else {
			// 			Sect.add_tag(tag_name, data);
			// 		}
			// 	}
			// }

			bool find_section(const std::string &sect_name, Section &Sect) {
				auto search_str = prep_string(sect_name);

				for (auto i=0; i<nSect; i++) {
					if (sections[i].name == search_str) {
						Sect = sections[i];
						return true;
					}

					else {
						bool is_there = sections[i].find_section(search_str, Sect);
						if (is_there == true) return true;
					}
				}

				return false;
			}

			void add_section(const std::string &sect_to_add,
						     const std::string &sect_to_change="") {

				// If we add to root
				if (sect_to_change == "") {
					auto Sect = create_section(sect_to_add);
					sections.push_back(Sect);
					nSect++;
				}

				// Add to another section.
				else {
					Section Sect("");
					bool is_there = find_section(sect_to_change, Sect);
					if (not is_there) {
						std::cerr << "No section named: '" << sect_to_change << "'.\n\n";
					} 

					else {
						std::cout << Sect.sections.size() << std::endl;
						Sect.add_section(sect_to_add);
						std::cout << Sect.sections.size() << std::endl;
					}
				}
			}

			std::string construct_file() {
				std::string file_txt = "";
				for (auto sect: sections){
					file_txt = file_txt + sect.construct_file() + "\n";
				}
				
				for (auto tag: tags)
					file_txt = file_txt + tag.print_tag() + "\n";

				return file_txt;
			}
	};

}




#endif