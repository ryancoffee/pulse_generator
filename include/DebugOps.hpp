#include <iostream>

namespace DebugOps
{
	int return_early(const int r){
		std::cerr << "\n"
			<< "===========================================\n"
			<< "======== returning early ==================\n"
			<< "======== exit code " << r << "================\n"
			<< std::flush;
		return r;
	}

	void pushout(std::string str){
		std::cerr << str << std::endl << std::flush;
	}

	template <typename T>
	void pushout(std::string str, T n){
		std::cerr << str << std::to_string(n) << std::endl << std::flush;
	}
}
