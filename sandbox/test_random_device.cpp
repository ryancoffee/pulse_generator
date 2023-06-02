#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <vector>
#include <DataOps.hpp>

int main(void){
	std::random_device rng;
	std::cout << "OK, that 'std::random_device rng;' worked" << std::endl;
	std::vector<double> v = {-1.,3.,4.,10.,-100.,101.};
	for (size_t i = 0;i<v.size();++i){
	std::cout << v[i] << "\t";
	}
	std::cout << "\n" << std::flush;

	return 0;
}
