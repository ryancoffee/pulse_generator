#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <DataOps.hpp>
#include <cstdint> 
#include <limits>

using namespace DataOps;

int main(void){
        std::vector<double> v = {-1.,3.,4.,10.,-100.,101.};
	std::cout << v;
	auto bounds = std::minmax_element(v.begin(),v.end());
	std::cout << "min = " << *bounds.first << " and max = " << *bounds.second << "\n" << std::flush;
	double scale;
	long int max = std::numeric_limits<int8_t>::max();
	long int min = std::numeric_limits<int8_t>::min();
	if (*bounds.second>std::abs(*bounds.first)){
		scale = max/ *bounds.second;
	} else {
		scale = min/ *bounds.first;
	}
	std::vector<int> result(v.size());
	std::transform(v.begin(),v.end(),result.begin(),[scale](double x){return int(scale*x);});
	std::cout << result;
	std::cout << "max = " << max << " and min = " << min << std::endl;

	std::string filestr("../data_fs/reference/fromNikita/OUTPUT_YAG_Density_1000.dat");
	std::ifstream infile(filestr.c_str(),std::ios::in);
	std::vector< std::vector<double> > data;
	infile >> data;
	double tstep = data[1][0] - data[0][0];
	for (size_t i=0;i<20;++i){ //data.size();i+= 10){
		std::cout << data[i][0] << "\t" << data[i][5] << "\t" << (data[i+1][5] - data[i][5])/tstep << "\n";
	}
	std::cout << std::flush;

	std::vector<double> xData = {0,1,5,8,20,30};
	std::vector<double> yData = {0,0,50,80,120,130};
	std::vector<double> samples(1,0); 
	std::vector<double> values; 
	samples.reserve(100);
	values.reserve(100); 
	double step = .75;
	while (samples.size()<100){
		samples.push_back(samples.back()+step);
	}
	for (size_t i=0;i<samples.size();++i){
		std::cout << samples[i] << "\t" << interpolate(xData,yData,samples[i]) << "\n" << std::flush;
	}
	


        return 0;
}

