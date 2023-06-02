#include <Params.hpp>

ScanParams::ScanParams(void)
: dpdelay(double(0))
, usechirpnoise(false)
, userandphase(false)
{

	std::cerr << "Constructor ScanPrams()" << std::endl;

	chirpvec.resize(4,double(0));
	std::string filebase = std::string(getenv("filebase"));

	delays_distributionPtr = new std::normal_distribution<double>(
			double(atof(getenv("delays_mean"))),
			double(atof(getenv("delays_std"))));
	delays_unidistributionPtr = new std::uniform_real_distribution<double> (
			double(atof(getenv("delays_mean")))-double(atof(getenv("delays_std"))),
			double(atof(getenv("delays_mean")))+double(atof(getenv("delays_std"))));
	xray_pos_distributionPtr = new std::normal_distribution<double> (
			double(atof(getenv("xray_pos_distribution_normal_mean"))),
			double(atof(getenv("xray_pos_distribution_normal_std"))));
	laser_pos_distributionPtr = new std::normal_distribution<double> (
			double(atof(getenv("laser_pos_distribution_normal_mean"))),
			double(atof(getenv("laser_pos_distribution_normal_std"))));
	xray_distributionPtr = new std::lognormal_distribution<double> (
			double(atof(getenv("xray_distribution_lognormal_mean"))),
			double(atof(getenv("xray_distribution_lognormal_std"))));
	laser_distributionPtr = new std::lognormal_distribution<double> (
			double(atof(getenv("laser_distribution_lognormal_mean"))),
			double(atof(getenv("laser_distribution_lognormal_std"))));


}

ScanParams::~ScanParams(void)
{
	delete delays_distributionPtr;
	delete delays_unidistributionPtr;
	delete xray_pos_distributionPtr;
	delete laser_pos_distributionPtr;
	delete xray_distributionPtr;
	delete laser_distributionPtr;
	if (usechirpnoise){
		delete chirpnoiseDistPtr;
		delete TODnoiseDistPtr;
		delete FODnoiseDistPtr;
		delete fifthODnoiseDistPtr;
	}
}


/* =========== chirp interfaces ============= */

void ScanParams::chirp(double second, double third, double fourth = double(0), double fifth = double(0))
{
	chirpvec[0] = second;
	chirpvec[1] = third;
	chirpvec[2] = fourth;
	chirpvec[3] = fifth;
}

void ScanParams::initchirpnoise(double second,double third,double fourth = double(0),double fifth = double(0))
{
	chirpnoiseDistPtr = new std::normal_distribution<double>( chirpvec[0], second );
	TODnoiseDistPtr = new std::normal_distribution<double>( chirpvec[1], third );
	FODnoiseDistPtr = new std::normal_distribution<double>( chirpvec[2], fourth );
	fifthODnoiseDistPtr = new std::normal_distribution<double>( chirpvec[3], fifth );
}

std::vector<double> & ScanParams::getchirpnoise(void)
{
	std::vector<double> v(4,double(0));
	v[0] = (*chirpnoiseDistPtr)(rng);
	v[1] = (*TODnoiseDistPtr)(rng);
	v[2] = (*FODnoiseDistPtr)(rng);
	v[3] = (*fifthODnoiseDistPtr)(rng);
	return v;

}


