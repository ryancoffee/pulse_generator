#include <Params.hpp>

Params::Params(void):
	usechirpnoise(false)
{

	std::cerr << "Constructor Prams()" << std::endl;

	chirpvec.resize(4,float(0));
	std::string filebase = std::string(getenv("filebase"));

	delays_distributionPtr = new std::normal_distribution<float>(
			float(atof(getenv("delays_mean"))),
			float(atof(getenv("delays_std"))));
	laser_distributionPtr = new std::uniform_real_distribution<float> (
			float(atof(getenv("amp_mean"))) - float(atof(getenv("amp_std"))),
			float(atof(getenv("amp_mean"))) + float(atof(getenv("amp_std"))),
			);


}

Params::~Params(void)
{
	delete delays_distributionPtr;
	delete laser_distributionPtr;
	if (usechirpnoise){
		delete chirpnoiseDistPtr;
		delete TODnoiseDistPtr;
		delete FODnoiseDistPtr;
		delete fifthODnoiseDistPtr;
	}
}


/* =========== chirp interfaces ============= */

void Params::initchirp(float second, float third, float fourth = float(0), float fifth = float(0))
{
	chirpvec[0] = second;
	chirpvec[1] = third;
	chirpvec[2] = fourth;
	chirpvec[3] = fifth;
}

void Params::initchirpnoise(float second,float third,float fourth = float(0),float fifth = float(0))
{
	chirpnoiseDistPtr = new std::normal_distribution<float>( chirpvec[0], second );
	TODnoiseDistPtr = new std::normal_distribution<float>( chirpvec[1], third );
	FODnoiseDistPtr = new std::normal_distribution<float>( chirpvec[2], fourth );
	fifthODnoiseDistPtr = new std::normal_distribution<float>( chirpvec[3], fifth );
}

std::vector<float> & Params::getchirpnoise(void)
{
	std::vector<float> v(4,float(0));
	v[0] = (*chirpnoiseDistPtr)(rng);
	v[1] = (*TODnoiseDistPtr)(rng);
	v[2] = (*FODnoiseDistPtr)(rng);
	v[3] = (*fifthODnoiseDistPtr)(rng);
	return v;

}


