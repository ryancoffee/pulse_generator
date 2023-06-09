#include <Params.hpp>

Params::Params(std::string filebase,float delays_mean,float delays_std,float amp_mean,float amp_std)
	:fbase(filebase)
{
	//std::cerr << "Constructor Params()" << std::endl;

	delays_distributionPtr = new std::uniform_real_distribution<float>(
			delays_mean - delays_std,
			delays_mean + delays_std
			);
	laser_distributionPtr = new std::uniform_real_distribution<float> (
			amp_mean - amp_std,
			amp_mean + amp_std
			);
}


Params::~Params(void)
{
	delete delays_distributionPtr;
	delete laser_distributionPtr;
	delete chirpDistPtr;
	delete TODDistPtr;
	delete FODDistPtr;
	delete fifthODDistPtr;
	std::cerr << "Leaving Params::~Params()\n" << std::flush;
}

float Params::getDelay(void){return (*delays_distributionPtr)(rng);}
float Params::getAmp(void){return (*laser_distributionPtr)(rng);}

Params & Params::initChirp(float second,float third,float fourth,float fifth)
{
	chirpDistPtr = new std::uniform_real_distribution<float>( -second, second );
	TODDistPtr = new std::uniform_real_distribution<float>( -third, third );
	FODDistPtr = new std::uniform_real_distribution<float>( -fourth, fourth );
	fifthODDistPtr = new std::uniform_real_distribution<float>( -fifth, fifth );
	return *this;
}

std::vector<float> & Params::getChirp(void)
{
	std::vector<float> v(4,float(0));
	v[0] = (*chirpDistPtr)(rng);
	v[1] = (*TODDistPtr)(rng);
	v[2] = (*FODDistPtr)(rng);
	v[3] = (*fifthODDistPtr)(rng);
	return v;
}

Params & Params::setnulims(float low, float high)
{
	nu_low = low;
	nu_high = high;
	return *this;
}
