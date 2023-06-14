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
	delete amp0thDistPtr;
	delete amp1stDistPtr;
	delete amp2ndDistPtr;
	delete amp3rdDistPtr;
	delete amp4thDistPtr;
	delete amp5thDistPtr;
	delete chirpDistPtr;
	delete TODDistPtr;
	delete FODDistPtr;
	delete fifthODDistPtr;
	std::cerr << "Leaving Params::~Params()\n" << std::flush;
}

float Params::getDelay(void){return (*delays_distributionPtr)(rng);}
float Params::getAmp(void){return (*laser_distributionPtr)(rng);}

Params & Params::initAmp(float zeroth,float first,float second,float third,float fourth,float fifth)
{
	amp0thDistPtr = new std::uniform_real_distribution<float>( 0, zeroth );
	amp1stDistPtr = new std::uniform_real_distribution<float>( -first, first );
	amp2ndDistPtr = new std::uniform_real_distribution<float>( -second, second );
	amp3rdDistPtr = new std::uniform_real_distribution<float>( -third, third );
	amp4thDistPtr = new std::uniform_real_distribution<float>( -fourth, fourth );
	amp5thDistPtr = new std::uniform_real_distribution<float>( -fifth, fifth );
	return *this;
}
Params & Params::initChirp(float second,float third,float fourth,float fifth)
{
	chirpDistPtr = new std::uniform_real_distribution<float>( -second, second );
	TODDistPtr = new std::uniform_real_distribution<float>( -third, third );
	FODDistPtr = new std::uniform_real_distribution<float>( -fourth, fourth );
	fifthODDistPtr = new std::uniform_real_distribution<float>( -fifth, fifth );
	return *this;
}

std::vector<float> & Params::getAmpMod(std::vector<float> & v)
{
	assert(5==v.size());
	v[0] = (*amp0thDistPtr)(rng);
	v[1] = (*amp1stDistPtr)(rng);
	v[2] = (*amp2ndDistPtr)(rng);
	v[3] = (*amp3rdDistPtr)(rng);
	v[4] = (*amp4thDistPtr)(rng);
	v[5] = (*amp5thDistPtr)(rng);
	return v;
}
std::vector<float> & Params::getChirp(std::vector<float> & v)
{
	//std::vector<float> v(4,float(0));
	assert(4==v.size());
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
