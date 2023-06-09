#ifndef PARAMS_H

#define PARAMS_H

#include <Pulse.hpp>
#include <Constants.hpp>
#include <random>
#include <cmath>

using namespace Constants;


class Params 
{
	friend class PulseFreq;

	public:
	Params(std::string filebase = "../data",float delays_mean=0.,float delays_std=50.,float amp_mean=1.,float amp_std=1.);
	~Params();

	size_t setNpulses(size_t in){npulses=in;return npulses;}
	size_t getNpulses(void){return npulses;}
	std::string filebase(std::string in){fbase=in;return fbase;}
	std::string filebase(void){return fbase;}
	float setTspan( float in) { tspan = in; return tspan;} 
	float getTspan( void) { return tspan;} 

	float lambda_0(float in){ lam0 = in; return lam0;} 
	float lambda_0(void){return lam0;} 
	float lambda_width(float in){ lamwidth = in; return lamwidth;}
	float lambda_width(void){return lamwidth;}
	float lambda_onoff(float in ){ lamonoff = in; return lamonoff;} 
	float lambda_onoff(void){return lamonoff;} 
	float omega_low(void){return twopi<float>()/(lam0+lamwidth/float(2))*C_nmPfs<float>()*fsPau<float>();}
	float omega_high(void){return twopi<float>()/(lam0-lamwidth/float(2))*C_nmPfs<float>()*fsPau<float>();}
	float omega_width(void){return omega_high()-omega_low();}
	float omega_onoff(void){return ( twopi<float>()/(lam0+lamwidth/float(2)-lamonoff)*C_nmPfs<float>()*fsPau<float>() - omega_low() );}
	float omega0(void){ return (omega_high()+omega_low())/float(2);}

	Params & initChirp(float second=0,float third=0,float fourth=0,float fifth=0);
	std::vector<float> & getChirp(void);
	float getAmp(void);
	float getDelay(void);
	Params & setnulims(float low, float high);
	Params & set_lamsamples(size_t in){lamsamples = in; return *this;}
	Params & set_gain(float in){ gain = in; return *this;} 
	Params & set_noisescale(float in){ noisescale = in; return *this;}
	Params & set_sampleinterval(size_t in){ sampleinterval = in; return *this;}
	Params & set_saturate(uint16_t in){saturate = in; return *this;}

	private:

	size_t npulses;
	float tspan;
	float lam0,lamwidth,lamonoff;
	float nu_low,nu_high;
	size_t lamsamples;
	float gain;
	size_t sampleinterval;
	uint16_t saturate;
	float noisescale;

	std::string fbase;

	/* =========== random members =========== */
	//std::default_random_engine rng;
	std::random_device rng;

	std::uniform_real_distribution<float>* delays_distributionPtr;
	std::uniform_real_distribution<float>* laser_distributionPtr; 
	std::uniform_real_distribution<float>* chirpDistPtr;
	std::uniform_real_distribution<float>* TODDistPtr;
	std::uniform_real_distribution<float>* FODDistPtr;
	std::uniform_real_distribution<float>* fifthODDistPtr;
	std::normal_distribution<float>* chirpnoiseDistPtr;
	std::normal_distribution<float>* TODnoiseDistPtr;
	std::normal_distribution<float>* FODnoiseDistPtr;
	std::normal_distribution<float>* fifthODnoiseDistPtr;

};
#endif

