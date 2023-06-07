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
	Params();
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

	/* =========== chirp interfaces ============= */
	void initchirp(float second, float third, float fourth, float fifth);
	std::vector<float> & getchirp(void){return chirpvec;}
	std::vector<float> & getchirpnoise(void);

	float chirp(float in){chirpvec[0] = in; return chirpvec[0];}
	float TOD(float in){chirpvec[1] = in; return chirpvec[1];}
	float FOD(float in){chirpvec[2] = in; return chirpvec[2];}
	float fifthOD(float in){chirpvec[3] = in; return chirpvec[3];}
	float chirp(void){return chirpvec[0];}
	float TOD(void){return chirpvec[1];}
	float FOD(void){return chirpvec[2];}
	float fifthOD(void){return chirpvec[3];}

	bool addchirpnoise(bool in){usechirpnoise = in; return usechirpnoise;}
	bool addchirpnoise(void){return usechirpnoise;}
	void initchirpnoise(float second,float third,float fourth,float fifth);

	/* =========== random interfaces =========== */
	float delays_rand(void){return (*delays_distributionPtr)(rng);}
	float laser_inten_rand(void){return (*laser_distributionPtr)(rng);}


	private:

	size_t npulses;
	float tspan;
	float lam0,lamwidth,lamonoff;
	bool usechirpnoise;
	bool userandphase;

	std::string fbase;
	std::string calfbase;

	std::vector<float> chirpvec;//(4,float(0));

	/* =========== random members =========== */
	//std::default_random_engine rng;
	std::random_device rng;

	std::normal_distribution<float> * delays_distributionPtr;
	std::uniform_real_distribution<float>* laser_distributionPtr; 
	std::normal_distribution<float>* chirpnoiseDistPtr;
	std::normal_distribution<float>* TODnoiseDistPtr;
	std::normal_distribution<float>* FODnoiseDistPtr;
	std::normal_distribution<float>* fifthODnoiseDistPtr;

};
#endif

