#ifndef SCANPARAMS

#define SCANPARAMS

#include <Pulse.hpp>
#include <Constants.hpp>
#include <random>
#include <cmath>

using namespace Constants;


class ScanParams 
{
	friend class PulseFreq;

	public:
	ScanParams();
	~ScanParams();

	size_t nimages(size_t in){nimgs=in;return nimgs;}
	inline size_t nimages(void){return nimgs;}
	double dalpha(double in){dalph=in;return dalph;}
	inline double dalpha(void){return dalph;}
	double doublepulsedelay(double in){dpdelay=in;return dpdelay;}
	inline double doublepulsedelay(void){return dpdelay;}
	bool doublepulse(const bool in){dpbool=in;return dpbool;}
	inline bool doublepulse(void){return dpbool;}
	std::string filebase(std::string in){fbase=in;return fbase;}
	inline std::string filebase(void){return fbase;}
	std::string calfilebase(std::string in){calfbase=in;return calfbase;}
	inline std::string calfilebase(void){return calfbase;}
	double tspan( double in) { tspn = in; return tspn;} 
	inline double tspan( void) { return tspn;} 

	double lambda_0(double in){ lam_0 = in; return lam_0;} 
	inline double lambda_0(void){return lam_0;} 
	double lambda_width(double in){ lam_w = in; return lam_w;}
	inline double lambda_width(void){return lam_w;}
	double lambda_onoff(double in ){ lam_onoff = in; return lam_onoff;} 
	inline double lambda_onoff(void){return lam_onoff;} 
	inline double omega_low(void){return twopi<double>()/(lam_0+lam_w/double(2))*C_nmPfs<double>()*fsPau<double>();}
	inline double omega_high(void){return twopi<double>()/(lam_0-lam_w/double(2))*C_nmPfs<double>()*fsPau<double>();}
	inline double omega_width(void){return omega_high()-omega_low();}
	inline double omega_onoff(void){return ( twopi<double>()/(lam_0+lam_w/double(2)-lam_onoff)*C_nmPfs<double>()*fsPau<double>() - omega_low() );}
	inline double omega0(void){ return (omega_high()+omega_low())/double(2);}

	size_t ngroupsteps(size_t in){ngrpsteps = in; return ngrpsteps;}
	inline size_t ngroupsteps(void){return ngrpsteps;}

	double groupdelay(double in){grpdelay = in; return grpdelay;}
	inline double groupdelay(void){return grpdelay;}

	double backdelay(double in){bkdelay = in; return bkdelay;}
	inline double backdelay(void){return bkdelay;}

	double groupstep(void) {return grpdelay/ngrpsteps ;}
	inline double backstep(void) {return bkdelay/ngrpsteps;}

	size_t netalon(size_t in){netln = in; return netln;}
	inline size_t netalon(void){return netln;}

	double etalondelay(double in){ etdelay = in; return etdelay;}
	inline double etalondelay(void){ return etdelay;}

	double etalonreflectance(double in){ etreflect = in; return etreflect;}
	inline double etalonreflectance(void){ return etreflect;}
	double interferedelay(double in){interferedly = in; return interferedly;}
	inline double interferedelay(void){return interferedly;}

	/* =========== chirp interfaces ============= */
	void chirp(double second, double third, double fourth, double fifth);
	std::vector<double> & getchirp(void){return chirpvec;}
	std::vector<double> & getchirpnoise(void);

	inline double chirp(double in){chirpvec[0] = in; return chirpvec[0];}
	inline double TOD(double in){chirpvec[1] = in; return chirpvec[1];}
	inline double FOD(double in){chirpvec[2] = in; return chirpvec[2];}
	inline double fifthOD(double in){chirpvec[3] = in; return chirpvec[3];}
	inline double chirp(void){return chirpvec[0];}
	inline double TOD(void){return chirpvec[1];}
	inline double FOD(void){return chirpvec[2];}
	inline double fifthOD(void){return chirpvec[3];}

	inline bool addchirpnoise(bool in){usechirpnoise = in; return usechirpnoise;}
	inline bool addchirpnoise(void){return usechirpnoise;}
	void initchirpnoise(double second,double third,double fourth,double fifth);

	inline bool addrandomphase(bool in) {userandphase = in; return userandphase;}
	inline bool addrandomphase(void) {return userandphase;}

	/* =========== random interfaces =========== */
	inline double delays_uniform(void){return (*delays_distributionPtr)(rng);}
	inline double delays_normal(void){return (*delays_unidistributionPtr)(rng);}
	inline double xray_pos_rand(void){return (*xray_pos_distributionPtr)(rng);}
	inline double laser_pos_rand(void){return (*laser_pos_distributionPtr)(rng);}
	inline double xray_inten_rand(void){return (*xray_distributionPtr)(rng);}
	inline double laser_inten_rand(void){return (*laser_distributionPtr)(rng);}


	private:

	size_t nimgs,ngrpsteps,netln;
	double interferedly;
	double grpdelay,bkdelay,etdelay,etreflect;
	double dalph;
	double tspn;
	double dpdelay;
	double lam_0,lam_w,lam_onoff;
	bool dpbool;
	bool usechirpnoise;
	bool userandphase;

	std::string fbase;
	std::string calfbase;

	std::vector<double> chirpvec;//(4,double(0));

	/* =========== random members =========== */
	//std::default_random_engine rng;
	std::random_device rng;

	std::normal_distribution<double> * delays_distributionPtr;
	std::uniform_real_distribution<double>* delays_unidistributionPtr;
	std::normal_distribution<double> *xray_pos_distributionPtr;
	std::normal_distribution<double>* laser_pos_distributionPtr;
	std::lognormal_distribution<double>* xray_distributionPtr;
	std::lognormal_distribution<double>* laser_distributionPtr; 
	std::normal_distribution<double>* chirpnoiseDistPtr;
	std::normal_distribution<double>* TODnoiseDistPtr;
	std::normal_distribution<double>* FODnoiseDistPtr;
	std::normal_distribution<double>* fifthODnoiseDistPtr;

};
#endif

