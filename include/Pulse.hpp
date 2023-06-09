// Pulse class definition

#ifndef PULSE_H
#define PULSE_H


// standard includes
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <random>
#include <memory>
#include <complex>
#include <fftw3.h>

// my headers
#include <Constants.hpp>
#include <DataOps.hpp>
#include <Params.hpp>

// my definitions
using namespace Constants;

class PulseTime {

	public:
		PulseTime(float strength_in = 1e-3 * 0.696, float width_in = 50, float t0_in = 0.0) : 
			strength(strength_in * auenergy<float>()/Eh<float>() * std::pow(aufor10PW<float>(),int(2))), 
			Ctau(width_in * root_pi<float>() / fsPau<float>() / 2.0),
			t0(t0_in / fsPau<float>())
	{
		//    std::clog << "Creating Pulse " << this << std::endl;
	}
		~PulseTime()
		{
			//    std::clog << "Destroying Pulse " << this << std::endl;
		}

		void setstrength(const float in);
		void setwidth(const float in);
		void sett0(const float in);

		float getstrength() { return strength; }
		float getCtau() { return Ctau; }
		float gett0() { return t0; }

		bool getenvelope(const float t,float *FF,float *dFFdt) 
		{
			if ( inpulse(t) ){
				*FF = strength * ( std::pow( cos(half_pi<float>()*(t-t0)/Ctau) , int(2)) );
				*dFFdt = -strength/2 * ( pi<float>()/Ctau * sin(pi<float>()*(t-t0)/Ctau));
				return true;
			} else {
				*FF = 0.0;
				*dFFdt = 0.0;
				return false;
			}
		}
		bool getenvelope(const float t,float *FF) 
		{
			if ( inpulse(t) ){
				*FF = strength * ( std::pow( cos(half_pi<float>()*(t-t0)/Ctau), int(2) ) );
				return true;
			} else {
				*FF = 0.0;
				return false;
			}
		}

	private:
		float strength, Ctau, t0;

		bool inpulse(const float t) 
		{
			if (t >= -Ctau && t <= Ctau){
				return true;
			} else {
				return false;
			}
		}
};


class Params;

class PulseFreq 
{

	public:
		PulseFreq & operator=(const PulseFreq & rhs); // assignment
		PulseFreq & operator+=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
		PulseFreq & operator-=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
		PulseFreq & operator*=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
		PulseFreq & operator/=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
		PulseFreq & operator*=(const float s);
		PulseFreq & diffamps(const PulseFreq &rhs);
		PulseFreq & normamps(const PulseFreq &rhs);
		PulseFreq & interfere(const PulseFreq &rhs);

		void setplans(const PulseFreq & rhs);
		void setmasterplans(fftw_plan * const forward,fftw_plan * const backward);
		void setancillaryplans(fftw_plan * const r2hc,fftw_plan * const hc2r,fftw_plan * const r2hc_2x,fftw_plan * const hc2r_2x);

	public:

		PulseFreq(Params &params);
		//PulseFreq(const float omcenter_in,const float omwidth_in,const float omonoff_in, float tspan_in);
		PulseFreq(const PulseFreq &rhs); // copy constructor
		~PulseFreq(void);

		bool addrandomphase(void);

		inline void scale(const float in){
			DataOps::mul(cvec,in,samples);
			cvec2rhophi();
		}
		inline float maxsignal(void){
			return std::pow(*std::max_element(rhovec.begin(),rhovec.end()),int(2));
		}
		inline unsigned get_samples(void) {return getsamples();}
		inline unsigned getsamples(void) {return samples;}
		inline unsigned getdt(void) {return dtime;}
		void fft_totime(void) {
			if (intime)
				std::cerr << "died here at fft_totime():\t domain() = " << domain() << "\n" << std::flush;
			assert (infreq);
			fftw_execute_dft(*FTplan_backwardPtr.get(),(fftw_complex*)cvec,(fftw_complex*)cvec);
			DataOps::mul(cvec,(float)(1./std::sqrt(samples)),samples);
			cvec2rhophi();
			infreq = false;
			intime = true;
		}
		void fft_tofreq(void) {
			if (infreq)
				std::cerr << "died here at fft_tofreq():\t domain() = " << domain() << "\n" << std::flush;
			assert (intime);
			fftw_execute_dft(*FTplan_forwardPtr.get(),(fftw_complex*)cvec,(fftw_complex*)cvec);
			DataOps::mul(cvec,(float)(1./std::sqrt(samples)),samples);
			cvec2rhophi();
			infreq=true;
			intime=false;
		}

		inline bool is_intime(void){return intime;}
		inline bool is_infreq(void){return infreq;}
		std::string domain(void){return (intime ? std::string("Time") : std::string("Frequency"));}

		int addchirp(float chirp_in) {
			assert (infreq);
			phase_GDD = chirp_in;
			addGDDtoindex(0,1);
			addGDDtoindex(samples/2,-1);
			for (unsigned i = 1; i<samples/2;i++){
				addGDDtoindex(i,1);
				addGDDtoindex(samples-i,-1);
			}
			return 0;
		}

		int addchirp(float* chirp_in) {
			assert (infreq);
			phase_GDD = chirp_in[0];
			phase_TOD = chirp_in[1];
			phase_4th = chirp_in[2];
			phase_5th = chirp_in[3];
			addGDDtoindex(0,1);
			addGDDtoindex(samples/2,-1);
			addTODtoindex(0,1);
			addTODtoindex(samples/2,-1);
			add4thtoindex(0,1);
			add4thtoindex(samples/2,-1);
			add5thtoindex(0,1);
			add5thtoindex(samples/2,-1);
			for (unsigned i = 1; i<samples/2;i++){
				addGDDtoindex(i,1);
				addGDDtoindex(samples-i,-1);
				addTODtoindex(i,1);
				addTODtoindex(samples-i,-1);
				add4thtoindex(i,1);
				add4thtoindex(samples-i,-1);
				add5thtoindex(i,1);
				add5thtoindex(samples-i,-1);
			}
			return 0;
		}

		int addchirp(std::vector<float> & chirp_in) {
			if (intime){
				std::cerr << "whoops, trying to add phase in the time domain\n" << std::flush;
				return 1;
			}
			assert(chirp_in.size()==4);
			phase_GDD = chirp_in[0];
			phase_TOD = chirp_in[1];
			phase_4th = chirp_in[2];
			phase_5th = chirp_in[3];
			addGDDtoindex(0,1);
			addGDDtoindex(samples/2,-1);
			addTODtoindex(0,1);
			addTODtoindex(samples/2,-1);
			add4thtoindex(0,1);
			add4thtoindex(samples/2,-1);
			add5thtoindex(0,1);
			add5thtoindex(samples/2,-1);
			for (unsigned i = 1; i<samples/2;i++){
				addGDDtoindex(i,1);
				addGDDtoindex(samples-i,-1);
				addTODtoindex(i,1);
				addTODtoindex(samples-i,-1);
				add4thtoindex(i,1);
				add4thtoindex(samples-i,-1);
				add5thtoindex(i,1);
				add5thtoindex(samples-i,-1);
			}
			return 0;
		}


		void attenuate(float attenfactor);
		void delay(float delayin); // expects delay in fs
		void phase(float phasein); // expects delay in units of pi , i.e. 1.0 = pi phase flip 

		int modulateamp_time(const std::vector<float> & modulation) {
			if (infreq){
				std::cerr << "whoops, trying time modulation but in frequency domain\n" << std::flush;
				return 1; }
			if (modulation.size() < samples){
				std::cerr << "size mismatch, out of range in modulateamp_time()" << std::endl;
				return 2;
			}
			modampatindx(0,modulation);
			modampatindx(samples/2,modulation);
			for (unsigned i = 1;i<samples/2;i++){
				modampatindx(i,modulation);
				modampatindx(samples-i,modulation);
			}
			return 0;
		}
		int modulateamp_time(void) {
			if (infreq){
				std::cerr << "whoops, trying time modulation but in frequency domain\t" << std::flush;
				return 1; }
			modampatindx(0);
			modampatindx(samples/2);
			for (unsigned i = 1;i<samples/2;i++){
				modampatindx(i);
				modampatindx(samples-i);
			}
			return 0;
		}

		int modulatephase_time(const std::vector<float> & modulation) {
			if (infreq){
				std::cerr << "whoops, trying time modulation but in frequency domain\n" << std::flush;
				return 1;
			}
			if (modulation.size() < samples){
				std::cerr << "size mismatch, out of range in modulateamp_time()" << std::endl;
				return 2;
			}
			modphaseatindx(0,modulation);
			modphaseatindx(samples/2,modulation);
			for (unsigned i = 1;i<samples/2;i++){
				modphaseatindx(i,modulation);
				modphaseatindx(samples-i,modulation);
			}
			return 0;
		}
		int modulatephase_time(void) {
			if (infreq){
				std::cerr << "whoops, trying time modulation but in frequency domain\n" << std::flush;
				return 1;
			}
			modphaseatindx(0);
			modphaseatindx(samples/2);
			for (unsigned i = 1;i<samples/2;i++){
				modphaseatindx(i);
				modphaseatindx(samples-i);
			}
			return 0;
		}

		void print_amp(std::ofstream & outfile);
		void print_phase(std::ofstream & outfile);
		void print_phase_powerspectrum(std::ofstream & outfile);
		void printwavelengthbins(std::ofstream * outfile);
		void appendwavelength(std::ofstream * outfile);
		void appendwavelength_deriv(std::ofstream * outfile);
		void appendwavelength_bin(std::ofstream * outfile);
		void appendfrequency(std::ofstream * outfile);
		void appendnoisy(std::ofstream * outfile);
		void appendfrequencybins(std::ofstream * outfile);
		void printfrequencybins(std::ofstream * outfile);
		void printfrequency(std::ofstream * outfile);
		void printfrequencydelay(std::ofstream * outfile, const float *delay);
		void printfrequencydelaychirp(std::ofstream * outfile, const float *delay,const float *chirp);
		void printtime(std::ofstream * outfile);
		void printwavelength(std::ofstream * outfile,const float *delay);
		float gettime(unsigned ind){return (time[ind]*fsPau<float>());}


	private:

		float m_noisescale;
		size_t m_sampleinterval;
		size_t m_lamsamples;
		unsigned long m_gain;
		unsigned m_saturate;

		size_t sampleround;

		bool intime,infreq;
		unsigned samples;
		unsigned startind,stopind,onwidth,offwidth;
		float tspan;
		float domega,lambda_center,lambda_width,omega_center,omega_width,omega_high;
		float omega_onwidth,omega_offwidth;
		float phase_GDD,phase_TOD,phase_4th,phase_5th;
		float dtime,time_center,time_wdith;

		// FFTW variables //
		// fftw defining the plans before first instantiation, this allows only one forward and backward plan to be created 
		std::shared_ptr<fftw_plan> FTplan_forwardPtr;
		std::shared_ptr<fftw_plan> FTplan_backwardPtr;
		std::shared_ptr<fftw_plan> FTplan_r2hcPtr;
		std::shared_ptr<fftw_plan> FTplan_hc2rPtr;
		std::shared_ptr<fftw_plan> FTplan_r2hc_2xPtr;
		std::shared_ptr<fftw_plan> FTplan_hc2r_2xPtr;

		std::complex<float> * cvec; // this is still fftw_malloc() for sake of fftw memory alignment optimization
		std::int32_t * ovec; // this is still fftw_malloc() for sake of fftw memory alignment optimization
		double * r_vec; // this will get fftw_malloc() for sake of fftw memory alignment 
		double * hc_vecFT; // this will get fftw_malloc() for sake of fftw memory alignment 
		double * r_vec_2x; // this will get fftw_malloc() for sake of fftw memory alignment 
		double * hc_vec_2xFT; // this will get fftw_malloc() for sake of fftw memory alignment 

		std::vector<float> rhovec;
		std::vector<float> modamp;
		std::vector<float> omega;
		std::vector<float> time;
		std::vector<float> modphase;
		std::vector<float> phivec;

		float nu0;

		unsigned i_low, i_high;

		void rhophi2cvec(const size_t i){ cvec[i] = std::polar(rhovec[i],phivec[i]);}
		void cvec2rhophi(const size_t i){ rhovec[i] = std::abs(cvec[i]); phivec[i] = std::arg(cvec[i]); }
		void rhophi2cvec(void);
		void cvec2rhophi(void);

		std::random_device rng;

		void buildvectors(const size_t s);
		void factorization(void);
		void killvectors(void);

		void addGDDtoindex(const unsigned indx,const int omega_sign) {
			phivec[indx] += omega_sign*phase_GDD*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(2));
			rhophi2cvec(indx);
		}	
		void addTODtoindex(const unsigned indx,const int omega_sign) {
			phivec[indx] += omega_sign*phase_TOD*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(3));
			rhophi2cvec(indx);
		}	
		void add4thtoindex(const unsigned indx,const int omega_sign) {
			phivec[indx] += omega_sign*phase_4th*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(4));
			rhophi2cvec(indx);
		}	
		void add5thtoindex(const unsigned indx,const int omega_sign) {
			phivec[indx] += omega_sign*phase_5th*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(5));
			rhophi2cvec(indx);
		}	

		void modampatindx(const unsigned indx,const std::vector<float> & modvec) {
			rhovec[indx] *= modvec[indx];
			rhophi2cvec(indx);
		}
		void modampatindx(const unsigned indx) {
			rhovec[indx] *= modamp[indx];
			rhophi2cvec(indx);
		}
		/*
		   E(t) = E0exp(iwt)=E0 exp(i w(t) (t-t0))
		   = E0 exp(i (w0+GDD*t) (t-t0) )
		   = E0 exp(i ((w0 * t) - (GDD*t*t0) + (GDD*t**2) - (w0*t0)) so adding a delay means adding an arbitrary phase of w0*t0 and subtracting GDD*t*t0 ( a linear term in phase which we know only dials the pulse position in z, or likewise phase in time... unless t0 is actually t0(t), then you have to include this one as well as w0*t0(t).
		   = E0 exp(i ((w0 * t) + (GDD*t**2) - (w0*t0) - (GDD*t*t0) ) )
		   = E0 exp(i w(t)*t) exp(i -w(t)*t0(t))
		   = E(t) exp(i -(w0 + 1/GDD*t + 1/TOD*t**2 + 1/FOD*t**3)*t0(t))
		 */
		void modphaseatindx(const unsigned indx,const std::vector<float> & modvec) {
			if (modvec[indx]!=0){
				phivec[indx] += modvec[indx];
				rhophi2cvec(indx);
			}
		}
		void modphaseatindx(const unsigned indx) {
			if (modphase[indx] != 0){
				phivec[indx] += modphase[indx];
				rhophi2cvec(indx);
			}
		}
		float rising(const unsigned indx) {
			return std::pow(sin(float(Constants::half_pi<float>()*(indx-startind)/onwidth)),int(2));
		}
		float falling(const unsigned indx) {
			return std::pow(sin(float(Constants::half_pi<float>()*(stopind-indx)/offwidth)),int(2));
		}
};

#endif
