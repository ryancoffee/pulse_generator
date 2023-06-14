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

		PulseTime & setstrength(const float in);
		PulseTime & setwidth(const float in);
		PulseTime & sett0(const float in);

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

		PulseFreq & setplans(const PulseFreq & rhs);
		PulseFreq & setancillaryplans(const PulseFreq & rhs);
		PulseFreq & setmasterplans(fftw_plan * const forward,fftw_plan * const backward);
		PulseFreq & setmasterancillaryplans(fftw_plan * const r2hc,fftw_plan * const hc2r,fftw_plan * const r2hc_2x,fftw_plan * const hc2r_2x);

	public:

		PulseFreq(Params &params);
		//PulseFreq(const float omcenter_in,const float omwidth_in,const float omonoff_in, float tspan_in);
		PulseFreq(const PulseFreq &rhs); // copy constructor
		~PulseFreq(void);

		bool addrandomphase(void);

		inline PulseFreq & scale(const float in){
			DataOps::mul(cvec,(double)in,samples);
			cvec2rhophi();
			return *this;
		}
		inline float maxsignal(void){
			return std::pow(*std::max_element(rhovec.begin(),rhovec.end()),int(2));
		}
		inline unsigned get_samples(void) {return getsamples();}
		inline unsigned getsamples(void) {return samples;}
		inline unsigned getdt(void) {return dtime;}
		PulseFreq & fft_totime(void) {
			if (intime)
				std::cerr << "died here at fft_totime():\t domain() = " << domain() << "\n" << std::flush;
			assert (infreq);
			fftw_execute_dft(*FTplan_backwardPtr,(fftw_complex*)cvec,(fftw_complex*)cvec);
			DataOps::mul(cvec,(double)(1./std::sqrt(samples)),samples);
			cvec2rhophi();
			infreq = false;
			intime = true;
			return *this;
		}
		PulseFreq & fft_tofreq(void) {
			if (infreq)
				std::cerr << "died here at fft_tofreq():\t domain() = " << domain() << "\n" << std::flush;
			assert (intime);
			fftw_execute_dft(*FTplan_forwardPtr,(fftw_complex*)cvec,(fftw_complex*)cvec);
			DataOps::mul(cvec,(double)(1./std::sqrt(samples)),samples);
			cvec2rhophi();
			infreq=true;
			intime=false;
			return *this;
		}

		PulseFreq & fillphase(std::vector < float > & data);
		PulseFreq & fillspect(std::vector < float > & data);
		PulseFreq & fillfreq(std::vector < float > & data);
		PulseFreq & filltime(std::vector < float > & data);
		PulseFreq & filltime_envelope(std::vector < float > & data);

		inline bool is_intime(void){return intime;}
		inline bool is_infreq(void){return infreq;}
		std::string domain(void){return (intime ? std::string("Time") : std::string("Frequency"));}

		PulseFreq & addchirp(float chirp_in) {
			assert (infreq);
			phase_GDD = chirp_in;
			addGDDtoindex(0,1);
			addGDDtoindex(samples/2,-1);
			for (unsigned i = 1; i<samples/2;i++){
				addGDDtoindex(i,1);
				addGDDtoindex(samples-i,-1);
			}
			return *this;
		}

		PulseFreq & addchirp(float* chirp_in) {
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
			return *this;
		}

		PulseFreq & setchirp(std::vector<float> & chirp_in) {
			if (intime){
				std::cerr << "whoops, trying to add phase in the time domain\n" << std::flush;
				return *this;
			}
			assert(chirp_in.size()==4);
			phase_GDD = chirp_in[0];
			phase_TOD = chirp_in[1];
			phase_4th = chirp_in[2];
			phase_5th = chirp_in[3];
			setGDDtoindex(0,1);
			setGDDtoindex(samples/2,-1);
			addTODtoindex(0,1);
			addTODtoindex(samples/2,-1);
			add4thtoindex(0,1);
			add4thtoindex(samples/2,-1);
			add5thtoindex(0,1);
			add5thtoindex(samples/2,-1);
			for (unsigned i = 1; i<samples/2;i++){
				setGDDtoindex(i,1);
				setGDDtoindex(samples-i,-1);
				addTODtoindex(i,1);
				addTODtoindex(samples-i,-1);
				add4thtoindex(i,1);
				add4thtoindex(samples-i,-1);
				add5thtoindex(i,1);
				add5thtoindex(samples-i,-1);
			}
			return *this;
		}
		PulseFreq & addchirp(std::vector<float> & chirp_in) {
			if (intime){
				std::cerr << "whoops, trying to add phase in the time domain\n" << std::flush;
				return *this;
			}
			assert(chirp_in.size()==4);
			phase_GDD = chirp_in[0];
			phase_TOD = chirp_in[1];
			phase_4th = chirp_in[2];
			phase_5th = chirp_in[3];
			addAllPhase();
			return *this;
		}
		PulseFreq * mulamp(){
			HERE
		}
		PulseFreq & mulamp(std::vector<float> & chirp_in);
		PulseFreq & addAllPhase(void);
		PulseFreq & mulAllAmp(void);

		PulseFreq & attenuate(float attenfactor);
		PulseFreq & delay(float delayin); // expects delay in fs
		PulseFreq & setdelay(float delayin); // expects delay in fs
		PulseFreq & phase(float phasein); // expects phase in units of pi , i.e. 1.0 = pi phase flip 
		PulseFreq & rebuildvectors(float newgain); // rebuilds to new amplitude and (eventually) center omega and rollonoff and spectral envelope.

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

		PulseFreq & modulateamp_freq(const std::vector<double> & modulation);
		PulseFreq & modulatephase_freq(const std::vector<double> & modulation);

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

		void print_amp(void);
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
		float m_gain;
		unsigned m_saturate;

		size_t sampleround;

		bool intime,infreq;
		unsigned samples;
		unsigned startind,stopind,onwidth,offwidth;
		float tspan;
		float domega,lambda_center,lambda_width,omega_center,omega_width,omega_high;
		float omega_onwidth,omega_offwidth;
		float phase_GDD,phase_TOD,phase_4th,phase_5th;
		float amp_0th,amp_1st,amp_2nd,amp_3rd,amp_4th,amp_5th;
		float dtime,time_center,time_wdith;

		// FFTW variables //
		fftw_plan * FTplan_forwardPtr;
		fftw_plan * FTplan_backwardPtr;
		fftw_plan * FTplan_r2hcPtr;
		fftw_plan * FTplan_hc2rPtr;
		fftw_plan * FTplan_r2hc_2xPtr;
		fftw_plan * FTplan_hc2r_2xPtr;

		std::complex<double> * cvec; // this is still fftw_malloc() for sake of fftw memory alignment optimization
		std::int32_t * ovec; // this is still fftw_malloc() for sake of fftw memory alignment optimization
		double * r_vec; // this will get fftw_malloc() for sake of fftw memory alignment 
		double * hc_vecFT; // this will get fftw_malloc() for sake of fftw memory alignment 
		double * r_vec_2x; // this will get fftw_malloc() for sake of fftw memory alignment 
		double * hc_vec_2xFT; // this will get fftw_malloc() for sake of fftw memory alignment 

		std::vector<double> rhovec;
		std::vector<double> modamp;
		std::vector<double> omega;
		std::vector<double> time;
		std::vector<double> modphase;
		std::vector<double> phivec;

		float nu0;

		unsigned i_low, i_high;

		PulseFreq & rhophi2cvec(const size_t i){ cvec[i] = std::polar(rhovec[i],phivec[i]);return *this;}
		PulseFreq & cvec2rhophi(const size_t i){ rhovec[i] = std::abs(cvec[i]); phivec[i] = std::arg(cvec[i]);return *this; }
		PulseFreq & rhophi2cvec(void);
		PulseFreq & cvec2rhophi(void);

		std::random_device rng;

		PulseFreq & buildvectors(const size_t s);
		PulseFreq & factorization(void);
		PulseFreq & killvectors(void);

		PulseFreq & setGDDtoindex(const unsigned indx,const int omega_sign);
		PulseFreq & setTODtoindex(const unsigned indx,const int omega_sign);
		PulseFreq & set4thtoindex(const unsigned indx,const int omega_sign);
		PulseFreq & set5thtoindex(const unsigned indx,const int omega_sign);
		PulseFreq & addGDDtoindex(const unsigned indx,const int omega_sign);
		PulseFreq & addTODtoindex(const unsigned indx,const int omega_sign);
		PulseFreq & add4thtoindex(const unsigned indx,const int omega_sign);
		PulseFreq & add5thtoindex(const unsigned indx,const int omega_sign);


		PulseFreq &  modampatindx(const unsigned indx,const std::vector<float> & modvec) {
			rhovec[indx] *= modvec[indx];
			rhophi2cvec(indx);
			return *this;
		}
		PulseFreq & modampatindx(const unsigned indx) {
			rhovec[indx] *= modamp[indx];
			rhophi2cvec(indx);
			return *this;
		}
		/*
		   E(t) = E0exp(iwt)=E0 exp(i w(t) (t-t0))
		   = E0 exp(i (w0+GDD*t) (t-t0) )
		   = E0 exp(i ((w0 * t) - (GDD*t*t0) + (GDD*t**2) - (w0*t0)) so adding a delay means adding an arbitrary phase of w0*t0 and subtracting GDD*t*t0 ( a linear term in phase which we know only dials the pulse position in z, or likewise phase in time... unless t0 is actually t0(t), then you have to include this one as well as w0*t0(t).
		   = E0 exp(i ((w0 * t) + (GDD*t**2) - (w0*t0) - (GDD*t*t0) ) )
		   = E0 exp(i w(t)*t) exp(i -w(t)*t0(t))
		   = E(t) exp(i -(w0 + 1/GDD*t + 1/TOD*t**2 + 1/FOD*t**3)*t0(t))
		 */
		PulseFreq &  modphaseatindx(const unsigned indx,const std::vector<float> & modvec) {
			if (modvec[indx]!=0){
				phivec[indx] += modvec[indx];
				rhophi2cvec(indx);
			}
			return *this;
		}
		PulseFreq & modphaseatindx(const unsigned indx) {
			if (modphase[indx] != 0){
				phivec[indx] += modphase[indx];
				rhophi2cvec(indx);
			}
			return *this;
		}
		float rising(const unsigned indx) {
			return std::pow(sin(float(Constants::half_pi<float>()*(indx-startind)/onwidth)),int(2));
		}
		float falling(const unsigned indx) {
			return std::pow(sin(float(Constants::half_pi<float>()*(stopind-indx)/offwidth)),int(2));
		}
};

#endif
