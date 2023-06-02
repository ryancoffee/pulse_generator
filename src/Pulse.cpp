// Pulse clas implimentation
// standard includes
#include <cmath>
#include <iterator>


// my headers
#include <Pulse.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <fftw3.h>
#include <algorithm>
#include <Constants.hpp>
#include <DataOps.hpp>
#include <random>
#include <cassert>
#include <cstdint> // for appendwavelegth() and int32_t
#include <limits> // for std::numeric_limits<short int>::max() and min()

using namespace Constants;
using namespace DataOps;

PulseFreq::PulseFreq(const double omcenter_in=(0.55*fsPau<double>()),const double omwidth_in=(0.15*fsPau<double>()),const double omonoff_in=(0.1*fsPau<double>()), double tspan_in=(10000.0/fsPau<double>())):
	omega_center(omcenter_in),
	omega_width(omwidth_in ),
	omega_high( std::max(4.0*(omcenter_in + omwidth_in),10.0*omcenter_in) ),
	domega( 2.0*pi<double>()/tspan_in),
	intime(false),
	infreq(true),
	m_noisescale(1e-3),
	m_sampleinterval(2),
	m_saturate(1<<12),
	m_gain(1000000),
	m_lamsamples(1<<10),
	sampleround(1000),
	cvec(NULL),
	r_vec(NULL),
	hc_vecFT(NULL),
	r_vec_2x(NULL),
	hc_vec_2xFT(NULL)
{
	std::cout << "In constructor PulseFreq()" << std::endl;
	i_low =  (unsigned)(double( atof( getenv("nu_low") ) )* twopi<double>()*fsPau<double>()/domega);
	i_high =  (unsigned)(double( atof( getenv("nu_high") ) )* twopi<double>()*fsPau<double>()/domega);
	samples = (( (unsigned)(2.0 * omega_high / domega))/sampleround + 1 ) *sampleround;// dt ~ .1fs, Dt ~ 10000fs, dom = 2*pi/1e4, omhigh = 2*pi/.1fs, samples = omhigh/dom*2
	dtime = tspan_in/double(samples);
	// no longer forcing gaussian sin2
	// omega_onwidth = omega_offwidth = omega_width/2.0; // forcing sin2 gaussian spectrum
	omega_onwidth = omonoff_in;
	omega_offwidth = omonoff_in;

	buildvectors(samples);
	nu0=omcenter_in/(2.0*pi<double>())*fsPau<double>();
	phase_GDD=phase_TOD=phase_4th=phase_5th=0.0;
	m_lamsamples = (size_t)atoi(getenv("lamsamples"));
	m_gain = (float)atof(getenv("gain"));
	m_noisescale = (double)atof(getenv("noisescale") ) ;
	m_sampleinterval = (size_t)atoi(getenv("sampleinterval"));
	m_saturate = uint16_t( atoi( getenv("saturate")));
	std::cout << "exiting constructor PulseFreq()" << std::endl;
}


PulseFreq::PulseFreq(const PulseFreq &rhs) // deep-ish copy constructor
	:samples(rhs.samples)
	,omega_center(rhs.omega_center)
	,omega_width(rhs.omega_width)
	,omega_high(rhs.omega_high)
	,omega_onwidth(rhs.omega_onwidth)
	,omega_offwidth(rhs.omega_offwidth)
	,domega(rhs.domega)
	,intime(rhs.intime)
	,infreq(rhs.infreq)
	,i_low(rhs.i_low) 
	,i_high(rhs.i_high)
	,m_noisescale(rhs.m_noisescale)
	,m_sampleinterval(rhs.m_sampleinterval)
	,m_saturate(rhs.m_saturate)
	,m_gain(rhs.m_gain)
	,m_lamsamples(rhs.m_lamsamples)
	,sampleround(1000)
{
	//std::cerr << "\t\t\t+++++  Copy constructor of PulseFreq::PulseFreq(PulseFreq &rhs)\t\tsamples = " << samples << "\n" << std::flush;
	DataOps::clone(omega,rhs.omega);
	DataOps::clone(time,rhs.time);

	startind = rhs.startind;stopind=rhs.stopind;onwidth=rhs.onwidth;offwidth=rhs.offwidth;
	tspan = rhs.tspan;
	lambda_center=rhs.lambda_center;lambda_width=rhs.lambda_width;
	phase_GDD=rhs.phase_GDD;phase_TOD=rhs.phase_TOD;phase_4th=rhs.phase_4th;phase_5th=rhs.phase_5th;

	dtime = rhs.dtime;time_center=rhs.time_center;time_wdith=rhs.time_wdith;


	nu0=rhs.nu0;
	FTplan_forwardPtr = rhs.FTplan_forwardPtr; 
	FTplan_backwardPtr = rhs.FTplan_backwardPtr; 
	// HERE HERE HERE HERE  is the seg fault //
	//std::cerr << "\tEntering buildvectors\t" << std::flush;
	buildvectors(samples);
	//std::cerr << "\tleft buildvectors\t" << std::flush;

	DataOps::clone(rhovec,rhs.rhovec);
	DataOps::clone(phivec,rhs.phivec);
	DataOps::clone(cvec,rhs.cvec,samples);
	DataOps::clone(r_vec,rhs.r_vec,samples);
	DataOps::clone(hc_vecFT,rhs.hc_vecFT,samples);
	DataOps::clone(r_vec_2x,rhs.r_vec_2x,2*samples);
	DataOps::clone(hc_vec_2xFT,rhs.hc_vec_2xFT,2*samples);
	DataOps::clone(modamp,rhs.modamp);
	DataOps::clone(modphase,rhs.modphase);
	//std::cerr << "\t leaving Copy constructor of PulseFreq::PulseFreq(PulseFreq &rhs)\n" << std::flush;
}

PulseFreq & PulseFreq::operator=(const PulseFreq & rhs) // shallow-ish assignment
{
	//std::cerr << "\n\n\t\t########### !!!!!!! copying into a PulseFreq, make sure same enumber of samples !!!!!! ##########\n\n" << std::flush;
	//std::cerr << "\t\t\t+++++  assignment copy of PulseFreq::operator= with " << samples << " samples and " << rhs.getsamples() << " on rhs \n" << std::flush;
	samples=rhs.samples;
	omega_center=rhs.omega_center;
	omega_width=rhs.omega_width;
	omega_high=rhs.omega_high;
	omega_onwidth=rhs.omega_onwidth;
	omega_offwidth=rhs.omega_offwidth;
	domega=rhs.domega;
	intime=rhs.intime;
	infreq=rhs.infreq;
	omega=rhs.omega; // these are static vectors... i'm trying not to make copies just yet
	time=rhs.time; // these are static vectors... i'm trying not to make copies just yet
	i_low=rhs.i_low; 
	i_high=rhs.i_high; 
	m_noisescale=rhs.m_noisescale;
	m_sampleinterval=rhs.m_sampleinterval;
	m_saturate=rhs.m_saturate;
	m_gain=rhs.m_gain;
	m_lamsamples=rhs.m_lamsamples;

	startind = rhs.startind;stopind=rhs.stopind;onwidth=rhs.onwidth;offwidth=rhs.offwidth;
	tspan = rhs.tspan;
	lambda_center=rhs.lambda_center;lambda_width=rhs.lambda_width;
	phase_GDD=rhs.phase_GDD;phase_TOD=rhs.phase_TOD;phase_4th=rhs.phase_4th;phase_5th=rhs.phase_5th;

	dtime = rhs.dtime;time_center=rhs.time_center;time_wdith=rhs.time_wdith;

	nu0=rhs.nu0;
	FTplan_forwardPtr = rhs.FTplan_forwardPtr; 
	FTplan_backwardPtr = rhs.FTplan_backwardPtr; 

	DataOps::clone(rhovec,rhs.rhovec);
	DataOps::clone(phivec,rhs.phivec);
	DataOps::clone(cvec,rhs.cvec,samples);
	DataOps::clone(r_vec,rhs.r_vec,samples);
	DataOps::clone(hc_vecFT,rhs.hc_vecFT,samples);
	DataOps::clone(r_vec_2x,rhs.r_vec_2x,2*samples);
	DataOps::clone(hc_vec_2xFT,rhs.hc_vec_2xFT,2*samples);

	DataOps::clone(modamp,rhs.modamp);
	DataOps::clone(modphase,rhs.modphase);
	return *this;
}

PulseFreq::~PulseFreq(void){
	killvectors();
}


void PulseFreq::rhophi2cvec(void)
{
	for (size_t i=0;i<samples;i++){
		cvec[i] = std::polar(rhovec[i],phivec[i]);
	}
}
void PulseFreq::cvec2rhophi(void)
{
	for (size_t i=0;i<samples;i++){
		rhovec[i] = std::abs(cvec[i]);
		phivec[i] = std::arg(cvec[i]);
	}
}

PulseFreq & PulseFreq::operator+=(const PulseFreq &rhs){
	DataOps::sum(cvec,rhs.cvec,samples);
	cvec2rhophi();
	return *this;
}

PulseFreq & PulseFreq::operator-=(const PulseFreq &rhs){
	DataOps::diff(cvec,rhs.cvec,samples);
	cvec2rhophi();
        return *this;
}

PulseFreq & PulseFreq::operator*=(const PulseFreq &rhs){
	DataOps::mul(cvec,rhs.cvec,samples);
	cvec2rhophi();
	return *this;
}

PulseFreq & PulseFreq::operator*=(const double s){
	DataOps::mul(cvec,s,samples);
	cvec2rhophi();
	return *this;
}

PulseFreq & PulseFreq::operator/=(const PulseFreq &rhs){
	DataOps::div(cvec,rhs.cvec,samples);
	cvec2rhophi();
        return *this;
}

PulseFreq & PulseFreq::interfere(const PulseFreq &rhs){
	*this += rhs;
        return *this;
}

PulseFreq & PulseFreq::diffamps(const PulseFreq &rhs){
	rhovec -= rhs.rhovec;
	std::fill(phivec.begin(),phivec.end(),0);
	rhophi2cvec();
        return *this;
}

PulseFreq & PulseFreq::normamps(const PulseFreq &rhs){
	rhovec /= rhs.rhovec;
	std::fill(phivec.begin(),phivec.end(),0);
	rhophi2cvec();
        return *this;
}

void PulseFreq::print_amp(std::ofstream & outfile)
{
	outfile << "# amp\n";
	outfile << rhovec << std::endl;
}
void PulseFreq::print_phase(std::ofstream & outfile)
{
	outfile << "# phase\n";
	outfile << phivec << std::endl;
}
void PulseFreq::print_phase_powerspectrum(std::ofstream & outfile)
{
/*
	double * phase = (double *) fftw_malloc(sizeof(double) * samples);
	double * phaseFT = (double *) fftw_malloc(sizeof(double) * samples);
	fftw_plan plan_r2hc = fftw_plan_r2r_1d(samples,
			phase,
			phaseFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
*/
	std::copy(phivec.begin(),phivec.end(),r_vec);
	fftw_execute_r2r(*FTplan_r2hcPtr.get(),r_vec,hc_vecFT);
	

	outfile << "# power spectrum of the Fourier phase\n";
	outfile << std::pow(hc_vecFT[0],int(2)) << "\n";
	for (size_t i = 1; i<samples/2;++i){
		outfile << std::pow(hc_vecFT[i],int(2)) + std::pow(hc_vecFT[samples-i],int(2)) << "\n";
	}
	outfile << std::pow(hc_vecFT[samples/2],int(2)) << std::endl;
	outfile << std::endl;
}

bool PulseFreq::addrandomphase(void)
{
	if (!infreq){
		std::cerr << "died here at addrandomphase()" << std::endl;
		return false;
	}
	size_t sz = samples*2; // doubling the vector to mirror it so that DFT handles the phase well

	double * randphase = (double *) fftw_malloc(sizeof(double) * sz);
	double * randphaseFT = (double *) fftw_malloc(sizeof(double) * sz);
/*

	fftw_plan plan_r2hc = fftw_plan_r2r_1d(sz,
			randphase,
			randphaseFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
	fftw_plan plan_hc2r = fftw_plan_r2r_1d(sz,
			randphaseFT,
			randphase,
			FFTW_HC2R,
			FFTW_MEASURE
			);

*/


	std::uniform_real_distribution<double> distribution(
		(double(atof(getenv("randphase_mean")))-double(atof(getenv("randphase_std"))))*Constants::pi<double>(),
		(double(atof(getenv("randphase_mean")))+double(atof(getenv("randphase_std"))))*Constants::pi<double>()
		);

	double phase = distribution(rng);
	randphase[0] = phase;
	randphase[samples/2] = -phase;
	for (size_t i = 1; i<samples/2;i++){
		phase = distribution(rng);
		randphase[i] = phase;
		randphase[samples-i] = -phase;
	}
	for (size_t i=sz-1;i>sz/2-1;--i){
		randphase[i] = randphase[sz-i];
	}

	size_t lowpass = (size_t)(atoi(getenv("phaseNoiseLowpass")));
	std::cerr << "\n======== lowpass is " << lowpass << " =======\n" << std::flush;

	fftw_execute_r2r(*FTplan_r2hc_2xPtr.get(),randphase,randphaseFT);
	std::fill(randphaseFT+lowpass,randphaseFT+sz-lowpass,0.);
	for (size_t i=1;i<lowpass;++i){
		double filter = std::pow(std::cos(double(i)/(double(lowpass)) * Constants::half_pi<double>() ),int(2));
		randphaseFT[i] *= filter;
		randphaseFT[sz-i] *= filter;
	}
	randphaseFT[sz/2] = 0.;
	fftw_execute_r2r(*FTplan_hc2r_2xPtr.get(),randphaseFT,randphase);

	for (size_t i=0;i<samples;++i){
		phivec[i] += randphase[i]/samples;
	}
	rhophi2cvec();
	return true;
}


void PulseTime::setstrength(const double in)
{
  strength = in * auenergy<double>()/Eh<double>() * std::pow(aufor10PW<double>(),int(2));
}

void PulseTime::setwidth(const double in)
{
  Ctau = in * Constants::root_pi<double>()/ Constants::fsPau<double>() / 2.0;
}

void PulseTime::sett0(const double in)
{
  t0 = in / fsPau<double>();
}

void PulseFreq::attenuate(double attenfactor){
	rhovec *= attenfactor;
	rhophi2cvec();
}
void PulseFreq::phase(double phasein){ // expects delay in units of pi , i.e. 1.0 = pi phase flip 
	if(intime){
		fft_tofreq();
	}
	phivec += phasein*Constants::pi<double>();
	rhophi2cvec();
	if(intime){
		fft_totime();
	}
}
void PulseFreq::delay(double delayin){ // expects delay in fs
	if(intime){
		fft_tofreq();
	}
	for (unsigned i=0;i<samples;i++){
		phivec[i] += omega[i]*delayin/fsPau<double>();
	}
	rhophi2cvec();
	if(intime){
		fft_totime();
	}
}


void PulseFreq::printfrequency(std::ofstream * outfile){
	double nu,lambda;
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lambda = C_nmPfs<double>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << std::pow(rhovec[i],int(2)) << "\t" << phivec[i] << "\n";
	}
}
void PulseFreq::printwavelengthbins(std::ofstream * outfile)
{
	std::vector<double> x(2);
	x.front() = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_low];
	x.back() = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_high-1];
	double dlam = (x.front()-x.back())/double(m_lamsamples);
        for (size_t i = 0;i<m_lamsamples;++i){
		(*outfile) << x.back() + i*dlam << "\t";
        }
	(*outfile) << "\n";
	return;
}
void PulseFreq::appendwavelength(std::ofstream * outfile)
{
	std::vector<double> x(i_high-i_low);
	std::vector<double> y(i_high-i_low);	
	for (size_t i=0;i<y.size();++i){
		x[i] = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_low+i];
		//y[i] = std::pow(rhovec[i_low+i],int(2)) * 200000000000;
		y[i] = std::min(std::pow(rhovec[i_low+i],int(2)) * m_gain,double(m_saturate));
	}
	double dlam = (x.front()-x.back())/double(m_lamsamples);
	boost::math::interpolators::barycentric_rational<double> interpolant(x.data(), y.data(), y.size());
	for (size_t i=0;i<m_lamsamples;++i){
		(*outfile) << uint16_t(interpolant(x.back()+i*dlam)) << "\t";
	}
	(*outfile) << std::endl;
	return;
}
void PulseFreq::appendwavelength_deriv(std::ofstream * outfile)
{
	std::vector<double> x(i_high-i_low);
	std::vector<double> y(i_high-i_low);	
	std::vector<int16_t> resultvec(m_lamsamples);
	std::vector<double> diffvec(m_lamsamples);

	for (size_t i=0;i<y.size();++i){
		x[i] = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_low+i];
		y[i] = std::pow(rhovec[i_low+i],int(2)) * m_gain;
	}
	double dlam = (x.front()-x.back())/double(m_lamsamples);
	boost::math::interpolators::barycentric_rational<double> interpolant(x.data(), y.data(), y.size());
	for (size_t i=0;i<m_lamsamples;++i){
		diffvec[i] = interpolant(x.back()+(i+1)*dlam) - interpolant(x.back()+i*dlam);
		//diffvec[i] -= interpolant(x.back()+(i+10)*dlam) - interpolant(x.back()+(i+11)*dlam);
	}
        double scale;
        int16_t max = std::numeric_limits<int16_t>::max();
        int16_t min = std::numeric_limits<int16_t>::min();
	auto bounds = std::minmax_element(diffvec.begin(),diffvec.end());
        if (*bounds.second>std::abs(*bounds.first)){
                scale = max/ *bounds.second;
        } else {
                scale = min/ *bounds.first;
        }
        std::transform(diffvec.begin(),diffvec.end(),resultvec.begin(),[scale](double x){return int16_t(scale*x);});
        *(outfile) << resultvec;

	return;
}
void PulseFreq::appendwavelength_bin(std::ofstream * outfile)
{
	std::vector<double> x(i_high-i_low);
	std::vector<double> y(i_high-i_low);	
	for (size_t i=0;i<y.size();++i){
		x[i] = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_low+i];
		//y[i] = std::pow(rhovec[i_low+i],int(2)) * 200000000000;
		y[i] = std::min(std::pow(rhovec[i_low+i],int(2)) * m_gain,double(m_saturate));
	}
	double dlam = (x.front()-x.back())/double(m_lamsamples);
	boost::math::interpolators::barycentric_rational<double> interpolant(x.data(), y.data(), y.size());
	for (size_t i=0;i<m_lamsamples;++i){
		(*outfile) << int32_t(interpolant(x.back()+i*dlam));
	}
	return;
}
void PulseFreq::appendfrequency(std::ofstream * outfile){
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		uint16_t val = std::min(uint16_t(rhovec[i] * m_gain),uint16_t(m_saturate));
       		(*outfile) << std::pow(val,int(2)) << "\t";
        }
	(*outfile) << std::endl;
}

void PulseFreq::appendnoisy(std::ofstream * outfile){
	std::normal_distribution<double> norm_dist( 0.0, m_noisescale);
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		double outval = std::pow(rhovec[i],int(2)) + norm_dist(rng);
       		(*outfile) << outval << "\t" ;
        }
	(*outfile) << std::endl;
}

void PulseFreq::printfrequencybins(std::ofstream * outfile){
	double nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lam = C_nmPfs<double>()/nu;
       		(*outfile) << nu << "\t" << lam << "\n";
        }
	(*outfile) << "\n";
}
void PulseFreq::appendfrequencybins(std::ofstream * outfile){
	double nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
       		(*outfile) << nu << "\t";
        }
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelay(std::ofstream * outfile, const double *delay){
	double nu,lambda,thisdelay;
	thisdelay = (*delay)*fsPau<double>();
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lambda = C_nmPfs<double>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec[i] << "\t" << phivec[i] << "\t" << thisdelay << "\n";
	}
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelaychirp(std::ofstream * outfile, const double *delay,const double *chirp){
	double nu,lambda,thisdelay,thischirp,reltime;
	thisdelay = (*delay)*fsPau<double>();
	thischirp = (*chirp)*std::pow(fsPau<double>(),int(2));
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\tchirp[fs^2]\treldelays[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lambda = C_nmPfs<double>()/nu;
		reltime = (nu-nu0)*thischirp*10.0;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec[i] << "\t" << phivec[i] << "\t" << thisdelay << "\t" << thischirp << "\t" << reltime << "\n";
	}
	(*outfile) << "\n";
}

void PulseFreq::printwavelength(std::ofstream * outfile,const double *delay){
	double nu,lambda,lamlast;
	lamlast=0;
        for (unsigned i = i_high;i>i_low;i-=(unsigned)(atoi(getenv("sampleinterval")))){
                nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
                lambda = (double)((int)((C_nmPfs<double>()/nu)*10.0))/10.0;
		if (lambda>lamlast & (int)(lambda*10)%5 == 0){
                	(*outfile) << lambda << "\t" << (*delay) << "\t" << std::pow(rhovec[i],int(2)) << "\n";
			lamlast = lambda;
		}
        }
	(*outfile) << "\n";
}

void PulseFreq::printtime(std::ofstream * outfile){
	(*outfile) << ("#time\treal\timag\n");
	for (unsigned i = samples/2;i<samples;i++){
		(*outfile) << (time[i]*fsPau<double>()) << "\t" << cvec[i].real() << "\t" << cvec[i].imag() << "\n";
	}
	for (unsigned i = 0;i<samples/2;i++){
		(*outfile) << (time[i]*fsPau<double>()) << "\t" << cvec[i].real() << "\t" << cvec[i].imag() << "\n";
	}

} 


void PulseFreq::buildvectors(const size_t s){
	//std::cerr << "allocating with fftw_malloc with samples = " << samples << std::endl << std::flush;
	cvec = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * size_t(s));
        std::fill(cvec,cvec + samples,std::complex<double>(0));
	r_vec = (double *) fftw_malloc(sizeof(double) * size_t(s));
        std::fill(r_vec,r_vec + s,double(0));
	//std::cerr << "\t\t...allocated with fftw_malloc with samples = " << samples << std::endl << std::flush;
	hc_vecFT = (double *) fftw_malloc(sizeof(double) * size_t(s));
        std::fill(hc_vecFT,hc_vecFT + s,double(0));
	//std::cerr << "allocating with fftw_malloc with samples = " << (2*samples) << std::endl << std::flush;
	r_vec_2x = (double *) fftw_malloc(sizeof(double) * size_t(s) * 2);
        //std::fill(r_vec_2x,r_vec_2x + 2*samples,double(0));
	hc_vec_2xFT = (double *) fftw_malloc(sizeof(double) * size_t(s) * 2);
        std::fill(hc_vec_2xFT,hc_vec_2xFT + 2*s,double(0));
	//std::cerr << "\t\t...allocated with fftw_malloc with 2*samples = " << (2*samples) << std::endl << std::flush;

	rhovec.resize(s,0.0);
	phivec.resize(s,0.0);
	modamp.resize(s,1.0);
	modphase.resize(s,0.0);
	omega.resize(s);
	time.resize(s);

	omega[0] = 0.0;
	time[0] = 0.0;
	omega[s/2] = -(double)(s/2)*domega;
	time[s/2] = -(double)(s/2)*dtime;

	startind = (unsigned)((omega_center-(omega_width/2.0))/domega);
	stopind = (unsigned)((omega_center+(omega_width/2.0))/domega);
	onwidth = (unsigned)(omega_onwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	offwidth = (unsigned)(omega_offwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	/*
	std::vector<unsigned> temp(4);
	temp = {startind,onwidth,stopind,offwidth};
	std::cerr << "(startind,onwidth,stopind,offwidth) = ( " << temp << " )\n" << std::flush;
	*/

	for (unsigned i = 1; i<startind;i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
	}
	for (unsigned i = startind;i<startind+onwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = rising(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[s-i] = -domega*(double)i;
		time[s-i] = -dtime*(double)i;
		rhovec[s-i] = rising(i);
		cvec[s-i] = std::polar(rhovec[s-i],phivec[s-i]);
	}
	//std::cerr << "Made it to HERE HERE HERE with startind+onwidth = "  << (startind+onwidth) << " and stopind-offwidth = " << (stopind-offwidth) << "\n" << std::flush;
	for (unsigned i = startind+onwidth;i<stopind-offwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = 1.0;
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[s-i] = -domega*(double)i;
		time[s-i] = -dtime*(double)i;
		rhovec[s-i] = 1.0;
		cvec[s-i] = std::polar(rhovec[s-i],phivec[s-i]);
	}
	//std::cerr << "Trying to make it to HERE NOW\n" << std::flush;
	for (unsigned i = stopind-offwidth;i<stopind; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = falling(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[s-i] = -domega*i;
		time[s-i] = -dtime*i;
		rhovec[s-i] = falling(i);
		cvec[s- i] = std::polar(rhovec[s- i],phivec[s- i]);
	}
	for (unsigned i = stopind;i<s/2; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[s-i] = -domega*i;
		time[s-i] = -dtime*i;
	}
	
}
void PulseFreq::killvectors(void){
	fftw_free(cvec);
	fftw_free(r_vec);
	fftw_free(hc_vecFT);
	fftw_free(r_vec_2x);
	fftw_free(hc_vec_2xFT);
	cvec = NULL;
	r_vec = hc_vecFT = r_vec_2x = hc_vec_2xFT = NULL;
}

void PulseFreq::setplans(const PulseFreq & rhs)
{
	FTplan_forwardPtr = rhs.FTplan_forwardPtr;
	FTplan_backwardPtr = rhs.FTplan_backwardPtr;
}
void PulseFreq::setmasterplans(fftw_plan * const forward,fftw_plan * const backward)
{
	assert(FTplan_forwardPtr.use_count()==0 && FTplan_backwardPtr.use_count()==0);
	*forward = fftw_plan_dft_1d(samples, 
			reinterpret_cast<fftw_complex*>(cvec),
			reinterpret_cast<fftw_complex*>(cvec), 
			FFTW_FORWARD, FFTW_ESTIMATE);
	*backward = fftw_plan_dft_1d(samples, 
			reinterpret_cast<fftw_complex*>(cvec), 
			reinterpret_cast<fftw_complex*>(cvec), 
			FFTW_BACKWARD, FFTW_ESTIMATE);
	FTplan_forwardPtr = std::make_shared<fftw_plan> (*forward);
	FTplan_backwardPtr = std::make_shared<fftw_plan> (*backward);
}
void PulseFreq::setancillaryplans(fftw_plan * const r2hc,fftw_plan * const hc2r,fftw_plan * const r2hc_2x,fftw_plan * const hc2r_2x)
{

	assert(FTplan_r2hcPtr.use_count()==0
			&& FTplan_hc2rPtr.use_count()==0
			&& FTplan_r2hc_2xPtr.use_count()==0
			&& FTplan_hc2r_2xPtr.use_count()==0);
	*r2hc = fftw_plan_r2r_1d(samples,
			r_vec,
			hc_vecFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
	*hc2r = fftw_plan_r2r_1d(samples,
			hc_vecFT,
			r_vec,
			FFTW_HC2R,
			FFTW_MEASURE
			);
	*r2hc_2x = fftw_plan_r2r_1d(2*samples,
			r_vec_2x,
			hc_vec_2xFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
	*hc2r_2x = fftw_plan_r2r_1d(2*samples,
			hc_vec_2xFT,
			r_vec_2x,
			FFTW_HC2R,
			FFTW_MEASURE
			);
	FTplan_r2hcPtr = std::make_shared<fftw_plan> (*r2hc);
	FTplan_hc2rPtr = std::make_shared<fftw_plan> (*hc2r);
	FTplan_r2hc_2xPtr = std::make_shared<fftw_plan> (*r2hc_2x);
	FTplan_hc2r_2xPtr = std::make_shared<fftw_plan> (*hc2r_2x);

}
