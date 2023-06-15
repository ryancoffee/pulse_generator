// Pulse clas implimentation
// standard includes
#include <cmath>
#include <iterator>

#include <boost/math/interpolators/barycentric_rational.hpp>
#include <algorithm>
#include <random>
#include <cassert>
#include <cstdint> // for appendwavelegth() and int32_t
#include <limits> // for std::numeric_limits<short int>::max() and min()
#include <functional>

// my headers
#include <Pulse.hpp>

using namespace Constants;
using namespace DataOps;


/*
 *
 *
 *
 *  PulseTime
 *
 *
 *
 */


PulseTime & PulseTime::setstrength(const float in)
{
  strength = in * auenergy<float>()/Eh<float>() * std::pow(aufor10PW<float>(),int(2));
  return *this;
}

PulseTime & PulseTime::setwidth(const float in)
{
  Ctau = in * Constants::root_pi<float>()/ Constants::fsPau<float>() / 2.0;
  return *this;
}

PulseTime & PulseTime::sett0(const float in)
{
  t0 = in / Constants::fsPau<float>();
  return *this;
}





/*
 *
 *
 *
 *
 * PulseFreq
 *
 *
 *
 *
 */

PulseFreq::PulseFreq(Params &params)
:omega_center(params.omega0())
,omega_width(params.omega_width())
,omega_high( std::max(4.0*(params.omega0() + params.omega_width()),10.0*params.omega0()) )
,domega( 2.0*pi<float>()/params.getTspan())
,intime(false)
,infreq(true)
,m_noisescale(1e-3)
,m_sampleinterval(2)
,m_saturate(1<<12)
,m_gain(1)
,m_lamsamples(1<<10)
,sampleround(1000)
,cvec(NULL)
,r_vec(NULL)
,hc_vecFT(NULL)
,r_vec_2x(NULL)
,hc_vec_2xFT(NULL)
{
	i_low =  (size_t)(params.nu_low * twopi<float>()*fsPau<float>()/domega);
	i_high =  (size_t)(params.nu_high * twopi<float>()*fsPau<float>()/domega);
	samples = (( (size_t)(2.0 * omega_high / domega))/sampleround + 1 ) *sampleround;// dt ~ .1fs, Dt ~ 10000fs, dom = 2*pi/1e4, omhigh = 2*pi/.1fs, samples = omhigh/dom*2
	dtime = params.getTspan()/float(samples);
	// no longer forcing gaussian sin2
	// omega_onwidth = omega_offwidth = omega_width/2.0; // forcing sin2 gaussian spectrum
	omega_onwidth = params.omega_onoff();
	omega_offwidth = params.omega_onoff();

	buildvectors(samples);
	nu0=omega_center/(2.0*pi<float>())*fsPau<float>();
	phase_GDD=phase_TOD=phase_4th=phase_5th=0.0;
	amp_0th = 1.0;
	amp_1st = amp_2nd = amp_3rd = amp_4th = amp_5th = 0.0;
	m_lamsamples = params.lamsamples;
	m_gain = params.gain;
	m_noisescale = params.noisescale;
	m_sampleinterval = params.sampleinterval;
	m_saturate = params.saturate;
	//std::cout << "exiting constructor PulseFreq()" << std::endl;
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
	std::cerr << "\t\t\t+++++  Copy constructor of PulseFreq::PulseFreq(PulseFreq &rhs)\t\tsamples = " << samples << "\n" << std::flush;
	DataOps::clone(omega,rhs.omega);
	DataOps::clone(time,rhs.time);

	startind = rhs.startind;stopind=rhs.stopind;onwidth=rhs.onwidth;offwidth=rhs.offwidth;
	tspan = rhs.tspan;
	lambda_center=rhs.lambda_center;lambda_width=rhs.lambda_width;
	phase_GDD=rhs.phase_GDD;phase_TOD=rhs.phase_TOD;phase_4th=rhs.phase_4th;phase_5th=rhs.phase_5th;
	amp_0th = rhs.amp_0th;amp_1st = rhs.amp_1st;amp_2nd = rhs.amp_2nd;amp_3rd = rhs.amp_3rd;amp_4th = rhs.amp_4th;amp_5th = rhs.amp_5th;
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

	return *this;
}

PulseFreq::~PulseFreq(void){
	killvectors();
	std::cerr << "Leaving PulseFreq::~PulseFreq()\n" << std::flush;
}


PulseFreq & PulseFreq::rhophi2cvec(void)
{
	for (size_t i=0;i<samples;i++){
		cvec[i] = std::polar(rhovec[i],phivec[i]);
	}
	return *this;
}
PulseFreq & PulseFreq::cvec2rhophi(void)
{
	for (size_t i=0;i<samples;i++){
		rhovec[i] = std::abs(cvec[i]);
		phivec[i] = std::arg(cvec[i]);
	}
	return *this;
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

PulseFreq & PulseFreq::operator*=(const float s){
	DataOps::mul(cvec,(double)s,samples);
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

void PulseFreq::print_amp(void)
{
	fft_totime();
	std::cout << "# amp**2\n";
	for (size_t i=0;i<rhovec.size();i++)
		std::cout << std::real(cvec[i]) << std::endl;
	std::cout << "# done\n\n" << std::flush;
	fft_tofreq();
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
	std::copy(phivec.begin(),phivec.end(),r_vec);
	fftw_execute_r2r(*FTplan_r2hcPtr,r_vec,hc_vecFT);
	

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

	std::uniform_real_distribution<float> distribution(
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

	fftw_execute_r2r(*FTplan_r2hc_2xPtr,randphase,randphaseFT);
	std::fill(randphaseFT+lowpass,randphaseFT+sz-lowpass,0.);
	for (size_t i=1;i<lowpass;++i){
		double filter = std::pow(std::cos(double(i)/(double(lowpass)) * Constants::half_pi<double>() ),int(2));
		randphaseFT[i] *= filter;
		randphaseFT[sz-i] *= filter;
	}
	randphaseFT[sz/2] = 0.;
	fftw_execute_r2r(*FTplan_hc2r_2xPtr,randphaseFT,randphase);

	for (size_t i=0;i<samples;++i){
		phivec[i] += randphase[i]/samples;
	}
	rhophi2cvec();
	return true;
}



PulseFreq & PulseFreq::attenuate(float attenfactor){
	rhovec *= (double)attenfactor;
	rhophi2cvec();
	return *this;
}
PulseFreq & PulseFreq::phase(float phasein){ // expects delay in units of pi , i.e. 1.0 = pi phase flip 
	if(intime){
		fft_tofreq();
	}
	phivec += phasein*Constants::pi<double>();
	rhophi2cvec();
	if(intime){
		fft_totime();
	}
	return *this;
}
PulseFreq & PulseFreq::setdelay(float delayin){ // expects delay in fs
	if(intime){
		fft_tofreq();
	}
	for (unsigned i=0;i<samples;i++){
		phivec[i] = omega[i]*(double)delayin/fsPau<double>();
	}
	rhophi2cvec();
	if(intime){
		fft_totime();
	}
	return *this;
}
PulseFreq & PulseFreq::delay(float delayin){ // expects delay in fs
	if(intime){
		fft_tofreq();
	}
	for (unsigned i=0;i<samples;i++){
		phivec[i] += omega[i]*(double)delayin/fsPau<double>();
	}
	rhophi2cvec();
	if(intime){
		fft_totime();
	}
	return *this;
}


void PulseFreq::printfrequency(std::ofstream * outfile){
	float nu,lambda;
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<float>())/fsPau<float>());
		lambda = C_nmPfs<float>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << std::pow(rhovec[i],int(2)) << "\t" << phivec[i] << "\n";
	}
}
void PulseFreq::printwavelengthbins(std::ofstream * outfile)
{
	std::vector<float> x(2);
	x.front() = C_nmPfs<float>()*2.0*pi<float>()*fsPau<float>()/omega[i_low];
	x.back() = C_nmPfs<float>()*2.0*pi<float>()*fsPau<float>()/omega[i_high-1];
	float dlam = (x.front()-x.back())/float(m_lamsamples);
        for (size_t i = 0;i<m_lamsamples;++i){
		(*outfile) << x.back() + i*dlam << "\t";
        }
	(*outfile) << "\n";
	return;
}
void PulseFreq::appendwavelength(std::ofstream * outfile)
{
	std::vector<float> x(i_high-i_low);
	std::vector<float> y(i_high-i_low);	
	for (size_t i=0;i<y.size();++i){
		x[i] = C_nmPfs<float>()*2.0*pi<float>()*fsPau<float>()/omega[i_low+i];
		//y[i] = std::pow(rhovec[i_low+i],int(2)) * 200000000000;
		y[i] = std::min((float)(std::pow(rhovec[i_low+i],int(2))),float(m_saturate));
	}
	float dlam = (x.front()-x.back())/float(m_lamsamples);
	boost::math::interpolators::barycentric_rational<float> interpolant(x.data(), y.data(), y.size());
	for (size_t i=0;i<m_lamsamples;++i){
		(*outfile) << uint16_t(interpolant(x.back()+i*dlam)) << "\t";
	}
	(*outfile) << std::endl;
	return;
}
void PulseFreq::appendwavelength_deriv(std::ofstream * outfile)
{
	std::vector<float> x(i_high-i_low);
	std::vector<float> y(i_high-i_low);	
	std::vector<int16_t> resultvec(m_lamsamples);
	std::vector<float> diffvec(m_lamsamples);

	for (size_t i=0;i<y.size();++i){
		x[i] = C_nmPfs<float>()*2.0*pi<float>()*fsPau<float>()/omega[i_low+i];
		y[i] = std::pow(rhovec[i_low+i],int(2));
	}
	float dlam = (x.front()-x.back())/float(m_lamsamples);
	boost::math::interpolators::barycentric_rational<float> interpolant(x.data(), y.data(), y.size());
	for (size_t i=0;i<m_lamsamples;++i){
		diffvec[i] = interpolant(x.back()+(i+1)*dlam) - interpolant(x.back()+i*dlam);
		//diffvec[i] -= interpolant(x.back()+(i+10)*dlam) - interpolant(x.back()+(i+11)*dlam);
	}
        float scale;
        int16_t max = std::numeric_limits<int16_t>::max();
        int16_t min = std::numeric_limits<int16_t>::min();
	auto bounds = std::minmax_element(diffvec.begin(),diffvec.end());
        if (*bounds.second>std::abs(*bounds.first)){
                scale = max/ *bounds.second;
        } else {
                scale = min/ *bounds.first;
        }
        std::transform(diffvec.begin(),diffvec.end(),resultvec.begin(),[scale](float x){return int16_t(scale*x);});
        *(outfile) << resultvec;

	return;
}
void PulseFreq::appendwavelength_bin(std::ofstream * outfile)
{
	std::vector<float> x(i_high-i_low);
	std::vector<float> y(i_high-i_low);	
	for (size_t i=0;i<y.size();++i){
		x[i] = C_nmPfs<float>()*2.0*pi<float>()*fsPau<float>()/omega[i_low+i];
		//y[i] = std::pow(rhovec[i_low+i],int(2)) * 200000000000;
		y[i] = std::min((float)std::pow(rhovec[i_low+i],int(2)),float(m_saturate));
	}
	float dlam = (x.front()-x.back())/float(m_lamsamples);
	boost::math::interpolators::barycentric_rational<float> interpolant(x.data(), y.data(), y.size());
	for (size_t i=0;i<m_lamsamples;++i){
		(*outfile) << int32_t(interpolant(x.back()+i*dlam));
	}
	return;
}
void PulseFreq::appendfrequency(std::ofstream * outfile){
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		uint16_t val = std::min(uint16_t(rhovec[i]),uint16_t(m_saturate));
       		(*outfile) << std::pow(val,int(2)) << "\t";
        }
	(*outfile) << std::endl;
}

void PulseFreq::appendnoisy(std::ofstream * outfile){
	std::normal_distribution<float> norm_dist( 0.0, m_noisescale);
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		float outval = std::pow(rhovec[i],int(2)) + norm_dist(rng);
       		(*outfile) << outval << "\t" ;
        }
	(*outfile) << std::endl;
}

void PulseFreq::printfrequencybins(std::ofstream * outfile){
	float nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<float>())/fsPau<float>());
		lam = C_nmPfs<float>()/nu;
       		(*outfile) << nu << "\t" << lam << "\n";
        }
	(*outfile) << "\n";
}
void PulseFreq::appendfrequencybins(std::ofstream * outfile){
	float nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<float>())/fsPau<float>());
       		(*outfile) << nu << "\t";
        }
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelay(std::ofstream * outfile, const float *delay){
	float nu,lambda,thisdelay;
	thisdelay = (*delay)*fsPau<float>();
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<float>())/fsPau<float>());
		lambda = C_nmPfs<float>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec[i] << "\t" << phivec[i] << "\t" << thisdelay << "\n";
	}
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelaychirp(std::ofstream * outfile, const float *delay,const float *chirp){
	float nu,lambda,thisdelay,thischirp,reltime;
	thisdelay = (*delay)*fsPau<float>();
	thischirp = (*chirp)*std::pow(fsPau<float>(),int(2));
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\tchirp[fs^2]\treldelays[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<float>())/fsPau<float>());
		lambda = C_nmPfs<float>()/nu;
		reltime = (nu-nu0)*thischirp*10.0;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec[i] << "\t" << phivec[i] << "\t" << thisdelay << "\t" << thischirp << "\t" << reltime << "\n";
	}
	(*outfile) << "\n";
}

void PulseFreq::printwavelength(std::ofstream * outfile,const float *delay){
	float nu,lambda,lamlast;
	lamlast=0;
        for (unsigned i = i_high;i>i_low;i-=(unsigned)(atoi(getenv("sampleinterval")))){
                nu = (omega[i]/(2.0*pi<float>())/fsPau<float>());
                lambda = (float)((int)((C_nmPfs<float>()/nu)*10.0))/10.0;
		if (lambda>lamlast & (int)(lambda*10)%5 == 0){
                	(*outfile) << lambda << "\t" << (*delay) << "\t" << std::pow(rhovec[i],int(2)) << "\n";
			lamlast = lambda;
		}
        }
	(*outfile) << "\n";
}



PulseFreq & PulseFreq::modulateamp_time(const std::vector<double> & modulation)
{
	if (infreq){
		std::cerr << "whoops, trying time modulation but in frequency domain\n" << std::flush;
		return *this;
	}
	assert(modulation.size() ==samples);
	std::transform(modulation.begin(),modulation.end(),rhovec.begin(),rhovec.begin(),[](const double &x,double &y){return (x>0?y*std::tanh(x):double(0));});
	return *this;
}
PulseFreq & PulseFreq::modulatephase_time(const std::vector<double> & modulation)
{
	if (infreq){
		std::cerr << "whoops, trying time modulation but in frequency domain\n" << std::flush;
		return *this;
	}
	assert(modulation.size() ==samples);
	std::transform(modulation.begin(),modulation.end(),phivec.begin(),phivec.begin(),[](const double &x,double &y){return (y + x);});
	return *this;
}

PulseFreq & PulseFreq::modulateamp_freq(const std::vector<double> & modulation) 
{
	if (intime){
		std::cerr << "whoops, trying time modulation but in time domain\n" << std::flush;
		return *this;
	}
	assert(modulation.size() == samples/2);
	std::transform(modulation.begin(),modulation.end(),rhovec.begin(),rhovec.begin(),[](const double &x,double &y){return (x>0?y*std::tanh(x):double(0));});
	std::transform(modulation.begin(),modulation.end(),rhovec.rbegin(),rhovec.rbegin(),[](const double &x,double &y){return (x>0?y*std::tanh(x):double(0));});
	rhophi2cvec();
	return *this;
}

PulseFreq & PulseFreq::modulatephase_freq(const std::vector<double> & modulation) 
{
	if (intime){
		std::cerr << "whoops, trying time modulation but in time domain\n" << std::flush;
		return *this;
	}
	assert(modulation.size() == samples/2);
	std::transform(modulation.begin(),modulation.end(),phivec.begin(),phivec.begin(),[](const double &x,double &y){return (y + x);});
	std::transform(modulation.begin(),modulation.end(),phivec.rbegin(),phivec.rbegin(),[](const double &x,double &y){return (y - x);});
	rhophi2cvec();
	return *this;
}

PulseFreq & PulseFreq::fillphase(std::vector < float > & data)
{
	if (intime)
		fft_tofreq();
	assert(samples==data.size());
	for (size_t i=0;i<samples;i++){
		data[i] = float(std::arg(cvec[i]));
	}
	return *this;
}
PulseFreq & PulseFreq::fillspect(std::vector < float > & data)
{
	if (intime)
		fft_tofreq();
	assert(samples==data.size());
	for (size_t i=0;i<samples;i++){
		data[i] = float(std::norm(cvec[i]));
	}
	return *this;
}
PulseFreq & PulseFreq::fillfreq(std::vector < float > & data)
{
	if (intime)
		fft_tofreq();
	assert(samples==data.size());
	for (size_t i=0;i<samples;i++){
		data[i] = float(std::real(cvec[i]));
	}
	return *this;
}
PulseFreq & PulseFreq::filltime_envelope(std::vector < float > & data)
{
	if (infreq)
		fft_totime();
	assert(samples==data.size());
	for (size_t i=0;i<samples/2;i++){
		data[i] = float(std::norm(cvec[i+samples/2]));
		data[i+samples/2] = float(std::norm(cvec[i]));
	}
	return *this;
}
PulseFreq & PulseFreq::filltime(std::vector < float > & data)
{
	if (infreq)
		fft_totime();	
	assert(data.size()==samples);
	for (size_t i=0;i<samples/2;i++){
		data[i] = float(cvec[i+samples/2].real());
		data[i+samples/2] = float(cvec[i].real());
	}
	return *this;
}

void PulseFreq::printtime(std::ofstream * outfile){
	(*outfile) << ("#time\treal\timag\n");
	for (unsigned i = samples/2;i<samples;i++){
		(*outfile) << (time[i]*fsPau<float>()) << "\t" << cvec[i].real() << "\t" << cvec[i].imag() << "\n";
	}
	for (unsigned i = 0;i<samples/2;i++){
		(*outfile) << (time[i]*fsPau<float>()) << "\t" << cvec[i].real() << "\t" << cvec[i].imag() << "\n";
	}

} 


PulseFreq & PulseFreq::rebuildvectors(float newgain){
	m_gain = newgain;

	startind = (unsigned)((omega_center-(omega_width/2.0))/domega);
	stopind = (unsigned)((omega_center+(omega_width/2.0))/domega);
	onwidth = (unsigned)(omega_onwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	offwidth = (unsigned)(omega_offwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2

	for (unsigned i = 1; i<startind;i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
	}
	for (unsigned i = startind;i<startind+onwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = m_gain*rising(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[samples-i] = -domega*(float)i;
		time[samples-i] = -dtime*(float)i;
		rhovec[samples-i] = m_gain*rising(i);
		cvec[samples-i] = std::polar(rhovec[samples-i],phivec[samples-i]);
	}
	for (unsigned i = startind+onwidth;i<stopind-offwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = m_gain; 
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[samples-i] = -domega*(float)i;
		time[samples-i] = -dtime*(float)i;
		rhovec[samples-i] = m_gain;
		cvec[samples-i] = std::polar(rhovec[samples-i],phivec[samples-i]);
	}
	for (unsigned i = stopind-offwidth;i<stopind; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = m_gain*falling(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
		rhovec[samples-i] = m_gain*falling(i);
		cvec[samples - i] = std::polar(rhovec[samples - i],phivec[samples - i]);
	}
	for (unsigned i = stopind;i<samples/2; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
	}
	return *this;
}

PulseFreq & PulseFreq::buildvectors(const size_t s){
	//std::cerr << "allocating with fftw_malloc with samples = " << samples << std::endl << std::flush;
	cvec = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * size_t(s));
        std::fill(cvec,cvec + samples,std::complex<float>(0));
	r_vec = (double *) fftw_malloc(sizeof(double) * s);
        std::fill(r_vec,r_vec + s,double(0));
	//std::cerr << "\t\t...allocated with fftw_malloc with samples = " << samples << std::endl << std::flush;
	hc_vecFT = (double *) fftw_malloc(sizeof(double) * size_t(s));
        std::fill(hc_vecFT,hc_vecFT + s,double(0));
	//std::cerr << "allocating with fftw_malloc with samples = " << (2*samples) << std::endl << std::flush;
	r_vec_2x = (double *) fftw_malloc(sizeof(double) * s * 2);
        //std::fill(r_vec_2x,r_vec_2x + 2*samples,float(0));
	hc_vec_2xFT = (double *) fftw_malloc(sizeof(double) * s * 2);
        std::fill(hc_vec_2xFT,hc_vec_2xFT + 2*s,double(0));
	//std::cerr << "\t\t...allocated with fftw_malloc with 2*samples = " << (2*samples) << std::endl << std::flush;

	rhovec.resize(s,0.0);
	phivec.resize(s,0.0);
	omega.resize(s);
	time.resize(s);
	samples = s;

	omega[0] = 0.0;
	time[0] = 0.0;
	omega[s/2] = -(float)(s/2)*domega;
	time[s/2] = -(float)(s/2)*dtime;

	startind = (unsigned)((omega_center-(omega_width/2.0))/domega);
	stopind = (unsigned)((omega_center+(omega_width/2.0))/domega);
	onwidth = (unsigned)(omega_onwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	offwidth = (unsigned)(omega_offwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2

	for (unsigned i = 1; i<startind;i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
	}
	for (unsigned i = startind;i<startind+onwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = m_gain*rising(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[s-i] = -domega*(float)i;
		time[s-i] = -dtime*(float)i;
		rhovec[s-i] = m_gain*rising(i);
		cvec[s-i] = std::polar(rhovec[s-i],phivec[s-i]);
	}
	for (unsigned i = startind+onwidth;i<stopind-offwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = m_gain; 
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[s-i] = -domega*(float)i;
		time[s-i] = -dtime*(float)i;
		rhovec[s-i] = m_gain;
		cvec[s-i] = std::polar(rhovec[s-i],phivec[s-i]);
	}
	for (unsigned i = stopind-offwidth;i<stopind; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = m_gain*falling(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[s-i] = -domega*i;
		time[s-i] = -dtime*i;
		rhovec[s-i] = m_gain*falling(i);
		cvec[s- i] = std::polar(rhovec[s- i],phivec[s- i]);
	}
	for (unsigned i = stopind;i<s/2; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[s-i] = -domega*i;
		time[s-i] = -dtime*i;
	}
	return *this;	
}
PulseFreq & PulseFreq::killvectors(void){
	if (cvec != NULL)
		fftw_free(cvec);
	if (r_vec != NULL)
		fftw_free(r_vec);
	if (hc_vecFT != NULL)
		fftw_free(hc_vecFT);
	if (r_vec_2x != NULL)
		fftw_free(r_vec_2x);
	if (hc_vec_2xFT != NULL)
		fftw_free(hc_vec_2xFT);
	cvec = NULL;
	r_vec = hc_vecFT = r_vec_2x = hc_vec_2xFT = NULL;
	return *this;
}

PulseFreq & PulseFreq::setplans(const PulseFreq & rhs)
{
	//std::cerr << "in PulseFreq::setplans() FTplan_forwardPtr.use_count() = " << FTplan_forwardPtr.use_count() << std::endl << std::flush;
	//assert(FTplan_forwardPtr.use_count()>0 && FTplan_backwardPtr.use_count()>0);
	FTplan_forwardPtr = rhs.FTplan_forwardPtr;
	FTplan_backwardPtr = rhs.FTplan_backwardPtr;
	return *this;
}
PulseFreq & PulseFreq::setmasterplans(fftw_plan * forward,fftw_plan * backward)
{
	//assert(FTplan_forwardPtr.use_count()==0 && FTplan_backwardPtr.use_count()==0);
	*forward = fftw_plan_dft_1d(samples, 
			reinterpret_cast<fftw_complex*>(cvec),
			reinterpret_cast<fftw_complex*>(cvec), 
			FFTW_FORWARD, FFTW_ESTIMATE);
	*backward = fftw_plan_dft_1d(samples, 
			reinterpret_cast<fftw_complex*>(cvec), 
			reinterpret_cast<fftw_complex*>(cvec), 
			FFTW_BACKWARD, FFTW_ESTIMATE);
	FTplan_forwardPtr = forward;
	FTplan_backwardPtr = backward;
	return *this;
}
PulseFreq & PulseFreq::setancillaryplans(const PulseFreq & rhs)
{
	/*
	 assert(FTplan_r2hcPtr.use_count()>0
			&& FTplan_hc2rPtr.use_count()>0
			&& FTplan_r2hc_2xPtr.use_count()>0
			&& FTplan_hc2r_2xPtr.use_count()>0);
	*/
	FTplan_r2hcPtr = rhs.FTplan_r2hcPtr;
	FTplan_hc2rPtr = rhs.FTplan_hc2rPtr;
	FTplan_r2hc_2xPtr = rhs.FTplan_r2hc_2xPtr;
	FTplan_hc2r_2xPtr = rhs.FTplan_hc2r_2xPtr;
	return *this;
}
PulseFreq & PulseFreq::setmasterancillaryplans(fftw_plan * r2hc,fftw_plan * hc2r,fftw_plan * r2hc_2x,fftw_plan * hc2r_2x)
{

	/*assert(FTplan_r2hcPtr.use_count()==0
			&& FTplan_hc2rPtr.use_count()==0
			&& FTplan_r2hc_2xPtr.use_count()==0
			&& FTplan_hc2r_2xPtr.use_count()==0);
			*/
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
	FTplan_r2hcPtr = r2hc;
	FTplan_hc2rPtr = hc2r;
	FTplan_r2hc_2xPtr = r2hc_2x;
	FTplan_hc2r_2xPtr = hc2r_2x;
	return *this;
}


PulseFreq & PulseFreq::mulAllAmp(void)
{
	std::vector<double> modvec(samples/2,amp_0th);
	for (size_t i=0;i<modvec.size();i++){
		modvec[i] += std::pow(Constants::fsPau<double>(),int(-1))*amp_1st * ((double)omega[i]-omega_center);
		modvec[i] += std::pow(Constants::fsPau<double>(),int(-2))*amp_2nd*std::pow((double)omega[i]-omega_center,int(2));
		modvec[i] += std::pow(Constants::fsPau<double>(),int(-3))*amp_3rd*std::pow((double)omega[i]-omega_center,int(3));
		modvec[i] += std::pow(Constants::fsPau<double>(),int(-4))*amp_4th*std::pow((double)omega[i]-omega_center,int(4));
		modvec[i] += std::pow(Constants::fsPau<double>(),int(-5))*amp_5th*std::pow((double)omega[i]-omega_center,int(5));
	}
	modulateamp_freq(modvec);
	return *this;
}
PulseFreq & PulseFreq::addAllPhase(void)
{
	std::vector<double> modvec(samples/2,0.);
	//std::cerr << "Constants::fsPau<double>() = " << Constants::fsPau<double>() << std::endl;
	//std::cerr << "domega = " << domega << std::endl;
	for (size_t i=0;i<modvec.size();i++){
		modvec[i] = std::pow(Constants::fsPau<double>(),int(-2))*phase_GDD*std::pow((double)omega[i]-omega_center,int(2));
		modvec[i] += std::pow(Constants::fsPau<double>(),int(-3))*phase_TOD*std::pow((double)omega[i]-omega_center,int(3));
		modvec[i] += std::pow(Constants::fsPau<double>(),int(-4))*phase_4th*std::pow((double)omega[i]-omega_center,int(4));
		modvec[i] += std::pow(Constants::fsPau<double>(),int(-5))*phase_5th*std::pow((double)omega[i]-omega_center,int(5));
	}
	modulatephase_freq(modvec);
	return *this;
}


PulseFreq & PulseFreq::setGDDtoindex(const unsigned indx,const int omega_sign) {
	phivec[indx] = omega_sign*phase_GDD*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(2));
	rhophi2cvec(indx);
	return *this;
}	
PulseFreq & PulseFreq::setTODtoindex(const unsigned indx,const int omega_sign) {
	phivec[indx] = omega_sign*phase_TOD*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(3));
	rhophi2cvec(indx);
	return *this;
}	
PulseFreq &  PulseFreq::set4thtoindex(const unsigned indx,const int omega_sign) {
	phivec[indx] = omega_sign*phase_4th*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(4));
	rhophi2cvec(indx);
	return *this;
}	
PulseFreq &  PulseFreq::set5thtoindex(const unsigned indx,const int omega_sign) {
	phivec[indx] = omega_sign*phase_5th*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(5));
	rhophi2cvec(indx);
	return *this;
}	
PulseFreq & PulseFreq::addGDDtoindex(const unsigned indx,const int omega_sign) {
	phivec[indx] += omega_sign*phase_GDD*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(2));
	rhophi2cvec(indx);
	return *this;
}	
PulseFreq & PulseFreq::addTODtoindex(const unsigned indx,const int omega_sign) {
	phivec[indx] += omega_sign*phase_TOD*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(3));
	rhophi2cvec(indx);
	return *this;
}	
PulseFreq &  PulseFreq::add4thtoindex(const unsigned indx,const int omega_sign) {
	phivec[indx] += omega_sign*phase_4th*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(4));
	rhophi2cvec(indx);
	return *this;
}	
PulseFreq &  PulseFreq::add5thtoindex(const unsigned indx,const int omega_sign) {
	phivec[indx] += omega_sign*phase_5th*std::pow(omega[indx]-(float(omega_sign)*omega_center),int(5));
	rhophi2cvec(indx);
	return *this;
}	
