#include "generate_dataset.hpp"
#include <cstdlib>
#include <algorithm>
#include <complex>
#include <memory>
#include <fftw3.h>

#include <vector>
#include <random>
#include <chrono>

#include <DebugOps.hpp>

#include <cstddef> // for nullptr
#include <cstdint> // for writing as an int32_t and int16_t

using namespace Constants;

/* Here's main */
int main(int argc, char* argv[])
{
	std::time_t tstart = std::time(nullptr);
	std::cout << "\t\t======================================\n"
		<< "\t\t======= scan_material started ========\n"
		<< "\t\t===== " << std::asctime(std::localtime(&tstart)) 
		<< "\t\t===== on host " << getenv("HOSTNAME") << "\n"
		<< "\t\t===== using " << getenv("nthreads") << " threads\n"
		<< "\t\t======================================\n" << std::flush;

	unsigned nthreads = (unsigned)atoi( getenv("nthreads") );

	Params params;
	params.npulses(size_t(atoi(getenv("npulses"))));
	params.filebase(std::string(getenv("filebase")));

	std::vector<float> runtimes(params.npulses(),0); // for benchmarking the processors

	params.dalpha((atof(getenv("drifting_alpha")))*pi<double>()/params.npulses());

	if (params.doublepulse(atoi(getenv("doublepulse")) > 0 )){
		params.doublepulsedelay(atof( getenv("doublepulsedelay") ) ) ; // this one gets used directly in atomic units
	}
	params.lambda_0(atof( getenv("lambda0") ));
	params.lambda_width( atof( getenv("lambda_width") ));
	params.lambda_onoff( atof( getenv("lambda_onoff") ));
	params.tspan((atof( getenv("tspan") ) )/fsPau<double>());

	params.ngroupsteps(atoi( getenv("ngroupsteps") ));
	params.groupdelay(atof(getenv("groupdelay")));
	params.backdelay(atof(getenv("backdelay")));
	params.netalon(atoi(getenv("netalon")));


	params.etalonreflectance(atof(getenv("etalon")));
	params.etalondelay(atof(getenv("etalondelay")));
	params.interferedelay((double)atof(getenv("interferedelay")));

	params.chirp(
			( atof( getenv("chirp") ) ) / std::pow(fsPau<float>(),int(2)), // the difference in slopes at omega_low versus omega_high must equal tspan
			( atof( getenv("TOD") ) ) / std::pow(fsPau<float>(),int(3)),
			( atof( getenv("FOD") ) ) / std::pow(fsPau<float>(),int(4)),
			( atof( getenv("fifthOD") ) ) / std::pow(fsPau<float>(),int(5))
			);


	if (params.addchirpnoise(atoi(getenv("usechirpnoise"))>0)){
		params.initchirpnoise( 
				( atof( getenv("chirpnoise") ) ) / std::pow(fsPau<float>(),int(2)), 
				( atof( getenv("TODnoise") ) ) / std::pow(fsPau<float>(),int(3)),
				( atof( getenv("FODnoise") ) ) / std::pow(fsPau<float>(),int(4)),
				( atof( getenv("fifthODnoise") ) ) / std::pow(fsPau<float>(),int(5))
				);
	}




	FiberBundle masterbundle(boost::lexical_cast<size_t>(atoi(getenv("nfibers"))));
	masterbundle.fiberdiameter(boost::lexical_cast<float>(atof(getenv("fiberdiam"))));
	masterbundle.laserdiameter(boost::lexical_cast<float>(atof(getenv("laserdiam"))));
	masterbundle.xraydiameter(boost::lexical_cast<float>(atof(getenv("xraydiam"))));
	masterbundle.set_fsPmm(boost::lexical_cast<float>(atof(getenv("bundle_fsPmm"))));
	masterbundle.scalePolarCoords();

	std::cout << "\t\tshuffle fibers?\t";
	if (getenv("shuffle_fibers"))
	{
		std::cout << "yes\n";masterbundle.shuffle_output();
	} else {
		std::cout << "no\n";
	}

	masterbundle.Ixray(float(1.));
	masterbundle.Ilaser(float(1.));
	std::string filename = params.filebase() + "fibermap.out";
	std::cout << "fibermap file = " << filename << std::endl << std::flush;
	std::ofstream mapfile(filename.c_str(),std::ios::out);
	//masterbundle.print_mapping(mapfile,double(0.0));
	mapfile.close();

	// file for delay bins
	ofstream outbins(std::string(params.filebase() + "delaybins.out").c_str(),ios::out); 
	for (size_t f=0;f<masterbundle.get_nfibers();++f){
		outbins << masterbundle.delay(f) << "\n";
	}
	outbins.close();

	MatResponse masterresponse(
			0,															// stepdelay
			(double)( atof( getenv("stepwidth") ) ),								// stepwidth
			((double)( atof( getenv("attenuation") ) ) - 1.0) * masterbundle.Ixray() / params.ngroupsteps() + 1.0,	// attenuation
			(double)( atof( getenv("phase") ) ) * masterbundle.Ixray() / params.ngroupsteps()				// phase
			);
	masterresponse.aalphabbeta(
			(double)( atof( getenv("a") ) ),		// a
			(double)( atof( getenv("alpha" ) ) ),		// alpha
			(double)( atof( getenv("b") ) ),		// b
			(double)( atof( getenv("beta") ) )		// beta
			);

	masterresponse.setreflectance(params.etalonreflectance());
	masterresponse.setetalondelay(params.etalondelay());

	std::cout << "initializing masterpulse and masterplans" << std::endl << std::flush;
	PulseFreq masterpulse(params.omega0(),params.omega_width(),params.omega_onoff(),params.tspan());

	fftw_plan forward;
	fftw_plan backward;
	fftw_plan plan_r2hc;
	fftw_plan plan_hc2r;
	fftw_plan plan_r2hc_2x;
	fftw_plan plan_hc2r_2x;
	masterpulse.setmasterplans(&forward,&backward);
	masterpulse.setancillaryplans(& plan_r2hc,& plan_hc2r,& plan_r2hc_2x,& plan_hc2r_2x);
	masterpulse.addchirp(params.getchirp());							// chirp that ref pulse

	double xrayphoton_energy = double(atof(getenv("xrayphoton_energy")));
	masterresponse.bandgap(double(atof(getenv("bandgap_eV")))); //
	if (getenv("usediamond")){
		masterresponse.fill_carriersvec(masterpulse,xrayphoton_energy);
	} else {
		std::string carriersfilename = getenv("carriersfile");
		std::cerr << "carriersfilename = " << carriersfilename << "\n" << std::flush;
		std::ifstream Nikita_file(carriersfilename.c_str(),std::ios::in);
		masterresponse.fill_carriersvec(masterpulse,Nikita_file);
	}


	std::time_t tstop = std::time(nullptr);
	std::cout << "\tIt has taken " << (tstop-tstart) << " s so far for initializing masterpulse and building fftw plans\n" << std::flush;

	CalibMat calibration(boost::lexical_cast<size_t>(atoi(getenv("ncalibdelays")))
			, boost::lexical_cast<double>(atof(getenv("fsWindow"))));
	if (!getenv("skipcalibration"))
	{
		std::cout << "\t\t############ entering calibration ###########\n" << std::flush;
		calibration.set_center(boost::lexical_cast<double>(atof(getenv("delays_mean"))));
		std::cout << "\t\t====== delays =======\n";
		for (size_t i = 0 ; i< calibration.get_ndelays(); ++i){
			std::cout << calibration.get_delay(i) << " ";
		}
		std::cout << std::endl << std::flush;
	}

	if (!getenv("skipcalibration"))
	{
		// Setup the shared pulse arrays
		std::vector< PulseFreq > calpulsearray(calibration.get_ndelays(),masterpulse);

#pragma omp parallel num_threads(nthreads) default(shared) shared(masterpulse)
		{ // begin parallel region 1
			size_t tid = omp_get_thread_num();

			// all non-shared objects must be created inside the parallel section for default is shared if defined outside
			// http://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-data.html
			PulseFreq threadpulse(masterpulse);

			calpulse.delay(params.randomdelay()); // expects this in fs 

#pragma omp barrier

#pragma omp master
			{
				std::cout << "|\t done with calibration delays\n" << std::flush;
				//std::cerr << "\t\t###### Made it here too #####\n\t\t##### should call only once in master #######\n" << std::flush;
				// print out the calibration as ascii for now //
				// print rows in order, eventually in tf_record or matrix or so. //
				std::string calfilename = params.calfilebase() + "interference.calibration";
				std::string derivfilename = params.calfilebase() + "interference.calibration.derivative";
				std::string calfilename_delays = params.calfilebase() + "interference.calibration.delays";
				std::string calfilename_wavelengths = params.calfilebase() + "interference.calibration.wavelengths";
				ofstream calibrationstream(calfilename.c_str(),ios::out); 
				ofstream derivstream(derivfilename.c_str(),ios::out); 
				ofstream calibrationstream_delays(calfilename_delays.c_str(),ios::out); 
				ofstream calibrationstream_wavelengths(calfilename_wavelengths.c_str(),ios::out); 
				/*
				std::string bin_calfilename = params.filebase() + "interference.calibration.bin";
				ofstream bin_calibrationstream(bin_calfilename.c_str(),ios::out | ios::binary); 
				std::cout << "\tcalibration filename out = " << calfilename << "\n\t and \t" << bin_calfilename << std::endl;
				*/
				calibrationstream << "# wavelengths\n#";
				derivstream << "# wavelengths\n#";
				calpulsearray[0].printwavelengthbins(&calibrationstream);
				calpulsearray[0].printwavelengthbins(&derivstream);
				calpulsearray[0].printwavelengthbins(&calibrationstream_wavelengths);
				calibrationstream << "# delays\n#";
				derivstream << "# delays\n#";
				calibrationstream_delays << "# delays\n";
				for (size_t i = 0 ; i< calibration.get_ndelays(); ++i){
					calibrationstream << calibration.get_delay(i) << "\t";
					derivstream << calibration.get_delay(i) << "\t";
					calibrationstream_delays << calibration.get_delay(i) << "\t";
				}
				calibrationstream << "\n";
				derivstream << "\n";
				calibrationstream_delays << "\n";

				for (size_t n=0;n<calpulsearray.size();++n){
					calpulsearray[n].appendwavelength(&calibrationstream);
					calpulsearray[n].appendwavelength_deriv(&derivstream);
					//	calpulsearray[n].appendwavelength_bin(&bin_calibrationstream);
				}

				calibrationstream.close();
				derivstream.close();
				calibrationstream_delays.close();
				calibrationstream_wavelengths.close();
				//bin_calibrationstream.close();
				std::cout << "Finished with the calibration image/matrix\n" << std::flush;


			}

#pragma omp master
			{
				std::cout << "\t\t############ ending parallel region 1 ###########\n" << std::flush;
			}
		} // end parallel region 1

	} // end if (!getenv("skipcalibration"))


	//############## Images section ##############

#pragma omp parallel num_threads(nthreads) default(shared) shared(masterpulse)
		{
	if (!getenv("skipimages"))
	{
		std::cout << "\t\t############ entering parallel/images ###########\n" << std::flush;
			size_t tid = omp_get_thread_num();

			PulseFreq crosspulse(masterpulse);
			PulseFreq etalonpulse(masterpulse);
			PulseFreq crossetalonpulse(masterpulse);
			std::vector< PulseFreq > pulsearray(nfibers,PulseFreq(masterpulse));
#pragma omp barrier

			if (params.addrandomphase(atoi(getenv("addrandomphase"))>0))
			{
				masterpulse.addrandomphase();
				std::string filename = params.filebase() + "spectralphaseFTpower.dat";
				std::ofstream outfile(filename.c_str(),std::ios::out);
				masterpulse.print_phase_powerspectrum(outfile);
				outfile.close();
				filename = params.filebase() + "spectralphase.dat";
				outfile.open(filename.c_str(),std::ios::out);
				masterpulse.print_phase(outfile);
				outfile.close();
				filename = params.filebase() + "spectralamp.dat";
				outfile.open(filename.c_str(),std::ios::out);
				masterpulse.print_amp(outfile);
				outfile.close();
			}


#pragma omp for schedule(dynamic)
			for (size_t n=0;n<params.npulses();++n)	{ // outermost loop for npulses to produce //
				PulseFreq pulse(masterpulse);
				//std::cerr << "\tinside the parallel region 2 for images loop n = " << n << " in thread " << tid << "\n" << std::flush;
				if (n==0 & tid==0) {
					std::cout << "========================================================================="
						<<   "\n\t\t ==== http://www.fftw.org/fftw3_doc/Advanced-Complex-DFTs.html ===="
						<<   "\n\t\t ====         use this for defining multiple fibers as         ===="
						<<   "\n\t\t ====         contiguous blocks for row-wise FFT as 2D         ===="
						<<   "\n\t\t ==================================================================\n" << std::flush;
				}

				std::time_t runstart = std::time(nullptr);

				double t0 = params.delays_uniform();
				double startdelay(0);

				parabundle = masterbundle;


				parabundle.Ixray(params.xray_inten_rand());
				parabundle.Ilaser(params.laser_inten_rand());
				parabundle.delay_angle(params.dalpha()*double(n));
				parabundle.center_Ixray(params.xray_pos_rand(),params.xray_pos_rand());
				parabundle.center_Ilaser(params.laser_pos_rand(),params.laser_pos_rand());



				//DebugOps::pushout(std::string("Running image " + std::to_string(n) + " for t0 = " + std::to_string(t0) + " in threaded for loop, thread " + std::to_string(tid)));
				std::string mapfilename = params.filebase() + "fibermap.out." + std::to_string(n);
				//std::cout << "fibermap file = " << mapfilename << std::endl << std::flush;
				std::ofstream mapfile(mapfilename.c_str(),std::ios::out);
				parabundle.print_mapping(mapfile,t0);
				mapfile.close();


				for(size_t f = 0; f < parabundle.get_nfibers(); f++)
				{ // begin fibers loop
					pulse = masterpulse;
					crosspulse = masterpulse;
					startdelay = t0 + parabundle.delay(f);
					pulse.scale(parabundle.Ilaser(f));
					crosspulse.scale(parabundle.Ilaser(f));

					pararesponse = masterresponse;

					if (getenv("scale_fibers")){
						pararesponse.setscale(parabundle.Ixray(f));
						//std::cerr << "parabundle.Ixray(" << f << ") = " << parabundle.Ixray(f) << "\n" << std::flush;
					}

					if (params.addchirpnoise()){
						std::vector<double> noise(params.getchirpnoise());
						pulse.addchirp(noise); 
						crosspulse.addchirp(noise); 
					}

					crosspulse.delay(params.interferedelay()); // delay in the frequency domain
					pulse.fft_totime();
					crosspulse.fft_totime();

					for(size_t g=0;g<params.ngroupsteps();g++){ // begin groupsteps loop
						pararesponse.setdelay(startdelay - g*params.groupstep()); // forward propagating, x-rays advance on the optical
						pararesponse.setstepvec_amp(pulse);
						pararesponse.setstepvec_phase(pulse);
						pararesponse.setstepvec_amp(crosspulse);
						pararesponse.setstepvec_phase(crosspulse);
						if (params.doublepulse()){
							pararesponse.addstepvec_amp(pulse,params.doublepulsedelay());
							pararesponse.addstepvec_phase(pulse,params.doublepulsedelay());
							pararesponse.addstepvec_amp(crosspulse,params.doublepulsedelay());
							pararesponse.addstepvec_phase(crosspulse,params.doublepulsedelay());
						}
						// this pulls down the tail of the response so vector is periodic on nsamples	
						pararesponse.buffervectors(pulse); 
						pararesponse.buffervectors(crosspulse); 
						pulse.modulateamp_time();
						pulse.modulatephase_time();
						crosspulse.modulateamp_time();
						crosspulse.modulatephase_time();
					}// end groupsteps loop
					//std::cerr << "tid = " << tid << "\tpulse/crosspulse.domain() = " << pulse.domain() << "/" << crosspulse.domain() << "\n" << std::flush;


					for (size_t e=0;e<params.netalon();e++){ // begin etalon loop
						//std::cerr << "tid = " << tid << "\tpulse/crosspulse.domain() = " << pulse.domain() << "/" << crosspulse.domain() << "\n" << std::flush;
						//std::cerr << "\n\t\t ---- starting etalon at " << e << " ----\n" << std::flush;
						// back propagation step //
						double etalondelay = startdelay - double(e+1) * (pararesponse.getetalondelay()); 
						// at front surface, x-rays see counter-propagating light from one full etalon delay

						etalonpulse = pulse;
						crossetalonpulse = crosspulse;
						//std::cerr << "etalonpulse/crossetalonpulse.domain() = " << etalonpulse.domain() << "/" << crossetalonpulse.domain() << "\n" << std::flush;

						for(size_t g=0;g<params.ngroupsteps();g++){
							pararesponse.setdelay(etalondelay + g*params.backstep()); 
							// counterpropagating, x-rays work backwards through the optical

							pararesponse.setstepvec_amp(etalonpulse);
							pararesponse.setstepvec_phase(etalonpulse);
							pararesponse.setstepvec_amp(crossetalonpulse);
							pararesponse.setstepvec_phase(crossetalonpulse);
							if (params.doublepulse()){
								pararesponse.addstepvec_amp(etalonpulse,params.doublepulsedelay());
								pararesponse.addstepvec_phase(etalonpulse,params.doublepulsedelay());
								pararesponse.addstepvec_amp(crossetalonpulse,params.doublepulsedelay());
								pararesponse.addstepvec_phase(crossetalonpulse,params.doublepulsedelay());
							}
							pararesponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							pararesponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							etalonpulse.modulateamp_time();
							etalonpulse.modulatephase_time();
							crossetalonpulse.modulateamp_time();
							crossetalonpulse.modulatephase_time();
						}
						// forward propagation //
						//std::cerr << "\t\t\t ########### // forward propagation // #############\n" << std::flush;
						for(size_t g=0;g<params.ngroupsteps();g++){
							pararesponse.setdelay(startdelay - g*params.groupstep()); // forward propagating, x-rays advance on the optical
							pararesponse.setstepvec_amp(etalonpulse);
							pararesponse.setstepvec_phase(etalonpulse);
							pararesponse.setstepvec_amp(crossetalonpulse);
							pararesponse.setstepvec_phase(crossetalonpulse);
							if (params.doublepulse()){
								pararesponse.addstepvec_amp(etalonpulse,params.doublepulsedelay());
								pararesponse.addstepvec_phase(etalonpulse,params.doublepulsedelay());
								pararesponse.addstepvec_amp(crossetalonpulse,params.doublepulsedelay());
								pararesponse.addstepvec_phase(crossetalonpulse,params.doublepulsedelay());

							}
							pararesponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							pararesponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							etalonpulse.modulateamp_time();
							etalonpulse.modulatephase_time();
							crossetalonpulse.modulateamp_time();
							crossetalonpulse.modulatephase_time();
						}
						//std::cerr << "etalonpulse/crossetalonpulse.domain() = " << etalonpulse.domain() << "/" << crossetalonpulse.domain() << "\n" << std::flush;
						etalonpulse.fft_tofreq();
						crossetalonpulse.fft_tofreq();
						etalonpulse.delay(pararesponse.getetalondelay()); // delay and attenuate in frequency domain
						etalonpulse.attenuate(pow(pararesponse.getreflectance(),(int)2));
						crossetalonpulse.delay(pararesponse.getetalondelay()); // delay and attenuate in frequency domain
						crossetalonpulse.attenuate(pow(pararesponse.getreflectance(),(int)2));
						etalonpulse.fft_totime();
						crossetalonpulse.fft_totime();
						pulse += etalonpulse;
						crosspulse += crossetalonpulse;
					} // end etalon loop


					pulse.fft_tofreq();
					crosspulse.fft_tofreq();
					pulse.delay(params.interferedelay()); // expects this in fs // time this back up to the crosspulse
					pulse -= crosspulse;
					// std::cerr << "\n\n\t\t\t\t============== testing... just before the push_back() ==============\n\n" << std::flush;
					pulsearray[f] = pulse;
				} // end nfibers loop


				std::string filename = params.filebase() + "interference.out." + std::to_string(n);// + ".tid" + std::to_string(tid);
				//std::cerr << "testing: filename = " << filename << "\n" << std::flush;
				ofstream interferestream(filename.c_str(),ios::out); // use app to append delays to same file.

				//std::cout << "tid = " << tid << ": interfere filename out = " << filename << std::endl;
				std::complex<double> z_laser = parabundle.center_Ilaser();
				std::complex<double> z_xray = parabundle.center_Ixray();
				interferestream << "#delay for image = \t" << t0 
					<< "\n#Ilaser = \t" << parabundle.Ilaser()
					<< "\n#Ixray = \t" << parabundle.Ixray()
					<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
					<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
					<< "\n#alpha = \t" << parabundle.delay_angle() 
					<< std::endl;
				interferestream << "#";
				pulsearray[0].printwavelengthbins(&interferestream);
				for (size_t f=0;f<pulsearray.size();f++){
					pulsearray[f].scale(parabundle.Ilaser(f)); 
					pulsearray[f].appendwavelength(&interferestream);
				}
				if (tid % 10 < 2){
					for (size_t f=0;f<parabundle.get_nfibers();f++){
						int max = boost::lexical_cast<double>(getenv("gain")) * pulsearray[f].maxsignal();
						for (size_t i=0;i<std::log(max);++i){
							std::cout << '.';
						}
						std::cout << "|";
					}
					std::cout << "\timg = " << n << " in tid = " << tid << "\n" << std::flush;
				}
				interferestream.close();

				std::time_t runstop = std::time(nullptr);
				runtimes[n] = float(runstop - runstart);

			} // outermost loop for npulses to produce //

#pragma omp master
			{
				std::cout << "\t\t############ ending parallel region 2 ###########\n" << std::flush;
			}

#pragma omp barrier

			//std::cerr << "\n\t... trying to leave parallel region 2" << std::endl;
	} // end if (!getenv("skipimages")
		} // end parallel region

	//std::cout << "\n ---- just left parallel region ----" << std::endl;
	std::cout << "params lambda_0: " << params.lambda_0() << std::endl;
	fftw_destroy_plan(forward);
	fftw_destroy_plan(backward);

	tstop = std::time(nullptr);
	tstop -= tstart;
	std::cout << "\t\t======================================\n"
		<< "\t\t======== scan_material stopped =======\n"
		<< "\t\t===== " << std::asctime(std::localtime(&tstop)) 
		<< "\t\t===== in " << tstop << " s \n"
		<< "\t\t======================================\n" << std::flush;
	std::string timesfilename = params.filebase() + "runtimes.log";
	std::ofstream timesout(timesfilename.c_str(),std::ios::app);
	timesout << "#########################################\n" << std::flush;
	timesout << "# HOSTNAME:\t" << getenv("HOSTNAME") << "\n" << std::flush;
	timesout << "# total time (seconds):\t" << tstop << "\n" << std::flush;
	timesout << "# npulses :\t" << params.npulses() << "\n" << std::flush;
	timesout << "# nthreads :\t" << nthreads << "\n" << std::flush;
	timesout << "# mean time / run:\t" << DataOps::mean(runtimes) 
	timesout << "\n" << std::flush;
	timesout.close();

	return 0;
}


