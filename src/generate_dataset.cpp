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
	size_t numpulses = 1<<10;



	Params params;
	params.lambda_0(atof( getenv("lambda0") ));
	params.lambda_width( atof( getenv("lambda_width") ));
	params.lambda_onoff( atof( getenv("lambda_onoff") ));
	params.tspan((atof( getenv("tspan") ) )/fsPau<double>());
	params.filebase(std::string(getenv("filebase")));


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

	std::time_t tstop = std::time(nullptr);
	std::cout << "\tIt has taken " << (tstop-tstart) << " s so far for initializing masterpulse and building fftw plans\n" << std::flush;

	// Setup the shared pulse arrays
	std::vector< PulseFreq > pulsearray(params.getNpulses(),masterpulse);
	std::vector<float> runtimes(params.getNpulses(),0); // for benchmarking the processors

#pragma omp parallel num_threads(nthreads) shared(masterpulse) private(params)
	{ // begin parallel region 1

		// all non-shared objects must be created inside the parallel section for default is shared if defined outside
		// http://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-data.html

		size_t tid = omp_get_thread_num();

#pragma omp for schedule(dynamic)
		for (size_t n=0;n<params.getNpulses();++n)	{ // outermost loop for npulses to produce //
			std::time_t runstart = std::time(nullptr);
			std::cerr << "\tinside the parallel region 1 for pulses loop n = " << n << " in thread " << tid << "\n" << std::flush;

			PulseFreq pulse(masterpulse);
			params.randomchirp();
			params.randomchirpnoise();
			pulse.addchirp(params.getchirp());						
			pulse.scale(params.getAmp_rand());

			float t0 = params.delays_rand();

			//DebugOps::pushout(std::string("Running pulse " + std::to_string(n) + " for t0 = " + std::to_string(t0) + " in threaded for loop, thread " + std::to_string(tid)));

			pulsearray[f] = pulse;

			std::time_t runstop = std::time(nullptr);

			runtimes[n] = float(runstop - runstart);

			} // outermost loop for npulses to produce //

#pragma omp master
			{
				std::cout << "\t\t############ ending parallel region 1 ###########\n" << std::flush;
			}

#pragma omp barrier

	} // end parallel region 1

	std::cout << "\n ---- just left parallel region ----" << std::endl;
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


