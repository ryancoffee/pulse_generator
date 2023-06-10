#include "generate_dataset.hpp"
#include <cstdlib>
#include <algorithm>
#include <complex>
#include <memory>
#include <fftw3.h>
#include <H5Cpp.h>

#include <vector>
#include <random>
#include <chrono>

#include <DebugOps.hpp>

#include <cstddef> // for nullptr
#include <cstdint> // for writing as an int32_t and int16_t

using namespace Constants;
H5::FloatType h5float( H5::PredType::NATIVE_FLOAT );

/* Here's main */
int main(int argc, char* argv[])
{

	std::time_t tstart = std::time(nullptr);
	std::cout << "\t\t======================================\n"
		<< "\t\t======= gen pulses started ========\n"
		<< "\t\t===== " << std::asctime(std::localtime(&tstart)) 
		<< "\t\t===== on host " << std::getenv("HOSTNAME") << "\n"
		<< "\t\t===== using " << std::getenv("nthreads") << " threads\n"
		<< "\t\t======================================\n" << std::flush;

	unsigned nthreads = (unsigned)atoi( std::getenv("nthreads") );
	size_t numpulses = 1<<10;

	std::string filebase(std::getenv("filebase"));
	float delays_mean = float(atof(std::getenv("delays_mean")));
	float delays_std = float(atof(std::getenv("delays_std")));
	float amp_mean = float(atof(std::getenv("amp_mean")));
	float amp_std = float(atof(std::getenv("amp_std")));

	std::cout << "filebase in main() is " << filebase << std::endl << std::flush;
	Params params(filebase,delays_mean,delays_std,amp_mean,amp_std);
	params.lambda_0(atof( std::getenv("lambda_0") ));
	params.lambda_width( atof( std::getenv("lambda_width") ));
	params.lambda_onoff( atof( std::getenv("lambda_onoff") ));
	params.setTspan((atof( std::getenv("tspan") ) )/fsPau<float>());

	float second = float( atof( getenv("chirp") ) );
	float third = float( atof( getenv("TOD") ) );
	float fourth = float( atof( getenv("FOD") ) );
	float fifth = float( atof( getenv("fifthOD") ) );
	params.initChirp(second,third,fourth,fifth).setnulims(float( atof( std::getenv("nu_low") ) ) , float( atof( std::getenv("nu_high") ) ));


	params.set_lamsamples((size_t)atoi(getenv("lamsamples")))
		.set_gain((float)atoi(getenv("gain")))
		.set_noisescale((float)atof(getenv("noisescale") ))
		.set_sampleinterval((size_t)atoi(getenv("sampleinterval")))
		.set_saturate(uint16_t( atoi( getenv("saturate"))));



	std::cout << "initializing masterpulse and masterplans" << std::endl << std::flush;
	PulseFreq masterpulse(params);

	fftw_plan forward;
	fftw_plan backward;
	fftw_plan plan_r2hc;
	fftw_plan plan_hc2r;
	fftw_plan plan_r2hc_2x;
	fftw_plan plan_hc2r_2x;
	masterpulse.setmasterplans(&forward,&backward).setmasterancillaryplans(& plan_r2hc,& plan_hc2r,& plan_r2hc_2x,& plan_hc2r_2x);

	std::time_t tstop = std::time(nullptr);
	std::cout << "\tIt has taken " << (tstop-tstart) << " s so far for initializing masterpulse and building fftw plans\n" << std::flush;

	// Setup the shared pulse arrays
	std::vector<float> runtimes(params.getNpulses(),0); // for benchmarking the processors

	std::vector< std::vector< float > > waves; 
	for (size_t i=0;i<params.getNpulses();i++){
		waves.push_back(std::vector< float >(masterpulse.getsamples(),0.));
	}

#pragma omp parallel num_threads(nthreads) shared(masterpulse,waves) private(params)
	{ // begin parallel region 1

		// all non-shared objects must be created inside the parallel section for default is shared if defined outside
		// http://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-data.html

		size_t tid = omp_get_thread_num();
#pragma omp master
			{
				std::cout << "\t\t############ ending parallel region 1 ###########\n" << std::flush;
			}
		PulseFreq pulse(params);
		pulse.setplans(masterpulse)
			.setancillaryplans(masterpulse);

#pragma omp for schedule(dynamic)
		for (size_t n=0;n<params.getNpulses();++n)	{ // outermost loop for npulses to produce //
			std::time_t runstart = std::time(nullptr);
			std::cerr << "\tinside the parallel region 1 for pulses loop n = " << n << " in thread " << (int)tid << "\n" << std::flush;

			pulse.addchirp(params.getChirp())
				.scale(params.getAmp())
				.delay(params.getDelay());

			pulse.fft_totime().filltime(n,waves);

			std::time_t runstop = std::time(nullptr);
			runtimes[n] = float(runstop - runstart);
			} // outermost loop for npulses to produce //
#pragma omp barrier
#pragma omp master
			{
				std::cout << "\t\t############ ending parallel region 1 ###########\n" << std::flush;
			}
	} // end parallel region 1


	std::cout << "\n ---- just left parallel region ----" << std::endl;
	std::cout << "\n ---------- runtimes are -----------" << std::endl;
	for (size_t i=0;i<runtimes.size();i++)
		std::cout << runtimes[i] << "\t";
	std::cout << std::endl;

	std::cout << "params lambda_0: " << params.lambda_0() << std::endl;
	std::cerr << "Destroying plans" << std::endl << std::flush;
	fftw_destroy_plan(forward);
	fftw_destroy_plan(backward);
	forward = backward = NULL;



	std::cout << "\n ---- Writing H5 file ----" << std::endl;

	H5::IntType h5uint16( H5::PredType::NATIVE_USHORT );
	h5uint16.setOrder( H5T_ORDER_LE );
	H5::IntType h5int16( H5::PredType::NATIVE_SHORT );
	h5int16.setOrder( H5T_ORDER_LE );
	H5::IntType h5uint32( H5::PredType::NATIVE_UINT );
	h5uint32.setOrder( H5T_ORDER_LE );
	H5::IntType h5int32( H5::PredType::NATIVE_INT );
	h5int32.setOrder( H5T_ORDER_LE );
	H5::StrType h5string(0, H5T_VARIABLE);


	std::tm * local_time = std::localtime(nullptr);

	std::string fname(params.filebase());
        std::stringstream tstartstream,tstopstream;
        std::stringstream ss;
        ss      << "-" << local_time->tm_year + 1900
                << "-" << local_time->tm_mon + 1
                << "-" << local_time->tm_mday
                << "-h"<< local_time->tm_hour
                << "-m"<< local_time->tm_min
                << ".h5";
        fname += ss.str();
        std::cout << "outfile = " << fname << std::endl;
       	H5::H5File * hfilePtr = new H5::H5File ( fname , H5F_ACC_TRUNC );
        H5::Group * sansPtr = new H5::Group( hfilePtr->createGroup( "sans" )); // sans noise
        H5::Group * avecPtr = new H5::Group( hfilePtr->createGroup( "avec" )); // avec noise

	const uint8_t rank(1);
	hsize_t dims[1] = {waves.back().size()};
	for (size_t n=0;n<waves.size();n++){
		std::string pulsename = "pulse_" + std::to_string((int)n) + "_real";
		H5::DataSpace * dataspace = new H5::DataSpace( rank , dims );
		H5::DataSet * datasetPtr = new H5::DataSet( sansPtr->createDataSet( pulsename, h5float, *dataspace ) );
		datasetPtr->write( waves[n].data(), h5float);

        	delete datasetPtr;
        	delete dataspace;
	}
        delete sansPtr;
        delete avecPtr;
        delete hfilePtr;

	std::time_t allstop = std::time(nullptr);
	std::cout << "\t\t======================================\n"
		<< "\t\t======== generate pulses stopped =======\n"
		<< "\t\t===== " << std::asctime(std::localtime(&allstop)) 
		<< "\t\t===== in " << float(allstop-tstart) << " s ====\n"
		<< "\t\t======================================\n" << std::flush;
	/*
	std::string timesfilename = params.filebase() + "runtimes.log";
	std::ofstream timesout(timesfilename.c_str(),std::ios::app);
	timesout << "#########################################\n" << std::flush;
	timesout << "# HOSTNAME:\t" << getenv("HOSTNAME") << "#####\n" << std::flush;
	timesout << "# total time (seconds):\t" << tstop << "#######\n" << std::flush;
	timesout << "# npulses :\t" << params.getNpulses() << "#####\n" << std::flush;
	timesout << "# nthreads :\t" << nthreads << "#####\n" << std::flush;
	timesout << "# mean time / run:\t" << DataOps::mean(runtimes) << "#####\n" << std::flush;
	timesout << "#########################################\n" << std::flush;
	timesout.close();
	*/

	return 0;
}


