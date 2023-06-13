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

/* Here's main */
//int main(int argc, char* argv[])
int main(void)
{
	std::time_t tstart = std::time(nullptr);
	std::cout << "\t\t======================================\n"
		<< "\t\t======= gen pulses started ========\n"
		<< "\t\t===== " << std::asctime(std::localtime(&tstart)) 
		<< "\t\t===== on host " << std::getenv("HOSTNAME") << "\n"
		<< "\t\t===== using " << std::getenv("nthreads") << " threads\n"
		<< "\t\t======================================\n" << std::flush;

	unsigned nthreads = (unsigned)atoi( std::getenv("nthreads") );

	std::string filebase(std::getenv("filebase"));
	float delays_mean = float(atof(std::getenv("delays_mean")));
	float delays_std = float(atof(std::getenv("delays_std")));
	float amp_mean = float(atof(std::getenv("amp_mean")));
	float amp_std = float(atof(std::getenv("amp_std")));

	std::cout << "filebase in main() is " << filebase << std::endl << std::flush;

	float lam0 = float(atof( std::getenv("lambda_0") ));
	float lamw = float(atof( std::getenv("lambda_width") ));
	float lamonoff = float(atof( std::getenv("lambda_onoff") ));
	float tspn = float(atof( std::getenv("tspan") ));

	float second = float( atof( getenv("chirp") ) );
	float third = float( atof( getenv("TOD") ) );
	float fourth = float( atof( getenv("FOD") ) );
	float fifth = float( atof( getenv("fifthOD") ) );
	float nulow = float( atof( std::getenv("nu_low") ) );
	float nuhigh = float( atof( std::getenv("nu_high") ) );



	size_t npulses = (size_t)atoi( std::getenv("npulses"));
	size_t lamsamples = (size_t)atoi(getenv("lamsamples"));
	size_t gain = (size_t)atoi(getenv("gain"));
	float noisescale = (float)atof(getenv("noisescale") );
	size_t sampleinterval = (size_t)atoi(getenv("sampleinterval"));
	uint32_t saturate = uint16_t( atoi( getenv("saturate")));

	std::vector<float> runtimes(nthreads,0); // for benchmarking the processors

	Params masterparams(filebase,delays_mean,delays_std,amp_mean,amp_std);
	masterparams.lambda_0(lam0)
		.lambda_width(lamw)
		.lambda_onoff(lamonoff)
		.setTspan(tspn/fsPau<float>())
		.setNpulses(npulses);

	masterparams.initChirp(second,third,fourth,fifth).setnulims(nulow,nuhigh);

	PulseFreq masterpulse(masterparams);

	std::vector< std::vector< float > > waves; 
	for (size_t i=0;i<nthreads*masterparams.getNpulses();i++){
		waves.push_back(std::vector< float >(masterpulse.getsamples(),0.));
	}
	std::cerr << "waves.size()\t" << (int)(waves.size()) << "\twaves.front().size()\t" << (int)(waves.front().size()) << std::endl << std::flush;

#pragma omp parallel num_threads(nthreads) shared(waves)
	{ // begin parallel region 1

		// all non-shared objects must be created inside the parallel section for default is shared if defined outside
		// http://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-data.html

		size_t tid = omp_get_thread_num();
		std::cerr << "\tentered thread " << (int)tid << "\n" << std::flush;

		Params params(filebase,delays_mean,delays_std,amp_mean,amp_std);
		params.lambda_0(lam0)
			.lambda_width(lamw)
			.lambda_onoff(lamonoff)
			.setTspan(tspn/fsPau<float>())
			.setNpulses(npulses);

		params.initChirp(second,third,fourth,fifth).setnulims(nulow,nuhigh);
		params.set_lamsamples(lamsamples)
			.set_gain(gain)
			.set_noisescale(noisescale)
			.set_sampleinterval(sampleinterval)
			.set_saturate(saturate);

#pragma omp barrier
		std::cout << "initializing pulse and plans" << std::endl << std::flush;
		PulseFreq pulse(params);

		fftw_plan forward;
		fftw_plan backward;
		fftw_plan plan_r2hc;
		fftw_plan plan_hc2r;
		fftw_plan plan_r2hc_2x;
		fftw_plan plan_hc2r_2x;

		std::cerr << "Here" << std::endl << std::flush;

		pulse.setmasterplans(&forward,&backward).setmasterancillaryplans(& plan_r2hc,& plan_hc2r,& plan_r2hc_2x,& plan_hc2r_2x);
#pragma omp barrier

		std::cerr << "Here too" << std::endl << std::flush;


		std::time_t tstop = std::time(nullptr);
		std::cout << "\tIt has taken " << (tstop-tstart) << " s for initializing pulse and fftw plans in tid " << int(tid) << std::endl << std::flush;

		std::time_t runstart = std::time(nullptr);
#pragma omp barrier
//#pragma omp for 
		for (size_t n=0;n<params.getNpulses();++n)	{ // outermost loop for npulses to produce //
		//for (size_t n=0;n<10;++n)	{ // outermost loop for npulses to produce //
			std::cerr << "\tloop n = " << (int)n << " in thread " << (int)tid << "\n" << std::flush;

			pulse.scale(params.getAmp());
			pulse.delay(params.getDelay());
			std::cerr << "Here in main" << std::endl << std::flush;
			std::vector<float> chirpvec(4,0.);
			pulse.addchirp(params.getChirp(chirpvec));


			std::cerr << "Here too in main" << std::endl << std::flush;
			pulse.fft_totime().filltime(waves[tid*nthreads + n]);
			std::cerr << "tid = " << (int)tid << std::endl;

			} // outermost loop for npulses to produce //
		std::time_t runstop = std::time(nullptr);
		runtimes[tid] = float(runstop - runstart);

#pragma omp barrier
		
		std::cout << "\t\t############ ending parallel region 1 ###########\n" << std::flush;
		std::cout << "\t\t################## thread id "<< int(tid) << " ##################\n" << std::flush;
		std::cout << "\t\t##############params lambda_0: " << params.lambda_0() << "##############\n" << std::flush;

		std::cerr << "\t\t############## Destroying plans ################\n" << std::flush;
		fftw_destroy_plan(forward);
		fftw_destroy_plan(backward);
		forward = backward = NULL;

	} // end parallel region 1


	std::cout << "\n ---- just left parallel region ----" << std::endl;
	std::cout << "\n ---------- runtimes are -----------" << std::endl;
	for (size_t i=0;i<runtimes.size();i++)
		std::cout << runtimes[i] << "\t";
	std::cout << std::endl;



	std::cout << "\n ---- Writing H5 file ----" << std::endl;

	H5::FloatType h5float( H5::PredType::NATIVE_FLOAT );
	H5::IntType h5uint16( H5::PredType::NATIVE_USHORT );
	h5uint16.setOrder( H5T_ORDER_LE );
	H5::IntType h5int16( H5::PredType::NATIVE_SHORT );
	h5int16.setOrder( H5T_ORDER_LE );
	H5::IntType h5uint32( H5::PredType::NATIVE_UINT );
	h5uint32.setOrder( H5T_ORDER_LE );
	H5::IntType h5int32( H5::PredType::NATIVE_INT );
	h5int32.setOrder( H5T_ORDER_LE );
	H5::StrType h5string(0, H5T_VARIABLE);


	std::time_t local = std::time(nullptr);
	std::tm * local_time = std::localtime(&local);

        std::stringstream ss;
        ss      << "-" << local_time->tm_year + 1900
                << "-" << local_time->tm_mon + 1
                << "-" << local_time->tm_mday
                << "-h"<< local_time->tm_hour
                << "-m"<< local_time->tm_min
                << ".h5";
	std::string fname = filebase + ss.str();
        std::cout << "outfile = " << fname << std::endl;
       	H5::H5File * hfilePtr = new H5::H5File ( fname , H5F_ACC_TRUNC );
        H5::Group * sansPtr = new H5::Group( hfilePtr->createGroup( "sans" )); // sans noise
        H5::Group * avecPtr = new H5::Group( hfilePtr->createGroup( "avec" )); // avec noise

	const uint8_t rank(1);
	hsize_t dims[1];
       	dims[0] = waves.back().size();
	for (size_t n=0;n<waves.size();n++){
		std::string pulsename = "/sans/pulse_" + std::to_string((int)n) + "_real";
		// std::cerr << "n = " << (int)n << "\tpulsename = " << pulsename << std::endl << std::flush;
		H5::DataSpace * dataspace = new H5::DataSpace( rank , dims );
		H5::DataSet * datasetPtr = new H5::DataSet( sansPtr->createDataSet( pulsename, h5float, *dataspace ) );


		std::vector<float>::iterator maxelement;
		maxelement = std::max_element(waves[n].begin(),waves[n].end());
		std::cerr << "address " << std::distance(waves[n].begin(),maxelement) << ", value = " << *(maxelement) << "\n";
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


