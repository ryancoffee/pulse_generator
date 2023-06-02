#ifndef DATAOPS_H
#define DATAOPS_H

#include <map>
#include <iterator>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
#include <numeric>
#include <math.h>
#include <cmath>
#include <complex>
#include <boost/multi_array.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <typeinfo>

namespace DataOps 
{
	template <typename T>
		std::vector<T>& sample_every(std::vector<T>& lhs, const std::vector<T>& rhs, const size_t nskip = 10, const size_t offset = 0)
		{
			if (lhs.size() < 1){
				lhs.resize(int(rhs.size()/nskip),T(0));
			}
			if (lhs.size() > (rhs.size()+offset)/nskip ){
				lhs.resize(size_t((rhs.size()+offset)/nskip));
			}
			for (size_t i = 0; i<lhs.size(); ++i){
				lhs[i] = T(rhs[nskip*i+offset]);
			}
			return lhs;
		}

	template <typename T>
		T* clone(T* lhs,const T* rhs,const size_t n)
		{
			std::copy(rhs, rhs + n, lhs);
			return lhs;
		}

	template <typename T>
		std::vector<T>& clone(std::vector<T> & lhs, const std::vector<T> & rhs)
		{
			lhs.resize(rhs.size());
			std::copy(rhs.begin(),rhs.end(),lhs.begin()); 
			return lhs;
		}

	template <typename T>
		T& gauss(const T & xin, const T& x0, const T& w,T& a =T(1),T& y0 = T(0)){
			T x = xin-x0;
			return (T) (y0+a*std::exp(- std::pow(x/w,int(2)))) ;
		}

	template <typename X_t,typename Y_t> 
		Y_t polynomial(const X_t xin,const std::vector<Y_t> c){
			assert(c.size()>1);
			Y_t x = Y_t(xin) - c.front();
			Y_t y = c[1];
			for (unsigned i=1; i<c.size();++i){
				y += c[i] * std::pow(x,(int)i-1);
			}
			return y;
		}

	template <typename T>
		void fixedfilter(T * vec, const size_t sz){
			T noise(T(.1));
			vec[0] *= T(1);
			vec[sz/2] *= T(1)/(T(1) + noise/10*std::pow(T(sz/2),2));
			for( size_t i=1;i<sz/2;++i){
				vec[i] *= T(1)/(T(1) + noise/T(10)*std::pow(T(i),2));
				vec[sz-i] *= T(1)/(T(1) + noise/T(10)*std::pow(T(i),2));
			}
			return;
		}

	template <typename T>
		bool inwin(const T val, const std::vector<T> win){
			assert(win.size()==2);
			return (val >= win[0] && val < win[1]);
		}

	template <typename T>
		void endstozero(T * vec, const size_t sz, const size_t steps)
		{
			T mean1(0);
			T mean2(0);
			T x1(0);
			T x2(0);
			for (size_t i=0;i<steps;++i){
				mean1 += vec[i];
				x1 += T(i);
				mean2 += vec[sz-i-1];
				x2 += T(i);
			}
			for (size_t i=0; i<sz; ++i){
				vec[i] -= (mean2-mean1)/(x2-x1)*(T(i)-x1) + mean1/T(steps);
			}
		}
	
	template <typename T>
		void sin2roll(T * vec, const size_t sz, const T center, const T width)
		{
			for (size_t i = 0;i<sz;++i){
				T x = T(M_PI)*(T(i)-center)/(T(2)*width);
				if (std::abs(x) > T(1)){
					vec[i] = T(0);
				} else {
					vec[i] *= std::pow(std::cos(x),int(2));
				}
			}
			return;
		}

	template <typename T>
		void gaussroll(T * vec, const size_t sz, const T center, const T width)
		{
			for (size_t i = 0;i<sz;++i){
				T x = (T(i)-center)/width;
				vec[i] *= std::exp(-1.*std::pow(x,2));
			}
			return;
		}

	template <typename T,typename V>
		V interpolate(const std::map<T,V> &data, T x)
		{
			typedef typename std::map<T,V>::const_iterator MapIterator;
			MapIterator i = data.upper_bound(x);
			MapIterator l = i;
			double slope;
			double span = (double) (x - l->first);
			if (i==data.begin())
			{
				MapIterator u=i;
				++u;
				slope = (double)(u->second - i->second)/(double)(u->first - i->first);
				return (V)( i->second + (slope * span/(u->first - i->first) ) );
			}
			--l;
			if(i==data.end())
			{
				MapIterator ll=l;
				--ll;
				slope = (double)(l->second - ll->second)/(double)(l->first - ll->first);
				return (V)( l->second + (slope * span/(l->first - ll->first) ) ) ;
			}
			slope = (double)( i->second - l->second ) / (double)( i->first - l->first );
			return (V)(l->second + (slope * span/(i->first - l->first)) );
		}

	template <typename T>
		void condense(std::vector< T > & invec,const unsigned newlen)
		{
			double ratio = ((double)invec.size()/(double)newlen);

			if (ratio <= 1){std::cerr << "Cannot condense2(): size mismatch" << std::endl;return;}
			for (unsigned j=1;j<(int)ratio;++j){
				invec[0] += invec[j];
			}
			invec[0] /= ratio;
			for (unsigned i=1;i<newlen;++i){
				invec[i] = invec[(int)(i*ratio)];
				for (unsigned j=1;j<ratio;++j){
					invec[i] += invec[(int)(i*ratio + j)];
				}
				invec[i] /= (double)ratio;
			}
			invec.resize(newlen);

			return;
		}

	template <typename T>
		void condense(std::vector< std::vector < T > > & inmat, const unsigned newlen)
		{
			for (unsigned i=0;i<inmat.size();++i)
				condense(inmat[i],newlen);
			return;
		}

	template <typename T>
		std::istream& operator >> (std::istream & ins, std::vector<T> & record)
		{
			record.clear();
			std::string line;
			getline( ins, line );
			const char head='#';
			if (line.find(head) != std::string::npos){
				return ins;
			}
			std::istringstream iss( (std::string)line );
			T value;
			while (iss >> value){
				record.push_back(value);
			}
			return ins;
		}

	template <typename T>
		std::istream& operator >> (std::istream & ins, std::vector< std::vector<T> > & matrix)
		{
			matrix.clear();
			std::vector<T> record;
			while (ins >> record)
			{
				if (record.size() > 0)
					matrix.push_back( record );
			}
			return ins;
		}
	/*
	   template <typename T, int dim>
	   std::ostream& operator << (std::ostream & outs, boost::multi_array<T,dim> & record)
	   {
	   for (unsigned i=0;i<record.shape()[0];++i){
	   outs << record[i] << "\t";
	   }
	   outs << std::flush;
	   return outs;
	   }
	 */

	template <typename T,typename V>
		std::ostream& operator << (std::ostream & outs, std::pair<T,V> outpair)
		{
			outs << "(" << outpair.first << "," << outpair.second << ")\t" << std::flush;
			return outs;
		}

	template <typename T>
		std::ostream& operator << (std::ostream & outs, std::vector< T > & vec)
		{
			for (unsigned i=0;i<vec.size();++i){
				outs << vec[i] << "\t";
			}
			outs << "\n" << std::flush;
			return outs;
		}

	template <typename T>
		std::ostream& operator << (std::ostream & outs, std::vector< std::vector <T> > & matrix)
		{
			for (unsigned i=0;i<matrix.size();++i){
				for (unsigned j=0;j<matrix[i].size();++j){
					outs << matrix[i][j] << "\t";
				}
				outs << "\n";
			}
			outs << std::flush;
			return outs;
		}

	template <typename T>
		std::complex<T>*& sum(std::complex<T>*& lhs,std::complex<T>* const & rhs,const size_t n)
		{
			for (size_t i = 0;i<n;++i){
				lhs[i] += rhs[i];
			}
			/*
			std::transform(lhs,lhs+n,rhs,[](std::complex<T>* d_lhs, std::complex<T>* d_rhs){
					return *d_lhs + *d_rhs;
					});
			*/
			return lhs;
		}

	template <typename T>
		std::complex<T>* diff(std::complex<T>*& lhs,std::complex<T>* const & rhs,const size_t n)
		{
			for (size_t i = 0;i<n;++i){
				lhs[i] -= rhs[i];
			}
			/*
			std::transform(lhs,lhs+n,rhs,[](std::complex<T>* d_lhs, std::complex<T>* d_rhs){
					return *d_lhs - *d_rhs;
					});
			*/
			return lhs;
		}

	template <typename T>
		std::complex<T>* mul(std::complex<T>*& lhs,std::complex<T>* const & rhs,const size_t n)
		{
			for (size_t i = 0;i<n;++i){
				lhs[i] *= rhs[i];
			}
			/*
			std::transform(lhs,lhs+n,rhs,lhs,[](std::complex<T>* d_lhs, std::complex<T>* d_rhs){
					return (*d_lhs) * (*d_rhs);
					});
			*/
			return lhs;
		}

	template <typename T>
		std::complex<T>* mul(std::complex<T>*& lhs,const T& scale,const size_t n)
		{
			std::transform(lhs, lhs+n, lhs, std::bind2nd(std::multiplies<std::complex<T> >(),scale));
			//std::cerr << "using this complex vec scaling function in DataOps::" << std::endl << std::flush;
			return lhs;
		}
	template <typename T>
		std::complex<T>* mul(std::complex<T>*& lhs,std::complex<T>& scale,const size_t n)
		{
			std::transform(lhs, lhs+n, lhs, std::bind2nd(std::multiplies<std::complex<T> >(),scale));
			return lhs;
		}


	template <typename T>
		std::complex<T>* div(std::complex<T>*& lhs,std::complex<T>* const & rhs,const size_t n)
		{
			for (size_t i=0;i<n;++i){
				lhs[i] /= rhs[i];
			}
			/*
			std::transform(lhs,lhs+n,rhs,lhs,[](std::complex<T>* d_lhs, std::complex<T>* d_rhs){
					return *d_lhs / *d_rhs;
					});
			*/
			return lhs;
		}

	template <typename T>
		inline bool inwin(std::vector<T> win,T val)
		{
			std::sort(win.begin(),win.end());
			return (val >= win.front() && val < win.back());
		}

	template <typename T>
		std::vector<T> & expscaled(const T scale, std::vector<T> & x)
		{
			x *= scale;
			std::transform(x.begin(), x.end(), x.begin(), [&](T xval){return std::exp(xval);});
			return x;
		}

	template <typename T>
		std::vector<std::complex<T> > & operator *= (std::vector<std::complex<T> > & vec,std::complex<T> val)
		{
			std::transform(vec.begin(), vec.end(), vec.begin(), std::bind2nd(std::multiplies<std::complex<T> >(),val));
			return vec;
		}

	template <typename T>
		std::vector<std::complex<T> > & operator *= (std::vector<std::complex<T> > & vec,T val)
		{
			std::transform(vec.begin(), vec.end(), vec.begin(), std::bind2nd(std::multiplies<std::complex<T> >(),val));
			return vec;
		}

	template <typename T>
		std::vector<T> & operator *= (std::vector<T>& vec,T const val)
		{
			std::transform(vec.begin(), vec.end(), vec.begin(), std::bind2nd(std::multiplies<T>(),val));
			return vec;
		}

	//================ HERE HERE HERE HERE ================//
	template <typename T>
		std::vector<T> & operator /= (std::vector<T>& vec,const T val)
		{
			assert(std::abs(val) != T(0) );
			std::transform(vec.begin(), vec.end(), vec.begin(), std::bind2nd(std::divides<T>(),val));
			return vec;
		}

	template <typename T>
		std::vector<T> & operator += (std::vector<T>& vec, const T val)
		{
			std::transform(vec.begin(), vec.end(), vec.begin(), std::bind2nd(std::plus<T>(),val));
			return vec;
		}

	template <typename T,typename T2>
		std::vector<T> & operator -= (std::vector<T>& vec, const T2 val)
		{
			std::transform(vec.begin(), vec.end(), vec.begin(), std::bind2nd(std::minus<T>(),val));
			return vec;
		}

	template <typename T>
		std::vector<T>& operator += (std::vector<T>& resvec,const std::vector<T> & srcvec)
		{
			assert(resvec.size() == srcvec.size());
			for (unsigned i=0;i<resvec.size();++i){
				resvec[i] += srcvec[i];
			}
			return resvec;
		}

	template <typename T>
		std::vector<T>& operator -= (std::vector<T>& resvec,const std::vector<T> & srcvec)
		{
			assert(resvec.size() == srcvec.size());
			for (unsigned i=0;i<resvec.size();++i){
				resvec[i] -= srcvec[i];
			}
			return resvec;
		}
	template <typename T>
		std::vector<T>& operator *= (std::vector<T>& resvec,const std::vector<T> & srcvec)
		{
			assert(resvec.size() == srcvec.size());
			for (unsigned i=0;i<resvec.size();++i){
				resvec[i] *= srcvec[i];
			}
			return resvec;
		}
	template <typename T>
		std::vector<T>& operator /= (std::vector<T>& resvec,const std::vector<T> & srcvec)
		{
			assert(resvec.size() == srcvec.size());
			for (unsigned i=0;i<resvec.size();++i){
				resvec[i] /= srcvec[i];
			}
			return resvec;
		}



	template <typename T>
		inline T projection(std::vector<T> &in1,std::vector<T> in2,std::vector<bool> & mask){
			assert(in1.size()==in2.size());
			assert(in1.size() == mask.size());
			T ip(0);
			for (unsigned i=0;i<in1.size();++i){
				if (mask[i]){
					ip += in1[i] * in2[i];
				}
			}
			return ip;	
		}

	inline size_t sum(std::vector<bool> & mask){
		size_t sum = 0;
		for (size_t i = 0 ; i<mask.size();++i){
			if (mask[i]==true) sum++;
		}
		return sum;
	}

	template <typename T>
		inline T projection(std::vector<T> &in1,std::vector<T> in2){
			assert(in1.size()==in2.size());
			return std::inner_product(in1.begin(),in1.end(),in2.begin(),0.);	
		}

	template <typename T>
		inline void safe_normalize(std::vector<T> & in)
		{
			T ip(0);
			for (unsigned i=0;i<in.size();++i){
				ip += (in[i] * in[i]);
			}
			assert(ip!=T(0));
			T scale = T(1)/std::sqrt(ip);
			std::transform(in.begin(), in.end(), in.begin(), std::bind2nd(std::multiplies<T>(),scale) );
		}
	
	template <typename T>
		inline void sqr_normalize(std::vector<T> & in,std::vector<bool>& mask)
		{
			assert(in.size() == mask.size());
			T ip = T(0);
			for (unsigned i=0;i<in.size();++i)
				if (mask[i]) {ip += (in[i] * in[i]);}
			assert(ip!=T(0));
			T scale = T(1)/std::sqrt(ip);
			std::transform(in.begin(), in.end(), in.begin(), std::bind2nd(std::multiplies<T>(),scale) );
		}
	
	template <typename T>
		inline void sqr_normalize(std::vector<T> & in) {
			T scale = sqrt( std::inner_product(in.begin(), in.end(), in.begin(),T(0)) );
			assert(scale != T(0));
			std::transform(in.begin(), in.end(), in.begin(), std::bind2nd(std::divides<T>(),scale) );
		};

	template <typename T>
		inline T removemean(T* in, size_t sz) {
			T sum = std::accumulate(in, in + sz, T(0));
			T mean = sum / T(sz);
			std::transform(in, in + sz, in, std::bind2nd(std::minus<T>(),mean));
			return mean;
		};

	template <typename T>
		inline T removemean(std::vector<T> & in,std::vector<bool> & mask) {
			assert(in.size() == mask.size());
			T sum(0);
			unsigned nvals(0);
			for (unsigned i=0;i<in.size();++i){
				if (mask[i]) { 
					sum += in[i];
					++nvals;
				}
			}
			T mean = sum / T(nvals);
			std::transform(in.begin(), in.end(), in.begin(), std::bind2nd(std::minus<T>(),mean));
			return mean;
		};
	
	template <typename T>
		inline T removemean(std::vector<T> & in) {
			T sum = std::accumulate(in.begin(), in.end(), T(0));
			T mean = sum / T(in.size());
			std::transform(in.begin(), in.end(), in.begin(), std::bind2nd(std::minus<T>(),mean));
			return mean;
		};
	
	template <typename T>
		inline T mean(std::vector<T> &in){
			T mean = std::accumulate(in.begin(), in.end(), T(0));
			mean /= T(in.size());
			return mean;
		};




	// WHoah  logarithmic stuff  //

	template <typename T>
		inline void meanstdlog(std::vector<T> & in, T& mean, T& std)
		{
			// I want the log of the values, then sum, then divide by size for mean
			// Diff is the log(values)-mean(log valuse)
			mean = std::log( std::accumulate(in.begin(),in.end(),T(1),std::multiplies<T>()) );
			mean /= T(in.size());
			std::vector<T> diff(in.size());
			for (unsigned i=0;i<diff.size();++i){
				diff[i] = std::log(in[i]) - mean;
			}
			T sq_sum = std::inner_product(diff.begin(),diff.end(),diff.begin(),T(0));
			std = std::sqrt(sq_sum / T(diff.size()));
		}

	template <typename T>
		inline T meanlog(std::vector<T> & in,T offset){
			T mean = offset;
			mean += std::log( std::accumulate(in.begin(),in.end(),T(1),std::multiplies<T>()) );
			mean /= T(in.size());
			return mean;
		};
	
	template <typename T>
		inline T stdlog(std::vector<T> & in, T mean){
			std::vector<T> diff(in.size());
			for (unsigned i=0;i<diff.size();++i){
				diff[i] = std::log(in[i]) - mean;
			}
			T sq_sum = std::inner_product(diff.begin(),diff.end(),diff.begin(),T(0));
			T stdev = std::sqrt(sq_sum / T(diff.size()));
			return stdev;
		};

	template <typename T>
		void detrend(T * vec, const size_t sz)
		{
			T num(0);
			T den(0);
			T xm = T(sz-1)/T(2);
			T ym = removemean(vec,sz);
			for (size_t i=0;i<sz;++i){
				num += (i-xm)*vec[i];
				den += std::pow((i-xm),int(2));
			}
			T beta = num/den;
			for (size_t i=0;i<sz;++i){
				vec[i] -= beta*T(i) - beta*xm;
			}
			return;
		}

	template <typename T>
		T  interpolate(const std::vector<T> &xData, const std::vector<T> & yData,T x)
		{

			if(x < xData.front()){return yData.front();}
			if(x > xData.back()){return yData.back();}
			size_t i = 0;
			while (x>xData[i]){++i;}
			size_t l=i-1;
			T delta=(x - xData[l])/(xData[i] - xData[l]);
			return delta*yData[i] + (1-delta) *yData[l];
		}
}

#endif
