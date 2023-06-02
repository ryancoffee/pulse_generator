// Noise class definition

#ifndef NOISE_H
#define NOISE_H


// standard includes
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <random>

// my headers

// my definitions

template <typename T>
class Noise {

public:
	Noise(T & meanin, T & stdevin) 
	: mean(meanin)
	, stdev(stdevin)
	{
		norm_dist = std::normal_distribution<T>(mean,stdev);
	}
	~Noise(){}

	inline void addnoise(std::vector<T> & vec){
		for (size_t i=0;i<vec.size();++i){
			vec[i] += norm_dist(rng);
		}
	}
	inline T& addnoise(T & value){
		value += norm_dist(rng);
		return value;
	}

private:
	std::random_device rng;
	std::normal_distribution<T> norm_dist;
	T mean;
	T stdev;
};

#endif 
