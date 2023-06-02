#include <stdio.h>
#include <tensorflow/c/c_api.h>
#include <random>
//#include <boost/python/numpy.hpp>
#include <hdf5/H5DataType.h>

int main() {
	  printf("Hello from TensorFlow C library version %s\n", TF_Version());
	    return 0;
}
