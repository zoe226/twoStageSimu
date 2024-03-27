#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <complex>

#include "param.hpp"
using namespace std;



// 1. read_bin
// 2. adc or 1dfft（real or complex） vector<double>(3) or data[size]
// 3. 1dfft() -> 点云选取

// input:  binary file （parasys and input data and compensate_mat （ adc or
// 1dfft）） output: r, v, theta, snr, point_num

// 根据frametype选择
// coarse
// fine

// libs: fft, blas






int main(int argc, char* argv[]) {
	// 1. get info
	// coarse frame data file parse code
	std::string filename = "ADCData_Frame_000000_20240314_17_22_19_2876_ShareMemory.bin";
	BinFile bin_file;
	Virtual_array virtual_array;
	vector<complex<float>> compensate;

	get_info(filename, bin_file, virtual_array);
	
	compensate.reserve(bin_file.para_sys.TxNum * bin_file.para_sys.RxNum);
	for (uint16_t i = 0; i < bin_file.para_sys.TxNum * bin_file.para_sys.RxNum; i++) {
		compensate.emplace_back(bin_file.compensate_mat[2*i], bin_file.compensate_mat[2 * i + 1]);
	}

	// 2. deal
	// Result fun(Binfile &info);

	return 0;
}