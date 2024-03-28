#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <complex>
#include <stdexcept>

#include "param.hpp"
using namespace std;

// 1. read_bin √
// 2. adc or 1dfft（real or complex） vector<double>(3) or data[size]
// 3. 1dfft() -> 点云选取

// input:  binary file （parasys and input data and compensate_mat （ adc or
// 1dfft）） output: r, v, theta, snr, point_num

// 根据frametype选择
// coarse
// fine

// libs: fft, blas
void func_signal_process_coarse(vector<vector<float>>& TOI, vector<vector<float>>& point_info, const string &filename,ParaSys &para_sys, Virtual_array &virtual_array, unique_ptr<int16_t[]> &radarInputdata, vector<complex<float>> &compensate_mat) {

}

void func_signal_process_fine(vector<vector<float>>& result_data, const string &filename, ParaSys& para_sys, Virtual_array& virtual_array, unique_ptr<int16_t[]>& radarInputdata, vector<complex<float>>& compensate_mat, vector<vector<float>>& TOI) {
	// result_data is detect_target_data

}

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

	vector<vector<float>> point_info;
	vector<vector<float>> result_data;
	uint16_t TOI_max = 5000;
	uint16_t point_max = 5000;
	uint16_t result_max = 15000;
	// 2. deal
	// Result fun(Binfile &info);
	if (bin_file.para_sys.frame_type == 0) {
		TOI.resize(TOI_max,vector<float>(4));
		if (bin_file.para_sys.work_mode == 0) {
			point_info.resize(point_max, vector<float>(11));
		}
		else if (bin_file.para_sys.work_mode == 1) {
			point_info.resize(point_max, vector<float>(6));
		}
		else
		{
			throw invalid_argument("Input value is illegal.");
		}

		func_signal_process_coarse(TOI, point_info, filename, bin_file.para_sys, virtual_array, bin_file.input_data, compensate);

	}
	else if (bin_file.para_sys.frame_type == 1) {
		result_data.resize(result_max, vector<float>(9));
		// analysis_step = 3
		if (bin_file.para_sys.analysis_step == 3) {
			func_signal_process_fine(result_data, filename, bin_file.para_sys, virtual_array, bin_file.input_data, compensate, TOI);
		}
		
	}
	else
	{
		throw  invalid_argument("Input value is illegal.");
	}
	return 0;
}