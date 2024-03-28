#include "param.hpp"
#include "coarseFrameSignalProcess.hpp"

void func_winR_process(vector<vector<vector<float>>> &WinRout_RangeSampleNum_ChirpNum_RxNum, vector<uint8_t> &fft_win_type, unique_ptr<int16_t[]> &radarInputdata) {
	// nowin,hann,hamm
}

void FFTD_SpatialFFT_CFAR_CoarseFrame(vector<vector<float>>& point_info, vector<vector<float>>& TOI, ParaSys& para_sys, Virtual_array& virtual_array, vector<complex<float>>& compensate_mat, vector<vector<vector<complex<float>>>>& CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum) {

	// signal process

	// cfar 

}

void func_signal_process_coarse(vector<vector<float>>& TOI, vector<vector<float>>& point_info, const string& filename, ParaSys& para_sys, Virtual_array& virtual_array, unique_ptr<int16_t[]>& radarInputdata, vector<complex<float>>& compensate_mat) {

	vector<vector<vector<complex<float>>>> CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum(para_sys.CoarseRangeNum, vector<vector<complex<float>>>(para_sys.VelocityNum, vector<complex<float>>(para_sys.RxNum)));

	if (para_sys.InputHasDonePreProcess == 0) {
		// winr and fftr
		vector<vector<vector<float>>> WinRout_RangeSampleNum_ChirpNum_RxNum(para_sys.RangeSampleNum, vector<vector<float>>(para_sys.VelocityNum, vector<float>(para_sys.RxNum)));
		func_winR_process(WinRout_RangeSampleNum_ChirpNum_RxNum, para_sys.fft_win_type, radarInputdata);
	}
	else if (para_sys.InputHasDonePreProcess == 1) {
		// reshape 
	}
	else {
		throw invalid_argument("Input value is illegal.");
	}

	FFTD_SpatialFFT_CFAR_CoarseFrame(point_info, TOI, para_sys, virtual_array, compensate_mat, CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum);
}