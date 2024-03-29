#define _USE_MATH_DEFINES
#include "param.hpp"
#include "coarseFrameSignalProcess.hpp"
#include <math.h>
#include <cmath>


void rectWin(vector<float>& win_rect, uint16_t N) {
	for (int i = 0; i < N; i++) {
		win_rect[i] = 1.0;
	}
}

void hanningWin(vector<float> &win_hanning, uint16_t N) {
	for (int i = 0; i < N; i++) {
		win_hanning[i] = 0.5 * (1 - cos(2 * M_PI * i / (N - 1)));
	}
}

void hammingWin(vector<float>& win_hamming, uint16_t N) {
	for (int i = 0; i < N; i++) {
		win_hamming[i] = 0.54 - 0.46 * cos(2 * M_PI * i / (N - 1));
	}
}

void func_winR_process(vector<vector<vector<float>>>& WinRout_RangeSampleNum_ChirpNum_RxNum, vector<uint8_t>& fft_win_type, uint16_t win_len, uint16_t chirpNum, uint16_t rxNum, unique_ptr<int16_t[]>& radarInputdata) {
	// paraWin Generate
	uint8_t winr_type = fft_win_type[0];
	vector<float> win_Coef(win_len);
	switch (winr_type)
	{
	case 0:
	{
		rectWin(win_Coef, win_len);
		break;
	}
	case 1:
	{
		hanningWin(win_Coef, win_len);
		break;
	}
	case 2:
	{
		hammingWin(win_Coef, win_len);
		break;
	}
	default:
		break;
	}
	// add win
	for (int chirpIdx = 0; chirpIdx < chirpNum; chirpIdx++) 
	{
		for (int rxIdx = 0; rxIdx < rxNum; rxIdx++)
		{
			for (int sampIdx = 0; sampIdx < win_len; sampIdx++)
			{
				WinRout_RangeSampleNum_ChirpNum_RxNum[sampIdx][chirpIdx][rxIdx] = win_Coef[sampIdx] * radarInputdata[sampIdx + rxIdx * win_len + chirpIdx * win_len * rxNum];
			}
		}
	}
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
		func_winR_process(WinRout_RangeSampleNum_ChirpNum_RxNum, para_sys.fft_win_type, para_sys.RangeSampleNum, para_sys.VelocityNum, para_sys.RxNum, radarInputdata);

	}
	else if (para_sys.InputHasDonePreProcess == 1) {
		// reshape 
	}
	else {
		throw invalid_argument("Input value is illegal.");
	}

	FFTD_SpatialFFT_CFAR_CoarseFrame(point_info, TOI, para_sys, virtual_array, compensate_mat, CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum);
}