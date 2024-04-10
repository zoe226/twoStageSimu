#define _USE_MATH_DEFINES
#include "param.hpp"
#include "coarseFrameSignalProcess.hpp"
#include <math.h>
#include <cmath>
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <unordered_set>

void func_winD_process(vector<vector<vector<vector<complex<float>>>>>& WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum,vector<uint8_t>& fft_win_type,uint16_t win_len, uint16_t rangeSampleNum, uint16_t MIMONum, uint16_t waveLocNum, vector<vector<vector<vector<complex<float>>>>>& CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum) {
	uint8_t wind_type = fft_win_type[1];
	vector<float> win_Coef(win_len);
	switch (wind_type)
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
	for (uint16_t waveLocIdx = 0; waveLocIdx < waveLocNum; waveLocIdx++)
	{
		for (uint16_t mimoIdx = 0; mimoIdx < MIMONum; mimoIdx++)
		{
			for (uint16_t sampleIdx = 0; sampleIdx < rangeSampleNum; sampleIdx++)
			{
				for (uint16_t chirpIdx = 0; chirpIdx < win_len; chirpIdx++)
				{
					WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum[chirpIdx][sampleIdx][mimoIdx][waveLocIdx] = win_Coef[chirpIdx] * CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum[chirpIdx][sampleIdx][mimoIdx][waveLocIdx];
				}
			}
		}
	}
}

void findDuplicates(vector<uint16_t>& output, const std::vector<uint16_t>& nums1, const std::vector<uint16_t>& nums2) {
	std::unordered_set<uint16_t> numSet(nums1.begin(), nums1.end());
	for (uint16_t num : nums2) {
		if (numSet.find(num) != numSet.end()) {
			output.push_back(num);
		}
	}
}

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
	for (uint16_t chirpIdx = 0; chirpIdx < chirpNum; chirpIdx++) 
	{
		for (uint16_t rxIdx = 0; rxIdx < rxNum; rxIdx++)
		{
			for (uint16_t sampIdx = 0; sampIdx < win_len; sampIdx++)
			{
				WinRout_RangeSampleNum_ChirpNum_RxNum[sampIdx][chirpIdx][rxIdx] = win_Coef[sampIdx] * radarInputdata[sampIdx + rxIdx * win_len + chirpIdx * win_len * rxNum];
			}
		}
	}
}

void FFTD_SpatialFFT_CFAR_CoarseFrame(vector<vector<float>>& point_info, vector<vector<float>>& TOI, ParaSys& para_sys, Virtual_array& virtual_array, vector<complex<float>>& compensate_mat, vector<vector<vector<complex<float>>>>& CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum) {
	uint16_t MIMONum = para_sys.TxNum * para_sys.RxNum;
	uint16_t VirtArrHorGridLen = 0;
	uint16_t VirtArrVertGridLen = 0;
	VirtArrHorGridLen = virtual_array.pos_in_mat[0][0];
	for (int i = 1; i < MIMONum; i++) {
		VirtArrHorGridLen = std::max(VirtArrHorGridLen, virtual_array.pos_in_mat[i][0]);
	}
	VirtArrHorGridLen = VirtArrHorGridLen + 1;
	VirtArrVertGridLen = virtual_array.pos_in_mat[0][1];
	for (int i = 1; i < MIMONum; i++) {
		VirtArrVertGridLen = std::max(VirtArrVertGridLen, virtual_array.pos_in_mat[i][1]);
	}
	VirtArrVertGridLen = VirtArrVertGridLen + 1;
	vector<vector<vector<complex<float>>>> CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum(para_sys.VelocityNum, vector<vector<complex<float>>>(para_sys.CoarseRangeNum, vector<complex<float>>(MIMONum,0.0)));
	vector<vector<vector<vector<complex<float>>>>> CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum(para_sys.VelocityNum,vector<vector<vector<complex<float>>>>(para_sys.CoarseRangeNum,vector<vector<complex<float>>>(MIMONum,vector<complex<float>>(para_sys.waveLocNum,0.0))));
	vector<vector<vector<vector<complex<float>>>>> WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum(para_sys.VelocityNum, vector<vector<vector<complex<float>>>>(para_sys.CoarseRangeNum, vector<vector<complex<float>>>(MIMONum, vector<complex<float>>(para_sys.waveLocNum, 0.0))));
	vector<vector<vector<vector<complex<float>>>>> FFT2D_VeloFFTNum_CoarseRangeBinNum_MIMONum_WaveLocNum(para_sys.VelocityNum, vector<vector<vector<complex<float>>>>(para_sys.CoarseRangeNum, vector<vector<complex<float>>>(MIMONum, vector<complex<float>>(para_sys.waveLocNum,0.0))));
	vector<vector<vector<vector<complex<float>>>>> SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum(para_sys.AngleHorNum,vector<vector<vector<complex<float>>>>(VirtArrHorGridLen,vector<vector<complex<float>>>(para_sys.VelocityNum,vector<complex<float>>(para_sys.CoarseRangeNum,0.0))));
	vector<vector<vector<vector<complex<float>>>>> SpatialFFT_AngleVertNum_AngleHorNum_VeloFFTNum_RangeBinNum(para_sys.AngleVertNum,vector<vector<vector<complex<float>>>>(para_sys.AngleHorNum,vector<vector<complex<float>>>(para_sys.VelocityNum,vector<complex<float>>(para_sys.CoarseRangeNum,0.0))));
	
	// signal process
	uint8_t CoHerAccu_en = 1;
	uint8_t DDMA_EN = 0;
	if (para_sys.waveform_option == 8 || para_sys.waveform_option == 9) 
	{
		DDMA_EN = 1;
	}
	else
	{
		DDMA_EN = 0;
	}
	if (CoarseFrame_CFARdim > 1)
	{
		// step1:数据重排
		for (uint16_t coarseRangeIdx = 0; coarseRangeIdx < para_sys.CoarseRangeNum; coarseRangeIdx++)
		{
			for (uint16_t mimo_index = 0; mimo_index < MIMONum; mimo_index++)
			{
				uint8_t tx_index = 0;
				uint8_t rx_index = 0;
				tx_index = floor(mimo_index / para_sys.RxNum);
				rx_index = fmod(mimo_index, para_sys.RxNum);
				for (uint16_t seq_index = 0; seq_index < para_sys.TxReuseNum; seq_index++)
				{
					uint16_t chirp_index = para_sys.all_tx_seq_pos[tx_index][seq_index] - 1;
					CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum[chirp_index][coarseRangeIdx][mimo_index] = CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[coarseRangeIdx][chirp_index][rx_index];
				}
			}
		}
		if (CoHerAccu_en == 1)
		{
			for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum; chirpIdx++)
			{
				for (uint16_t coarseRangeIdx = 0; coarseRangeIdx < para_sys.CoarseRangeNum; coarseRangeIdx++)
				{
					for (uint16_t MIMOIdx = 0; MIMOIdx < MIMONum; MIMOIdx++)
					{
						CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum[chirpIdx][coarseRangeIdx][MIMOIdx][0] = CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum[chirpIdx][coarseRangeIdx][MIMOIdx];
					}
				}
			}
		}
		else
		{
			//vector<vector<vector<uint16_t>>> tempWaveLoc(para_sys.VelocityNum/para_sys.TxNum/para_sys.waveLocNum,vector<vector<uint16_t>>(para_sys.waveLocNum,vector<uint16_t>(para_sys.TxNum)));
			for (uint16_t tx_index = 0; tx_index < para_sys.TxNum; tx_index++)
			{
				vector<uint16_t> tempIdx(para_sys.TxReuseNum);
				for (uint16_t seq_index = 0; seq_index < para_sys.TxReuseNum; seq_index++)
				{
					tempIdx[seq_index] = para_sys.all_tx_seq_pos[tx_index][seq_index];
				}
				for (uint16_t waveLocIdx = 0; waveLocIdx < para_sys.waveLocNum; waveLocIdx++)
				{
					vector<uint16_t> chirpSeqPerWave(para_sys.TxReuseNum * para_sys.TxGroupNum / para_sys.waveLocNum);
					for (uint16_t i = 0; i < para_sys.TxReuseNum * para_sys.TxGroupNum / para_sys.waveLocNum; i++)
					{
						chirpSeqPerWave[i] = para_sys.waveLocChirpSeq[waveLocIdx][i];
					}
					vector<uint16_t> chirpSeqPerWaveTx(para_sys.VelocityNum / para_sys.TxNum / para_sys.waveLocNum);
					findDuplicates(chirpSeqPerWaveTx, chirpSeqPerWave, tempIdx);
					/*for (uint16_t j = 0; j < para_sys.VelocityNum / para_sys.TxNum / para_sys.waveLocNum; j++)
					{
						tempWaveLoc[j][waveLocIdx][tx_index] = chirpSeqPerWaveTx[j];
					}*/
					for (uint16_t rx_index = 0; rx_index < para_sys.RxNum; rx_index++)
					{
						uint16_t mimoNum = tx_index * para_sys.RxNum + rx_index;
						for (uint16_t m = 0; m < para_sys.CoarseRangeNum; m++)
						{
							for (uint16_t n = 0; n < para_sys.VelocityNum / para_sys.TxNum / para_sys.waveLocNum; n++)
							{
								CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum[chirpSeqPerWaveTx[n] - 1][m][mimoNum][waveLocIdx] = CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum[chirpSeqPerWaveTx[n] - 1][m][mimoNum];
							}
						}
					}
				}
			}


		}
		// step2:速度维加窗
		func_winD_process(WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum, para_sys.fft_win_type, para_sys.VelocityNum, para_sys.CoarseRangeNum, MIMONum, para_sys.waveLocNum, CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum);
		// test code
		vector<complex<float>>  windChirp(para_sys.CoarseRangeNum);
		for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
		{
			windChirp[sampleIdx] = WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum[5][sampleIdx][24][0];
		}
		// step3:速度维FFT
		for (uint16_t waveIdx = 0; waveIdx < para_sys.waveLocNum; waveIdx++)
		{
			for (uint16_t mimoIdx = 0; mimoIdx < MIMONum; mimoIdx++)
			{
				for (uint16_t rangeIdx = 0; rangeIdx < para_sys.CoarseRangeNum; rangeIdx++)
				{
					fftw_complex* in;
					fftw_complex* out;
					fftw_plan pfftd;
					in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * para_sys.VelocityNum);
					out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * para_sys.VelocityNum);
					for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum; chirpIdx++)
					{
						in[chirpIdx][0] = real(WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum[chirpIdx][rangeIdx][mimoIdx][waveIdx]);
						in[chirpIdx][1] = imag(WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum[chirpIdx][rangeIdx][mimoIdx][waveIdx]);
					}
					pfftd = fftw_plan_dft_1d(para_sys.VelocityNum, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
					fftw_execute(pfftd);
					for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum/2; chirpIdx++)
					{
						FFT2D_VeloFFTNum_CoarseRangeBinNum_MIMONum_WaveLocNum[chirpIdx + para_sys.VelocityNum / 2][rangeIdx][mimoIdx][waveIdx] = complex<float>(out[chirpIdx][0], out[chirpIdx][1]);
					}
					for (uint16_t chirpIdx = para_sys.VelocityNum / 2; chirpIdx < para_sys.VelocityNum; chirpIdx++)
					{
						FFT2D_VeloFFTNum_CoarseRangeBinNum_MIMONum_WaveLocNum[chirpIdx - para_sys.VelocityNum / 2][rangeIdx][mimoIdx][waveIdx] = complex<float>(out[chirpIdx][0], out[chirpIdx][1]);
					}
					fftw_destroy_plan(pfftd);
					fftw_free(in);
					fftw_free(out);
				}
			}
		}
		// test code
		vector<complex<float>> fftdChirp(para_sys.CoarseRangeNum);
		for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
		{
			fftdChirp[sampleIdx] = FFT2D_VeloFFTNum_CoarseRangeBinNum_MIMONum_WaveLocNum[5][sampleIdx][24][0];
		}
	}
	if (CoarseFrame_CFARdim > 2)
	{

	}
	if (CoarseFrame_CFARdim > 3)
	{

	}

	// cfar 

}

void func_signal_process_coarse(vector<vector<float>>& TOI, vector<vector<float>>& point_info, const string& filename, ParaSys& para_sys, Virtual_array& virtual_array, unique_ptr<int16_t[]>& radarInputdata, vector<complex<float>>& compensate_mat) {

	vector<vector<vector<complex<float>>>> CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum(para_sys.CoarseRangeNum, vector<vector<complex<float>>>(para_sys.VelocityNum, vector<complex<float>>(para_sys.RxNum)));

	if (para_sys.InputHasDonePreProcess == 0) {
		// winr 
		vector<vector<vector<float>>> WinRout_RangeSampleNum_ChirpNum_RxNum(para_sys.RangeSampleNum, vector<vector<float>>(para_sys.VelocityNum, vector<float>(para_sys.RxNum)));
		func_winR_process(WinRout_RangeSampleNum_ChirpNum_RxNum, para_sys.fft_win_type, para_sys.RangeSampleNum, para_sys.VelocityNum, para_sys.RxNum, radarInputdata);
		
		// test code
		vector<float> temp(para_sys.RangeSampleNum);
		for (int i = 0; i < para_sys.RangeSampleNum; i++) {
			temp[i] = WinRout_RangeSampleNum_ChirpNum_RxNum[i][45][23];
		}

		// fftr
		//vector<vector<vector<complex<float>>>> CoarseRangeFFT_SampleNum_ChirpNum_RxNum(para_sys.RangeSampleNum, vector<vector<complex<float>>>(para_sys.VelocityNum, vector<complex<float>>(para_sys.RxNum)));
		for (int i = 0; i < para_sys.VelocityNum; i++) {
			for (int j = 0; j < para_sys.RxNum; j++) {
				double* in;
				fftw_complex* out;
				fftw_plan pfft;
				in = (double*)fftw_malloc(sizeof(double) * para_sys.RangeSampleNum);
				out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * para_sys.RangeSampleNum);
				for (int k = 0; k < para_sys.RangeSampleNum; k++) {
					in[k] = WinRout_RangeSampleNum_ChirpNum_RxNum[k][i][j];
				}
				pfft = fftw_plan_dft_r2c_1d(para_sys.RangeSampleNum,in,out,FFTW_ESTIMATE);
				fftw_execute(pfft);
				for (int k = 0; k < para_sys.RangeSampleNum/2; k++) {
					CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[k][i][j] = complex<float>(out[k][0],out[k][1]);
				}
				fftw_destroy_plan(pfft);
				fftw_free(in);
				fftw_free(out);
			}
		}
		// test code
		vector<complex<float>> tempfftr(para_sys.RangeSampleNum/2);
		for (int i = 0; i < para_sys.RangeSampleNum / 2; i++) {
			tempfftr[i] = CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[i][45][23];
		}
	}
	else if (para_sys.InputHasDonePreProcess == 1) {
		// reshape 

	}
	else {
		throw invalid_argument("Input value is illegal.");
	}

	FFTD_SpatialFFT_CFAR_CoarseFrame(point_info, TOI, para_sys, virtual_array, compensate_mat, CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum);
}