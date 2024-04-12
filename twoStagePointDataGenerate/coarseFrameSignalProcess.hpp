#include "param.hpp"

void func_winA_process(vector<vector<vector<vector<vector<complex<float>>>>>>& winAout_xNum_yNum_VeloFFTNum_RangeNum_waveLocNum, vector<uint8_t>& fft_win_type, uint16_t win_len, uint16_t VirtArrVertGridLen, uint16_t VelocityNum, uint16_t RangeSampleNum, uint16_t waveLocNum, vector<vector<vector<vector<vector<complex<float>>>>>>& result_xNum_yNum_VeloFFTNum_RangeNum);
void func_Spatial_Reorder(vector<vector<vector<vector<vector<complex<float>>>>>>& result_xNum_yNum_VeloFFTNum_RangeNum, vector<vector<vector<vector<complex<float>>>>>& FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum, uint8_t Array_option, uint16_t MIMONum, Virtual_array& virtual_array, uint16_t VirtArrHorGridLen, uint16_t VirtArrVertGridLen, uint16_t VelocityNum, uint16_t CoarseRangeNum, uint16_t waveLocNum);
void func_winD_process(vector<vector<vector<vector<complex<float>>>>>& WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum, vector<uint8_t>& fft_win_type, uint16_t win_len, uint16_t rangeSampleNum, uint16_t MIMONum, uint16_t waveLocNum, vector<vector<vector<vector<complex<float>>>>>& CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum);
void findDuplicates(vector<uint16_t>& output, const std::vector<uint16_t>& nums1, const std::vector<uint16_t>& nums2);
void rectWin(vector<float>& win_rect, uint16_t N);
void hanningWin(vector<float>& win_hanning, uint16_t N);
void hammingWin(vector<float>& win_hamming, uint16_t N);
void func_winR_process(vector<vector<vector<float>>>& WinRout_RangeSampleNum_ChirpNum_RxNum, vector<uint8_t>& fft_win_type, uint16_t win_len, uint16_t chirpNum, uint16_t rxNum, unique_ptr<int16_t[]>& radarInputdata);
void FFTD_SpatialFFT_CFAR_CoarseFrame(vector<vector<float>>& point_info, vector<vector<float>>& TOI, ParaSys& para_sys, Virtual_array& virtual_array, vector<complex<float>>& compensate_mat, vector<vector<vector<complex<float>>>>& CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum);
void func_signal_process_coarse(vector<vector<float>>& TOI, vector<vector<float>>& point_info, const string& filename, ParaSys& para_sys, Virtual_array& virtual_array, unique_ptr<int16_t[]>& radarInputdata, vector<complex<float>>& compensate_mat);
