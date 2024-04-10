#include "param.hpp"

void findDuplicates(vector<uint16_t>& output, const std::vector<uint16_t>& nums1, const std::vector<uint16_t>& nums2);
void rectWin(vector<float>& win_rect, uint16_t N);
void hanningWin(vector<float>& win_hanning, uint16_t N);
void hammingWin(vector<float>& win_hamming, uint16_t N);
void func_winR_process(vector<vector<vector<float>>>& WinRout_RangeSampleNum_ChirpNum_RxNum, vector<uint8_t>& fft_win_type, uint16_t win_len, uint16_t chirpNum, uint16_t rxNum, unique_ptr<int16_t[]>& radarInputdata);
void FFTD_SpatialFFT_CFAR_CoarseFrame(vector<vector<float>>& point_info, vector<vector<float>>& TOI, ParaSys& para_sys, Virtual_array& virtual_array, vector<complex<float>>& compensate_mat, vector<vector<vector<complex<float>>>>& CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum);
void func_signal_process_coarse(vector<vector<float>>& TOI, vector<vector<float>>& point_info, const string& filename, ParaSys& para_sys, Virtual_array& virtual_array, unique_ptr<int16_t[]>& radarInputdata, vector<complex<float>>& compensate_mat);
