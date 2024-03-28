#include "param.hpp"

void ModifyRDFFT_SpatialFFT_CFAR_LessMemory_TOI_V2(vector<vector<float>>& result_data, ParaSys& para_sys, Virtual_array& virtual_array, vector<complex<float>>& compensate_mat, vector<vector<float>>& TOI, vector<vector<vector<complex<float>>>>& CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum);
void func_signal_process_fine(vector<vector<float>>& result_data, const string& filename, ParaSys& para_sys, Virtual_array& virtual_array, unique_ptr<int16_t[]>& radarInputdata, vector<complex<float>>& compensate_mat, vector<vector<float>>& TOI);
