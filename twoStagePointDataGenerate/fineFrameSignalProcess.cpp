#include "param.hpp"
#include "fineFrameSignalProcess.hpp"

void ModifyRDFFT_SpatialFFT_CFAR_LessMemory_TOI_V2(vector<vector<float>>& result_data, ParaSys& para_sys, Virtual_array& virtual_array, vector<complex<float>>& compensate_mat, vector<vector<float>>& TOI, vector<vector<vector<complex<float>>>>& CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum) {

}

void func_signal_process_fine(vector<vector<float>>& result_data, const string& filename, ParaSys& para_sys, Virtual_array& virtual_array, unique_ptr<int16_t[]>& radarInputdata, vector<complex<float>>& compensate_mat, vector<vector<float>>& TOI) {

	vector<vector<vector<complex<float>>>> CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum(para_sys.CoarseRangeNum, vector<vector<complex<float>>>(para_sys.VelocityNum, vector<complex<float>>(para_sys.RxNum)));

	// result_data is detect_target_data
	if (para_sys.InputHasDonePreProcess == 0)
	{
		// winr and fftr
	}
	else if (para_sys.InputHasDonePreProcess == 1)
	{
		// reshape
	}
	else
	{
		throw invalid_argument("Input value is illegal.");
	}
	ModifyRDFFT_SpatialFFT_CFAR_LessMemory_TOI_V2(result_data, para_sys, virtual_array, compensate_mat, TOI, CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum);
}
