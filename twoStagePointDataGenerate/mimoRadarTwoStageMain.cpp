#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
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
struct ConReg {
	int PowerDetEn;
	float SnrThre;
	int RextendNum;
	int VextendNum;
};
struct CFAR {
	int ProCellNum_R_u2;
	int ProCellNum_V_u2;
	int RefCellNum_1D_u5;
};
struct SpatialPara {
	float Vsnr_Thre_db;
	float HorProfThre_db;
	int PeakSen;
	int HorOffsetidx;
	float HorOffPeakThre_db;
	float Global_noise_Thre_db;
	float min_snr_Hor_db;
	float VerThreshold;
	int HorCfarGuardNum;
	int HorCfarReferNum;
	int Horminidx;
	int Hormaxidx;
	int Verminidx;
	int Vermaxidx;
};
struct Union {
	int CompareIdx;
	int MaxValIdx;
};

struct DetPara {
	ConReg con_reg;
	int ChToProcess_Num_u7;
	vector<char> Switch_DataCompress_u2; // (2); ‘11’bin
	int PassAbsData_u1;
	int PassChMeanData_u1;
	int PassPeakS_u1;
	int OneDEnable;
	int Reg_PGA_u15;
	int PeakS_Enable_u1;
	int RangeCellNum_u10;
	int ChirpNum_u11;
	int DetectCell_RIndex_Min_u10;
	int DetectCell_RIndex_Max_u10;
	int DetectCell_VIndex_Min_u11;
	int DetectCell_VIndex_Max_u11;
	vector<char> CFARTypeSwitch_u2; // (2);'11'
	int LogicTestFlag_u1;
	CFAR cfar_para;
	int Loc_OSCFAR_u5;
	int Index_Chirp_NotMove_OSCFAR_u11;
	int Threshold_RangeDim_For_2D_OSCFAR_u9;
	int asicsavedata_flag;
	int IndexCHIP_u3;
	SpatialPara spatial_para;
	Union union_para;
};

struct ParaSys {
	uint8_t vel_target_generated_by_moniqi;
	int need_spatial_cal;
	uint8_t InputHasDonePreProcess;
	uint16_t RangeSampleNum;
	uint16_t TxReuseNum;
	uint16_t TxNum;
	int TxGroupNum;
	uint16_t RxNum;
	float ChirpPRI;
	vector<vector<uint16_t>> all_tx_seq_pos; // (TxGroupNum, vector<int>(TxReuseNum));
	vector<float> delta_hop_freq_fstart_seq; // (TxGroupNum* TxReuseNum);
	float hop_freq_low_boundry;
	uint8_t Array_option;
	float MinElementSpaceHorRelToLamda;
	float MinElementSpaceVertRelToLamda;
	int txNum_in_group;
	vector<vector<int>> TxGroup; // (TxGroupNum, vector<int>(txNum_in_group));
	int analysis_step;
	int waveLocNum;
	vector<int> waveLocChirpSeq; // (TxGroupNum* TxReuseNum);
	int ValidCoarseRangeBin_StartIndex;
	int ValidCoarseRangeBinNum;
	int DopplerBandNum;
	vector<int> fft_win_type; // (3);
	vector<float> VelocityAxis; //(TxGroupNum* TxReuseNum)
	vector<float> CoarseRangeAxis; // (RangeSampleNum/2)
	vector<float> FineRangeAxis;
	vector<float> AngleHorAxis; //(128)
	vector<float> AngleVertAxis; // (32)
	vector<int> snr_dB_different_dim; //(4)
	int detect_option; // nouse
	float static_vel_thres;  //nouse
	float blind_zone;
	vector<int> cfar_include_order; //(4)
	vector<int> cfar_exclude_order; //(4)
	float cfar_2d_tf;
	DetPara det_para;
	float lightSpeed;
	int plot_option;
	int save_process_data;
	int work_mode;
	int frame_type;  // different with matlab, add to parasys
};

struct BinFile {
	ParaSys para_sys;
	float* input_data;
	float* compensate_mat;
};

struct Virtual_array {
	vector<vector<int>> pos;  // size need compute
	int x_max_pos;
	int x_min_pos;
	int y_max_pos;
	int y_min_pos;
	int x_space;
	int y_space;
	vector<vector<int>> pos_in_mat;
};

BinFile get_info(std::string file) {
	BinFile bin_file;
	ParaSys sys;
	bin_file.para_sys = sys;
	bin_file.input_data = new float(1);
	bin_file.compensate_mat = new float(1);
	return bin_file;
}

//Virtual_array get_virtual_array(const vector<int>& rx_x_pos_renorm, const vector<int>& rx_y_pos_renorm, const vector<int>& tx_x_pos_renorm,const vector<int>& tx_y_pos_renorm) {
//	Virtual_array virtual_array;
//	int Tx_num = size(tx_x_pos_renorm);
//	int Rx_num = size(rx_x_pos_renorm);
//	vector<int> Relative_Tx1_pos_x(Tx_num);
//	vector<int> Relative_Tx1_pos_y(Tx_num);
//	for (int i = 0; i < Tx_num; i++){
//		Relative_Tx1_pos_x[i] = tx_x_pos_renorm[i] - tx_x_pos_renorm[0];
//		Relative_Tx1_pos_y[i] = tx_y_pos_renorm[i] - tx_y_pos_renorm[0];
//	}
//
//	int repsize_rx = Tx_num * Rx_num;
//	vector<int> Virtual_Rx_pos_x(repsize_rx);
//	vector<int> Virtual_Rx_pos_y(repsize_rx);
//	for (int i = 0; i < repsize_rx; i++) {
//		Virtual_Rx_pos_x[i] = rx_x_pos_renorm[int(i % 48)] + Relative_Tx1_pos_x[int(i/48)];
//		Virtual_Rx_pos_y[i] = rx_y_pos_renorm[int(i % 48)] + Relative_Tx1_pos_y[int(i/48)];
//	}
//
//	vector<vector<int>> pos(repsize_rx, vector<int>(2));
//	for (int i = 0; i < repsize_rx; i++) {
//		pos[i][0] = Virtual_Rx_pos_x[i];
//		pos[i][1] = Virtual_Rx_pos_y[i];
//	}
//	virtual_array.pos = pos;
//
//	int max_Rx_x;
//	int min_Rx_x;
//	int max_Rx_y;
//	int min_Rx_y;
//	max_Rx_x = *max_element(Virtual_Rx_pos_x.begin(), Virtual_Rx_pos_x.end());
//	min_Rx_x = *min_element(Virtual_Rx_pos_x.begin(), Virtual_Rx_pos_x.end());
//	max_Rx_y = *max_element(Virtual_Rx_pos_y.begin(), Virtual_Rx_pos_y.end());
//	min_Rx_y = *min_element(Virtual_Rx_pos_y.begin(), Virtual_Rx_pos_y.end());
//	virtual_array.x_max_pos = max_Rx_x;
//	virtual_array.x_min_pos = min_Rx_x;
//	virtual_array.y_max_pos = max_Rx_x;
//	virtual_array.y_min_pos = min_Rx_y;
//
//	sort(Virtual_Rx_pos_x.begin(), Virtual_Rx_pos_x.end());
//	sort(Virtual_Rx_pos_y.begin(), Virtual_Rx_pos_y.end());
//	vector<int> Virtual_Rx_pos_x_reorder_degrade(repsize_rx - 1);
//	for (int i = 0; i < repsize_rx - 1; i++) {
//		Virtual_Rx_pos_x_reorder_degrade[i] = Virtual_Rx_pos_x[i + 1] - Virtual_Rx_pos_x[i];
//	}
//	int x_no_zero_num = 0;
//	vector<int> Virtual_Rx_pos_x_reorder_degrade_no_zero;
//	for (int i = 0; i < repsize_rx - 1; i++) {
//		if (Virtual_Rx_pos_x_reorder_degrade[i] != 0) {
//			x_no_zero_num++;
//			Virtual_Rx_pos_x_reorder_degrade_no_zero.push_back(Virtual_Rx_pos_x_reorder_degrade[i]);
//		}
//	}
//	vector<int> Virtual_Rx_pos_y_reorder_degrade(repsize_rx - 1);
//	for (int i = 0; i < repsize_rx - 1; i++) {
//		Virtual_Rx_pos_y_reorder_degrade[i] = Virtual_Rx_pos_y[i + 1] - Virtual_Rx_pos_y[i];
//	}
//	int y_no_zero_num = 0;
//	vector<int> Virtual_Rx_pos_y_reorder_degrade_no_zero;
//	for (int i = 0; i < repsize_rx - 1; i++) {
//		if (Virtual_Rx_pos_y_reorder_degrade[i] != 0) {
//			y_no_zero_num++;
//			Virtual_Rx_pos_y_reorder_degrade_no_zero.push_back(Virtual_Rx_pos_y_reorder_degrade[i]);
//		}
//	}
//
//	if (x_no_zero_num == 0) {
//		virtual_array.x_space = 0;
//		virtual_array.y_max_pos = INT16_MAX;
//	}
//	else {
//		virtual_array.x_space = *min_element(Virtual_Rx_pos_x_reorder_degrade_no_zero.begin(), Virtual_Rx_pos_x_reorder_degrade_no_zero.end());
//		if (virtual_array.x_space < 1) {
//			// not int type?
//			virtual_array.x_space = 1;
//		}
//	}
//
//	if (y_no_zero_num == 0) {
//		virtual_array.y_space = 0;
//		virtual_array.y_max_pos = INT16_MAX;
//	}
//	else {
//		virtual_array.y_space = *min_element(Virtual_Rx_pos_y_reorder_degrade_no_zero.begin(), Virtual_Rx_pos_y_reorder_degrade_no_zero.end());
//		if (virtual_array.y_space < 1) {
//			virtual_array.y_space = 1;
//		}
//	}
//
//	vector<vector<int>> relative_pos(repsize_rx, vector<int>(2));
//	for (int i = 0; i < repsize_rx; i++) {
//		relative_pos[i][0] = virtual_array.pos[i][0] - virtual_array.x_min_pos;
//		relative_pos[i][1] = virtual_array.pos[i][1] - virtual_array.y_min_pos;
//	}
//
//	vector<vector<int>> pos_in_mat(repsize_rx, vector<int>(2));
//	if (virtual_array.x_space != 0) {
//		for (int i = 0; i < repsize_rx; i++) {
//			pos_in_mat[i][0] = relative_pos[i][0] / virtual_array.x_space;
//		}
//	}
//	if (virtual_array.y_space != 0) {
//		for (int i = 0; i < repsize_rx; i++) {
//			pos_in_mat[i][1] = relative_pos[i][1] / virtual_array.y_space;
//		}
//	}
//	virtual_array.pos_in_mat = pos_in_mat;  // matlab add one?
//}

//struct ParaArray {
//	vector<int> rx_x_pos_renorm;
//	vector<int> rx_y_pos_renorm;
//	vector<int> tx_x_pos_renorm;
//	vector<int> tx_y_pos_renorm;
//	int min_element_space_x_relto_semilamda;
//	int min_element_space_y_relto_semilamda;
//};
//ParaArray getArray()
//{
//	// only option = 3;
//	int tx_num = 48;
//	int rx_num = 48;
//	int min_element_space_x_relto_semilamda = 1;
//	int min_element_space_y_relto_semilamda = 3;
//	vector<int> Rx_pos_x{0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52};
//	vector<int> Rx_pos_y{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36};
//	vector<int> Tx_pos_x{0,0,0,0,0,0,0,0,0,0,0,0,0,3,3,3,3,3,3,3,3,3,3,3,3,48,48,48,48,48,48,48,48,48,48,48,48,51,51,51,51,51,51,51,51,51,51,51,51};
//	vector<int> Tx_pos_y{6,9,12,15,18,21,24,27,30,33,36,39,6,9,12,15,18,21,24,27,30,33,36,39,-3,0,3,6,9,12,15,18,21,24,27,30,-3,0,3,6,9,12,15,18,21,24,27,30};
//
//	ParaArray para_array;
//	para_array.tx_x_pos_renorm = Tx_pos_x;
//	para_array.tx_y_pos_renorm = Tx_pos_y;
//	para_array.rx_x_pos_renorm = Rx_pos_x;
//	para_array.rx_y_pos_renorm = Rx_pos_y;
//	para_array.min_element_space_x_relto_semilamda = min_element_space_x_relto_semilamda;
//	para_array.min_element_space_y_relto_semilamda = min_element_space_y_relto_semilamda;
//}

BinFile parseDataFile(vector<unsigned char> &inputBuffer) {
	BinFile bin_file;
	if (inputBuffer.size() < 83886080) {
		std::cerr << "Buffer does not contain enough data to parse." << endl;
			return bin_file;
	}
	memcpy(&bin_file.para_sys.InputHasDonePreProcess, inputBuffer.data() + 260,sizeof(unsigned char));
	memcpy(&bin_file.para_sys.vel_target_generated_by_moniqi, inputBuffer.data() + 261, sizeof(unsigned char));
	memcpy(&bin_file.para_sys.RangeSampleNum, inputBuffer.data() + 262, 2*sizeof(unsigned char));
	memcpy(&bin_file.para_sys.TxReuseNum, inputBuffer.data() + 264, 2*sizeof(unsigned char));
	memcpy(&bin_file.para_sys.TxNum, inputBuffer.data() + 266, 2*sizeof(unsigned char));
	memcpy(&bin_file.para_sys.RxNum, inputBuffer.data() + 268, 2*sizeof(unsigned char));
	vector<vector<uint16_t>> all_tx_seq_pos(bin_file.para_sys.TxNum, vector<uint16_t>(bin_file.para_sys.TxReuseNum, 0));
	bin_file.para_sys.all_tx_seq_pos = all_tx_seq_pos;
	for (int i = 0; i < bin_file.para_sys.TxNum; i++) {
		memcpy(&bin_file.para_sys.all_tx_seq_pos[i][0], inputBuffer.data() + 270 + i* 2 * bin_file.para_sys.TxReuseNum, 2 * bin_file.para_sys.TxReuseNum * sizeof(unsigned char));
	}
	vector<float> delta_hop_freq_fstart_seq(bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 0);
	bin_file.para_sys.delta_hop_freq_fstart_seq = delta_hop_freq_fstart_seq;
	memcpy(&bin_file.para_sys.delta_hop_freq_fstart_seq[0], inputBuffer.data() + 270 + 2 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 4 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.ChirpPRI, inputBuffer.data() + 270 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 4 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.hop_freq_low_boundry, inputBuffer.data() + 274 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 4 * sizeof(unsigned char));
	// memcpy(&bin_file.para_sys.Virtual_array_pos_Hor, inputBuffer.data() + 278 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 4 * sizeof(unsigned char));

	memcpy(&bin_file.para_sys.Array_option, inputBuffer.data() + 278 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, sizeof(unsigned char));
	memcpy(&bin_file.para_sys.MinElementSpaceHorRelToLamda, inputBuffer.data() + 279 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 4 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.MinElementSpaceVertRelToLamda, inputBuffer.data() + 283 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 4 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.analysis_step, inputBuffer.data() + 287 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 1 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.fft_win_type, inputBuffer.data() + 288 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 4 * sizeof(unsigned char));
	// memcpy(&bin_file.para_sys.CoarseRangeNum, inputBuffer.data() + 291 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 2 * sizeof(unsigned char));
	// memcpy(&bin_file.para_sys.FineRangeAxis, inputBuffer.data() + 295 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 4 * sizeof(unsigned char));


	//std::cout << static_cast<uint8_t>(bin_file.para_sys.InputHasDonePreProcess) << std::endl;
	return bin_file;
}

int main(int argc, char* argv[]) {
	// 1. get info
	auto bin_info = get_info("test");
	auto& para_sys = bin_info.para_sys;
	auto& input_data = bin_info.input_data;
	auto& compensate_mat = bin_info.compensate_mat;

	// coarse frame data file parse code
	std::string filename = "ADCData_Frame_000000_20240314_17_22_19_2876_ShareMemory.bin";

	const size_t bufferSize = 83886080;

	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open()) {
		std::cerr << "无法打开文件" << std::endl;
		return 1;
	}

	std::vector<unsigned char> buffer(bufferSize, 0);
	file.read(reinterpret_cast<char*>(buffer.data()), bufferSize);
	size_t bytesRead = file.gcount();

	bin_info = parseDataFile(buffer);

	// 2. deal
	// Result fun(Binfile &info);

	return 0;
}