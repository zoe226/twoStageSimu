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
	uint8_t PowerDetEn;
	float SnrThre;
	uint16_t RextendNum;
	uint16_t VextendNum;
};

struct CFAR {
	uint8_t ProCellNum_R_u2;
	uint8_t ProCellNum_V_u2;
	uint8_t RefCellNum_1D_u5;
};

struct SpatialPara {
	float Vsnr_Thre_db;
	float HorProfThre_db;
	uint8_t PeakSen;
	uint16_t HorOffsetidx;
	float HorOffPeakThre_db;
	float Global_noise_Thre_db;
	float min_snr_Hor_db;
	float VerThreshold;
	uint16_t HorCfarGuardNum;
	uint16_t HorCfarReferNum;
	uint16_t Horminidx;
	uint16_t Hormaxidx;
	uint16_t Verminidx;
	uint16_t Vermaxidx;
};

struct Union {
	uint16_t CompareIdx;
	uint16_t MaxValIdx;
};

struct DetPara {
	ConReg con_reg;
	uint16_t ChToProcess_Num_u7;
	vector<char> Switch_DataCompress_u2; // (2); ‘11’bin
	uint8_t PassAbsData_u1;
	uint8_t PassChMeanData_u1;
	uint8_t PassPeakS_u1;
	uint8_t OneDEnable;
	uint16_t Reg_PGA_u15;
	uint8_t PeakS_Enable_u1;
	uint16_t RangeCellNum_u10;
	uint16_t ChirpNum_u11;
	uint16_t DetectCell_RIndex_Min_u10;
	uint16_t DetectCell_RIndex_Max_u10;
	uint16_t DetectCell_VIndex_Min_u11;
	uint16_t DetectCell_VIndex_Max_u11;
	vector<char> CFARTypeSwitch_u2; // (2);'11'
	uint8_t LogicTestFlag_u1;
	CFAR cfar_para;
	uint8_t Loc_OSCFAR_u5;
	uint16_t Index_Chirp_NotMove_OSCFAR_u11;
	uint16_t Threshold_RangeDim_For_2D_OSCFAR_u9;
	uint8_t asicsavedata_flag;
	uint8_t IndexCHIP_u3;
	SpatialPara spatial_para;
	Union union_para;
};

struct ParaSys {
	uint8_t vel_target_generated_by_moniqi;
	uint8_t need_spatial_cal;
	uint8_t InputHasDonePreProcess;
	uint16_t RangeSampleNum;
	uint16_t TxReuseNum;
	uint16_t TxNum;
	uint16_t TxGroupNum;
	uint16_t RxNum;
	float ChirpPRI;
	vector<vector<uint16_t>> all_tx_seq_pos; // (TxGroupNum, vector<int>(TxReuseNum));
	vector<float> delta_hop_freq_fstart_seq; // (TxGroupNum* TxReuseNum);
	float hop_freq_low_boundry;
	uint8_t Array_option;
	float MinElementSpaceHorRelToLamda;
	float MinElementSpaceVertRelToLamda;
	uint16_t txNum_in_group;
	vector<vector<uint16_t>> TxGroup; // (TxGroupNum, vector<int>(txNum_in_group));
	uint8_t analysis_step;
	uint8_t waveLocNum;
	vector<vector<uint16_t>> waveLocChirpSeq; // (waveLocNum,totalChirpNum/waveLocNum);
	uint16_t ValidCoarseRangeBin_StartIndex;
	uint16_t ValidCoarseRangeBinNum;
	uint16_t DopplerBandNum;
	vector<uint8_t> fft_win_type; // (3);
	vector<float> VelocityAxis; //(TxGroupNum* TxReuseNum)
	vector<float> CoarseRangeAxis; // (RangeSampleNum/2)
	vector<float> FineRangeAxis;
	vector<float> AngleHorAxis; //(128)
	vector<float> AngleVertAxis; // (32)
	vector<float> snr_dB_different_dim; //(4)
	uint8_t detect_option; // nouse
	float static_vel_thres;  //nouse
	float blind_zone;
	vector<uint8_t> cfar_include_order; //(4)
	vector<uint8_t> cfar_exclude_order; //(4)
	float cfar_2d_tf;
	DetPara det_para;
	float lightSpeed;
	uint8_t plot_option;
	uint8_t save_process_data;
	uint8_t work_mode;
	uint8_t frame_type;  // different with matlab, add to parasys
	uint16_t CoarseRangeNum;  // raw frame head parameter,size of coarseRangeAxis
	uint16_t FineRangeNum;  // raw frame head parameter,size of coarseRangeAxis
	uint16_t VelocityNum;  // raw frame head parameter,size of coarseRangeAxis
	uint16_t AngleHorNum;  // raw frame head parameter,size of coarseRangeAxis
	uint16_t AngleVertNum;  // raw frame head parameter,size of coarseRangeAxis
};

struct Virtual_array {
	vector<vector<int16_t>> pos;  // size need compute
	uint16_t x_max_pos;
	uint16_t x_min_pos;
	uint16_t y_max_pos;
	uint16_t y_min_pos;
	float x_space;
	float y_space;
	vector<vector<uint16_t>> pos_in_mat;
};

void  get_virtual_array(Virtual_array &virtual_array, const vector<int16_t>& rx_x_pos_renorm, const vector<int16_t>& rx_y_pos_renorm, const vector<int16_t>& tx_x_pos_renorm,const vector<int16_t>& tx_y_pos_renorm) {
	uint16_t Tx_num = size(tx_x_pos_renorm);
	uint16_t Rx_num = size(rx_x_pos_renorm);
	vector<int16_t> Relative_Tx1_pos_x(Tx_num);
	vector<int16_t> Relative_Tx1_pos_y(Tx_num);
	for (int i = 0; i < Tx_num; i++){
		Relative_Tx1_pos_x[i] = tx_x_pos_renorm[i] - tx_x_pos_renorm[0];
		Relative_Tx1_pos_y[i] = tx_y_pos_renorm[i] - tx_y_pos_renorm[0];
	}

	uint16_t repsize_rx = Tx_num * Rx_num;
	vector<int16_t> Virtual_Rx_pos_x(repsize_rx);
	vector<int16_t> Virtual_Rx_pos_y(repsize_rx);
	for (int i = 0; i < repsize_rx; i++) {
		Virtual_Rx_pos_x[i] = rx_x_pos_renorm[int(i % 48)] + Relative_Tx1_pos_x[int(i/48)];
		Virtual_Rx_pos_y[i] = rx_y_pos_renorm[int(i % 48)] + Relative_Tx1_pos_y[int(i/48)];
	}

	vector<vector<int16_t>> pos(repsize_rx, vector<int16_t>(2));
	for (int i = 0; i < repsize_rx; i++) {
		pos[i][0] = Virtual_Rx_pos_x[i];
		pos[i][1] = Virtual_Rx_pos_y[i];
	}
	virtual_array.pos = pos;

	uint16_t max_Rx_x;
	uint16_t min_Rx_x;
	uint16_t max_Rx_y;
	uint16_t min_Rx_y;
	max_Rx_x = *max_element(Virtual_Rx_pos_x.begin(), Virtual_Rx_pos_x.end());
	min_Rx_x = *min_element(Virtual_Rx_pos_x.begin(), Virtual_Rx_pos_x.end());
	max_Rx_y = *max_element(Virtual_Rx_pos_y.begin(), Virtual_Rx_pos_y.end());
	min_Rx_y = *min_element(Virtual_Rx_pos_y.begin(), Virtual_Rx_pos_y.end());
	virtual_array.x_max_pos = max_Rx_x;
	virtual_array.x_min_pos = min_Rx_x;
	virtual_array.y_max_pos = max_Rx_x;
	virtual_array.y_min_pos = min_Rx_y;

	sort(Virtual_Rx_pos_x.begin(), Virtual_Rx_pos_x.end());
	sort(Virtual_Rx_pos_y.begin(), Virtual_Rx_pos_y.end());
	vector<int16_t> Virtual_Rx_pos_x_reorder_degrade(repsize_rx - 1);
	for (int i = 0; i < repsize_rx - 1; i++) {
		Virtual_Rx_pos_x_reorder_degrade[i] = Virtual_Rx_pos_x[i + 1] - Virtual_Rx_pos_x[i];
	}
	uint16_t x_no_zero_num = 0;
	vector<int16_t> Virtual_Rx_pos_x_reorder_degrade_no_zero;
	for (int i = 0; i < repsize_rx - 1; i++) {
		if (Virtual_Rx_pos_x_reorder_degrade[i] != 0) {
			x_no_zero_num++;
			Virtual_Rx_pos_x_reorder_degrade_no_zero.push_back(Virtual_Rx_pos_x_reorder_degrade[i]);
		}
	}
	vector<int16_t> Virtual_Rx_pos_y_reorder_degrade(repsize_rx - 1);
	for (int i = 0; i < repsize_rx - 1; i++) {
		Virtual_Rx_pos_y_reorder_degrade[i] = Virtual_Rx_pos_y[i + 1] - Virtual_Rx_pos_y[i];
	}
	uint16_t y_no_zero_num = 0;
	vector<int16_t> Virtual_Rx_pos_y_reorder_degrade_no_zero;
	for (int i = 0; i < repsize_rx - 1; i++) {
		if (Virtual_Rx_pos_y_reorder_degrade[i] != 0) {
			y_no_zero_num++;
			Virtual_Rx_pos_y_reorder_degrade_no_zero.push_back(Virtual_Rx_pos_y_reorder_degrade[i]);
		}
	}

	if (x_no_zero_num == 0) {
		virtual_array.x_space = 0;
		virtual_array.y_max_pos = INT16_MAX;
	}
	else {
		virtual_array.x_space = *min_element(Virtual_Rx_pos_x_reorder_degrade_no_zero.begin(), Virtual_Rx_pos_x_reorder_degrade_no_zero.end());
		if (virtual_array.x_space < 1) {
			// not int type?
			virtual_array.x_space = 1;
		}
	}

	if (y_no_zero_num == 0) {
		virtual_array.y_space = 0;
		virtual_array.y_max_pos = INT16_MAX;
	}
	else {
		virtual_array.y_space = *min_element(Virtual_Rx_pos_y_reorder_degrade_no_zero.begin(), Virtual_Rx_pos_y_reorder_degrade_no_zero.end());
		if (virtual_array.y_space < 1) {
			virtual_array.y_space = 1;
		}
	}

	vector<vector<int16_t>> relative_pos(repsize_rx, vector<int16_t>(2));
	for (int i = 0; i < repsize_rx; i++) {
		relative_pos[i][0] = virtual_array.pos[i][0] - virtual_array.x_min_pos;
		relative_pos[i][1] = virtual_array.pos[i][1] - virtual_array.y_min_pos;
	}

	vector<vector<uint16_t>> pos_in_mat(repsize_rx, vector<uint16_t>(2));
	if (virtual_array.x_space != 0) {
		for (int i = 0; i < repsize_rx; i++) {
			pos_in_mat[i][0] = relative_pos[i][0] / virtual_array.x_space;
		}
	}
	if (virtual_array.y_space != 0) {
		for (int i = 0; i < repsize_rx; i++) {
			pos_in_mat[i][1] = relative_pos[i][1] / virtual_array.y_space;
		}
	}
	virtual_array.pos_in_mat = pos_in_mat;  // matlab add one?
}

struct ParaArray {
	vector<int16_t> rx_x_pos_renorm;
	vector<int16_t> rx_y_pos_renorm;
	vector<int16_t> tx_x_pos_renorm;
	vector<int16_t> tx_y_pos_renorm;
	float min_element_space_x_relto_semilamda;
	float min_element_space_y_relto_semilamda;
};
 void getArray(ParaArray &para_array)
{
	// only option = 3;
	uint16_t tx_num = 48;
	uint16_t rx_num = 48;
	float min_element_space_x_relto_semilamda = 1;
	float min_element_space_y_relto_semilamda = 3;
	vector<int16_t> Rx_pos_x{0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52};
	vector<int16_t> Rx_pos_y{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36};
	vector<int16_t> Tx_pos_x{0,0,0,0,0,0,0,0,0,0,0,0,0,3,3,3,3,3,3,3,3,3,3,3,3,48,48,48,48,48,48,48,48,48,48,48,48,51,51,51,51,51,51,51,51,51,51,51,51};
	vector<int16_t> Tx_pos_y{6,9,12,15,18,21,24,27,30,33,36,39,6,9,12,15,18,21,24,27,30,33,36,39,-3,0,3,6,9,12,15,18,21,24,27,30,-3,0,3,6,9,12,15,18,21,24,27,30};

	para_array.tx_x_pos_renorm = Tx_pos_x;
	para_array.tx_y_pos_renorm = Tx_pos_y;
	para_array.rx_x_pos_renorm = Rx_pos_x;
	para_array.rx_y_pos_renorm = Rx_pos_y;
	para_array.min_element_space_x_relto_semilamda = min_element_space_x_relto_semilamda;
	para_array.min_element_space_y_relto_semilamda = min_element_space_y_relto_semilamda;
}

 struct BinFile {
	 ParaSys para_sys;
	 std::unique_ptr<int16_t[]> input_data;
	 std::unique_ptr<float[]> compensate_mat;
 };

void parseDataFile(vector<unsigned char> &inputBuffer, BinFile &bin_file) {
	if (inputBuffer.size() < 83886080) {
		std::cerr << "Buffer does not contain enough data to parse." << endl;
	}
	memcpy(&bin_file.para_sys.InputHasDonePreProcess, inputBuffer.data() + 260,sizeof(unsigned char));
	memcpy(&bin_file.para_sys.vel_target_generated_by_moniqi, inputBuffer.data() + 261, sizeof(unsigned char));
	memcpy(&bin_file.para_sys.RangeSampleNum, inputBuffer.data() + 262, 2*sizeof(unsigned char));
	memcpy(&bin_file.para_sys.TxReuseNum, inputBuffer.data() + 264, 2*sizeof(unsigned char));
	memcpy(&bin_file.para_sys.TxNum, inputBuffer.data() + 266, 2*sizeof(unsigned char));
	memcpy(&bin_file.para_sys.RxNum, inputBuffer.data() + 268, 2*sizeof(unsigned char));
	bin_file.para_sys.all_tx_seq_pos.resize(bin_file.para_sys.TxNum, vector<uint16_t>(bin_file.para_sys.TxReuseNum));
	for (int i = 0; i < bin_file.para_sys.TxNum; i++) {
		bin_file.para_sys.all_tx_seq_pos[i].assign(reinterpret_cast<uint16_t*>(inputBuffer.data() + 270 + i * 2 * bin_file.para_sys.TxReuseNum), reinterpret_cast<uint16_t*>(inputBuffer.data() + 270 + i * 2 * bin_file.para_sys.TxReuseNum+ 2 * bin_file.para_sys.TxReuseNum * sizeof(unsigned char)));
	}
	bin_file.para_sys.delta_hop_freq_fstart_seq.resize(bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum);
	bin_file.para_sys.delta_hop_freq_fstart_seq.assign(reinterpret_cast < float*>(inputBuffer.data()+270+ bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum*2), reinterpret_cast <float*>(inputBuffer.data()+270 + bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum * 2 + 4* bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum*sizeof(unsigned char)));
	memcpy(&bin_file.para_sys.ChirpPRI, inputBuffer.data() + 270 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 4 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.hop_freq_low_boundry, inputBuffer.data() + 274 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 4 * sizeof(unsigned char));
	// memcpy(&bin_file.para_sys.Virtual_array_pos_Hor, inputBuffer.data() + 278 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 4 * sizeof(unsigned char));
	// Virtual_array_pos_Vert
	// Virtual_array_pos_Hor_in_mat
	// Virtual_array_pos_Vert_in_mat

	memcpy(&bin_file.para_sys.Array_option, inputBuffer.data() + 278 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, sizeof(unsigned char));
	memcpy(&bin_file.para_sys.MinElementSpaceHorRelToLamda, inputBuffer.data() + 279 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 4 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.MinElementSpaceVertRelToLamda, inputBuffer.data() + 283 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 4 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.analysis_step, inputBuffer.data() + 287 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 1 * sizeof(unsigned char));
	bin_file.para_sys.fft_win_type.resize(3);
	bin_file.para_sys.fft_win_type.assign(inputBuffer.data() + 288 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, inputBuffer.data() + 288 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 3 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.CoarseRangeNum, inputBuffer.data() + 291 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum, 2 * sizeof(unsigned char));
	bin_file.para_sys.CoarseRangeAxis.resize(bin_file.para_sys.CoarseRangeNum);
	bin_file.para_sys.CoarseRangeAxis.assign(reinterpret_cast<float*>(inputBuffer.data() + 293 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum), reinterpret_cast<float*>(inputBuffer.data() + 293 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * sizeof(unsigned char) * bin_file.para_sys.CoarseRangeNum));
	memcpy(&bin_file.para_sys.FineRangeNum, inputBuffer.data() + 293 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum, 2 * sizeof(unsigned char));
	bin_file.para_sys.FineRangeAxis.resize(bin_file.para_sys.FineRangeNum);
	bin_file.para_sys.FineRangeAxis.assign(reinterpret_cast<float*>(inputBuffer.data() + 295 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum),reinterpret_cast<float*>(inputBuffer.data() + 295 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum*sizeof(unsigned char)));
	memcpy(&bin_file.para_sys.VelocityNum, inputBuffer.data() + 295 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum, 2 * sizeof(unsigned char));
	bin_file.para_sys.VelocityAxis.resize(bin_file.para_sys.VelocityNum);
	bin_file.para_sys.VelocityAxis.assign(reinterpret_cast<float*>(inputBuffer.data() + 297 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum), reinterpret_cast<float*>(inputBuffer.data() + 297 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * sizeof(unsigned char) * bin_file.para_sys.VelocityNum));
	memcpy(&bin_file.para_sys.AngleHorNum, inputBuffer.data() + 297 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum, 2 * sizeof(unsigned char));
	bin_file.para_sys.AngleHorAxis.resize(bin_file.para_sys.AngleHorNum);
	bin_file.para_sys.AngleHorAxis.assign(reinterpret_cast<float*>(inputBuffer.data() + 299 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum),reinterpret_cast<float*>(inputBuffer.data() + 299 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * sizeof(unsigned char) * bin_file.para_sys.AngleHorNum));
	memcpy(&bin_file.para_sys.AngleVertNum, inputBuffer.data() + 299 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum, 2 * sizeof(unsigned char));
	bin_file.para_sys.AngleVertAxis.resize(bin_file.para_sys.AngleVertNum);
	bin_file.para_sys.AngleVertAxis.assign(reinterpret_cast<float*>(inputBuffer.data() + 301 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum),reinterpret_cast<float*>(inputBuffer.data() + 301 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * sizeof(unsigned char)*bin_file.para_sys.AngleVertNum));
	memcpy(&bin_file.para_sys.ValidCoarseRangeBin_StartIndex, inputBuffer.data() + 301 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, 2 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.ValidCoarseRangeBinNum, inputBuffer.data() + 303 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, 2 * sizeof(unsigned char));
	bin_file.para_sys.snr_dB_different_dim.resize(4);
	bin_file.para_sys.snr_dB_different_dim.assign(reinterpret_cast<float*>(inputBuffer.data() + 305 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum), reinterpret_cast<float*>(inputBuffer.data() + 305 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum + 4 * sizeof(unsigned char) * 4));
	memcpy(&bin_file.para_sys.static_vel_thres, inputBuffer.data() + 322 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, 4 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.blind_zone, inputBuffer.data() + 326 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, 4 * sizeof(unsigned char));
	bin_file.para_sys.cfar_include_order.resize(4);
	bin_file.para_sys.cfar_include_order.assign(inputBuffer.data() + 330 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, inputBuffer.data() + 330 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum + 4 * sizeof(unsigned char));
	bin_file.para_sys.cfar_exclude_order.resize(4);
	bin_file.para_sys.cfar_exclude_order.assign(inputBuffer.data() + 334 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, inputBuffer.data() + 334 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum + 4 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.save_process_data, inputBuffer.data() + 338 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, 1 * sizeof(unsigned char));
	
	bin_file.compensate_mat = std::make_unique<float[]>(bin_file.para_sys.TxNum * bin_file.para_sys.RxNum * 2);
	bin_file.input_data = std::make_unique<int16_t[]>(bin_file.para_sys.RxNum*bin_file.para_sys.TxNum*bin_file.para_sys.TxReuseNum*bin_file.para_sys.RangeSampleNum * 2);
	memcpy(bin_file.compensate_mat.get(), inputBuffer.data() + 305 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, 4 * sizeof(unsigned char) * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum * 2);
	memcpy(bin_file.input_data.get(), inputBuffer.data() + 3328 + 340 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, 2 * sizeof(unsigned char) * bin_file.para_sys.RxNum * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum * bin_file.para_sys.RangeSampleNum * 2);
}
void getDetPara(DetPara &det_para) {
	;

}
void get_info(std::string filename, BinFile &bin_file) {
	// 输入文件名字，和bin_file
	// 输出解析文件赋值bin_file，帧头中没有的手动写入(注意ADC文件中有些参数两阶段平台是没有采用的，跳过处理)
	// 先打开文件确认能打开，然后解析文件的帧头来解析参数，最后将帧头中没有的参数部分手动配置(注意如TxGroupNum的参数要先于一些ADC中的文件来配置，类似这样的顺序要注意)
	// 帧头中包含virtual array的信息，确认是否将virtual array也作为输入的参数，对其部分值也进行配置
	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open()) {
		std::cerr << "无法打开文件!" << std::endl;
		return;
	}
	const size_t bufferSize = 83886080;
	std::vector<unsigned char> buffer(bufferSize, 0);
	file.read(reinterpret_cast<char*>(buffer.data()), bufferSize);
	size_t bytesRead = file.gcount();
	
	bin_file.para_sys.TxGroupNum = 1;
	bin_file.para_sys.need_spatial_cal = 1;
	bin_file.para_sys.txNum_in_group = 1;
	bin_file.para_sys.TxGroup.push_back({1});  // TxGroupNum*txNum_in_group
	bin_file.para_sys.waveLocNum = 1;
	// bin_file.para_sys.waveLocChirpSeq.resize(bin_file.para_sys.waveLocNum,vector<uint16_t>(cols));
	bin_file.para_sys.waveLocChirpSeq.push_back({0});  // waveLocNum * (totalChirpNum/waveLocNum)
	bin_file.para_sys.DopplerBandNum = 1;
	bin_file.para_sys.detect_option = 0;
	bin_file.para_sys.cfar_2d_tf = 1.8;
	//bin_file.para_sys.det_para;
	getDetPara(bin_file.para_sys.det_para);
	bin_file.para_sys.lightSpeed = 299792458.0f;
	bin_file.para_sys.plot_option = 0;
	bin_file.para_sys.work_mode = 0;
	bin_file.para_sys.frame_type = 0;

	parseDataFile(buffer, bin_file);


	/*ParaSys sys;
	bin_file.para_sys = sys;
	bin_file.input_data = std::make_unique<int16_t[]>(10);
	bin_file.compensate_mat = std::make_unique<float[]>(10);*/
}

int main(int argc, char* argv[]) {
	// 1. get info
	// coarse frame data file parse code
	std::string filename = "ADCData_Frame_000000_20240314_17_22_19_2876_ShareMemory.bin";

	BinFile bin_file;
	get_info(filename, bin_file);
	auto& para_sys = bin_file.para_sys;
	auto& input_data = bin_file.input_data;
	auto& compensate_mat = bin_file.compensate_mat;


	const size_t bufferSize = 83886080;

	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open()) {
		std::cerr << "无法打开文件" << std::endl;
		return 1;
	}

	std::vector<unsigned char> buffer(bufferSize, 0);
	file.read(reinterpret_cast<char*>(buffer.data()), bufferSize);
	size_t bytesRead = file.gcount();

	BinFile temp_bin_info;
	parseDataFile(buffer, temp_bin_info);

	// 2. deal
	// Result fun(Binfile &info);

	return 0;
}