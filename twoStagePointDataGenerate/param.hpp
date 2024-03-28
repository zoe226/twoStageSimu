#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <complex>

using namespace std;
extern const uint16_t CoarseFrame_CFARdim;
extern uint16_t FineFrame_CfarDim;
extern vector<vector<float>> TOI;

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
	vector<char> Switch_DataCompress_u2; // (2); ¡®11¡¯bin
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
	uint16_t RangeCell_left;
	uint16_t RangeCell_right;
	uint16_t VelocityCell_left;
	uint16_t VelocityCell_right;
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

struct ParaArray {
	vector<int16_t> rx_x_pos_renorm;
	vector<int16_t> rx_y_pos_renorm;
	vector<int16_t> tx_x_pos_renorm;
	vector<int16_t> tx_y_pos_renorm;
	float min_element_space_x_relto_semilamda;
	float min_element_space_y_relto_semilamda;
};

struct BinFile {
	ParaSys para_sys;
	std::unique_ptr<int16_t[]> input_data;
	std::unique_ptr<float[]> compensate_mat;
};

void get_virtual_array(Virtual_array& virtual_array, const vector<int16_t>& rx_x_pos_renorm, const vector<int16_t>& rx_y_pos_renorm, const vector<int16_t>& tx_x_pos_renorm, const vector<int16_t>& tx_y_pos_renorm);
void getArray(ParaArray& para_array);
void parseDataFile(vector<unsigned char>& inputBuffer, BinFile& bin_file);
void getDetPara(DetPara& det_para, ParaSys& para_sys);
void get_info(std::string filename, BinFile& bin_file, Virtual_array& virtual_array);
