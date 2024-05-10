#define _USE_MATH_DEFINES
#include "param.hpp"
#include "coarseFrameSignalProcess.hpp"
#include <math.h>
#include <cmath>
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <algorithm>

extern vector<float> R_Threshold1_u9;
extern vector<float> V_Threshold1_u9;
extern vector<float> PeakThreshold1_u30;

uint8_t func_DetermineEachElement(uint8_t LocComp_u5, uint8_t LeftLoc_u1, float Cell_ToCompare_u30, vector<float>& DataSet_RefCell_u30, uint16_t RefCellNum_u6) {
	uint8_t IsData_u1 = 0;
	uint8_t flag_Not_u1 = 0;
	uint8_t Data_Equal_u5 = 0;
	uint8_t Count_u5 = 0;
	for (size_t Loop2 = 0; Loop2 < RefCellNum_u6; Loop2++)
	{
		float Cell_Temp_u30 = DataSet_RefCell_u30[Loop2]; 
		if (LeftLoc_u1 == 1)
		{
			if (Cell_ToCompare_u30 > Cell_Temp_u30)
			{
				Count_u5 = Count_u5 + 1;
			}
		}
		else
		{
			if (Cell_ToCompare_u30 < Cell_Temp_u30)
			{
				Count_u5 = Count_u5 + 1;
			}
		}
		if (Cell_ToCompare_u30 == Cell_Temp_u30)
		{
			Data_Equal_u5 = Data_Equal_u5 + 1;
		}
		if (Count_u5 > LocComp_u5)
		{
			flag_Not_u1 = 1;
			break;
		}
	}
	if (flag_Not_u1 == 0)
	{
		if (Data_Equal_u5 == 1)
		{
			if (Count_u5 == LocComp_u5)
			{
				IsData_u1 = 1;
			}
		}
		else
		{
			uint8_t NumLess_u5 = Count_u5 + Data_Equal_u5;
			if (NumLess_u5 > LocComp_u5)
			{
				IsData_u1 = 1;
			}
		}
	}
	return IsData_u1;
}

float func_ObtainOSLocData(uint8_t Loc_OSCFAR_u5, vector<float>& DataSet_RefCell_u30, uint16_t RefCellNum_u6)
{
	float Data_OSLoc_u30;
	Data_OSLoc_u30 = DataSet_RefCell_u30[RefCellNum_u6 - 1];
	uint8_t LeftLoc_u1 = 1;
	uint8_t LocComp_u5 = Loc_OSCFAR_u5 - 1;
	uint8_t middleLoc = floor(RefCellNum_u6 / 2);
	if (Loc_OSCFAR_u5 > middleLoc)
	{
		LeftLoc_u1 = 0;
		LocComp_u5 = RefCellNum_u6 - Loc_OSCFAR_u5;
	}
	for (size_t Loop1 = 0; Loop1 < RefCellNum_u6-1; Loop1++)
	{
		float Cell_ToCompare_u30 = DataSet_RefCell_u30[Loop1];
		uint8_t IsData_u1 = 0;
		IsData_u1 = func_DetermineEachElement(LocComp_u5,LeftLoc_u1,Cell_ToCompare_u30,DataSet_RefCell_u30,RefCellNum_u6);
		if (IsData_u1 == 1)
		{
			Data_OSLoc_u30 = Cell_ToCompare_u30;
			break;
		}
	}
	return Data_OSLoc_u30;
}

void func_CFARChM_OS_1D_V(uint8_t& IsTarget_1D_V_u1, float& VSNR_u11, DetPara& det_para, float DataToDetect, uint16_t Index_V, uint16_t Index_R, MultiDimensionalVector<float, 2>& SpatialFFTVelSel_VeloNum_RangeNum)
{
	uint8_t Loc_OSCFAR_u5 = det_para.Loc_OSCFAR_u5;
	uint8_t V_Threshold_Num_u6 = 16;
	vector<float> DataSet_RefCell_u30(32, 0.0);
	IsTarget_1D_V_u1 = 0;
	VSNR_u11 = 0;
	float Data_OSLoc_u30;
	uint8_t Logic1_u1 = 0;
	uint8_t Logic2_u1 = 0;
	uint8_t Logic3_u1 = 0;
	uint8_t Logic4_u1 = 0;
	uint8_t Logic5_u1 = 0;
	uint16_t UpBoundary_u9 = 1 + det_para.cfar_para.ProCellNum_V_u2 + det_para.cfar_para.RefCellNum_1D_u5;
	uint16_t DownBoundary_u10 = det_para.ChirpNum_u11 - det_para.cfar_para.ProCellNum_R_u2 - det_para.cfar_para.RefCellNum_1D_u5;
	uint16_t RefCellNum_u6 = 2 * det_para.cfar_para.RefCellNum_1D_u5;
	if (UpBoundary_u9 >= DownBoundary_u10)
	{
		det_para.cfar_para.ProCellNum_V_u2 = 1;
		det_para.cfar_para.RefCellNum_1D_u5 = 1;
		uint16_t UpBoundary_u9 = 1 + det_para.cfar_para.ProCellNum_V_u2 + det_para.cfar_para.RefCellNum_1D_u5;
		uint16_t DownBoundary_u10 = det_para.ChirpNum_u11 - det_para.cfar_para.ProCellNum_R_u2 - det_para.cfar_para.RefCellNum_1D_u5;
	}
	if (Index_V >= 0  && Index_V <= 0 + det_para.cfar_para.ProCellNum_V_u2)
	{
		Logic1_u1 = 1;
	}
	if (Index_V > (0 + det_para.cfar_para.ProCellNum_V_u2) && Index_V < 0 + det_para.cfar_para.ProCellNum_V_u2 + det_para.cfar_para.RefCellNum_1D_u5)
	{
		Logic2_u1 = 1;
	}
	if (Index_V >= 0 + det_para.cfar_para.ProCellNum_V_u2 + det_para.cfar_para.RefCellNum_1D_u5 && Index_V <= det_para.ChirpNum_u11 - det_para.cfar_para.ProCellNum_V_u2 - det_para.cfar_para.RefCellNum_1D_u5 - 1)
	{
		Logic3_u1 = 1;
	}
	if (Index_V > det_para.ChirpNum_u11 - det_para.cfar_para.ProCellNum_V_u2 - det_para.cfar_para.RefCellNum_1D_u5 - 1 && Index_V < det_para.ChirpNum_u11 - det_para.cfar_para.ProCellNum_V_u2 - 1)
	{
		Logic4_u1 = 1;
	}
	if (Index_V >= det_para.ChirpNum_u11 - det_para.cfar_para.ProCellNum_V_u2 - 1)
	{
		Logic5_u1 = 1;
	}
	if (Logic1_u1 == 1)
	{
		for (size_t i = 0; i < RefCellNum_u6; i++)
		{
			if (i % 2 == 1)
			{
				uint16_t upIndex = det_para.ChirpNum_u11 - det_para.cfar_para.ProCellNum_V_u2 - det_para.cfar_para.RefCellNum_1D_u5 + Index_V + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({upIndex,Index_R});
			}
			else
			{
				uint16_t downIndex = Index_V + det_para.cfar_para.ProCellNum_V_u2 + 1 + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ downIndex,Index_R });
			}
		}
	}
	else if(Logic2_u1 == 1)
	{
		for (size_t i = 0; i < RefCellNum_u6; i++)
		{
			if (i % 2 == 1)
			{
				uint16_t upIndex = (det_para.ChirpNum_u11 - det_para.cfar_para.ProCellNum_V_u2 - det_para.cfar_para.RefCellNum_1D_u5 + Index_V + i / 2)%det_para.ChirpNum_u11;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ upIndex,Index_R });
			}
			else
			{
				uint16_t downIndex = Index_V + det_para.cfar_para.ProCellNum_V_u2 + 1 + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ downIndex,Index_R });
			}
		}
	}
	else if (Logic3_u1 == 1)
	{
		for (size_t i = 0; i < RefCellNum_u6; i++)
		{
			if (i % 2 == 1)
			{
				uint16_t upIndex = Index_V - det_para.cfar_para.ProCellNum_V_u2 - det_para.cfar_para.RefCellNum_1D_u5 + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ upIndex,Index_R });
			}
			else
			{
				uint16_t downIndex = Index_V + det_para.cfar_para.ProCellNum_V_u2 + 1 + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ downIndex,Index_R });
			}
		}
	}
	else if (Logic4_u1 == 1)
	{
		for (size_t i = 0; i < RefCellNum_u6; i++)
		{
			if (i % 2 == 1)
			{
				uint16_t upIndex = Index_V - det_para.cfar_para.ProCellNum_V_u2 - det_para.cfar_para.RefCellNum_1D_u5 + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ upIndex,Index_R });
			}
			else
			{
				uint16_t downIndex = (Index_V + det_para.cfar_para.ProCellNum_V_u2 + 1 + i / 2)%det_para.ChirpNum_u11;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ downIndex,Index_R });
			}
		}
	}
	else if (Logic5_u1 == 1)
	{
		for (size_t i = 0; i < RefCellNum_u6; i++)
		{
			if (i % 2 == 1)
			{
				uint16_t upIndex = Index_V - det_para.cfar_para.ProCellNum_V_u2 - det_para.cfar_para.RefCellNum_1D_u5 + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ upIndex,Index_R });
			}
			else
			{
				uint16_t downIndex = Index_V + det_para.cfar_para.ProCellNum_V_u2 + 1 - det_para.ChirpNum_u11 + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ downIndex,Index_R });
			}
		}
	}
	if (Loc_OSCFAR_u5 >= RefCellNum_u6)
	{
		Loc_OSCFAR_u5 = floor(RefCellNum_u6 / 2) + 1;
	}
	Data_OSLoc_u30 = func_ObtainOSLocData(Loc_OSCFAR_u5, DataSet_RefCell_u30, RefCellNum_u6);
	uint16_t thIdx = floor(Index_V / (det_para.ChirpNum_u11 / V_Threshold_Num_u6));
	float Threshold2_u39 = V_Threshold1_u9[thIdx] * Data_OSLoc_u30;
	float CellToDetect_u32 = DataToDetect * 2 * 2;
	if (CellToDetect_u32 > Threshold2_u39)
	{
		IsTarget_1D_V_u1 = 1;
		VSNR_u11 = DataToDetect / Data_OSLoc_u30;
	}
}

void func_CFARChM_OS_1D(uint8_t& IsTarget_1D_R_u1, float& RSNR_u11,DetPara& det_para, float DataToDetect, uint16_t Index_V,uint16_t Index_R, MultiDimensionalVector<float, 2>& SpatialFFTVelSel_VeloNum_RangeNum)
{
	uint8_t Loc_OSCFAR_u5 = det_para.Loc_OSCFAR_u5;
	uint8_t R_Threshold_Num_u6 = 32;
	vector<float> DataSet_RefCell_u30(32, 0.0);
	uint8_t Logic1_u1 = 0;
	uint8_t Logic2_u1 = 0;
	IsTarget_1D_R_u1 = 0;
	RSNR_u11 = 0;
	uint16_t LeftBoundary_u9 = det_para.cfar_para.ProCellNum_R_u2 + det_para.cfar_para.RefCellNum_1D_u5 - 1;
	uint16_t RightBoundary_u10 = det_para.RangeCellNum_u10 - det_para.cfar_para.ProCellNum_R_u2 - det_para.cfar_para.RefCellNum_1D_u5;
	uint8_t RefCellNum_u6 = 0;
	float Data_OSLoc_u30;

	if (LeftBoundary_u9 >= RightBoundary_u10)
	{
		det_para.cfar_para.ProCellNum_R_u2 = 1;
		det_para.cfar_para.RefCellNum_1D_u5 = 1;
		LeftBoundary_u9 = det_para.cfar_para.ProCellNum_R_u2 + det_para.cfar_para.RefCellNum_1D_u5 - 1;
		RightBoundary_u10 = det_para.RangeCellNum_u10 - det_para.cfar_para.ProCellNum_R_u2 - det_para.cfar_para.RefCellNum_1D_u5;
	}

	if (Index_R > LeftBoundary_u9)
	{
		Logic1_u1 = 1;
	}
	if (Index_R < RightBoundary_u10)
	{
		Logic2_u1 = 1;
	}
	if (Logic1_u1 == 1 && Logic2_u1 == 1)
	{
		RefCellNum_u6 = 2 * det_para.cfar_para.RefCellNum_1D_u5;
		for (size_t i = 0; i < RefCellNum_u6; i++)
		{	
			if ((i%2) == 1)
			{
				uint16_t leftRefIndex = Index_R - det_para.cfar_para.ProCellNum_R_u2 - det_para.cfar_para.RefCellNum_1D_u5 + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ Index_V, leftRefIndex });
			}
			else
			{
				uint16_t rightRefIndex = Index_R + det_para.cfar_para.ProCellNum_R_u2 + 1 + i / 2;
				DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ Index_V, rightRefIndex });
			}
		}
	}
	else if (Logic2_u1 == 1)
	{
		RefCellNum_u6 = det_para.cfar_para.RefCellNum_1D_u5;
		for (size_t i = 0; i < det_para.cfar_para.RefCellNum_1D_u5; i++)
		{
			uint16_t rightRefIndex = Index_R + det_para.cfar_para.ProCellNum_R_u2 + 1 + i;
			DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,rightRefIndex});
		}
		for (size_t i = det_para.cfar_para.RefCellNum_1D_u5; i < RefCellNum_u6; i++)
		{
			uint16_t leftRefIndex = i - det_para.cfar_para.RefCellNum_1D_u5 + 1;
			DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ Index_V,leftRefIndex });
		}
	}
	else
	{
		RefCellNum_u6 = det_para.cfar_para.RefCellNum_1D_u5;
		for (size_t i = 0; i < det_para.cfar_para.RefCellNum_1D_u5; i++)
		{
			uint16_t leftRefIndex = Index_R - det_para.cfar_para.ProCellNum_R_u2 - det_para.cfar_para.RefCellNum_1D_u5 + i;
			DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,leftRefIndex});
		}
		for (size_t i = det_para.cfar_para.RefCellNum_1D_u5; i < RefCellNum_u6; i++)
		{
			uint16_t rightRefIndex = Index_R + det_para.cfar_para.ProCellNum_R_u2 + 1 + i;
			DataSet_RefCell_u30[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({ Index_V,rightRefIndex });
		}
	}
	if (Loc_OSCFAR_u5 >= RefCellNum_u6 - 1)
	{
		Loc_OSCFAR_u5 = floor(RefCellNum_u6 / 2) + 1;
	}
	Data_OSLoc_u30 = func_ObtainOSLocData(Loc_OSCFAR_u5,DataSet_RefCell_u30,RefCellNum_u6);
	uint16_t thIdx = floor(Index_R / (det_para.RangeCellNum_u10 / R_Threshold_Num_u6));
	float Threshold2_u39 = R_Threshold1_u9[thIdx] * Data_OSLoc_u30;
	float CellToDetect_u32 = DataToDetect * 2 * 2;
	if (CellToDetect_u32 > Threshold2_u39)
	{
		IsTarget_1D_R_u1 = 1;
		RSNR_u11 = DataToDetect / Data_OSLoc_u30;
	}
}

void cfarOS_cal_single_point_2D_Cross(uint8_t& IsTarget_2D, vector<float>& SNR, DetPara& det_para, float DataToDetect, uint16_t Index_R, uint16_t Index_V, MultiDimensionalVector<float, 2>& SpatialFFTVelSel_VeloNum_RangeNum)
{
	uint8_t IsTarget_1D_R_u1 = 0;
	float RSNR_u11 = 0.0;
	func_CFARChM_OS_1D(IsTarget_1D_R_u1, RSNR_u11,det_para,DataToDetect,Index_V,Index_R,SpatialFFTVelSel_VeloNum_RangeNum);
	if (IsTarget_1D_R_u1 == 1 || det_para.LogicTestFlag_u1 == 1)
	{
		uint8_t IsTarget_1D_V_u1 = 0;
		float VSNR_u11 = 0.0;
		func_CFARChM_OS_1D_V(IsTarget_1D_V_u1,VSNR_u11,det_para, DataToDetect, Index_V, Index_R, SpatialFFTVelSel_VeloNum_RangeNum);
		if (IsTarget_1D_R_u1 == 1 && IsTarget_1D_V_u1 == 1)
		{
			IsTarget_2D = 1;
			SNR[0] = RSNR_u11;
			SNR[1] = VSNR_u11;

			float GateLimit_Stationary_u5 = std::max(3.0, floor(det_para.ChirpNum_u11 / pow(2,6)));
			float CellToDetect_u32 = DataToDetect * pow(2, 2);
			if (Index_R > 4 && Index_R < det_para.RangeCellNum_u10 - 3)
			{
				uint16_t ChirpDiff_s11 = Index_V - det_para.Index_Chirp_NotMove_OSCFAR_u11;
				if (det_para.Index_Chirp_NotMove_OSCFAR_u11 > GateLimit_Stationary_u5 && det_para.Index_Chirp_NotMove_OSCFAR_u11 <= det_para.ChirpNum_u11 - GateLimit_Stationary_u5)
				{
					if (ChirpDiff_s11 >= -GateLimit_Stationary_u5 && ChirpDiff_s11 <= GateLimit_Stationary_u5)
					{
						uint16_t IndexR_Left_u10 = Index_R - 2;
						uint16_t IndexR_Right_u10 = Index_R + 2;
						float SumChLeft_RIndex_u30 = SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,IndexR_Left_u10});
						float SumChRight_RIndex_u30 = SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,IndexR_Right_u10});
						float Threshold3_u30 = 0;
						if (SumChLeft_RIndex_u30 < SumChRight_RIndex_u30)
						{
							Threshold3_u30 = SumChLeft_RIndex_u30;
						}
						else
						{
							Threshold3_u30 = SumChRight_RIndex_u30;
						}
						float Threshold4_u39 = Threshold3_u30 * det_para.Threshold_RangeDim_For_2D_OSCFAR_u9;
						if (DataToDetect <= Threshold4_u39)
						{
							IsTarget_2D = 0;
							SNR[0] = 0;
							SNR[1] = 0;
						}
					}
					
				}
				else if (det_para.Index_Chirp_NotMove_OSCFAR_u11 <= GateLimit_Stationary_u5)
				{
					uint16_t LeftBoundary_u10 = det_para.ChirpNum_u11 - GateLimit_Stationary_u5 + det_para.Index_Chirp_NotMove_OSCFAR_u11;
					uint16_t RightBoundary_u11 = GateLimit_Stationary_u5 + det_para.Index_Chirp_NotMove_OSCFAR_u11;
					if (Index_V > LeftBoundary_u10 || Index_V <= RightBoundary_u11)
					{
						uint16_t IndexR_Left_u10 = Index_R - 2;
						uint16_t IndexR_Right_u10 = Index_R + 2;
						float SumChLeft_RIndex_u30 = SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,IndexR_Left_u10});
						float SumChRight_RIndex_u30 = SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,IndexR_Right_u10});
						float Threshold3_u30 = 0;
						if (SumChLeft_RIndex_u30 < SumChRight_RIndex_u30)
						{
							Threshold3_u30 = SumChLeft_RIndex_u30;
						}
						else
						{
							Threshold3_u30 = SumChRight_RIndex_u30;
						}
						float Threshold4_u39 = Threshold3_u30 * det_para.Threshold_RangeDim_For_2D_OSCFAR_u9;
						if (DataToDetect <= Threshold4_u39)
						{
							IsTarget_2D = 0; 
							SNR[0] = 0;
							SNR[1] = 0;
						}
					}
				}
				else if (det_para.Index_Chirp_NotMove_OSCFAR_u11 > det_para.ChirpNum_u11 - GateLimit_Stationary_u5)
				{
					uint16_t LeftBoundary_u10 = det_para.Index_Chirp_NotMove_OSCFAR_u11 - GateLimit_Stationary_u5;
					uint16_t RightBoundary_u11 = det_para.ChirpNum_u11 - det_para.Index_Chirp_NotMove_OSCFAR_u11;
					if (Index_V >= LeftBoundary_u10 || Index_V <= RightBoundary_u11)
					{
						uint16_t IndexR_Left_u10 = Index_R - 2;
						uint16_t IndexR_Right_u10 = Index_R + 2;
						float SumChLeft_RIndex_u30 = SpatialFFTVelSel_VeloNum_RangeNum.get({ Index_V,IndexR_Left_u10 });
						float SumChRight_RIndex_u30 = SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,IndexR_Right_u10});
						float Threshold3_u30 = 0;
						if (SumChLeft_RIndex_u30 < SumChRight_RIndex_u30)
						{
							Threshold3_u30 = SumChLeft_RIndex_u30;
						}
						else
						{
							Threshold3_u30 = SumChRight_RIndex_u30;
						}
						float Threshold4_u39 = Threshold3_u30 * det_para.Threshold_RangeDim_For_2D_OSCFAR_u9;
						if (DataToDetect <= Threshold4_u39)
						{
							IsTarget_2D = 0;
							SNR[0] = 0;
							SNR[1] = 0;
						}
					}
				}
			}
			else if(Index_R > det_para.RangeCellNum_u10 - 3)
			{
				 uint16_t ChirpDiff_s11 = Index_V - det_para.Index_Chirp_NotMove_OSCFAR_u11;
				 if (det_para.Index_Chirp_NotMove_OSCFAR_u11 > GateLimit_Stationary_u5 && det_para.Index_Chirp_NotMove_OSCFAR_u11 <= det_para.ChirpNum_u11 - GateLimit_Stationary_u5)
				 {
					 if (ChirpDiff_s11 >= -GateLimit_Stationary_u5 && ChirpDiff_s11 <= GateLimit_Stationary_u5)
					 {
						 uint16_t IndexR_Left_u10 = Index_R - 2;
						 float SumChLeft_RIndex_u30 = SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,IndexR_Left_u10});
						 float Threshold3_u30 = SumChLeft_RIndex_u30;
						 float Threshold4_u39 = Threshold3_u30 * det_para.Threshold_RangeDim_For_2D_OSCFAR_u9;
						 if (DataToDetect <= Threshold4_u39)
						 {
							 IsTarget_2D = 0;
							 SNR[0] = 0;
							 SNR[1] = 0;
						 }
					 }
				 }
				 else if (det_para.Index_Chirp_NotMove_OSCFAR_u11 <= GateLimit_Stationary_u5)
				 {
					 uint16_t LeftBoundary_u10 = det_para.ChirpNum_u11 - GateLimit_Stationary_u5 + det_para.Index_Chirp_NotMove_OSCFAR_u11;
					 uint16_t RightBoudnary_u11 = GateLimit_Stationary_u5 + det_para.Index_Chirp_NotMove_OSCFAR_u11;
					 if (Index_V > LeftBoundary_u10 || Index_V <= RightBoudnary_u11)
					 {
						 uint16_t IndexR_Left_u10 = Index_R - 2;
						 float SumChLeft_RIndex_u30 = SpatialFFTVelSel_VeloNum_RangeNum.get({ Index_V,IndexR_Left_u10 });
						 float Threshold3_u30 = SumChLeft_RIndex_u30;
						 float Threshold4_u39 = Threshold3_u30 * det_para.Threshold_RangeDim_For_2D_OSCFAR_u9;
						 if (DataToDetect <= Threshold4_u39)
						 {
							 IsTarget_2D = 0;
							 SNR[0] = 0;
							 SNR[1] = 0;
						 }
					 }
				 }
				 else if (det_para.Index_Chirp_NotMove_OSCFAR_u11 > det_para.ChirpNum_u11 - GateLimit_Stationary_u5)
				 {
					 uint16_t LeftBoundary_u10 = det_para.Index_Chirp_NotMove_OSCFAR_u11 - GateLimit_Stationary_u5;
					 uint16_t RightBoudnary_u11 = det_para.ChirpNum_u11 - det_para.Index_Chirp_NotMove_OSCFAR_u11;
					 if (Index_V >= LeftBoundary_u10 || Index_V <= RightBoudnary_u11)
					 {
						 uint16_t IndexR_Left_u10 = Index_R - 2;
						 float SumChLeft_RIndex_u30 = SpatialFFTVelSel_VeloNum_RangeNum.get({ Index_V,IndexR_Left_u10 });
						 float Threshold3_u30 = SumChLeft_RIndex_u30;
						 float Threshold4_u39 = Threshold3_u30 * det_para.Threshold_RangeDim_For_2D_OSCFAR_u9;
						 if (DataToDetect <= Threshold4_u39)
						 {
							 IsTarget_2D = 0;
							 SNR[0] = 0;
							 SNR[1] = 0;
						 }
					 }
				 }
			}
		}
		else
		{
			IsTarget_2D = 0;
			SNR[0] = 0;
			SNR[1] = 0;
		}
	}
	else
	{
		IsTarget_2D = 0;
		SNR[0] = 0;
		SNR[1] = 0;
	}
}

void func_PeakSearch_And_CFAR_2D_Cross(uint16_t& TarNum_Detected, vector<uint16_t>& peak_R, vector<uint16_t>& peak_V, vector<float>& peak_Val, vector<vector<float>>& peak_SNR, DetPara& det_para, MultiDimensionalVector<float,2>& SpatialFFTVelSel_VeloNum_RangeNum)
{
	uint16_t TarCountLimit_u9 = 200;
	uint8_t PeakSearchWin = 1;
	uint16_t ii = 0;
	for (uint16_t Index_R = det_para.DetectCell_RIndex_Min_u10-1; Index_R < det_para.DetectCell_RIndex_Max_u10; Index_R++)
	{
		for (uint16_t Index_V = det_para.DetectCell_VIndex_Min_u11-1; Index_V < det_para.DetectCell_VIndex_Max_u11; Index_V++)
		{
			float DataToDetect = SpatialFFTVelSel_VeloNum_RangeNum.get({ Index_V,Index_R });
			uint16_t V_Upboundary_u11 = 0;
			uint16_t V_Backboundary_u11 = 0;
			uint8_t Logic_PeakCondition = 0;
			if (Index_V < det_para.ChirpNum_u11-1)
			{
				V_Upboundary_u11 = Index_V + 1;
			}
			else
			{
				V_Upboundary_u11 = 0;
			}
			if (Index_V > 0)
			{
				V_Backboundary_u11 = Index_V - 1;
			}
			else
			{
				V_Backboundary_u11 = det_para.ChirpNum_u11 - 1;
			}
			if (det_para.PeakS_Enable_u1 == 1)
			{
				if (PeakSearchWin == 0)
				{
					if (Index_R == 0)
					{
						uint16_t R_front = Index_R + 1;
						if (DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({ Index_V,R_front}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({ V_Upboundary_u11,Index_R }) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,Index_R}))
						{
							Logic_PeakCondition = 1;
						}
					}
					else if (Index_R == det_para.RangeCellNum_u10 - 1)
					{
						uint16_t R_back = Index_R - 1;
						if (DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,R_back}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({ V_Upboundary_u11,Index_R }) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,Index_R}))
						{
							Logic_PeakCondition = 1;
						}
					}
					else
					{
						uint16_t R_front = Index_R + 1;
						uint16_t R_back = Index_R - 1;
						if (DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,R_back }) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,R_front }) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Upboundary_u11,Index_R}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,Index_R}))
						{
							Logic_PeakCondition = 1;
						}
					}
				}
				else
				{
					if (Index_R == 0)
					{
						uint16_t R_front = Index_R + 1;
						if (DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,R_front}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,R_front}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Upboundary_u11,R_front}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Upboundary_u11,Index_R}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,Index_R}))
						{
							Logic_PeakCondition = 1;
						}
					}
					else if (Index_R == det_para.RangeCellNum_u10-1)
					{
						uint16_t R_back = Index_R - 1;
						if (DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,R_back}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Upboundary_u11,R_back}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,R_back}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({ V_Upboundary_u11,Index_R }) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,Index_R}))
						{
							Logic_PeakCondition = 1;
						}
					}
					else
					{
						uint16_t R_back = Index_R - 1;
						uint16_t R_front = Index_R + 1;
						if (DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,R_back}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Upboundary_u11,R_back}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,R_back}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({Index_V,R_front}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Upboundary_u11,R_front}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,R_front}) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({ V_Upboundary_u11,Index_R }) && \
							DataToDetect > SpatialFFTVelSel_VeloNum_RangeNum.get({V_Backboundary_u11,Index_R}))
						{
							Logic_PeakCondition = 1;
						}

					}
				}
			}
			else
			{
				Logic_PeakCondition = 1;
			}
			uint8_t ptIdx = ceil(Index_R/(det_para.RangeCellNum_u10/32));
			if (Logic_PeakCondition == 1 && DataToDetect < PeakThreshold1_u30[ptIdx])
			{
				Logic_PeakCondition = 0;
			}
			if (Logic_PeakCondition == 1)
			{
				uint8_t typeSwitch = 0;
				uint8_t IsTarget_2D = 0;
				vector<float> SNR(2,0.0);
				if (det_para.CFARTypeSwitch_u2.compare("00") == 0)
				{
					typeSwitch = 0;
				}
				else if(det_para.CFARTypeSwitch_u2.compare("01") == 0)
				{
					typeSwitch = 1;
				}
				else if (det_para.CFARTypeSwitch_u2.compare("10") == 0)
				{
					typeSwitch = 2;
				}
				else if (det_para.CFARTypeSwitch_u2.compare("11") == 0)
				{
					typeSwitch = 3;
				}
				switch (typeSwitch)
				{
				case 0:
					break;
				case 1:
					break;
				case 2:
					break;
				case 3:
					cfarOS_cal_single_point_2D_Cross(IsTarget_2D,SNR,det_para,DataToDetect,Index_R,Index_V, SpatialFFTVelSel_VeloNum_RangeNum);
					break;
				default:
					break;
				}
				if (IsTarget_2D == 1)
				{
					peak_R.push_back(Index_R);
					peak_V.push_back(Index_V);
					peak_Val.push_back(DataToDetect);
					SNR[0] = 20 * log10(SNR[0]);
					SNR[1] = 20 * log10(SNR[1]);
					peak_SNR.push_back(SNR);
					ii = ii + 1;
				}
				if (ii > TarCountLimit_u9 - 1)
				{
					break;
				}
			}
		}
		if (ii > TarCountLimit_u9 - 1)
		{
			break;
		}
	}
	if (ii <= 0)
	{
		peak_R = { 0 };
		peak_V = { 0 };
		TarNum_Detected = 0;
		peak_Val = { 0 };
		peak_SNR[0] = {0,0};
	}
	else
	{
		TarNum_Detected = ii;
	}
}

void cfar3d_cal_across_ArbitaryDim(MultiDimensionalVector<float, 2>& SpatialFFTVelSel_VeloNum_RangeNum, MultiDimensionalVector<float, 3>& DataInput, uint8_t SqueezeDim, uint8_t cfar_include_order, uint8_t cfar_exclude_order, float min_snr, uint8_t Switch3DMode, uint16_t horNum, uint16_t chirpNum, uint16_t sampleNum)
{
	// 方位维取最大值
	for (uint16_t chirpIdx = 0; chirpIdx < chirpNum; chirpIdx++)
	{
		for (uint16_t sampleIdx = 0; sampleIdx < sampleNum; sampleIdx++)
		{
			float tempMax = INT32_MIN;
			for (uint16_t horIdx = 0; horIdx < horNum; horIdx++)
			{
				if (DataInput.get({ horIdx,chirpIdx,sampleIdx }) > tempMax)
				{
					tempMax = DataInput.get({ horIdx,chirpIdx,sampleIdx });
				}
			}
			SpatialFFTVelSel_VeloNum_RangeNum.set({ chirpIdx,sampleIdx }, tempMax);
		}
	}

	// ...
}

void func_winA_process(vector<vector<vector<vector<vector<complex<float>>>>>>& winAout_xNum_yNum_VeloFFTNum_RangeNum_waveLocNum, vector<uint8_t>& fft_win_type, uint16_t win_len, uint16_t VirtArrVertGridLen, uint16_t VelocityNum, uint16_t RangeSampleNum, uint16_t waveLocNum, vector<vector<vector<vector<vector<complex<float>>>>>>& result_xNum_yNum_VeloFFTNum_RangeNum)
{
	uint8_t wina_type = fft_win_type[2];
	vector<float> win_Coef(win_len);
	switch (wina_type)
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
		for (uint16_t sampleIdx = 0; sampleIdx < RangeSampleNum; sampleIdx++)
		{
			for (uint16_t chirpIdx = 0; chirpIdx < VelocityNum; chirpIdx++)
			{
				for (uint16_t yIdx = 0; yIdx < VirtArrVertGridLen; yIdx++)
				{
					for (uint16_t xIdx = 0; xIdx < win_len; xIdx++)
					{
						winAout_xNum_yNum_VeloFFTNum_RangeNum_waveLocNum[xIdx][yIdx][chirpIdx][sampleIdx][waveLocIdx] = win_Coef[xIdx] * result_xNum_yNum_VeloFFTNum_RangeNum[xIdx][yIdx][chirpIdx][sampleIdx][waveLocIdx];
					}
				}
			}
		}
	}
}

void func_Spatial_Reorder(vector<vector<vector<vector<vector<complex<float>>>>>>& result_xNum_yNum_VeloFFTNum_RangeNum, vector<vector<vector<vector<complex<float>>>>>& FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum, uint8_t Array_option, uint16_t MIMONum, vector<vector<uint16_t>>& pos_in_mat, uint16_t VirtArrHorGridLen, uint16_t VirtArrVertGridLen, uint16_t VelocityNum, uint16_t CoarseRangeNum,uint16_t waveLocNum)
{
	vector<vector<vector<vector<complex<float>>>>> ModifyRDFFTReshape_VirtArrTotalNum_VeloNum_RangeBinNum(VirtArrHorGridLen * VirtArrVertGridLen, vector<vector<vector<complex<float>>>>(VelocityNum, vector<vector<complex<float>>>(CoarseRangeNum, vector<complex<float>>(waveLocNum, 0.0))));
	switch (Array_option)
	{
	case 0:
		for (uint16_t mimoIdx = 0; mimoIdx < MIMONum; mimoIdx++)
		{
			for (uint16_t chirpIdx = 0; chirpIdx < VelocityNum; chirpIdx++) 
			{
				for (uint16_t sampleIdx = 0; sampleIdx < CoarseRangeNum; sampleIdx++)
				{
					for (uint16_t waveLocIdx = 0; waveLocIdx < waveLocNum; waveLocIdx++)
					{
						result_xNum_yNum_VeloFFTNum_RangeNum[mimoIdx][0][chirpIdx][sampleIdx][waveLocIdx] = FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum[mimoIdx][chirpIdx][sampleIdx][waveLocIdx];
					}
				}
			}
		}
	case 1:
		for (uint16_t chirpIdx = 0; chirpIdx < VelocityNum; chirpIdx++)
		{
			for(uint16_t sampleIdx = 0; sampleIdx < CoarseRangeNum; sampleIdx++)
			{
				for (uint16_t waveLocIdx = 0; waveLocIdx < waveLocNum; waveLocIdx++)
				{
					for(uint16_t mimoIdx = 0; mimoIdx < MIMONum; mimoIdx++)
					{
						uint16_t  posIdx = pos_in_mat[mimoIdx][1] + 1 + (pos_in_mat[mimoIdx][0] + 1 - 1) * VirtArrVertGridLen - 1;
						ModifyRDFFTReshape_VirtArrTotalNum_VeloNum_RangeBinNum[posIdx][chirpIdx][sampleIdx][waveLocIdx] = FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum[mimoIdx][chirpIdx][sampleIdx][waveLocIdx];
					}
				}
			}
		}
		// MATLAB中waveloc在这里没了,是有问题的
		for (uint16_t chirpIdx = 0; chirpIdx < VelocityNum; chirpIdx++)
		{
			for (uint16_t sampleIdx = 0; sampleIdx < CoarseRangeNum; sampleIdx++)
			{
				for (uint16_t waveLocIdx = 0; waveLocIdx < waveLocNum; waveLocIdx++)
				{
					for (uint16_t totalIdx = 0; totalIdx < VirtArrHorGridLen*VirtArrVertGridLen; totalIdx++)
					{
						uint16_t xIdx = totalIdx/VirtArrVertGridLen;
						uint16_t yIdx = fmod(totalIdx, VirtArrVertGridLen);
						result_xNum_yNum_VeloFFTNum_RangeNum[xIdx][yIdx][chirpIdx][sampleIdx][waveLocIdx] = ModifyRDFFTReshape_VirtArrTotalNum_VeloNum_RangeBinNum[totalIdx][chirpIdx][sampleIdx][waveLocIdx];
					}
				}
			}
		}
	case 2:
		for (uint16_t chirpIdx = 0; chirpIdx < VelocityNum; chirpIdx++)
		{
			for (uint16_t sampleIdx = 0; sampleIdx < CoarseRangeNum; sampleIdx++)
			{
				for (uint16_t waveLocIdx = 0; waveLocIdx < waveLocNum; waveLocIdx++)
				{
					for (uint16_t mimoIdx = 0; mimoIdx < MIMONum; mimoIdx++)
					{
						uint16_t  posIdx = pos_in_mat[mimoIdx][1] + 1 + (pos_in_mat[mimoIdx][0] + 1 - 1) * VirtArrVertGridLen - 1;
						ModifyRDFFTReshape_VirtArrTotalNum_VeloNum_RangeBinNum[posIdx][chirpIdx][sampleIdx][waveLocIdx] = FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum[mimoIdx][chirpIdx][sampleIdx][waveLocIdx];
					}
				}
			}
		}
		// MATLAB中waveloc在这里没了,是有问题的
		for (uint16_t chirpIdx = 0; chirpIdx < VelocityNum; chirpIdx++)
		{
			for (uint16_t sampleIdx = 0; sampleIdx < CoarseRangeNum; sampleIdx++)
			{
				for (uint16_t waveLocIdx = 0; waveLocIdx < waveLocNum; waveLocIdx++)
				{
					for (uint16_t totalIdx = 0; totalIdx < VirtArrHorGridLen * VirtArrVertGridLen; totalIdx++)
					{
						uint16_t xIdx = totalIdx / VirtArrVertGridLen;
						uint16_t yIdx = fmod(totalIdx, VirtArrVertGridLen);
						result_xNum_yNum_VeloFFTNum_RangeNum[xIdx][yIdx][chirpIdx][sampleIdx][waveLocIdx] = ModifyRDFFTReshape_VirtArrTotalNum_VeloNum_RangeBinNum[totalIdx][chirpIdx][sampleIdx][waveLocIdx];
					}
				}
			}
		}
	case 3:
		for (uint16_t chirpIdx = 0; chirpIdx < VelocityNum; chirpIdx++)
		{
			for (uint16_t sampleIdx = 0; sampleIdx < CoarseRangeNum; sampleIdx++)
			{
				for (uint16_t waveLocIdx = 0; waveLocIdx < waveLocNum; waveLocIdx++)
				{
					for (uint16_t mimoIdx = 0; mimoIdx < MIMONum; mimoIdx++)
					{
						uint16_t  posIdx = pos_in_mat[mimoIdx][1] + 1 + (pos_in_mat[mimoIdx][0] + 1 - 1) * VirtArrVertGridLen - 1;
						ModifyRDFFTReshape_VirtArrTotalNum_VeloNum_RangeBinNum[posIdx][chirpIdx][sampleIdx][waveLocIdx] = FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum[mimoIdx][chirpIdx][sampleIdx][waveLocIdx];
					}
				}
			}
		}
		// MATLAB中waveloc在这里没了,是有问题的
		for (uint16_t chirpIdx = 0; chirpIdx < VelocityNum; chirpIdx++)
		{
			for (uint16_t sampleIdx = 0; sampleIdx < CoarseRangeNum; sampleIdx++)
			{
				for (uint16_t waveLocIdx = 0; waveLocIdx < waveLocNum; waveLocIdx++)
				{
					for (uint16_t totalIdx = 0; totalIdx < VirtArrHorGridLen * VirtArrVertGridLen; totalIdx++)
					{
						uint16_t xIdx = totalIdx / VirtArrVertGridLen;
						uint16_t yIdx = fmod(totalIdx, VirtArrVertGridLen);
						result_xNum_yNum_VeloFFTNum_RangeNum[xIdx][yIdx][chirpIdx][sampleIdx][waveLocIdx] = ModifyRDFFTReshape_VirtArrTotalNum_VeloNum_RangeBinNum[totalIdx][chirpIdx][sampleIdx][waveLocIdx];
					}
				}
			}
		}
	case 4:
		;
	default:
		break;
	}
}


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
	__TIC__(DEFVEC)
	vector<vector<vector<vector<vector<complex<float>>>>>> SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum(para_sys.AngleHorNum,vector<vector<vector<vector<complex<float>>>>>(VirtArrVertGridLen,vector<vector<vector<complex<float>>>>(para_sys.VelocityNum,vector<vector<complex<float>>>(para_sys.CoarseRangeNum,vector<complex<float>>(para_sys.waveLocNum, 0.0)))));
	__TOC__(DEFVEC)
	//vector<vector<vector<vector<vector<complex<float>>>>>> SpatialFFT_AngleVertNum_AngleHorNum_VeloFFTNum_RangeBinNum(para_sys.AngleVertNum,vector<vector<vector<vector<complex<float>>>>>(para_sys.AngleHorNum,vector<vector<vector<complex<float>>>>(para_sys.VelocityNum,vector<vector<complex<float>>>(para_sys.CoarseRangeNum,vector<complex<float>>(para_sys.waveLocNum, complex<float>(0.0,0.0))))));
	/*__TIC__(DEFVEC)
	MultiDimensionalVector<complex<float>, 5> SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum({ para_sys.AngleHorNum, VirtArrVertGridLen, para_sys.VelocityNum, para_sys.CoarseRangeNum, para_sys.waveLocNum });
	__TOC__(DEFVEC)*/

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
		/*vector<complex<float>>  windChirp(para_sys.CoarseRangeNum);
		for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
		{
			windChirp[sampleIdx] = WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum[5][sampleIdx][24][0];
		}*/
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
		/*vector<complex<float>> fftdChirp(para_sys.CoarseRangeNum);
		for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
		{
			fftdChirp[sampleIdx] = FFT2D_VeloFFTNum_CoarseRangeBinNum_MIMONum_WaveLocNum[5][sampleIdx][24][0];
		}*/
		
		// save fftd as a file
	/*	string filename_2dout = "fftdOut.bin";
		std::ofstream file(filename_2dout, std::ios::out | std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "Error: Cannot open file for writing." << std::endl;
		}
		for (uint16_t i = 0; i < para_sys.waveLocNum; i++)
		{
			for (uint16_t j = 0; j < MIMONum; j++)
			{
				for (uint16_t m = 0; m < para_sys.CoarseRangeNum; m++)
				{
					for (size_t n = 0; n < para_sys.VelocityNum; n++)
					{
						float real_part = FFT2D_VeloFFTNum_CoarseRangeBinNum_MIMONum_WaveLocNum[n][m][j][i].real();
						float imag_part = FFT2D_VeloFFTNum_CoarseRangeBinNum_MIMONum_WaveLocNum[n][m][j][i].imag();
						file.write(reinterpret_cast<const char*>(&real_part), sizeof(float));
						file.write(reinterpret_cast<const char*>(&imag_part), sizeof(float));
					}
				}
			}
		}
		file.close();*/
	}
	if (CoarseFrame_CFARdim > 2)
	{
		vector<vector<vector<vector<complex<float>>>>> FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum(MIMONum, vector<vector<vector<complex<float>>>>(para_sys.VelocityNum, vector<vector<complex<float>>>(para_sys.CoarseRangeNum, vector<complex<float>>(para_sys.waveLocNum, 0.0))));
		vector<vector<vector<vector<vector<complex<float>>>>>> result_xNum_yNum_VeloFFTNum_RangeNum(VirtArrHorGridLen, vector<vector<vector<vector<complex<float>>>>>(VirtArrVertGridLen, vector<vector<vector<complex<float>>>>(para_sys.VelocityNum, vector<vector<complex<float>>>(para_sys.CoarseRangeNum, vector<complex<float>>(para_sys.waveLocNum,0.0)))));
		vector<vector<vector<vector<vector<complex<float>>>>>> winAout_xNum_yNum_VeloFFTNum_RangeNum_waveLocNum(VirtArrHorGridLen, vector<vector<vector<vector<complex<float>>>>>(VirtArrVertGridLen, vector<vector<vector<complex<float>>>>(para_sys.VelocityNum, vector<vector<complex<float>>>(para_sys.CoarseRangeNum, vector<complex<float>>(para_sys.waveLocNum, 0.0)))));
		if (DDMA_EN == 1)
		{
			vector<vector<uint16_t>> pos_in_mat_ddma(48, vector<uint16_t>(2));
			for (uint16_t virIdx = 0; virIdx < 48; virIdx++)
			{
				pos_in_mat_ddma[virIdx][0] = (virtual_array.pos_in_mat[virIdx][0] + 1) / 2;
				pos_in_mat_ddma[virIdx][1] = virtual_array.pos_in_mat[virIdx][1];
			}
			uint16_t VirtArrHorGridLen_ddma = 27;
			uint16_t VirtArrVertGridLen_ddma = 2;
			for (uint16_t mimoIdx = 0; mimoIdx < MIMONum; mimoIdx++)
			{
				for (uint16_t dopplerIdx = 0; dopplerIdx < para_sys.VelocityNum; dopplerIdx++)
				{
					for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
					{
						for (uint16_t waveLocIdx = 0; waveLocIdx < para_sys.waveLocNum; waveLocIdx++)
						{
							FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum[mimoIdx][dopplerIdx][sampleIdx][waveLocIdx] = FFT2D_VeloFFTNum_CoarseRangeBinNum_MIMONum_WaveLocNum[dopplerIdx][sampleIdx][mimoIdx][waveLocIdx];
						}
					}
				}
			}
			func_Spatial_Reorder(result_xNum_yNum_VeloFFTNum_RangeNum, FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum, para_sys.Array_option, MIMONum, pos_in_mat_ddma, VirtArrHorGridLen, VirtArrVertGridLen, para_sys.VelocityNum, para_sys.CoarseRangeNum, para_sys.waveLocNum);

		}
		else
		{
			for (uint16_t mimoIdx = 0; mimoIdx < MIMONum; mimoIdx++)
			{
				for (uint16_t dopplerIdx = 0; dopplerIdx < para_sys.VelocityNum; dopplerIdx ++) 
				{
					for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
					{
						for (uint16_t waveLocIdx = 0; waveLocIdx < para_sys.waveLocNum; waveLocIdx ++)
						{
							FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum[mimoIdx][dopplerIdx][sampleIdx][waveLocIdx] = FFT2D_VeloFFTNum_CoarseRangeBinNum_MIMONum_WaveLocNum[dopplerIdx][sampleIdx][mimoIdx][waveLocIdx];
						}
					}
				}
			}
			//test code
			/*vector<complex<float>> beforeWinOne(MIMONum);
			for (uint16_t mimoIdx = 0; mimoIdx < MIMONum; mimoIdx++)
			{
				beforeWinOne[mimoIdx] = FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum[mimoIdx][24][46][0];
			}*/
			func_Spatial_Reorder(result_xNum_yNum_VeloFFTNum_RangeNum, FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum,para_sys.Array_option,MIMONum,virtual_array.pos_in_mat,VirtArrHorGridLen,VirtArrVertGridLen,para_sys.VelocityNum,para_sys.CoarseRangeNum,para_sys.waveLocNum);
			//test code
			/*vector<complex<float>> afterReorderOne(VirtArrHorGridLen);
			for (uint16_t xIdx = 0; xIdx < VirtArrHorGridLen; xIdx++)
			{
				afterReorderOne[xIdx] = result_xNum_yNum_VeloFFTNum_RangeNum[xIdx][0][56][89][0];
			}*/
		}
		// winA (matlab需要加一个分支，不加窗的分支)
		func_winA_process(winAout_xNum_yNum_VeloFFTNum_RangeNum_waveLocNum,para_sys.fft_win_type , VirtArrHorGridLen, VirtArrVertGridLen, para_sys.VelocityNum, para_sys.CoarseRangeNum, para_sys.waveLocNum,result_xNum_yNum_VeloFFTNum_RangeNum);
		// test code
		/*vector<complex<float>> winAone(VirtArrHorGridLen);
		for (uint16_t xIdx = 0; xIdx < VirtArrHorGridLen; xIdx++)
		{
			winAone[xIdx] = winAout_xNum_yNum_VeloFFTNum_RangeNum_waveLocNum[xIdx][0][89][40][0];
		}*/
		// FFTA
		__TIC__(USEVEC)
		for (uint16_t waveLocIdx = 0; waveLocIdx < para_sys.waveLocNum; waveLocIdx++)
		{
			for (uint16_t rangeIdx = 0; rangeIdx < para_sys.CoarseRangeNum; rangeIdx++)
			{
				for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum; chirpIdx++)
				{
					for (uint16_t yIdx = 0; yIdx < VirtArrVertGridLen; yIdx++)
					{
						fftw_complex* in;
						fftw_complex* out;
						fftw_plan pffta;
						in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * para_sys.AngleHorNum);
						out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * para_sys.AngleHorNum);
						for (uint16_t xIdx = 0; xIdx < VirtArrHorGridLen; xIdx++)
						{
							in[xIdx][0] = real(winAout_xNum_yNum_VeloFFTNum_RangeNum_waveLocNum[xIdx][yIdx][chirpIdx][rangeIdx][waveLocIdx]);
							in[xIdx][1] = imag(winAout_xNum_yNum_VeloFFTNum_RangeNum_waveLocNum[xIdx][yIdx][chirpIdx][rangeIdx][waveLocIdx]);
						}
						for (uint16_t xIdx = VirtArrHorGridLen; xIdx < para_sys.AngleHorNum; xIdx++)
						{
							in[xIdx][0] = 0;
							in[xIdx][1] = 0;
						}
						pffta = fftw_plan_dft_1d(para_sys.AngleHorNum, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
						fftw_execute(pffta);
						for (uint16_t horIdx = 0; horIdx < para_sys.AngleHorNum / 2; horIdx++)
						{
							SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum[horIdx + para_sys.AngleHorNum / 2][yIdx][chirpIdx][rangeIdx][waveLocIdx] = complex<float>(out[horIdx][0], out[horIdx][1]);
							//SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum.set({ uint16_t(horIdx + para_sys.AngleHorNum / 2),yIdx, chirpIdx, rangeIdx, waveLocIdx }, complex<float>(out[horIdx][0], out[horIdx][1]));
						}
						for (uint16_t horIdx = para_sys.AngleHorNum / 2; horIdx < para_sys.AngleHorNum; horIdx++)
						{
							SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum[horIdx - para_sys.AngleHorNum / 2][yIdx][chirpIdx][rangeIdx][waveLocIdx] = complex<float>(out[horIdx][0], out[horIdx][1]);
							//SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum.set({ uint16_t(horIdx - para_sys.AngleHorNum / 2),yIdx,chirpIdx,rangeIdx,waveLocIdx} , complex<float>(out[horIdx][0], out[horIdx][1]));
						}
						fftw_destroy_plan(pffta);
						fftw_free(in);
						fftw_free(out);
					}
				}
			}		
		}
		__TOC__(USEVEC)
		// test code
		//vector<complex<float>> fftaOne(para_sys.AngleHorNum);
		//for (uint16_t horIdx = 0; horIdx < para_sys.AngleHorNum; horIdx++)
		//{
		//	fftaOne[horIdx] = SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum[horIdx][0][560][54][0];
		//	//fftaOne[horIdx] = SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum.get({horIdx,0,560,54,0});
		//}

	    // save ffta as a file
		/*string filename_3dout = "fftaOut.bin";
		std::ofstream file(filename_3dout, std::ios::out | std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "Error: Cannot open file for writing." << std::endl;
		}
		for (uint16_t i = 0; i < para_sys.waveLocNum; i++)
		{
			for (uint16_t j = 0; j < para_sys.CoarseRangeNum; j++)
			{
				for (uint16_t m = 0; m < para_sys.VelocityNum; m++)
				{
					for (size_t n = 0; n < VirtArrVertGridLen; n++)
					{
						for (size_t k = 0; k < para_sys.AngleHorNum; k++)
						{
							float real_part = SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum[k][n][m][j][i].real();
							float imag_part = SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum[k][n][m][j][i].imag();
							file.write(reinterpret_cast<const char*>(&real_part), sizeof(float));
							file.write(reinterpret_cast<const char*>(&imag_part), sizeof(float));
						}
					}
				}
			}
		}
		file.close();*/
	}
	if (CoarseFrame_CFARdim > 3)
	{
		// 数据量太大暂不实现
		;
	}

	// cfar 
	for (uint16_t i = 0; i < 32; i++)
	{
		R_Threshold1_u9[i] = 5;
		PeakThreshold1_u30[i] = 0;
	}
	for (uint16_t i = 0; i < 16; i++)
	{
		V_Threshold1_u9[i] = 8;
	}
	switch (CoarseFrame_CFARdim)
	{
	case 2:
		;
	case 3:
		// 质心提取
		if (DDMA_EN == 1)
		{
			MultiDimensionalVector<float, 5> SpatialFFTA_ABS_AngleHorNum_yNum_VeloFFTNum_RangeNum({ para_sys.AngleHorNum, VirtArrVertGridLen, para_sys.VelocityNum, para_sys.CoarseRangeNum, para_sys.waveLocNum });
			MultiDimensionalVector<float, 4> SpatialFFTA_ABSMean_AngleHorNum_VeloFFTNum_RangeNum({ para_sys.AngleHorNum, para_sys.VelocityNum, para_sys.CoarseRangeNum, para_sys.waveLocNum });
			MultiDimensionalVector<float, 3> SpatialFFTA_ABSMeanMax_VeloFFTNum_RangeNum({ para_sys.VelocityNum, para_sys.CoarseRangeNum, para_sys.waveLocNum });

			for (uint16_t waveIdx = 0; waveIdx < para_sys.waveLocNum; waveIdx++)
			{
				for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
				{
					for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum; chirpIdx++)
					{
						for (uint16_t yIdx = 0; yIdx < VirtArrVertGridLen; yIdx++)
						{
							for (uint16_t horIdx = 0; horIdx < para_sys.AngleHorNum; horIdx++)
							{
								SpatialFFTA_ABS_AngleHorNum_yNum_VeloFFTNum_RangeNum.set({ horIdx,yIdx,chirpIdx,sampleIdx,waveIdx }, abs(SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum[horIdx][yIdx][chirpIdx][sampleIdx][waveIdx]));
							}
						}
					}
				}
			}
			for (uint16_t waveIdx = 0; waveIdx < para_sys.waveLocNum; waveIdx++)
			{
				for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
				{
					for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum; chirpIdx++)
					{
						for (uint16_t horIdx = 0; horIdx < para_sys.AngleHorNum; horIdx++)
						{
							float tempMean = 0;
							for (uint16_t yIdx = 0; yIdx < VirtArrVertGridLen; yIdx++)
							{
								tempMean = tempMean + SpatialFFTA_ABS_AngleHorNum_yNum_VeloFFTNum_RangeNum.get({ horIdx,yIdx,chirpIdx,sampleIdx,waveIdx });
							}
							SpatialFFTA_ABSMean_AngleHorNum_VeloFFTNum_RangeNum.set({ horIdx,chirpIdx,sampleIdx,waveIdx }, tempMean/ VirtArrVertGridLen);
						}
					}
				}
			}
			for (uint16_t waveIdx = 0; waveIdx < para_sys.waveLocNum; waveIdx++)
			{
				for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
				{
					for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum; chirpIdx++)
					{
						float maxValue = INT32_MIN;
						for (uint16_t horIdx = 0; horIdx < para_sys.AngleHorNum; horIdx++)
						{
							if (SpatialFFTA_ABSMean_AngleHorNum_VeloFFTNum_RangeNum.get({horIdx,chirpIdx,sampleIdx,waveIdx}) > maxValue)
							{
								maxValue = SpatialFFTA_ABSMean_AngleHorNum_VeloFFTNum_RangeNum.get({ horIdx,chirpIdx,sampleIdx,waveIdx });
							}
						}
						SpatialFFTA_ABSMeanMax_VeloFFTNum_RangeNum.set({chirpIdx,sampleIdx,waveIdx},maxValue);
					}
				}
			}
			// test code
			vector<float> tempAbs(para_sys.CoarseRangeNum);
			/*for (size_t i = 0; i < para_sys.CoarseRangeNum; i++)
			{
				tempAbs[i] = SpatialFFTA_ABSMeanMax_VeloFFTNum_RangeNum.get({ 15, i, 0 });
			}*/
		}
		else
		{
			MultiDimensionalVector<float, 5> SpatialFFTA_ABS_AngleHorNum_yNum_VeloFFTNum_RangeNum({ para_sys.AngleHorNum, VirtArrVertGridLen, para_sys.VelocityNum, para_sys.CoarseRangeNum, para_sys.waveLocNum });
			MultiDimensionalVector<float, 4> SpatialFFTA_ABSMean_AngleHorNum_VeloFFTNum_RangeNum({ para_sys.AngleHorNum, para_sys.VelocityNum, para_sys.CoarseRangeNum, para_sys.waveLocNum });
			for (uint16_t waveIdx = 0; waveIdx < para_sys.waveLocNum; waveIdx++)
			{
				for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
				{
					for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum; chirpIdx++)
					{
						for (uint16_t yIdx = 0; yIdx < VirtArrVertGridLen; yIdx++)
						{
							for (uint16_t horIdx = 0; horIdx < para_sys.AngleHorNum; horIdx++)
							{
								SpatialFFTA_ABS_AngleHorNum_yNum_VeloFFTNum_RangeNum.set({ horIdx,yIdx,chirpIdx,sampleIdx,waveIdx }, abs(SpatialFFTA_AngleHorNum_yNum_VeloFFTNum_RangeNum[horIdx][yIdx][chirpIdx][sampleIdx][waveIdx]));
							}
						}
					}
				}
			}
			for (uint16_t waveIdx = 0; waveIdx < para_sys.waveLocNum; waveIdx++)
			{
				for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
				{
					for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum; chirpIdx++)
					{
						for (uint16_t horIdx = 0; horIdx < para_sys.AngleHorNum; horIdx++)
						{
							float tempMean = 0;
							for (uint16_t yIdx = 0; yIdx < VirtArrVertGridLen; yIdx++)
							{
								tempMean = tempMean + SpatialFFTA_ABS_AngleHorNum_yNum_VeloFFTNum_RangeNum.get({ horIdx,yIdx,chirpIdx,sampleIdx,waveIdx });
							}
							SpatialFFTA_ABSMean_AngleHorNum_VeloFFTNum_RangeNum.set({ horIdx,chirpIdx,sampleIdx,waveIdx }, tempMean / VirtArrVertGridLen);
						}
					}
				}
			}
			// test code
			/*vector<float> fft3dAbs(para_sys.CoarseRangeNum);
			for (size_t i = 0; i < para_sys.CoarseRangeNum; i++)
			{
				fft3dAbs[i] = SpatialFFTA_ABSMean_AngleHorNum_VeloFFTNum_RangeNum.get({32,15,i,0});
			}*/
			uint8_t Switch3DMode = 0;
			uint8_t SqueezeDim = 1;
			for (uint16_t Beamcnt = 0; Beamcnt < para_sys.waveLocNum; Beamcnt++)
			{
				MultiDimensionalVector<float, 2> SpatialFFTVelSel_VeloNum_RangeNum({ para_sys.VelocityNum, para_sys.CoarseRangeNum });
				MultiDimensionalVector<float, 3> DataInput({ para_sys.AngleHorNum, para_sys.VelocityNum, para_sys.CoarseRangeNum });
				for (uint16_t horIdx = 0; horIdx < para_sys.AngleHorNum; horIdx++)
				{
					for (uint16_t chirpIdx = 0; chirpIdx < para_sys.VelocityNum; chirpIdx++)
					{
						for (uint16_t sampleIdx = 0; sampleIdx < para_sys.CoarseRangeNum; sampleIdx++)
						{
							DataInput.set({ horIdx,chirpIdx,sampleIdx }, SpatialFFTA_ABSMean_AngleHorNum_VeloFFTNum_RangeNum.get({ horIdx,chirpIdx,sampleIdx,Beamcnt }));
						}
					}
				}
				cfar3d_cal_across_ArbitaryDim(SpatialFFTVelSel_VeloNum_RangeNum,DataInput,SqueezeDim,para_sys.cfar_include_order[2], para_sys.cfar_exclude_order[2], para_sys.snr_dB_different_dim[2], Switch3DMode, para_sys.AngleHorNum,para_sys.VelocityNum,para_sys.CoarseRangeNum);
				// test code
				/*vector<float> max3dOne(para_sys.CoarseRangeNum);
				for (uint16_t i = 0; i < para_sys.CoarseRangeNum; i++)
				{
					max3dOne[i] = SpatialFFTVelSel_VeloNum_RangeNum.get({56,i});
				}*/
				uint16_t TarNum_Detected = 0;
				vector<uint16_t> peak_R;
				vector<uint16_t> peak_V;
				vector<float> peak_Val;
				vector<vector<float>> peak_SNR;
				func_PeakSearch_And_CFAR_2D_Cross(TarNum_Detected,peak_R,peak_V,peak_Val,peak_SNR,para_sys.det_para, SpatialFFTVelSel_VeloNum_RangeNum);
			}
		}
		// 生成联通区域

	case 4:
		// 暂不实现
	default:
		break;
	}
}

void func_signal_process_coarse(vector<vector<float>>& TOI, vector<vector<float>>& point_info, const string& filename, ParaSys& para_sys, Virtual_array& virtual_array, unique_ptr<int16_t[]>& radarInputdata, vector<complex<float>>& compensate_mat) {

	vector<vector<vector<complex<float>>>> CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum(para_sys.CoarseRangeNum, vector<vector<complex<float>>>(para_sys.VelocityNum, vector<complex<float>>(para_sys.RxNum)));
	//MultiDimensionalVector<complex<float>, 3> CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum({ para_sys.CoarseRangeNum ,para_sys.VelocityNum ,para_sys.RxNum });
	if (para_sys.InputHasDonePreProcess == 0) {
		// winr 
		vector<vector<vector<float>>> WinRout_RangeSampleNum_ChirpNum_RxNum(para_sys.RangeSampleNum, vector<vector<float>>(para_sys.VelocityNum, vector<float>(para_sys.RxNum)));
		func_winR_process(WinRout_RangeSampleNum_ChirpNum_RxNum, para_sys.fft_win_type, para_sys.RangeSampleNum, para_sys.VelocityNum, para_sys.RxNum, radarInputdata);
		
		// test code
		/*ector<float> temp(para_sys.RangeSampleNum);
		for (int i = 0; i < para_sys.RangeSampleNum; i++) {
			temp[i] = WinRout_RangeSampleNum_ChirpNum_RxNum[i][45][23];
		}*/

		// fftr
		//vector<vector<vector<complex<float>>>> CoarseRangeFFT_SampleNum_ChirpNum_RxNum(para_sys.RangeSampleNum, vector<vector<complex<float>>>(para_sys.VelocityNum, vector<complex<float>>(para_sys.RxNum)));
		for (uint16_t i = 0; i < para_sys.VelocityNum; i++) {
			for (uint16_t j = 0; j < para_sys.RxNum; j++) {
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
				for (uint16_t k = 0; k < para_sys.RangeSampleNum/2; k++) {
					CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[k][i][j] = complex<float>(out[k][0],out[k][1]);
					//CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum.set({k,i,j}, complex<float>(out[k][0], out[k][1]));
				}
				fftw_destroy_plan(pfft);
				fftw_free(in);
				fftw_free(out);
			}
		}
		// test code
		/*vector<complex<float>> tempfftr(para_sys.RangeSampleNum/2);
		for (int i = 0; i < para_sys.RangeSampleNum / 2; i++) {
			tempfftr[i] = CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[i][45][23];
		}*/
	}
	else if (para_sys.InputHasDonePreProcess == 1) {
		// reshape，复数的实部和虚部顺序存储
		for (uint16_t i = 0; i < para_sys.VelocityNum; i++)
		{
			for (uint16_t j = 0; j < para_sys.RxNum; j++)
			{
				for (uint16_t k = 0; k < para_sys.CoarseRangeNum; k++)
				{
					CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[k][i][j].real(radarInputdata[2 * k + j * para_sys.CoarseRangeNum + i * para_sys.CoarseRangeNum * para_sys.RxNum]);
					CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[k][i][j].imag(radarInputdata[2 * k + 1 + j * para_sys.CoarseRangeNum + i * para_sys.CoarseRangeNum * para_sys.RxNum]);
				}

			}

		}
	}
	else {
		throw invalid_argument("Input value is illegal.");
	}
	// test code
	/*vector<complex<float>> fftrOutOne(para_sys.CoarseRangeNum);
	for (size_t i = 0; i < para_sys.CoarseRangeNum; i++)
	{
		fftrOutOne[i] = CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[i][0][0];
	}*/
	// save the result as a bin file
	//string filename_1dout = "fftrOut.bin";
	//std::ofstream file(filename_1dout, std::ios::out | std::ios::binary);
	//if (!file.is_open()) {
	//	std::cerr << "Error: Cannot open file for writing." << std::endl;
	//}
	//for (uint16_t i = 0; i < para_sys.RxNum; i++)
	//{
	//	for (uint16_t j = 0; j < para_sys.VelocityNum; j++)
	//	{
	//		for (size_t k = 0; k < para_sys.CoarseRangeNum; k++)
	//		{
	//			float real_part = CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[k][j][i].real();
	//			float imag_part = CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum[k][j][i].imag();
	//			file.write(reinterpret_cast<const char*>(&real_part),sizeof(float));
	//			file.write(reinterpret_cast<const char*>(&imag_part),sizeof(float));
	//		}

	//	}

	//}
	//file.close();
	FFTD_SpatialFFT_CFAR_CoarseFrame(point_info, TOI, para_sys, virtual_array, compensate_mat, CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum);
}