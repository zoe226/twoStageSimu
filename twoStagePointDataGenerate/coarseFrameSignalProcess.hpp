#include "param.hpp"
#include <cassert>
#include <chrono>
#include <iostream>

using Clock = std::chrono::steady_clock;

#define __TIC__(tag) auto __##tag##_start_time = Clock::now();

#define __TOC__(tag)                                                           \
  auto __##tag##_end_time = Clock::now();                                      \
  std::cout << #tag << " : "                                                   \
            << std::chrono::duration_cast<std::chrono::microseconds>(          \
                   __##tag##_end_time - __##tag##_start_time)                  \
                   .count()                                                    \
            << "us" << std::endl;


template<typename T, size_t N>
class MultiDimensionalVector {
public:
    MultiDimensionalVector(const std::vector<size_t>& dimensions, const T& default_value = T())
        : dimensions_(dimensions), data_(calculateTotalSize(dimensions), default_value) {}

    // 获取多维向量中的元素
    T get(const std::vector<size_t>& indices) const {
        assert(indices.size() == N);
        size_t index = calculateIndex(indices);
        return data_[index];
    }

    // 设置多维向量中的元素
    void set(const std::vector<size_t>& indices, const T& value) {
        assert(indices.size() == N);
        size_t index = calculateIndex(indices);
        data_[index] = value;
    }

    std::vector<size_t> dimensions() const {
        return dimensions_;
    }

private:
    size_t calculateTotalSize(const std::vector<size_t>& dimensions) {
        size_t totalSize = 1;
        for (size_t dim : dimensions) {
            totalSize *= dim;
        }
        return totalSize;
    }

    size_t calculateIndex(const std::vector<size_t>& indices) const {
        size_t index = 0;
        size_t multiplier = 1;
        for (int i = N - 1; i >= 0; --i) {
            assert(indices[i] < dimensions_[i]);
            index += indices[i] * multiplier;
            multiplier *= dimensions_[i];
        }
        return index;
    }

    std::vector<size_t> dimensions_;
    std::vector<T> data_;
};

struct TOI_CFAR_element
{
    float rIdx;
    float vIdx;
    uint16_t Beamcnt;
    float pVal;
    float rSNR;
    float vSNR;
};

struct TOI_ConRegion_element
{
    float rIdx;
    float vIdx;
    float pVal;
    uint16_t tarcnt;
};

void Func_Beam_Union_2(vector<TOI_ConRegion_element>& TOI_ConRegion, uint16_t tarcnt, TOI_ConRegion_element& single_tar_Info, uint16_t CompareIdx, uint16_t MaxValIdx);
void func_getConRegion(vector<TOI_ConRegion_element>& TOI_ConRegion, uint16_t TarNum_union, vector<TOI_CFAR_element>& TOI_CFAR, MultiDimensionalVector<float, 3>& FFT2D_abs, ConReg& con_reg, uint16_t CoarseRangeNum, uint16_t VelocityNum);
void Func_Beam_Union(vector<TOI_CFAR_element>& TOI_CFAR, uint16_t Beamcnt, TOI_CFAR_element& toiEle, Union& union_para);
void Spatial_2D_det_V2(vector<uint16_t>& SpatialIdx, vector<float>& SpatialSNR, vector<float>& SpatialPow, vector<float>& FFT3D_abs_single_tar, float Vsnr_db, SpatialPara& spatial_para);
uint8_t func_DetermineEachElement(uint8_t LocComp_u5, uint8_t LeftLoc_u1, float Cell_ToCompare_u30, vector<float>& DataSet_RefCell_u30, uint16_t RefCellNum_u6);
float func_ObtainOSLocData(uint8_t Loc_OSCFAR_u5, vector<float>& DataSet_RefCell_u30, uint16_t RefCellNum_u6);
void func_CFARChM_OS_1D_V(uint8_t& IsTarget_1D_V_u1, float& VSNR_u11, DetPara& det_para, float DataToDetect, uint16_t Index_V, uint16_t Index_R, MultiDimensionalVector<float, 2>& SpatialFFTVelSel_VeloNum_RangeNum);
void func_CFARChM_OS_1D(uint8_t& IsTarget_1D_R_u1, float& RSNR_u11, DetPara& det_para, float DataToDetect, uint16_t Index_V, uint16_t Index_R, MultiDimensionalVector<float, 2>& SpatialFFTVelSel_VeloNum_RangeNum);
void func_PeakSearch_And_CFAR_2D_Cross(uint16_t& TarNum_Detected, vector<uint16_t>& peak_R, vector<uint16_t>& peak_V, vector<float>& peak_Val, vector<vector<float>>& peak_SNR, DetPara& det_para, MultiDimensionalVector<float, 2>& SpatialFFTVelSel_VeloNum_RangeNum);
void cfar3d_cal_across_ArbitaryDim(MultiDimensionalVector<float, 2>& SpatialFFTVelSel_VeloNum_RangeNum, MultiDimensionalVector<float, 3>& DataInput, uint8_t SqueezeDim, uint8_t cfar_include_order, uint8_t cfar_exclude_order, float min_snr, uint8_t Switch3DMode, uint16_t horNum, uint16_t chirpNum, uint16_t sampleNum);
void func_winA_process(vector<vector<vector<vector<vector<complex<float>>>>>>& winAout_xNum_yNum_VeloFFTNum_RangeNum_waveLocNum, vector<uint8_t>& fft_win_type, uint16_t win_len, uint16_t VirtArrVertGridLen, uint16_t VelocityNum, uint16_t RangeSampleNum, uint16_t waveLocNum, vector<vector<vector<vector<vector<complex<float>>>>>>& result_xNum_yNum_VeloFFTNum_RangeNum);
void func_Spatial_Reorder(vector<vector<vector<vector<vector<complex<float>>>>>>& result_xNum_yNum_VeloFFTNum_RangeNum, vector<vector<vector<vector<complex<float>>>>>& FFT2D_MIMONum_VeloFFTNum_CoarseRangeBinNum_WaveLocNum, uint8_t Array_option, uint16_t MIMONum, vector<vector<uint16_t>>& pos_in_mat, uint16_t VirtArrHorGridLen, uint16_t VirtArrVertGridLen, uint16_t VelocityNum, uint16_t CoarseRangeNum, uint16_t waveLocNum);
void func_winD_process(vector<vector<vector<vector<complex<float>>>>>& WinDout_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum, vector<uint8_t>& fft_win_type, uint16_t win_len, uint16_t rangeSampleNum, uint16_t MIMONum, uint16_t waveLocNum, vector<vector<vector<vector<complex<float>>>>>& CoarseRangeFFT_ChirpNum_RangeSampleNum_MIMONum_WaveLocNum);
void findDuplicates(vector<uint16_t>& output, const std::vector<uint16_t>& nums1, const std::vector<uint16_t>& nums2);
void rectWin(vector<float>& win_rect, uint16_t N);
void hanningWin(vector<float>& win_hanning, uint16_t N);
void hammingWin(vector<float>& win_hamming, uint16_t N);
void func_winR_process(vector<vector<vector<float>>>& WinRout_RangeSampleNum_ChirpNum_RxNum, vector<uint8_t>& fft_win_type, uint16_t win_len, uint16_t chirpNum, uint16_t rxNum, unique_ptr<int16_t[]>& radarInputdata);
void FFTD_SpatialFFT_CFAR_CoarseFrame(vector<TOI_CFAR_element>& TOI_CFAR, vector<TOI_ConRegion_element>& TOI_ConRegion, ParaSys& para_sys, Virtual_array& virtual_array, vector<complex<float>>& compensate_mat, vector<vector<vector<complex<float>>>>& CoarseRangeFFT_ValidCoarseRangeBinNum_ChirpNum_RxNum);
void func_signal_process_coarse(vector<TOI_CFAR_element>& TOI_CFAR, vector<TOI_ConRegion_element>& TOI_ConRegion, const string& filename, ParaSys& para_sys, Virtual_array& virtual_array, unique_ptr<int16_t[]>& radarInputdata, vector<complex<float>>& compensate_mat);
