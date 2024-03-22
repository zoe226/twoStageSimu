%{
    two stage simulation platform datafile format
    add head information into adc raw data without head
    not modify yet
%}
clc
clear
close all

fileList = dir('./coarseFrameHead_1Tx/*.mat');
for i = 1:numel(fileList)
    currentFile = fileList(i).name;
    load(['./coarseFrameHead_1Tx/',currentFile]);
end 
compensate_mat_IQ = zeros(1*48*2,1);
compensate_mat_IQ(1:2:end-1) = real(compensate_mat);
compensate_mat_IQ(2:2:end) = imag(compensate_mat);

filename = 'ADCData_Frame_000000_20240314_17_22_19_2876.bin';
fid = fopen(filename,'r');
data = fread(fid,'int16');
head_len = 256;
sample_num = 256;
tx_repeat_num = 768;
rx_num = 48;
tx_num = 1;
chirp_num = tx_repeat_num * tx_num;
adc_data_raw = data((head_len + ceil(chirp_num/1024)*1024)/2+1:end);
% adc_data = reshape(adc_data_raw,sample_num,rx_num,chirp_num);

Frame_index = 0;
Running_time = 0;
cfar_method = 2;
Frame_name = char(repmat('0',1,251));
InputHasDonePreProcess = 0;
Vel_target_generated_by_moniqi = 0;
RangeSampleNum = 256;
TxReuseNum = 768;
TxNum = 1;
RxNum = 48;
% all_tx_seq_pos[48 * 48] = { 0 };
% delta_hop_freq_fstart_seq[48 * 48] = { 0 };
ChirpPRI = 10.24e-6 + 5.84e-6;
hop_freq_low_boundry = 76.1e9;
% Virtual_array_pos_Hor[48 * 48] = { 0 };
% Virtual_array_pos_Vert[48 * 48] = { 0 };
% Virtual_array_pos_Hor_in_mat[48 * 48] = { 0 };
% Virtual_array_pos_Vert_in_mat[48 * 48] = { 0 };
Array_option = 3;
MinElementSpaceHorRelToLamda = 1;
MinElementSpaceVertRelToLamda = 18;
analysis_step = 3;
fft_win_type = [1,1,1];
CoarseRangeNum = 128;
% CoarseRangeAxis[256] = { 0 };
FineRangeNum = 1;
FineRangeAxis = [0];
VelocityNum = 768;
% VelocityAxis[2304] = { 0 };
AngleHorNum = 128;
% AngleHorAxis[128] = { 0 };
AngleVertNum = 32;
% AngleVertAxis[32] = { 0 };
ValidCoarseRangeBin_StartIndex = 1;
ValidCoarseRangeBinNum = 128;
% compensate_mat[48 * 48] = { 0 };
snr_dB_different_dim = [-1,17,5,5];
ignore_static_obj = 0;
static_vel_thres = 0.2;
blind_zone = 2;
cfar_include_order = [5,99,7,7];
cfar_exclude_order = [3,9,3,3];
save_process_data = 0;
algorithm_selection = 1;
% RadarInputData[256 * 48 * 48 * 48] = { 0 };
shorthead = char(repmat('1',1,3328));

global filenameOut;
filenameOut = 'ADCData_Frame_000000_20240314_17_22_19_2876_ShareMemory.bin';
fid = fopen(filenameOut,'wb');
fclose(fid);
startAddress = 0;
appendToBinaryFile(Frame_index,'uint32',startAddress);
startAddress = 4;
appendToBinaryFile(Running_time,'float',startAddress);
startAddress = 8;
appendToBinaryFile(cfar_method,'uint8',startAddress);
startAddress = 9;
appendToBinaryFile(Frame_name,'char',startAddress);
startAddress = 260;
appendToBinaryFile(InputHasDonePreProcess,'uint8',startAddress);
startAddress = 261;
appendToBinaryFile(Vel_target_generated_by_moniqi,'uint8',startAddress);
startAddress = 262;
appendToBinaryFile(RangeSampleNum,'uint16',startAddress);
startAddress = 264;
appendToBinaryFile(TxReuseNum,'uint16',startAddress);
startAddress = 266;
appendToBinaryFile(TxNum,'uint16',startAddress);
startAddress = 268;
appendToBinaryFile(RxNum,'uint16',startAddress);
startAddress = 270;
appendToBinaryFile(all_tx_seq_pos,'uint16',startAddress);
startAddress = 270 + TxReuseNum*TxNum*2;  % 此处TxNum应该是TxGroupNum 单发这两个值相等
appendToBinaryFile(delta_hop_freq_fstart_seq,'float32',startAddress);
startAddress = 270 + TxReuseNum*TxNum*6;  % 此处TxNum应该式TxGroupNum
appendToBinaryFile(ChirpPRI,'float32',startAddress);
startAddress = 274 + TxReuseNum*TxNum*6;
appendToBinaryFile(hop_freq_low_boundry,'float32',startAddress);
startAddress = 278 + TxReuseNum*TxNum*6; % 此处TxNum应该式TxGroupNum
appendToBinaryFile(Virtual_array_pos_Hor,'float32',startAddress);
startAddress = 278 + TxReuseNum*TxNum*6 + TxNum*RxNum*4;
appendToBinaryFile(Virtual_array_pos_Vert,'float32',startAddress);
startAddress = 278 + TxReuseNum*TxNum*6 + TxNum*RxNum*8;
appendToBinaryFile(Virtual_array_pos_Hor_in_mat,'uint16',startAddress);
startAddress = 278 + TxReuseNum*TxNum*6 + TxNum*RxNum*10;
appendToBinaryFile(Virtual_array_pos_Vert_in_mat,'uint16',startAddress);
startAddress = 278 + TxReuseNum*TxNum*6 + TxNum*RxNum*12;
appendToBinaryFile(Array_option,'uint8',startAddress);
startAddress = 279 + TxReuseNum*TxNum*6 + TxNum*RxNum*12;
appendToBinaryFile(MinElementSpaceHorRelToLamda,'float32',startAddress);
startAddress = 283 + TxReuseNum*TxNum*6 + TxNum*RxNum*12;
appendToBinaryFile(MinElementSpaceVertRelToLamda,'float32',startAddress);
startAddress = 287 + TxReuseNum*TxNum*6 + TxNum*RxNum*12;
appendToBinaryFile(analysis_step,'uint8',startAddress);
startAddress = 288 + TxReuseNum*TxNum*6 + TxNum*RxNum*12;
appendToBinaryFile(fft_win_type,'uint8',startAddress);
startAddress = 291 + TxReuseNum*TxNum*6 + TxNum*RxNum*12;
appendToBinaryFile(CoarseRangeNum,'uint16',startAddress);
startAddress = 293 + TxReuseNum*TxNum*6 + TxNum*RxNum*12;
appendToBinaryFile(CoarseRangeAxis,'float32',startAddress);
startAddress = 293 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4;
appendToBinaryFile(FineRangeNum,'uint16',startAddress);
startAddress = 295 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4;
appendToBinaryFile(FineRangeAxis,'float32',startAddress);
startAddress = 295 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4 + FineRangeNum*4;
appendToBinaryFile(VelocityNum,'uint16',startAddress);
startAddress = 297 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4 + FineRangeNum*4;
appendToBinaryFile(VelocityAxis,'float32',startAddress);
startAddress = 297 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4;
appendToBinaryFile(AngleHorNum,'uint16',startAddress);
startAddress = 299 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4;
appendToBinaryFile(AngleHorAxis,'float32',startAddress);
startAddress = 299 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4;
appendToBinaryFile(AngleVertNum,'uint16',startAddress);
startAddress = 301 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4;
appendToBinaryFile(AngleVertAxis,'float32',startAddress);
startAddress = 301 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(ValidCoarseRangeBin_StartIndex,'uint16',startAddress);
startAddress = 303 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(ValidCoarseRangeBinNum,'uint16',startAddress);
startAddress = 305 + TxReuseNum*TxNum*6 + TxNum*RxNum*12 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(compensate_mat_IQ,'float32',startAddress);
startAddress = 305 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(snr_dB_different_dim,'float32',startAddress);
startAddress = 321 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(ignore_static_obj,'uint8',startAddress);
startAddress = 322 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(static_vel_thres,'float32',startAddress);
startAddress = 326 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(blind_zone,'float32',startAddress);
startAddress = 330 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(cfar_include_order,'uint8',startAddress);
startAddress = 334 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(cfar_exclude_order,'uint8',startAddress);
startAddress = 338 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(save_process_data,'uint8',startAddress);
startAddress = 339 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(algorithm_selection,'uint8',startAddress);
startAddress = 340 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(shorthead,'char',startAddress);
startAddress = 3328+340 + TxReuseNum*TxNum*6 + TxNum*RxNum*20 + CoarseRangeNum*4 + FineRangeNum*4 + VelocityNum*4 + AngleHorNum*4 + AngleVertNum*4;
appendToBinaryFile(adc_data_raw,'int16',startAddress);

fid = fopen(filenameOut,'r');
data_file = fread(fid,'uint8');
fclose(fid);

function appendToBinaryFile(data,dataType,startAddress)
global filenameOut
if isempty(filenameOut)
    error('文件名未设置');
end
fid = fopen(filenameOut,'r+');
if fid == -1
    error('无法打开文件')
end
fseek(fid,startAddress,'bof');
fwrite(fid,data,dataType);
fclose(fid);
end
