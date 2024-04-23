#include "param.hpp"

const uint16_t CoarseFrame_CFARdim = 3;
uint16_t FineFrame_CfarDim = 4;
vector<vector<float>> TOI;
vector<float> R_Threshold1_u9(32, 0.0);
vector<float> V_Threshold1_u9(16, 0.0);
vector<float> PeakThreshold1_u30(32, 0.0);

void  get_virtual_array(Virtual_array& virtual_array, const vector<int16_t>& rx_x_pos_renorm, const vector<int16_t>& rx_y_pos_renorm, const vector<int16_t>& tx_x_pos_renorm, const vector<int16_t>& tx_y_pos_renorm) {
	uint16_t Tx_num = size(tx_x_pos_renorm);
	uint16_t Rx_num = size(rx_x_pos_renorm);
	vector<int16_t> Relative_Tx1_pos_x(Tx_num);
	vector<int16_t> Relative_Tx1_pos_y(Tx_num);
	for (int i = 0; i < Tx_num; i++) {
		Relative_Tx1_pos_x[i] = tx_x_pos_renorm[i] - tx_x_pos_renorm[0];
		Relative_Tx1_pos_y[i] = tx_y_pos_renorm[i] - tx_y_pos_renorm[0];
	}

	uint16_t repsize_rx = Tx_num * Rx_num;
	vector<int16_t> Virtual_Rx_pos_x(repsize_rx);
	vector<int16_t> Virtual_Rx_pos_y(repsize_rx);
	for (int i = 0; i < repsize_rx; i++) {
		Virtual_Rx_pos_x[i] = rx_x_pos_renorm[int(i % 48)] + Relative_Tx1_pos_x[int(i / 48)];
		Virtual_Rx_pos_y[i] = rx_y_pos_renorm[int(i % 48)] + Relative_Tx1_pos_y[int(i / 48)];
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
	virtual_array.y_max_pos = max_Rx_y;
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
	virtual_array.pos_in_mat = pos_in_mat;
}

void getArray(ParaArray& para_array)
{
	// only option = 3;
	uint16_t tx_num = 48;
	uint16_t rx_num = 48;
	float min_element_space_x_relto_semilamda = 1;
	float min_element_space_y_relto_semilamda = 3;
	vector<int16_t> Rx_pos_x{0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52};
	vector<int16_t> Rx_pos_y{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36};
	vector<int16_t> Tx_pos_x{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51};
	vector<int16_t> Tx_pos_y{6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, -3, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, -3, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30};

	para_array.tx_x_pos_renorm = Tx_pos_x;
	para_array.tx_y_pos_renorm = Tx_pos_y;
	para_array.rx_x_pos_renorm = Rx_pos_x;
	para_array.rx_y_pos_renorm = Rx_pos_y;
	para_array.min_element_space_x_relto_semilamda = min_element_space_x_relto_semilamda;
	para_array.min_element_space_y_relto_semilamda = min_element_space_y_relto_semilamda;
}

void parseDataFile(vector<unsigned char>& inputBuffer, BinFile& bin_file) {
	if (inputBuffer.size() < 83886080) {
		std::cerr << "Buffer does not contain enough data to parse." << endl;
	}
	memcpy(&bin_file.para_sys.InputHasDonePreProcess, inputBuffer.data() + 260, sizeof(unsigned char));
	memcpy(&bin_file.para_sys.vel_target_generated_by_moniqi, inputBuffer.data() + 261, sizeof(unsigned char));
	memcpy(&bin_file.para_sys.RangeSampleNum, inputBuffer.data() + 262, 2 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.TxReuseNum, inputBuffer.data() + 264, 2 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.TxNum, inputBuffer.data() + 266, 2 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.RxNum, inputBuffer.data() + 268, 2 * sizeof(unsigned char));
	bin_file.para_sys.all_tx_seq_pos.resize(bin_file.para_sys.TxNum, vector<uint16_t>(bin_file.para_sys.TxReuseNum));
	for (int i = 0; i < bin_file.para_sys.TxNum; i++) {
		bin_file.para_sys.all_tx_seq_pos[i].assign(reinterpret_cast<uint16_t*>(inputBuffer.data() + 270 + i * 2 * bin_file.para_sys.TxReuseNum), reinterpret_cast<uint16_t*>(inputBuffer.data() + 270 + i * 2 * bin_file.para_sys.TxReuseNum + 2 * bin_file.para_sys.TxReuseNum * sizeof(unsigned char)));
	}
	bin_file.para_sys.delta_hop_freq_fstart_seq.resize(bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum);
	bin_file.para_sys.delta_hop_freq_fstart_seq.assign(reinterpret_cast <float*>(inputBuffer.data() + 270 + bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum * 2), reinterpret_cast <float*>(inputBuffer.data() + 270 + bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum * 2 + 4 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum * sizeof(unsigned char)));
	memcpy(&bin_file.para_sys.ChirpPRI, inputBuffer.data() + 270 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 4 * sizeof(unsigned char));
	memcpy(&bin_file.para_sys.hop_freq_low_boundry, inputBuffer.data() + 274 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum, 4 * sizeof(unsigned char));
	// Virtual_array_pos_Hor
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
	bin_file.para_sys.FineRangeAxis.assign(reinterpret_cast<float*>(inputBuffer.data() + 295 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum), reinterpret_cast<float*>(inputBuffer.data() + 295 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum * sizeof(unsigned char)));
	memcpy(&bin_file.para_sys.VelocityNum, inputBuffer.data() + 295 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum, 2 * sizeof(unsigned char));
	bin_file.para_sys.VelocityAxis.resize(bin_file.para_sys.VelocityNum);
	bin_file.para_sys.VelocityAxis.assign(reinterpret_cast<float*>(inputBuffer.data() + 297 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum), reinterpret_cast<float*>(inputBuffer.data() + 297 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * sizeof(unsigned char) * bin_file.para_sys.VelocityNum));
	memcpy(&bin_file.para_sys.AngleHorNum, inputBuffer.data() + 297 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum, 2 * sizeof(unsigned char));
	bin_file.para_sys.AngleHorAxis.resize(bin_file.para_sys.AngleHorNum);
	bin_file.para_sys.AngleHorAxis.assign(reinterpret_cast<float*>(inputBuffer.data() + 299 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum), reinterpret_cast<float*>(inputBuffer.data() + 299 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * sizeof(unsigned char) * bin_file.para_sys.AngleHorNum));
	memcpy(&bin_file.para_sys.AngleVertNum, inputBuffer.data() + 299 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum, 2 * sizeof(unsigned char));
	bin_file.para_sys.AngleVertAxis.resize(bin_file.para_sys.AngleVertNum);
	bin_file.para_sys.AngleVertAxis.assign(reinterpret_cast<float*>(inputBuffer.data() + 301 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum), reinterpret_cast<float*>(inputBuffer.data() + 301 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * sizeof(unsigned char) * bin_file.para_sys.AngleVertNum));
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
	memcpy(bin_file.compensate_mat.get(), inputBuffer.data() + 305 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 12 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, 4 * sizeof(unsigned char) * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum * 2);
	bin_file.input_data = std::make_unique<int16_t[]>(bin_file.para_sys.RxNum * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum * bin_file.para_sys.RangeSampleNum * 2);
	memcpy(bin_file.input_data.get(), inputBuffer.data() + 3328 + 340 + 6 * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum + 20 * bin_file.para_sys.TxNum * bin_file.para_sys.RxNum + 4 * bin_file.para_sys.CoarseRangeNum + 4 * bin_file.para_sys.FineRangeNum + 4 * bin_file.para_sys.VelocityNum + 4 * bin_file.para_sys.AngleHorNum + 4 * bin_file.para_sys.AngleVertNum, 2 * sizeof(unsigned char) * bin_file.para_sys.RxNum * bin_file.para_sys.TxNum * bin_file.para_sys.TxReuseNum * bin_file.para_sys.RangeSampleNum * 2);
}
void getDetPara(DetPara& det_para, ParaSys& para_sys) {
	det_para.con_reg.PowerDetEn = 0;
	det_para.con_reg.SnrThre = 8;
	det_para.con_reg.RextendNum = 5;
	det_para.con_reg.VextendNum = 7;
	if ((FineFrame_CfarDim == 2) || (FineFrame_CfarDim == 3)) {
		det_para.con_reg.RextendNum = 1;
		det_para.con_reg.VextendNum = 1;
	}
	if (para_sys.work_mode == 0) {
		if (para_sys.frame_type == 0) {
			det_para.con_reg.RextendNum = 1;
			det_para.con_reg.VextendNum = 1;
		}
		else if (para_sys.frame_type == 1)
		{
			TOI.resize(para_sys.CoarseRangeNum, vector<float>(2));
			for (int i = 0; i < para_sys.CoarseRangeNum; i++) {
				TOI[i][0] = para_sys.CoarseRangeAxis[i];
				TOI[i][1] = para_sys.VelocityAxis[0];
			}
			para_sys.RangeCell_left = 0;
			para_sys.RangeCell_right = 0;
			para_sys.VelocityCell_left = 0;
			para_sys.VelocityCell_right = para_sys.TxGroupNum * para_sys.TxReuseNum - 1;
		}
	}
	else if (para_sys.work_mode == 1) {
		if (para_sys.frame_type == 1) {
			if ((det_para.con_reg.RextendNum != 1) && (det_para.con_reg.VextendNum != 1)) {
				para_sys.RangeCell_left = 0;
				para_sys.RangeCell_right = 0;
				para_sys.VelocityCell_left = 1;
				para_sys.VelocityCell_right = 1;
				FineFrame_CfarDim = 4;
			}
			else {
				para_sys.RangeCell_left = 3;
				para_sys.RangeCell_right = 3;
				para_sys.VelocityCell_left = 10;
				para_sys.VelocityCell_right = 10;
			}
		}
	}
	det_para.ChToProcess_Num_u7 = 1;
	string str = "11";
	det_para.Switch_DataCompress_u2.assign(str.begin(), str.end());
	det_para.PassAbsData_u1 = 1;
	det_para.PassChMeanData_u1 = 1;
	det_para.PassPeakS_u1 = 1;
	det_para.OneDEnable = 0;
	det_para.Reg_PGA_u15 = 272;
	det_para.PeakS_Enable_u1 = 1;

	if (para_sys.frame_type == 0) {
		det_para.RangeCellNum_u10 = para_sys.CoarseRangeNum;
		det_para.ChirpNum_u11 = para_sys.VelocityNum;
	}
	else {
		det_para.RangeCellNum_u10 = para_sys.FineRangeNum * (para_sys.RangeCell_left + para_sys.RangeCell_right + 1);
		det_para.ChirpNum_u11 = para_sys.VelocityCell_left + para_sys.VelocityCell_right + 1;
	}
	det_para.DetectCell_RIndex_Min_u10 = 2;
	det_para.DetectCell_RIndex_Max_u10 = det_para.RangeCellNum_u10 - 1;
	det_para.DetectCell_VIndex_Min_u11 = 1;
	det_para.DetectCell_VIndex_Max_u11 = det_para.ChirpNum_u11;

	string str_cfar = "11";
	det_para.CFARTypeSwitch_u2 = str_cfar;
	det_para.LogicTestFlag_u1 = 1;
	det_para.cfar_para.ProCellNum_R_u2 = 3;
	det_para.cfar_para.ProCellNum_V_u2 = 3;
	det_para.cfar_para.RefCellNum_1D_u5 = 8;
	det_para.Loc_OSCFAR_u5 = 10;
	det_para.Index_Chirp_NotMove_OSCFAR_u11 = 1;
	det_para.Threshold_RangeDim_For_2D_OSCFAR_u9 = 13;
	det_para.asicsavedata_flag = 0;
	det_para.IndexCHIP_u3 = 1;

	det_para.spatial_para.Vsnr_Thre_db = 1;
	det_para.spatial_para.HorProfThre_db = 16;
	det_para.spatial_para.PeakSen = 1;
	det_para.spatial_para.HorOffsetidx = 6;
	det_para.spatial_para.HorOffPeakThre_db = 8;
	det_para.spatial_para.Global_noise_Thre_db = 20;
	det_para.spatial_para.min_snr_Hor_db = 15;
	det_para.spatial_para.VerThreshold = 16;
	det_para.spatial_para.HorCfarGuardNum = 5;
	det_para.spatial_para.HorCfarReferNum = 25;

	det_para.spatial_para.Horminidx = 1;
	det_para.spatial_para.Hormaxidx = para_sys.AngleHorNum;
	det_para.spatial_para.Verminidx = 1;
	det_para.spatial_para.Vermaxidx = para_sys.AngleVertNum;
	if (para_sys.work_mode == 1) {
		if (para_sys.frame_type == 1) {
			if (det_para.con_reg.RextendNum != 1 && det_para.con_reg.VextendNum != 1) {
				det_para.spatial_para.Horminidx = 1;
				det_para.spatial_para.Hormaxidx = para_sys.AngleHorNum;
				det_para.spatial_para.Verminidx = 1;
				det_para.spatial_para.Vermaxidx = para_sys.AngleVertNum * det_para.RangeCellNum_u10 * det_para.ChirpNum_u11;
			}
		}
	}

	if (para_sys.frame_type == 0) {
		if (para_sys.work_mode == 0) {
			det_para.union_para.CompareIdx = 3;
			det_para.union_para.MaxValIdx = 9;
		}
		else {
			det_para.union_para.CompareIdx = 2;
			det_para.union_para.MaxValIdx = 4;
		}
	}
}

void get_info(std::string filename, BinFile& bin_file, Virtual_array& virtual_array) {
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
	bin_file.para_sys.TxGroup.push_back({ 1 });  // TxGroupNum*txNum_in_group
	bin_file.para_sys.waveLocNum = 1;
	bin_file.para_sys.DopplerBandNum = 1;
	bin_file.para_sys.detect_option = 0;
	bin_file.para_sys.cfar_2d_tf = 1.8;
	bin_file.para_sys.lightSpeed = 299792458.0f;
	bin_file.para_sys.plot_option = 0;
	bin_file.para_sys.work_mode = 0;
	bin_file.para_sys.frame_type = 0;
	bin_file.para_sys.waveform_option = 1;

	parseDataFile(buffer, bin_file);

	uint16_t cols = bin_file.para_sys.TxReuseNum * bin_file.para_sys.TxGroupNum / bin_file.para_sys.waveLocNum;
	bin_file.para_sys.waveLocChirpSeq.resize(bin_file.para_sys.waveLocNum, vector<uint16_t>(cols));
	for (int i = 0; i < bin_file.para_sys.waveLocNum; i++) {
		for (int j = 0; j < cols; j++) {
			bin_file.para_sys.waveLocChirpSeq[i][j] = j + 1;  // waveLocNum * (totalChirpNum/waveLocNum)
		}
	}

	bin_file.para_sys.work_mode = 0;
	bin_file.para_sys.frame_type = 0;
	getDetPara(bin_file.para_sys.det_para, bin_file.para_sys);

	ParaArray para_array;
	getArray(para_array);
	if (bin_file.para_sys.frame_type == 0) {
		// tx pos renorm 
		vector<uint16_t> tx_x_pos_renorm(bin_file.para_sys.txNum_in_group * bin_file.para_sys.TxGroupNum);
		vector<uint16_t> tx_y_pos_renorm(bin_file.para_sys.txNum_in_group * bin_file.para_sys.TxGroupNum);
		uint16_t txRealNum = 0;
		for (int i = 0; i < bin_file.para_sys.TxGroupNum; i++) {
			for (int j = 0; j < bin_file.para_sys.txNum_in_group; j++) {
				tx_x_pos_renorm[txRealNum] = para_array.tx_x_pos_renorm[bin_file.para_sys.TxGroup[i][j] - 1];
				tx_y_pos_renorm[txRealNum] = para_array.tx_x_pos_renorm[bin_file.para_sys.TxGroup[i][j] - 1];
				txRealNum++;
			}
		}
		para_array.tx_x_pos_renorm.assign(tx_x_pos_renorm.begin(), tx_x_pos_renorm.end());
		para_array.tx_y_pos_renorm.assign(tx_y_pos_renorm.begin(), tx_y_pos_renorm.end());
	}
	get_virtual_array(virtual_array, para_array.rx_x_pos_renorm, para_array.rx_y_pos_renorm, para_array.tx_x_pos_renorm, para_array.tx_y_pos_renorm);
	if (bin_file.para_sys.frame_type == 0) {
		bin_file.para_sys.MinElementSpaceHorRelToLamda = virtual_array.x_space / 2;
		bin_file.para_sys.MinElementSpaceVertRelToLamda = virtual_array.y_space / 2;
	}
}