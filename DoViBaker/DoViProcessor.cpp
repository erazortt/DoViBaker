#include "DoViProcessor.h"
#include <array>
#include <algorithm>

DoViProcessor::DoViProcessor(const char* rpuPath, IScriptEnvironment* env)
	: max_content_light_level(1000), creationError(true)
{
	ycc_to_rgb_coef[0] = 8192;
	ycc_to_rgb_coef[1] = 0;
	ycc_to_rgb_coef[2] = 12900;
	ycc_to_rgb_coef[3] = 8192;
	ycc_to_rgb_coef[4] = -1534;
	ycc_to_rgb_coef[5] = -3836;
	ycc_to_rgb_coef[6] = 8192;
	ycc_to_rgb_coef[7] = 15201;
	ycc_to_rgb_coef[8] = 0;

	ycc_to_rgb_offset[0] = 0;
	ycc_to_rgb_offset[1] = (1 << (containerBitDepth - 1)) << ycc_to_rgb_offset_scale_shifts;
	ycc_to_rgb_offset[2] = (1 << (containerBitDepth - 1)) << ycc_to_rgb_offset_scale_shifts;

	doviLib = ::LoadLibrary(L"dovi.dll"); // delayed loading, original name
	if (doviLib == NULL) {
		showMessage("DoViBaker: Cannot load dovi.dll", env);
		return;
	}

	dovi_parse_rpu_bin_file = (f_dovi_parse_rpu_bin_file)GetProcAddress(doviLib, "dovi_parse_rpu_bin_file");
	if (dovi_parse_rpu_bin_file == NULL) {
		showMessage("DoViBaker: Cannot load function dovi_parse_rpu_bin_file", env);
		return;
	}
	dovi_rpu_list_free = (f_dovi_rpu_list_free)GetProcAddress(doviLib, "dovi_rpu_list_free");
	dovi_rpu_get_header = (f_dovi_rpu_get_header)GetProcAddress(doviLib, "dovi_rpu_get_header");
	dovi_rpu_free_header = (f_dovi_rpu_free_header)GetProcAddress(doviLib, "dovi_rpu_free_header");
	dovi_rpu_get_data_nlq = (f_dovi_rpu_get_data_nlq)GetProcAddress(doviLib, "dovi_rpu_get_data_nlq");
	dovi_rpu_free_data_nlq = (f_dovi_rpu_free_data_nlq)GetProcAddress(doviLib, "dovi_rpu_free_data_nlq");
	dovi_rpu_get_vdr_dm_data = (f_dovi_rpu_get_vdr_dm_data)GetProcAddress(doviLib, "dovi_rpu_get_vdr_dm_data");
	dovi_rpu_free_vdr_dm_data = (f_dovi_rpu_free_vdr_dm_data)GetProcAddress(doviLib, "dovi_rpu_free_vdr_dm_data");
	dovi_rpu_get_data_mapping = (f_dovi_rpu_get_data_mapping)GetProcAddress(doviLib, "dovi_rpu_get_data_mapping");
	dovi_rpu_free_data_mapping = (f_dovi_rpu_free_data_mapping)GetProcAddress(doviLib, "dovi_rpu_free_data_mapping");

	rpus = dovi_parse_rpu_bin_file(rpuPath);
	if (rpus->error) {
		showMessage((std::string("DoViBaker: ")+rpus->error).c_str(), env);
		return;
	}

	pivot_value.resize(3);

	mapping_idc.resize(3);
	poly_order.resize(3);
	fp_poly_coef.resize(3);

	mmr_order.resize(3);
	fp_mmr_const.resize(3);
	fp_mmr_coef.resize(3);

	creationError = false;
}

DoViProcessor::~DoViProcessor()
{
	dovi_rpu_list_free(rpus);
}

void DoViProcessor::showMessage(const char* message, IScriptEnvironment* env)
{
	if (env)
		env->ThrowError(message);
	else
		printf(message);
}

bool DoViProcessor::intializeFrame(int frame, IScriptEnvironment* env) {
	DoviRpuOpaque* rpu = rpus->list[frame];
	const DoviRpuDataHeader* header = dovi_rpu_get_header(rpu);
	if (!header) {
		const char* error = dovi_rpu_get_error(rpu);
		showMessage((std::string("DoViBaker: ") + error).c_str(), env);
		return false;
	}

	if (header->guessed_profile != 7) {
		showMessage("DoViBaker: Expecting profile 7 rpu data.", env);
		return false;
	}
	
	std::string subprofile(header->subprofile);
	std::transform(subprofile.begin(), subprofile.end(), subprofile.begin(),
		[](unsigned char c) { return std::toupper(c); });
	is_fel = (subprofile.compare("FEL")==0);

	auto num_pivots_minus2 = header->num_pivots_minus_2;
	auto pred_pivot_value = header->pred_pivot_value;
	for (int cmp = 0; cmp < 3; cmp++) {
		num_pivots_minus1[cmp] = num_pivots_minus2[cmp] + 1;
		pivot_value[cmp] = std::vector<uint16_t>(num_pivots_minus1[cmp] + 1);
		pivot_value[cmp][0] = pred_pivot_value[cmp].data[0];
		for (int pivot_idx = 1; pivot_idx < num_pivots_minus1[cmp] + 1; pivot_idx++) {
			pivot_value[cmp][pivot_idx] = pivot_value[cmp][pivot_idx - 1] + pred_pivot_value[cmp].data[pivot_idx];
		}
	}

	out_bit_depth = header->vdr_bit_depth_minus_8 + 8;
	bl_bit_depth = header->bl_bit_depth_minus8 + 8;
	el_bit_depth = header->el_bit_depth_minus8 + 8;
	coeff_log2_denom = header->coefficient_log2_denom;
	disable_residual_flag = header->disable_residual_flag;

	if (header->nlq_method_idc != 0) {
		//https://ffmpeg.org/doxygen/trunk/dovi__rpu_8c_source.html
		showMessage("DoViBaker: Only method NLQ_LINEAR_DZ can be applied, NLQ_MU_LAW is not documented.", env);
		return false;
		//alternativlely we could just gracefully disable the nlq processing with disable_residual_flag=true
	}
	if (header->nlq_num_pivots_minus2 != 0) {
		showMessage("DoViBaker: Expecting nlq_num_pivots_minus2 to be 0.", env);
		return false;
		//alternativlely we could just gracefully disable the nlq processing with disable_residual_flag=true
	}

	const DoviRpuDataMapping* mapping_data = dovi_rpu_get_data_mapping(rpu);
	auto poly_order_minus1 = mapping_data->poly_order_minus1;
	auto poly_coef_int = mapping_data->poly_coef_int;
	auto poly_coef = mapping_data->poly_coef;
	for (int cmp = 0; cmp < 3; cmp++) {
		mapping_idc[cmp] = std::vector<uint8_t>(num_pivots_minus1[cmp]);
		poly_order[cmp] = std::vector<uint8_t>(num_pivots_minus1[cmp]);
		fp_poly_coef[cmp] = std::vector<std::vector<int32_t>>(num_pivots_minus1[cmp]);
		for (int pivot_idx = 0; pivot_idx < num_pivots_minus1[cmp]; pivot_idx++) {
			mapping_idc[cmp][pivot_idx] = mapping_data->mapping_idc[cmp].data[0];
			if (mapping_idc[cmp][pivot_idx] != 0) continue;
			poly_order[cmp][pivot_idx] = poly_order_minus1[cmp].data[pivot_idx] + 1; 
			fp_poly_coef[cmp][pivot_idx] = std::vector<int32_t>(poly_order[cmp][pivot_idx] + 1); // an order n equation has n+1 coefficients, thus +1!
			for (int coeff = 0; coeff < poly_order[cmp][pivot_idx] + 1; coeff++) {  // an order n equation has n+1 coefficients, thus +1!
				auto port_int = poly_coef_int[cmp].list[pivot_idx]->data[coeff];
				auto port_frac = poly_coef[cmp].list[pivot_idx]->data[coeff];
				fp_poly_coef[cmp][pivot_idx][coeff] = (port_int << coeff_log2_denom) + port_frac;
			}
		}
	}

	auto mmr_order_minus1 = mapping_data->mmr_order_minus1;
	auto mmr_constant_int = mapping_data->mmr_constant_int;
	auto mmr_constant = mapping_data->mmr_constant;
	auto mmr_coef_int = mapping_data->mmr_coef_int;
	auto mmr_coef = mapping_data->mmr_coef;

	for (int cmp = 0; cmp < 3; cmp++) {
		mmr_order[cmp] = std::vector<uint8_t>(num_pivots_minus1[cmp]);
		fp_mmr_const[cmp] = std::vector<int32_t>(num_pivots_minus1[cmp]);
		fp_mmr_coef[cmp] = std::vector<std::vector<std::vector<int32_t>>>(num_pivots_minus1[cmp]);
		for (int pivot_idx = 0; pivot_idx < num_pivots_minus1[cmp]; pivot_idx++) {
			if (mapping_idc[cmp][pivot_idx] != 1) continue;
			mmr_order[cmp][pivot_idx] = mmr_order_minus1[cmp].data[pivot_idx] + 1;
			auto constant_int = mmr_constant_int[cmp].data[pivot_idx];
			auto constant = mmr_constant[cmp].data[pivot_idx];
			fp_mmr_const[cmp][pivot_idx] = (constant_int << coeff_log2_denom) + constant;
			fp_mmr_coef[cmp][pivot_idx] = std::vector<std::vector<int32_t>>(mmr_order[cmp][pivot_idx] + 1); // an order n equation has n+1 coefficients, thus +1!
			for (int i = 0; i < mmr_order[cmp][pivot_idx] + 1; i++) { // an order n equation has n+1 coefficients, thus +1!
				fp_mmr_coef[cmp][pivot_idx][i] = std::vector<int32_t>(7);
				for (int j = 0; j < 7; j++) {
					auto port_int = mmr_coef_int[cmp].list[pivot_idx]->list[i]->data[j];
					auto port_frac = mmr_coef[cmp].list[pivot_idx]->list[i]->data[j];
					fp_mmr_coef[cmp][pivot_idx][i][j] = (port_int << coeff_log2_denom) + port_frac;
				}
			}
		}
	}

	const DoviRpuDataNlq* nlq_data = dovi_rpu_get_data_nlq(rpu);
	auto nlq_offsets = nlq_data->nlq_offset.list[0];
	auto vdr_in_max_int = nlq_data->vdr_in_max_int.list[0];
	auto vdr_in_max = nlq_data->vdr_in_max.list[0];
	auto linear_deadzone_slope_int = nlq_data->linear_deadzone_slope_int.list[0];
	auto linear_deadzone_slope = nlq_data->linear_deadzone_slope.list[0];
	auto linear_deadzone_threshold_int = nlq_data->linear_deadzone_threshold_int.list[0];
	auto linear_deadzone_threshold = nlq_data->linear_deadzone_threshold.list[0];

	for (int cmp = 0; cmp < 3; cmp++) {
		nlq_offset[cmp] = nlq_offsets->data[cmp];
		fp_hdr_in_max[cmp] = (vdr_in_max_int->data[cmp] << coeff_log2_denom) + vdr_in_max->data[cmp];
		fp_linear_deadzone_slope[cmp] = (linear_deadzone_slope_int->data[cmp] << coeff_log2_denom) + linear_deadzone_slope->data[cmp];
		fp_linear_deadzone_threshold[cmp] = (linear_deadzone_threshold_int->data[cmp] << coeff_log2_denom) + linear_deadzone_threshold->data[cmp];
	}

	if (header->vdr_dm_metadata_present_flag) {
		const DoviVdrDmData* vdr_dm_data = dovi_rpu_get_vdr_dm_data(rpu);

		max_pq = vdr_dm_data->dm_data.level1->max_pq;
		//max_content_light_level = pq2nits(vdr_dm_data->source_max_pq);

		//max_content_light_level = vdr_dm_data->dm_data.level6->max_content_light_level;
		max_content_light_level = pq2nits(max_pq);

		ycc_to_rgb_coef[0] = vdr_dm_data->ycc_to_rgb_coef0;
		ycc_to_rgb_coef[1] = vdr_dm_data->ycc_to_rgb_coef1;
		ycc_to_rgb_coef[2] = vdr_dm_data->ycc_to_rgb_coef2;
		ycc_to_rgb_coef[3] = vdr_dm_data->ycc_to_rgb_coef3;
		ycc_to_rgb_coef[4] = vdr_dm_data->ycc_to_rgb_coef4;
		ycc_to_rgb_coef[5] = vdr_dm_data->ycc_to_rgb_coef5;
		ycc_to_rgb_coef[6] = vdr_dm_data->ycc_to_rgb_coef6;
		ycc_to_rgb_coef[7] = vdr_dm_data->ycc_to_rgb_coef7;
		ycc_to_rgb_coef[8] = vdr_dm_data->ycc_to_rgb_coef8;

		ycc_to_rgb_offset[0] = vdr_dm_data->ycc_to_rgb_offset0 >> ycc_to_rgb_offset_scale_shifts;
		ycc_to_rgb_offset[1] = vdr_dm_data->ycc_to_rgb_offset1 >> ycc_to_rgb_offset_scale_shifts;
		ycc_to_rgb_offset[2] = vdr_dm_data->ycc_to_rgb_offset2 >> ycc_to_rgb_offset_scale_shifts;

		dovi_rpu_free_vdr_dm_data(vdr_dm_data);
	}

	dovi_rpu_free_data_mapping(mapping_data);
	dovi_rpu_free_data_nlq(nlq_data);
	dovi_rpu_free_header(header);
	return true;
}

uint16_t DoViProcessor::processSample(int cmp, uint16_t bl, uint16_t el, uint16_t mmrBlY, uint16_t mmrBlU, uint16_t mmrBlV) const {
	bl >>= (containerBitDepth - bl_bit_depth);
	int pivot_idx = getPivotIndex(cmp, bl);
	int v;
	if (cmp == 0 || mapping_idc[cmp][pivot_idx] == 0) {
		v = polynompialMapping(cmp, pivot_idx, bl);
	}
	else {
		mmrBlY >>= (containerBitDepth - bl_bit_depth);
		mmrBlU >>= (containerBitDepth - bl_bit_depth);
		mmrBlV >>= (containerBitDepth - bl_bit_depth);
		v = mmrMapping(cmp, pivot_idx, mmrBlY, mmrBlU, mmrBlV);
	}
	int r = 0;
	if (!disable_residual_flag) {
		el >>= (containerBitDepth - el_bit_depth);
		r = nonLinearInverseQuantization(cmp, el);
	}
	uint16_t h = signalReconstruction(v, r);
	h <<= (containerBitDepth - out_bit_depth);
	return h;
}

int DoViProcessor::getPivotIndex(int cmp, uint16_t s) const {
	int pivot_idx = num_pivots_minus1[cmp];
	for (int idx = 0; idx < num_pivots_minus1[cmp]; idx++) {
		if (s < pivot_value[cmp][idx + 1]) {
			pivot_idx = idx;
			break;
		}
	}
	return pivot_idx;
}

uint16_t DoViProcessor::polynompialMapping(int cmp, int pivot_idx, uint16_t s) const {
	if (s < pivot_value[cmp][0])
		s = pivot_value[cmp][0];
	if (s > pivot_value[cmp][num_pivots_minus1[cmp]])
		s = pivot_value[cmp][num_pivots_minus1[cmp]];
	// compute polynom at s in fixed point arithmetic
	int64_t ss = 1;
	int64_t shift = 20; // 2*(maximum BL_bit_depth)
	int64_t vv = 0;
	for (int i = 0; i <= poly_order[cmp][pivot_idx]; i++)
	{
		vv += fp_poly_coef[cmp][pivot_idx][i] * (ss << shift);
		ss *= s;
		shift -= bl_bit_depth;
	}
	vv = (vv < 0) ? 0 : vv;
	int64_t v = vv >> (4 + coeff_log2_denom);
	v = (v > 0xffff) ? 0xffff : v;
	return v;
}

uint16_t DoViProcessor::mmrMapping(int cmp, int pivot_idx, uint16_t s0, uint16_t s1, uint16_t s2) const {
	if (s0 < pivot_value[0][0])
		s0 = pivot_value[0][0];
	if (s0 > pivot_value[0][num_pivots_minus1[0]])
		s0 = pivot_value[0][num_pivots_minus1[0]];
	if (s1 < pivot_value[1][0])
		s1 = pivot_value[1][0];
	if (s1 > pivot_value[1][num_pivots_minus1[1]])
		s1 = pivot_value[1][num_pivots_minus1[1]];
	if (s2 < pivot_value[2][0])
		s2 = pivot_value[2][0];
	if (s2 > pivot_value[2][num_pivots_minus1[2]])
		s2 = pivot_value[2][num_pivots_minus1[2]];
	// constant
	int64_t tt[22];
	tt[0] = 1 << 20;
	//num_coeff = 1;
	// first order
	if (mmr_order[cmp][pivot_idx] >= 1) {
		tt[1] = s0 << (20 - bl_bit_depth);
		tt[2] = s1 << (20 - bl_bit_depth);
		tt[3] = s2 << (20 - bl_bit_depth);
		tt[4] = (s0 * s1) << (20 - 2 * bl_bit_depth);
		tt[5] = (s0 * s2) << (20 - 2 * bl_bit_depth);
		tt[6] = (s1 * s2) << (20 - 2 * bl_bit_depth);
		tt[7] = (tt[4] * tt[3]) >> 20;
	}
	// second order
	if (mmr_order[cmp][pivot_idx] >= 2) {
		tt[8] = (s0 * s0) << (20 - 2 * bl_bit_depth);
		tt[9] = (s1 * s1) << (20 - 2 * bl_bit_depth);
		tt[10] = (s2 * s2) << (20 - 2 * bl_bit_depth);
		tt[11] = (tt[4] * tt[4]) >> 20;
		tt[12] = (tt[5] * tt[5]) >> 20;
		tt[13] = (tt[6] * tt[6]) >> 20;
		tt[14] = (tt[7] * tt[7]) >> 20;
	}
	// third order
	if (mmr_order[cmp][pivot_idx] >= 3) {
		tt[15] = (tt[1] * tt[8]) >> 20;
		tt[16] = (tt[2] * tt[9]) >> 20;
		tt[17] = (tt[3] * tt[10]) >> 20;
		tt[18] = (tt[4] * tt[11]) >> 20;
		tt[19] = (tt[5] * tt[12]) >> 20;
		tt[20] = (tt[6] * tt[13]) >> 20;
		tt[21] = (tt[7] * tt[14]) >> 20;
	}
	int64_t rr = fp_mmr_const[cmp][pivot_idx] * tt[0];
	int cnt = 1;
	for (int i = 1; i <= mmr_order[cmp][pivot_idx]; i++) {
		for (int j = 0; j < 7; j++) {
			rr += fp_mmr_coef[cmp][pivot_idx][i][j] * tt[cnt];
			cnt++;
		}
	}
	rr = rr < 0 ? 0 : rr;
	int64_t v = (rr >> (4 + coeff_log2_denom));
	v = v > 0xffff ? 0xffff : v;
	return v;
}

int16_t DoViProcessor::nonLinearInverseQuantization(int cmp, uint16_t e) const {
	// coefficients
	int T = fp_linear_deadzone_threshold[cmp];
	int S = fp_linear_deadzone_slope[cmp];
	int R = fp_hdr_in_max[cmp];
	// input data
	int64_t rr = e - nlq_offset[cmp];
	int64_t r;
	if (rr == 0) {
		r = 0;
	}
	else {
		int sign = rr < 0 ? -1 : 1;
		rr <<= 1;
		rr -= sign;
		rr <<= (10 - el_bit_depth);
		// output data
		int64_t dq = rr * S;
		int64_t TT = (T << (10 - el_bit_depth + 1)) * sign;
		dq += TT;
		int64_t RR = (R << (10 - el_bit_depth + 1));
		if (dq > RR)
			dq = RR;
		else if (dq < -RR)
			dq = -RR;
		r = (dq >> (coeff_log2_denom - 5 - el_bit_depth));
	}
	return r;
}

uint16_t DoViProcessor::signalReconstruction(uint16_t v, int16_t r) const {
	int MAXOUT = (1 << out_bit_depth) - 1;
	int h = v;
	if (!disable_residual_flag)
		h += r;
	h += (1 << (15 - out_bit_depth));
	h >>= (16 - out_bit_depth);
	h = h < 0 ? 0 : h;
	h = h > MAXOUT ? MAXOUT : h;
	return h;
}
