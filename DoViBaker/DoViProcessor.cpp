#include <algorithm>
#include <array>
#include <string>

#include "DoViProcessor.h"


DoViProcessor::DoViProcessor(const char* rpuPath, IScriptEnvironment* env)
	: successfulCreation(false), rgbProof(false), nlqProof(false), desiredTrimPq(0), max_content_light_level(1000)
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

	successfulCreation = true;
}

DoViProcessor::~DoViProcessor()
{
	dovi_rpu_list_free(rpus);
	::FreeLibrary(doviLib);
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

	const DoviRpuDataMapping* mapping_data = dovi_rpu_get_data_mapping(rpu);
	if (!mapping_data) {
		const char* error = dovi_rpu_get_error(rpu);
		showMessage((std::string("DoViBaker: ") + error).c_str(), env);
		return false;
	}

	out_bit_depth = header->vdr_bit_depth_minus8 + 8;
	bl_bit_depth = header->bl_bit_depth_minus8 + 8;
	el_bit_depth = header->el_bit_depth_minus8 + 8;
	coeff_log2_denom = header->coefficient_log2_denom;
	disable_residual_flag = header->disable_residual_flag;

	for (int cmp = 0; cmp < 3; cmp++) {
        const DoviReshapingCurve curve = mapping_data->curves[cmp];
		num_pivots_minus1[cmp] = curve.num_pivots_minus2 + 1;

		pivot_value[cmp] = std::vector<uint16_t>(curve.pivots.len);
		pivot_value[cmp][0] = curve.pivots.data[0];

		for (int pivot_idx = 1; pivot_idx < curve.pivots.len; pivot_idx++) {
			pivot_value[cmp][pivot_idx] = pivot_value[cmp][pivot_idx - 1] + curve.pivots.data[pivot_idx];
		}

		mapping_idc[cmp] = std::vector<uint8_t>(num_pivots_minus1[cmp]);
		poly_order[cmp] = std::vector<uint8_t>(num_pivots_minus1[cmp]);
		fp_poly_coef[cmp] = std::vector<std::vector<int32_t>>(num_pivots_minus1[cmp]);

		mmr_order[cmp] = std::vector<uint8_t>(num_pivots_minus1[cmp]);
		fp_mmr_const[cmp] = std::vector<int32_t>(num_pivots_minus1[cmp]);
		fp_mmr_coef[cmp] = std::vector<std::vector<std::vector<int32_t>>>(num_pivots_minus1[cmp]);

		for (int pivot_idx = 0; pivot_idx < num_pivots_minus1[cmp]; pivot_idx++) {
			mapping_idc[cmp][pivot_idx] = curve.mapping_idc;

			if (curve.polynomial) {
				const DoviPolynomialCurve *poly_curve = curve.polynomial;
				
				auto poly_order_minus1 = poly_curve->poly_order_minus1;
				auto poly_coef_int = poly_curve->poly_coef_int;
				auto poly_coef = poly_curve->poly_coef;

				poly_order[cmp][pivot_idx] = poly_order_minus1.data[pivot_idx] + 1; 
				fp_poly_coef[cmp][pivot_idx] = std::vector<int32_t>(poly_order[cmp][pivot_idx] + 1); // an order n equation has n+1 coefficients, thus +1!
				for (int coeff = 0; coeff < poly_order[cmp][pivot_idx] + 1; coeff++) {  // an order n equation has n+1 coefficients, thus +1!
					auto port_int = poly_coef_int.list[pivot_idx]->data[coeff];
					auto port_frac = poly_coef.list[pivot_idx]->data[coeff];
					fp_poly_coef[cmp][pivot_idx][coeff] = (port_int << coeff_log2_denom) + port_frac;
				}
			} else if (curve.mmr) {
				const DoviMMRCurve *mmr_curve = curve.mmr;

				auto mmr_order_minus1 = mmr_curve->mmr_order_minus1;
				auto mmr_constant_int = mmr_curve->mmr_constant_int;
				auto mmr_constant = mmr_curve->mmr_constant;
				auto mmr_coef_int = mmr_curve->mmr_coef_int;
				auto mmr_coef = mmr_curve->mmr_coef;

				mmr_order[cmp][pivot_idx] = mmr_order_minus1.data[pivot_idx] + 1;
				auto constant_int = mmr_constant_int.data[pivot_idx];
				auto constant = mmr_constant.data[pivot_idx];
				fp_mmr_const[cmp][pivot_idx] = (constant_int << coeff_log2_denom) + constant;
				fp_mmr_coef[cmp][pivot_idx] = std::vector<std::vector<int32_t>>(mmr_order[cmp][pivot_idx]);

				for (int i = 0; i < mmr_order[cmp][pivot_idx]; i++) {
					fp_mmr_coef[cmp][pivot_idx][i] = std::vector<int32_t>(7);
					for (int j = 0; j < 7; j++) {
						auto port_int = mmr_coef_int.list[pivot_idx]->list[i]->data[j];
						auto port_frac = mmr_coef.list[pivot_idx]->list[i]->data[j];
						fp_mmr_coef[cmp][pivot_idx][i][j] = (port_int << coeff_log2_denom) + port_frac;
					}
				}
			}
		}
	}

	if (header->vdr_dm_metadata_present_flag) {
		const DoviVdrDmData* vdr_dm_data = dovi_rpu_get_vdr_dm_data(rpu);
		if (!vdr_dm_data) {
			const char* error = dovi_rpu_get_error(rpu);
			showMessage((std::string("DoViBaker: ") + error).c_str(), env);
			return false;
		}

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

		if (rgbProof) {
			ycc_to_rgb_coef[0] *= 2;
		}

		scene_refresh_flag = vdr_dm_data->scene_refresh_flag;
		signal_full_range_flag = vdr_dm_data->signal_full_range_flag;
		/*if (!signal_full_range_flag) {
			showMessage("DoViBaker: Limited range output signals have not been tested, no idea if that works.", env);
			return false;
		}*/

		min_pq = vdr_dm_data->dm_data.level1->min_pq;
		max_pq = vdr_dm_data->dm_data.level1->max_pq;
		//max_content_light_level = pq2nits(vdr_dm_data->source_max_pq);
		//max_content_light_level = vdr_dm_data->dm_data.level6->max_content_light_level;
		max_content_light_level = pq2nits(max_pq);

		skipTrim = true;
		if (desiredTrimPq) {
			skipTrim = false;
			avg_pq = vdr_dm_data->dm_data.level1->avg_pq;
			auto lvl2 = vdr_dm_data->dm_data.level2;
			availableTrimPqs = std::vector<uint16_t>(lvl2.len);
			trimInfoMissing = true;
			for (int i = 0; i < lvl2.len; i++) {
				availableTrimPqs[i] = lvl2.list[i]->target_max_pq;
				if (desiredTrimPq != lvl2.list[i]->target_max_pq) continue;
				trimInfoMissing = false;
				trim.slope = lvl2.list[i]->trim_slope;
				trim.offset = lvl2.list[i]->trim_offset;
				trim.power = lvl2.list[i]->trim_power;
				trim.chroma_weight = lvl2.list[i]->trim_chroma_weight;
				trim.saturation_gain = lvl2.list[i]->trim_saturation_gain;
				trim.tone_detail = lvl2.list[i]->ms_weight;
			}
			prepareTrimCoef();
		}

		dovi_rpu_free_vdr_dm_data(vdr_dm_data);
	}
	else {
		showMessage("DoViBaker: No DM Metadata avilable.", env);
		return false;
	}

	if (header->guessed_profile != 7) {
		dovi_rpu_free_data_mapping(mapping_data);
		dovi_rpu_free_header(header);
		return successfulCreation;
	}

	const DoviRpuDataNlq* nlq_data = mapping_data->nlq;
	if (!nlq_data) {
		const char* error = dovi_rpu_get_error(rpu);
		showMessage((std::string("DoViBaker: ") + error).c_str(), env);
		return false;
	}

	if (mapping_data->nlq_method_idc != 0) {
		//https://ffmpeg.org/doxygen/trunk/dovi__rpu_8c_source.html
		showMessage("DoViBaker: Only method NLQ_LINEAR_DZ can be applied, NLQ_MU_LAW is not documented.", env);
		return false;
		//alternativlely we could just gracefully disable the nlq processing with disable_residual_flag=true
	}
	if (mapping_data->nlq_num_pivots_minus2 != 0) {
		showMessage("DoViBaker: Expecting nlq_num_pivots_minus2 to be 0.", env);
		return false;
		//alternativlely we could just gracefully disable the nlq processing with disable_residual_flag=true
	}

	std::string el_type(header->el_type);
	std::transform(el_type.begin(), el_type.end(), el_type.begin(),
		[](unsigned char c) { return std::toupper(c); });
	is_fel = (el_type.compare("FEL")==0);

	auto nlq_offsets = nlq_data->nlq_offset;
	auto vdr_in_max_int = nlq_data->vdr_in_max_int;
	auto vdr_in_max = nlq_data->vdr_in_max;
	auto linear_deadzone_slope_int = nlq_data->linear_deadzone_slope_int;
	auto linear_deadzone_slope = nlq_data->linear_deadzone_slope;
	auto linear_deadzone_threshold_int = nlq_data->linear_deadzone_threshold_int;
	auto linear_deadzone_threshold = nlq_data->linear_deadzone_threshold;

	for (int cmp = 0; cmp < 3; cmp++) {
		nlq_offset[cmp] = nlq_offsets[cmp];
		fp_hdr_in_max[cmp] = (vdr_in_max_int[cmp] << coeff_log2_denom) + vdr_in_max[cmp];
		fp_linear_deadzone_slope[cmp] = (linear_deadzone_slope_int[cmp] << coeff_log2_denom) + linear_deadzone_slope[cmp];
		fp_linear_deadzone_threshold[cmp] = (linear_deadzone_threshold_int[cmp] << coeff_log2_denom) + linear_deadzone_threshold[cmp];
	}
	if (nlqProof) {
		fp_linear_deadzone_slope[0] *= 4;
	}

	dovi_rpu_free_data_mapping(mapping_data);
	dovi_rpu_free_header(header);
	return successfulCreation;
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
	for (int i = 0; i < mmr_order[cmp][pivot_idx]; i++) {
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

void DoViProcessor::prepareTrimCoef() {
	float x1 = trim.minNits = pq2nits(min_pq);
	float x2 = pq2nits(avg_pq);
	float x3 = trim.maxNits = pq2nits(max_pq);

	float y1 = targetMinNits;
	float y2 = sqrtf(x2 * sqrtf(targetMaxNits * targetMinNits));
	float y3 = targetMaxNits;

	float m[10];
	m[9] = x3 * y3 * (x1 - x2) + x2 * y2 * (x3 - x1) + x1 * y1 * (x2 - x3);
	m[0] = x2 * x3 * (y2 - y3); m[1] = x1 * x3 * (y3 - y1); m[2] = x1 * x2 * (y1 - y2);
	m[3] = x3 * y3 - x2 * y2; m[4] = x1 * y1 - x3 * y3; m[5] = x2 * y2 - x1 * y1;
	m[6] = x3 - x2; m[7] = x1 - x3; m[8] = x2 - x1;

	trim.ccc[0] = (m[0] * y1 + m[1] * y2 + m[2] * y3) / m[9];
	trim.ccc[1] = (m[3] * y1 + m[4] * y2 + m[5] * y3) / m[9];
	trim.ccc[2] = (m[6] * y1 + m[7] * y2 + m[8] * y3) / m[9];

	if (!trimInfoMissing) {
		trim.goP[0] = trim.slope / 4096.0 + 0.5;
		trim.goP[1] = trim.offset / 4096.0 - 0.5;
		trim.goP[2] = trim.power / 4096.0 + 0.5;
		trim.cS[0] = trim.chroma_weight / 4096.0 - 0.5;
		trim.cS[1] = trim.saturation_gain / 4096.0 - 0.5;
	}
}

void DoViProcessor::processTrim(uint16_t& ro, uint16_t& go, uint16_t& bo, const uint16_t& ri, const uint16_t& gi, const uint16_t& bi) const  {
	float dr = pq2nits(ri >> (containerBitDepth - out_bit_depth));
	float dg = pq2nits(gi >> (containerBitDepth - out_bit_depth));
	float db = pq2nits(bi >> (containerBitDepth - out_bit_depth));

	float er = (trim.ccc[0] + dr * trim.ccc[1]) / (1 + dr * trim.ccc[2]);
	float eg = (trim.ccc[0] + dg * trim.ccc[1]) / (1 + dg * trim.ccc[2]);
	float eb = (trim.ccc[0] + db * trim.ccc[1]) / (1 + db * trim.ccc[2]);

	if (trimInfoMissing) {
		ro = nits2pq(er) << (containerBitDepth - out_bit_depth);
		go = nits2pq(eg) << (containerBitDepth - out_bit_depth);
		bo = nits2pq(eb) << (containerBitDepth - out_bit_depth);
	}	else {
		float y3 = targetMaxNits;
		float fr = powf((std::clamp(((er / y3) * trim.goP[0]) + trim.goP[1], 0.0f, 1.0f)), trim.goP[2]) * y3;
		float fg = powf((std::clamp(((eg / y3) * trim.goP[0]) + trim.goP[1], 0.0f, 1.0f)), trim.goP[2]) * y3;
		float fb = powf((std::clamp(((eb / y3) * trim.goP[0]) + trim.goP[1], 0.0f, 1.0f)), trim.goP[2]) * y3;

		float Y = 0.22897 * fr + 0.69174 * fg + fb * 0.07929;
		float gr = fr * powf((1 + trim.cS[0]) * fr / Y, trim.cS[1]);
		float gg = fg * powf((1 + trim.cS[0]) * fg / Y, trim.cS[1]);
		float gb = fb * powf((1 + trim.cS[0]) * fb / Y, trim.cS[1]);

		ro = nits2pq(gr) << (containerBitDepth - out_bit_depth);
		go = nits2pq(gg) << (containerBitDepth - out_bit_depth);
		bo = nits2pq(gb) << (containerBitDepth - out_bit_depth);
	}
}
