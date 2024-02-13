#include "def.h"
#include <fstream>

double EOTFpq(double ep)
{
  static const double m1 = 2610.0 / 4096 / 4;
  static const double m2 = 2523.0 / 4096 * 128;
  static const double c3 = 2392.0 / 4096 * 32;
  static const double c2 = 2413.0 / 4096 * 32;
  static const double c1 = c3 - c2 + 1;

  const double epower = pow(ep, 1 / m2);
  const double num = std::max(epower - c1, 0.0);
  const double denom = c2 - c3 * epower;
  return pow(num / denom, 1 / m1);
}

double OOTFhlgInv(double fd, double yd)
{
  if (yd > 0) {
    static const double power = (1 - 1.2) / 1.2;
    return fd * pow(yd, power);
  }
  return 0;
}

double OETFhlg(double e)
{
  static const double a = 0.17883277;
  static const double b = 1 - 4 * a;
  static const double c = 0.5 - a * std::log(4 * a);

  if (e * 12 < 1) {
    return sqrt(3 * e);
  }
  return a * std::log(12 * e - b) + c;
}

double matchHlg2Sdr(double x) {
  return 2.8215 * x * x - 2.7758 * x + 1.6826;
}

double matchHlg2SdrD(double x) {
  return 2 * 2.8215 * x - 2.7758;
}

double hlg2sdr(double x, double kS) {
  if (x < 0.5) {
    // hlg and sdr are itentical below 203 nits
    return x;
  }
  if (x < kS) {
    // this function maps the hlg signal above 203 nits to the sdr signal range until it runs out of sdr signal range at 267.6nits hlg
    // thus when the sdr is viewed it appears idential to the hlg as long as the hlg stays below these 267.6 nits
    return x * matchHlg2Sdr(x);
  }
  // since the above function will run out of sdr signal space, we need to moderate it so that it does not clip
  // this flattening is done by a hermite spline (just like in hdr tonemapping)
  // the start point is defined by kS, p0 and m0 and the end point by kE, p1 and m1.
  double p0 = kS * matchHlg2Sdr(kS);
  double p1 = 1;
  double m0 = (kS * matchHlg2SdrD(kS) + 1 * matchHlg2Sdr(kS));
  double m1 = std::pow(1 / m0, 3.6231); //the power was chosen such that the curve has a monotonus falling derivative, in the whole validity range of kS=[0.5, 0.7]
  m0 *= (1 - kS);
  m1 *= (1 - kS);
  double kE = 1;
  float t = (x - kS) / (kE - kS);
  float p = ((2 * t - 3) * t * t + 1) * p0 + (((t - 2) * t + 1) * m0 + (-2 * t + 3) * t) * t * p1 + t * t * (t - 1) * m1;
  return p;
}

void showBestLutSizes() {
  printf("For unnormalized PQ inputs only the following LUT sizes should be used:\n");
  for (int s = 2; s < 139; s++) {
    if (s <= 69 && (s % 4) != 1) continue;
    if (s >= 70 && (s % 4) != 2) continue;
    if (s > 138) continue;
    printf("%i ", s);
  }
  printf("\n");
}

void warnAboutLutSize() {
  printf("\n");
  printf("WARNING:\n");
  showBestLutSizes();
  printf("\n");
}

#ifdef DOVI_LUTGEN
int main(int argc, char* argv[])
{
  if (argc < 4) {
    printf("usage:\n");
    printf("<OUTPUT_CUBE_FILE> <LUT_SIZE> <NORMALIZED_INPUT> <SDR_KNEE_START> \n");
    printf("\n");
    printf("examples for PQ to HLG conversions: \n");
    printf("DoViLutGen.exe pq2hlg.cube 65 0\n");
    printf("DoViLutGen.exe pq2hlg_normalizedInput.cube 50 1\n");
    printf("\n");
    printf("examples for PQ to SDR conversions: \n");
    printf("DoViLutGen.exe pq2sdr.cube 65 0 0.6 \n");
    printf("DoViLutGen.exe pq2sdr_normalizedInput.cube 50 1 0.6 \n");
    printf("\n");
    showBestLutSizes();
    return 1;
  }
  std::ofstream fsCube;
  fsCube.open(argv[1], std::ifstream::out);
  if (!fsCube.is_open()) {
    fsCube.close();
    printf("Unable to open output file\n");
    return 2;
  }

  int lutSize = std::atoi(argv[2]);
  bool normalizedInput = std::atoi(argv[3]);

  float kneeStart = 0;
  if (argc > 4) {
    kneeStart = std::atof(argv[4]);
    if (!(kneeStart > 0.5)) {
      printf("\n");
      printf("Only SDR knee values above 0.5 are allowed.\n");
      printf("To create a HLG LUT remove the SDR knee parameter!\n");
      return 2;
    }
  }
  bool sdr = false;
  std::string mode = "HLG";
  if (kneeStart > 0.5) {
    mode = "SDR";
    sdr = true;
    if (kneeStart > 0.7) {
      printf("The knee start cannot be above 0.7\n");
      return 2;
    }
  }

  if (!normalizedInput) {
    if (lutSize <= 69 && (lutSize % 4) != 1) {
      warnAboutLutSize();
    }
    if (lutSize >= 70 && (lutSize % 4) != 2) {
      warnAboutLutSize();
    }
    if (lutSize > 138) {
      warnAboutLutSize();
    }
  }

  if (sdr) {
    printf("Producing a LUT for PQ -> SDR conversions of size %i with flattending starting at %f\n", lutSize, kneeStart);
    fsCube << "# LUT for conversions from BT.2100 HDR PQ to BT.2020 SDR with flattening starting at " << kneeStart << std::endl;
  }
  else {
    printf("Producing LUT for PQ -> HLG conversions of size %i\n", lutSize);
    fsCube << "# LUT for conversions from BT.2100 HDR PQ to BT.2100 HDR HLG following the BBC specs with 75% reference white" << std::endl;
  }

  double inputScale = 1;
  if (normalizedInput) {
    printf("The LUT will expect the PQ input to be re-normalized to 1000 nits\n");
    inputScale = 0.7518271;
    fsCube << "# ATTENTION: This special LUT expects the PQ input to be re-normalized to 1000nits max brightness!" << std::endl;
    fsCube << "TITLE \"PQ_renorm_1000nits_to_" << mode << "\"" << std::endl;
  }
  else {
    printf("The LUT will expect usual PQ input\n");
    fsCube << "TITLE \"PQ_to_" << mode << "\"" << std::endl;
  }
  fsCube << "LUT_3D_SIZE " << lutSize << std::endl;

  // PQ to HLG conversion based on BT.2408-7 in conjunction with BT.2100-2
  for (int bi = 0; bi < lutSize; bi++) {
    double ebp = double(bi) / (lutSize - 1) * inputScale;
    double bd = EOTFpq(ebp) * 10000 / 1000;
    bd = (bd > 1) ? 1 : bd;
    for (int gi = 0; gi < lutSize; gi++) {
      double egp = double(gi) / (lutSize - 1) * inputScale;
      double gd = EOTFpq(egp) * 10000 / 1000;
      gd = (gd > 1) ? 1 : gd;
      for (int ri = 0; ri < lutSize; ri++) {
        double erp = double(ri) / (lutSize - 1) * inputScale;
        double rd = EOTFpq(erp) * 10000 / 1000;
        rd = (rd > 1) ? 1 : rd;

        double yd = 0.2627 * rd + 0.6780 * gd + 0.0593 * bd;
        double bs = OOTFhlgInv(bd, yd);
        double gs = OOTFhlgInv(gd, yd);
        double rs = OOTFhlgInv(rd, yd);
        double bg = OETFhlg(bs);
        double gg = OETFhlg(gs);
        double rg = OETFhlg(rs);

        if (sdr) {
          bg = hlg2sdr(bg, kneeStart);
          gg = hlg2sdr(gg, kneeStart);
          rg = hlg2sdr(rg, kneeStart);
        }

        bg = (bg > 1) ? 1 : bg;
        gg = (gg > 1) ? 1 : gg;
        rg = (rg > 1) ? 1 : rg;

        fsCube << rg << " " << gg << " " << bg << std::endl;
      }
    }
  }
  fsCube.close();
}
#endif // DOVI_LUTGEN