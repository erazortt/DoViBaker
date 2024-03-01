#include "def.h"
#include <fstream>
#include <algorithm>

double EOTFpq(double ep)
{
  static constexpr double m1 = 2610.0 / 4096 / 4;
  static constexpr double m2 = 2523.0 / 4096 * 128;
  static constexpr double c3 = 2392.0 / 4096 * 32;
  static constexpr double c2 = 2413.0 / 4096 * 32;
  static constexpr double c1 = c3 - c2 + 1;

  const double epower = pow(ep, 1 / m2);
  const double num = std::max(epower - c1, 0.0);
  const double denom = c2 - c3 * epower;
  return pow(num / denom, 1 / m1);
}

double OOTFhlgInv(double fd, double yd)
{
  if (yd > 0) {
    static constexpr double power = (1 - 1.2) / 1.2;
    return fd * pow(yd, power);
  }
  return 0;
}

double OETFhlg(double e)
{
  static constexpr double a = 0.17883277;
  static constexpr double b = 1 - 4 * a;
  static const double c = 0.5 - a * std::log(4 * a);

  if (e * 12 < 1) {
    return sqrt(3 * e);
  }
  return a * std::log(12 * e - b) + c;
}

double OETFsdr(double l) {
  static constexpr double a = 1.099296827;
  static constexpr double b = 0.018053969;
  static constexpr double c = a - b - 1;

  if (l < b)
    return 4.5 * l;
  return a * std::pow(l, 0.45) - (a - 1);
}

double EOTFsdr(double g) {
  static constexpr double a = 1.099296827;
  static constexpr double b = 0.018053969;
  static constexpr double c = a - b - 1;

  if (g < c)
    return g/4.5;
  return std::pow((g + (a - 1))/a, 1.0/0.45);
}


// this function maps the hlg signal above 203 nits to the sdr signal range until it runs out of sdr signal range at 267.6nits hlg
// thus when the sdr is viewed it appears idential to the hlg as long as the hlg stays below these 267.6 nits
double matchHlg2Sdr(double x) {
  if (x > 0.5) {
    return 2.8188 * x * x - 2.7718 * x + 1.6812;
  }
  // hlg and sdr are itentical below 203 nits
  return 1;
}
// derivative of the matching functions
double matchHlg2SdrD(double x) {
  if (x > 0.5) {
    return 2 * 2.8188 * x - 2.7718;
  }
  return 0;
}

// since the above function will run out of sdr signal space, we need to moderate it so that it does not clip
// this flattening is done by a hermite spline (just like in hdr tonemapping)
// the start point is defined by kS, p0 and m0 and the end point by kE, p1 and m1.
// p(t)=p0*h00(t)+m0*h10(t)+p1*h01(t)+m1*h11(t)
// with t(x)=(x-kS)/(kE-kS), h00(t)=(2*t-3)*t*t+1, h10(t)=(((t-2)*t+1)*t, h10=(-2*t+3)*t*t, h11(t)=t*t*(t-1)
// d^2p(t)/dt^2=p_2(t)=p0*h00_2(t)+m0*h10_2(t)+p1*h01_2(t)+m1*h11_2(t)
// with h00_2(t)=12*t-6, h10_2(t)=6*t-4, h01_2(t)=-12*t+6, h11_2(t)=6*t-2
double hlg2sdr(double x, double kS, double m1Factor) {
  if (x > kS) {
    double p0 = kS * matchHlg2Sdr(kS);
    double m0 = kS * matchHlg2SdrD(kS) + 1 * matchHlg2Sdr(kS);
    m0 *= (1 - kS); // must be scaled since we are not on the interval x=[0,1] but on the affine t=[0,1]
    double p1 = 1;

    // calculate the maximal m1 such that the curvature of p(t) never becomes positive for any t, especially not for t=1
    // p0*h00_2(1)+m0*h10_2(1)+p1*h01_2(1)+m1_max*h11_2(1)=0 <=> p0*h00_2(1)+m0*h10_2(1)+p1*h01_2(1)=-m1_max*h11_2(1) 
    // => p0*6+m0*2-p1*6=-4*m1_max <=> m1_max=-(p0*6+m0*2-p1*6)/4
    double m1max = -(p0*6+m0*2-p1*6)/4;
    double m1 = m1Factor * m1max;
    
    double kE = 1;
    double t = (x - kS) / (kE - kS);
    double p = ((2*t-3)*t*t+1)*p0 + ((t-2)*t+1)*t*m0 + (-2*t+3)*t*t*p1 + (t-1)*t*t*m1;
    return p;
  }
  return x * matchHlg2Sdr(x);
}

bool outOfRange(double r, double g, double b) {
  return b < 0 || b > 1 || g < 0 || g > 1 || r < 0 || r > 1;
}

void xyzFromRgb2020(double& x, double& y, double& z, double rl, double gl, double bl) {
  x = 0.636958048 * rl + 0.144616904 * gl + 0.168880975 * bl;
  y = 0.262700212 * rl + 0.677998072 * gl + 0.059301716 * bl;
  z = 0.000000000 * rl + 0.028072693 * gl + 1.060985058 * bl;
}
void rgb2020FromXYZ(double& rl, double& gl, double& bl, double x, double y, double z) {
  rl = 1.716651188 * x - 0.355670784 * y - 0.253366281 * z;
  gl = -0.666684352 * x + 1.616481237 * y + 0.015768546 * z;
  bl = 0.017639857 * x - 0.042770613 * y + 0.942103121 * z;
}

void xyzFromRgb709(double& x, double& y, double& z, double rl, double gl, double bl) {
  x = 0.412390799 * rl + 0.357584339 * gl + 0.180480788 * bl;
  y = 0.212639006 * rl + 0.715168679 * gl + 0.072192315 * bl;
  z = 0.019330819 * rl + 0.119194780 * gl + 0.950532152 * bl;
}
void rgb709FromXYZ(double& rl, double& gl, double& bl, double x, double y, double z) {
  rl = 3.240969942 * x - 1.537383178 * y - 0.498610760 * z;
  gl = -0.969243636 * x + 1.875967502 * y + 0.041555057 * z;
  bl = 0.055630080 * x - 0.203976959 * y + 1.056971514 * z;
}

void lmsFromXYZ(double& l, double& m, double& s, double x, double y, double z) {
  l = 0.8189330101 * x + 0.3618667424 * y - 0.1288597137 * z;
  m = 0.0329845436 * x + 0.9293118715 * y + 0.0361456387 * z;
  s = 0.0482003018 * x + 0.2643662691 * y + 0.6338517070 * z;
}
void XYZfromLms(double& x, double& y, double& z, double l, double m, double s) {
  x = 1.227013851 * l - 0.557799981 * m + 0.281256149 * s;
  y = -0.040580178 * l + 1.11225687 * m - 0.071676679 * s;
  z = -0.076381285 * l - 0.421481978 * m + 1.58616322 * s;
}

void labFromLms(double& L, double& a, double& b, double l, double m, double s) {
  l = std::pow(l, 0.33333333);
  m = std::pow(m, 0.33333333);
  s = std::pow(s, 0.33333333);
  L = 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s;
  a = 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s;
  b = 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s;
}
void lmsFromLab(double& l, double& m, double& s, double L, double a, double b) {
  l = 1.000000000 * L + 0.396337792 * a + 0.215803758 * b;
  m = 1.000000000 * L - 0.105561342 * a - 0.063854175 * b;
  s = 1.000000000 * L - 0.089484182 * a - 1.291485538 * b;
  l = l*l*l;
  m = m*m*m;
  s = s*s*s;
}

void labFromRgb2020(double& L, double& a, double& b, double rl, double gl, double bl) {
  double x, y, z;
  xyzFromRgb2020(x, y, z, rl, gl, bl);
  double l, m, s;
  lmsFromXYZ(l, m, s, x, y, z);
  labFromLms(L, a, b, l, m, s);
}
void labFromRgb709(double& L, double& a, double& b, double rl, double gl, double bl) {
  double x, y, z;
  xyzFromRgb709(x, y, z, rl, gl, bl);
  double l, m, s;
  lmsFromXYZ(l, m, s, x, y, z);
  labFromLms(L, a, b, l, m, s);
}
void rgb2020fromLab(double& rl, double& gl, double& bl, double L, double a, double b) {
  double l, m, s;
  lmsFromLab(l, m, s, L, a, b);
  double x, y, z;
  XYZfromLms(x, y, z, l, m, s);
  rgb2020FromXYZ(rl, gl, bl, x, y, z);
}
void rgb709fromLab(double& rl, double& gl, double& bl, double L, double a, double b) {
  double l, m, s;
  lmsFromLab(l, m, s, L, a, b);
  double x, y, z;
  XYZfromLms(x, y, z, l, m, s);
  rgb709FromXYZ(rl, gl, bl, x, y, z);
}

void convert2020To709(double& rl, double& gl, double& bl) {
  double x, y, z;
  xyzFromRgb2020(x, y, z, rl, gl, bl);
  rgb709FromXYZ(rl, gl, bl, x, y, z);
}

void clip2020Chroma(double& rl, double& gl, double& bl) {
  double L, a, b;
  labFromRgb2020(L, a, b, rl, gl, bl);
  double f=0.5;
  for (int d = 4; d < 16385; d *= 2) {
    rgb2020fromLab(rl, gl, bl, L, a*f, b*f);
    if (outOfRange(rl, gl, bl)) {
      f -= 1.0 / d;
    }
    else {
      f += 1.0 / d;
    }
  }
}

void clip709Chroma(double& rl, double& gl, double& bl) {
  double L, a, b;
  labFromRgb709(L, a, b, rl, gl, bl);
  double f = 0.5;
  for (int d = 4; d < 8193; d *= 2) {
    rgb709fromLab(rl, gl, bl, L, a * f, b * f);
    if (outOfRange(rl, gl, bl)) {
      f -= 1.0 / d;
    }
    else {
      f += 1.0 / d;
    }
  }
}

void reduceChroma(double& rl, double& gl, double& bl) {
  double L, a, b;
  labFromRgb2020(L, a, b, rl, gl, bl);
  double f = 0.9418;
  rgb2020fromLab(rl, gl, bl, L, a * f, b * f);
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
  if (argc < 3) {
    printf("usage:\n");
    printf("<OUTPUT_FILE> <LUT_SIZE> (<NORMALIZED_INPUT>) (<SDR>) (<GAIN>) (<COMPRESSION>)\n");
    printf("\n");
    printf("examples for PQ to HLG conversions:\n");
    printf("DoViLutGen.exe pq2hlg.cube 65\n");
    printf("DoViLutGen.exe pq2hlg_normalizedInput.cube 50 1\n");
    printf("\n");
    printf("examples for PQ to BT.2020 SDR conversions:\n");
    printf("DoViLutGen.exe pq2sdr2020.cube 65 0 1\n");
    printf("DoViLutGen.exe pq2sdr2020_normalizedInput.cube 50 1 1\n");
    printf("\n");
    printf("examples for PQ to BT.709 SDR conversions:\n");
    printf("DoViLutGen.exe pq2sdr709.cube 65 0 2\n");
    printf("DoViLutGen.exe pq2sdr709_normalizedInput.cube 50 1 2\n");
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
  bool normalizedInput = false;
  if (argc > 3) {
    normalizedInput = std::atoi(argv[3]);
  }

  std::string mode = "HLG";
  int sdr = false;
  double sdrGain = 0.0;
  double sdrCompression = 0.0;
  if (argc > 4) {
    sdr = std::atoi(argv[4]);
    if (sdr) {
      if (sdr > 1) {
        mode = "BT.709 SDR";
      }
      else {
        mode = "BT.2020 SDR";
      }
    }
    if (argc > 5) {
      sdrGain = std::atof(argv[5]);
      if (sdrGain < 0) {
        printf("gain cannot be negative");
        return 2;
      }
      if (argc > 6) {
        sdrCompression = std::atof(argv[6]);
        if (sdrCompression > 1) {
          printf("compression cannot be above 1");
          return 2;
        }
      }
    }
  }
  double sdrKneeStart = std::sqrt(sdrGain) * 0.21 + 0.5;
  double sdrKneeEndTangentFactor = (1 - std::sqrt(1 - sdrCompression));

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

  fsCube << "# Generated by DoViLutGen" << std::endl;
  fsCube << "# re-normalized input: " << normalizedInput << std::endl;
  fsCube << "# mode: " << mode << std::endl;
  if (sdr) {
    fsCube << "# midtone gain: " << sdrGain << " (kS=" << sdrKneeStart << ")" << std::endl;
    fsCube << "# highlight compression: " << sdrCompression << " (m1Factor=" << sdrKneeEndTangentFactor << ")" << std::endl;
    printf("Producing a LUT for PQ -> SDR conversions of size %i\n", lutSize);
    fsCube << "# LUT for conversions from BT.2100 HDR PQ to " << mode << std::endl;
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
    fsCube << "TITLE \"PQ renorm 1000nits to " << mode << "\"" << std::endl;
  }
  else {
    printf("The LUT will expect usual PQ input\n");
    fsCube << "TITLE \"PQ to " << mode << "\"" << std::endl;
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
        double rs = OOTFhlgInv(rd, yd);
        double gs = OOTFhlgInv(gd, yd);
        double bs = OOTFhlgInv(bd, yd);

        /*if (outOfRange(rs, gs, bs)) {
          clip2020Chroma(rs, gs, bs);
        }*/

        double rg = OETFhlg(rs);
        double gg = OETFhlg(gs);
        double bg = OETFhlg(bs);

        bg = std::clamp(bg, 0.0, 1.0);
        gg = std::clamp(gg, 0.0, 1.0);
        rg = std::clamp(rg, 0.0, 1.0);

        if (sdr) {
          rg = hlg2sdr(rg, sdrKneeStart, sdrKneeEndTangentFactor);
          gg = hlg2sdr(gg, sdrKneeStart, sdrKneeEndTangentFactor);
          bg = hlg2sdr(bg, sdrKneeStart, sdrKneeEndTangentFactor);

          if (sdr>1) {
            double rl = EOTFsdr(rg);
            double gl = EOTFsdr(gg);
            double bl = EOTFsdr(bg);
            reduceChroma(rl, gl, bl);
            convert2020To709(rl, gl, bl);
            /*if (outOfRange(rl, gl, bl)) {
              clip709Chroma(rl, gl, bl);
            }*/
            rg = OETFsdr(rl);
            gg = OETFsdr(gl);
            bg = OETFsdr(bl);
          }

          bg = std::clamp(bg, 0.0, 1.0);
          gg = std::clamp(gg, 0.0, 1.0);
          rg = std::clamp(rg, 0.0, 1.0);
        }

        fsCube << rg << " " << gg << " " << bg << std::endl;
      }
    }
  }
  fsCube.close();
}
#endif // DOVI_LUTGEN