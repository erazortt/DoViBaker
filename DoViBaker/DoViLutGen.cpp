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

  const double epower = std::pow(ep, 1 / m2);
  const double num = std::max(epower - c1, 0.0);
  const double denom = c2 - c3 * epower;
  return std::pow(num / denom, 1 / m1);
}

double OOTFhlgInv(double fd, double yd)
{
  if (yd > 0) {
    static constexpr double power = (1 - 1.2) / 1.2;
    return fd * std::pow(yd, power);
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
    // p_2(t)=0 <=> p0*h00_2(1)+m0*h10_2(1)+p1*h01_2(1)+m1_max*h11_2(1)=0 <=> p0*h00_2(1)+m0*h10_2(1)+p1*h01_2(1)=-m1_max*h11_2(1) 
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

void XYZfromRgb2020(double& x, double& y, double& z, double rl, double gl, double bl) {
  x = 0.636958048 * rl + 0.144616904 * gl + 0.168880975 * bl;
  y = 0.262700212 * rl + 0.677998072 * gl + 0.059301716 * bl;
  z = 0.000000000 * rl + 0.028072693 * gl + 1.060985058 * bl;
}
void rgb2020FromXYZ(double& rl, double& gl, double& bl, double x, double y, double z) {
  rl = 1.716651188 * x - 0.355670784 * y - 0.253366281 * z;
  gl = -0.666684352 * x + 1.616481237 * y + 0.015768546 * z;
  bl = 0.017639857 * x - 0.042770613 * y + 0.942103121 * z;
}

void XYZfromRgb709(double& x, double& y, double& z, double rl, double gl, double bl) {
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
  l = std::pow(l, 0.3333333333);
  m = std::pow(m, 0.3333333333);
  s = std::pow(s, 0.3333333333);
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
  XYZfromRgb2020(x, y, z, rl, gl, bl);
  double l, m, s;
  lmsFromXYZ(l, m, s, x, y, z);
  labFromLms(L, a, b, l, m, s);
}
void labFromRgb709(double& L, double& a, double& b, double rl, double gl, double bl) {
  double x, y, z;
  XYZfromRgb709(x, y, z, rl, gl, bl);
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
  XYZfromRgb2020(x, y, z, rl, gl, bl);
  rgb709FromXYZ(rl, gl, bl, x, y, z);
}

// slightly decreases the chroma such that the shrunken BT.2020 gamut tangentialy hits the BT.709 gamut
// this first happens along the blue edge of BT.709 (0,0,b)
void reduceChroma(double& rl, double& gl, double& bl) {
  double L, a, b;
  labFromRgb2020(L, a, b, rl, gl, bl);
  double f = 0.9418;
  rgb2020fromLab(rl, gl, bl, L, a * f, b * f);
}

typedef void (*rgbFromLab)(double& rl, double& gl, double& bl, double L, double a, double b);
typedef void (*labFromRgb)(double& L, double& a, double& b, double rl, double gl, double bl);

bool outOfRange(double r, double g, double b) {
  return b > 1. || g > 1. || r > 1.;
}

void findInRangeLab(double& Li, double& Ai, double& Bi, double ri, double gi, double bi, double sensL, double sensA, double sensB, rgbFromLab rgbFromLabFunc) {
  double f = 1;
  for (int d = 2; d < (1 << 29) + 1; d *= 2) {
    if (outOfRange(ri, gi, bi)) {
      f -= 1.0 / d;
    }
    else {
      f += 1.0 / d;
    }
    double facL = 1 - sensL * (1 - f);
    double facA = 1 - sensA * (1 - f);
    double facB = 1 - sensB * (1 - f);
    rgbFromLabFunc(ri, gi, bi, Li * facL, Ai * facA, Bi * facB);
  }

  double facL = 1 - sensL * (1 - f);
  double facA = 1 - sensA * (1 - f);
  double facB = 1 - sensB * (1 - f);
  Li *= facL;
  Ai *= facA;
  Bi *= facB;
}

void clipToPoint(double& ri, double& gi, double& bi, rgbFromLab rgbFromLabFunc, labFromRgb labFromRgbFunc) {
  // changes the chroma of an out-of-range color such that it becomes in-range, while leaving lightness and hue constant.
  // unfortunately this does not produce good results, the approch varying also lightness the way to go. See function below.

  double Li, Ai, Bi;
  labFromRgbFunc(Li, Ai, Bi, ri, gi, bi);

  double Lo = Li;
  double Ao = Ai;
  double Bo = Bi;
  findInRangeLab(Lo, Ao, Bo, ri, gi, bi, 0, 1, 1, rgbFromLabFunc);

  rgbFromLabFunc(ri, gi, bi, Lo, Ao, Bo);
}

void clipToTangentLine(double& ri, double& gi, double& bi, rgbFromLab rgbFromLabFunc, labFromRgb labFromRgbFunc) {
  // finds the most similar looking in-gamut color for an out-of-range color, while leaving hue constant.

  double Li, Ai, Bi;
  labFromRgbFunc(Li, Ai, Bi, ri, gi, bi);

  //first in-range point on the surface of the target gamut
  double L1 = Li;
  double A1 = Ai;
  double B1 = Bi;
  findInRangeLab(L1, A1, B1, ri, gi, bi, 0, 1, 1, rgbFromLabFunc);

  //second in-range point on the surface of the target gamut
  double L2 = Li;
  double A2 = Ai;
  double B2 = Bi;
  findInRangeLab(L2, A2, B2, ri, gi, bi, 1, 1, 1, rgbFromLabFunc);

  //calculate normal vector from given out-of-range point to the line connecting the two in-range points
  //the intersection point of that line to the normal is the nearest in-range point
  double r = (A1 * A1 + A2 * Ai - A1 * (A2 + Ai) + (B1 - B2) * (B1 - Bi) + (L1 - L2) * (L1 - Li)) / ((A1 - A2) * (A1 - A2) + (B1 - B2) * (B1 - B2) + (L1 - L2) * (L1 - L2));
  double Lo = L1 + (L2 - L1) * r;
  double Ao = A1 + (A2 - A1) * r;
  double Bo = B1 + (B2 - B1) * r;

  rgbFromLabFunc(ri, gi, bi, Lo, Ao, Bo);
}

void softClip709(double& rl, double& gl, double& bl) {
  clipToTangentLine(rl, gl, bl, rgb709fromLab, labFromRgb709);
}
void softClip2020(double& rl, double& gl, double& bl) {
  clipToTangentLine(rl, gl, bl, rgb2020fromLab, labFromRgb2020);
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

void selfTest() {
  {
    double r = 1.25;
    double g = 0;
    double b = 0;
    softClip2020(r, g, b);
  }
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
          softClip2020(rs, gs, bs);
          softClip2020(rs, gs, bs);
        }*/

        if (rs < 0)rs = 0;
        if (gs < 0)gs = 0;
        if (bs < 0)bs = 0;

        double rg = OETFhlg(rs);
        double gg = OETFhlg(gs);
        double bg = OETFhlg(bs);

        if (rg > 1)rg = 1;
        if (gg > 1)gg = 1;
        if (bg > 1)bg = 1;

        if (sdr) {
          rg = hlg2sdr(rg, sdrKneeStart, sdrKneeEndTangentFactor);
          gg = hlg2sdr(gg, sdrKneeStart, sdrKneeEndTangentFactor);
          bg = hlg2sdr(bg, sdrKneeStart, sdrKneeEndTangentFactor);

          if (sdr>1) {
            double rl = EOTFsdr(rg);
            double gl = EOTFsdr(gg);
            double bl = EOTFsdr(bg);
            //reduceChroma(rl, gl, bl);
            convert2020To709(rl, gl, bl);

            if (rl < 0)rl = 0;
            if (gl < 0)gl = 0;
            if (bl < 0)bl = 0;

            /*if (outOfRange(rl, gl, bl)) {
              softClip709(rl, gl, bl);
              softClip709(rl, gl, bl);
              if (rl < 0)rl = 0;
              if (gl < 0)gl = 0;
              if (bl < 0)bl = 0;
            }*/

            rg = OETFsdr(rl);
            gg = OETFsdr(gl);
            bg = OETFsdr(bl);
          }

          if (rg > 1)rg = 1;
          if (gg > 1)gg = 1;
          if (bg > 1)bg = 1;
        }

        fsCube << rg << " " << gg << " " << bg << std::endl;
      }
    }
  }
  fsCube.close();
}
#endif // DOVI_LUTGEN