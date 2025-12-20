#include "def.h"
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <filesystem>

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

/*double OETFsdr(double l) {
  static constexpr double a = 1.099296827;
  static constexpr double b = 0.018053969;
  static constexpr double c = a - b - 1;
  if (l < b)
    return 4.5 * l;
  return a * std::pow(l, 0.45) - (a - 1);
}
double OETFsdrInv(double g) {
  static constexpr double a = 1.099296827;
  static constexpr double b = 0.018053969;
  static constexpr double c = a - b - 1;
  if (g < c)
    return g/4.5;
  return std::pow((g + (a - 1))/a, 1.0/0.45);
}*/
double EOTFsdrInv(double l) {
  return std::pow(l, 1/2.4);
}
double EOTFsdr(double g) {
  return std::pow(g, 2.4);
}

static constexpr double unitHermiteSpline(double t, double y0, double m0, double y1, double m1) {
  // implementation of the hermite spline on the unit interval:
  // p(t)=p0*h00(t)+m0*h10(t)+p1*h01(t)+m1*h11(t)
  // with t(x)=(x-x0)/(x1-x0), h00(t)=(2*t-3)*t*t+1, h10(t)=(((t-2)*t+1)*t, h10=(-2*t+3)*t*t, h11(t)=t*t*(t-1)
  // variables used are also known as: x0=kS, x1=kE, y0=p0, y1=p1
  // d^2p(t)/dt^2=p_2(t)=p0*h00_2(t)+m0*h10_2(t)+p1*h01_2(t)+m1*h11_2(t)
  // with h00_2(t)=12*t-6, h10_2(t)=6*t-4, h01_2(t)=-12*t+6, h11_2(t)=6*t-2
  double p = ((2 * t - 3) * t * t + 1) * y0 + ((t - 2) * t + 1) * t * m0 + (-2 * t + 3) * t * t * y1 + (t - 1) * t * t * m1;
  return p;
}

static constexpr double hermiteSpline(double x, double x0, double y0, double m0, double x1, double y1, double m1) {
  // spline scaled from the unit interval to [x0,x1]
  if (x < x0) return y0;
  if (x > x1) return y1;
  m0 *= (x1 - x0);
  m1 *= (x1 - x0);
  double t = (x - x0) / (x1 - x0);
  double p = unitHermiteSpline(t, y0, m0, y1, m1);
  return p;
}

// this function maps the hlg signal to the sdr signal until it runs out of sdr signal range at 267.6nits hlg
// thus when the sdr is viewed it appears idential to the hlg as long as the hlg stays below these 267.6 nits
double matchHlg2Sdr(double x) {
  if (x > 0.5) {
    return 2.8188 * x * x - 2.7718 * x + 1.6812;
  }
  // hlg and sdr are itentical below 0.5 signal strength
  return 1;
}
// derivative of the matching functions
double matchHlg2SdrD(double x) {
  if (x > 0.5) {
    return 2 * 2.8188 * x - 2.7718;
  }
  return 0;
}

// since the function above will run out of sdr signal space, we need to moderate it so that it does not clip
// this flattening is done by a hermite spline (just like in hdr tonemapping)
// the start point is defined by kS, p0 and m0 and the end point by kE, p1 and m1.
double hlg2sdr(double x, double kS, double m1Factor) {
  if (x > kS) {
    double p0 = kS * matchHlg2Sdr(kS);
    double m0 = kS * matchHlg2SdrD(kS) + 1 * matchHlg2Sdr(kS);
    m0 *= (1 - kS); // must be scaled since we are not on the interval x=[0,1] but on the affine t=[0,1]
    double p1 = 1;

    // calculate the maximal m1 such that the curvature of the hermite spline p(t) does not becomes positive near t=1
    // d^2p(t)/dt^2=p_2(t) <=> p0*h00_2(1)+m0*h10_2(1)+p1*h01_2(1)+m1_max*h11_2(1)=0 <=> p0*h00_2(1)+m0*h10_2(1)+p1*h01_2(1)=-m1_max*h11_2(1) 
    // => p0*6+m0*2-p1*6=-4*m1_max <=> m1_max=-(p0*6+m0*2-p1*6)/4
    double m1max = -(p0*6+m0*2-p1*6)/4;
    double m1 = m1Factor * m1max;
    
    double kE = 1;
    double t = (x - kS) / (kE - kS);
    double p = unitHermiteSpline(t, p0, m0, p1, m1);
    return p;
  }
  return x * matchHlg2Sdr(x);
}

void xyFromXYZ(double& x, double &y, double X, double Y, double Z) {
  double norm = (X + Y + Z);
  x = X / norm;
  y = Y / norm;
}

void XYZfromRgb2020(double& X, double& Y, double& Z, double rl, double gl, double bl) {
  X = 0.636958048 * rl + 0.144616904 * gl + 0.168880975 * bl;
  Y = 0.262700212 * rl + 0.677998072 * gl + 0.059301716 * bl;
  Z = 0.000000000 * rl + 0.028072693 * gl + 1.060985058 * bl;
}
void rgb2020fromXYZ(double& rl, double& gl, double& bl, double X, double Y, double Z) {
  rl = 1.716651189 * X - 0.3556707849 * Y - 0.2533662813 * Z;
  gl =-0.6666843518 * X + 1.616481236 * Y + 0.01576854646 * Z;
  bl = 0.01763985741 * X - 0.04277061315 * Y + 0.942103121 * Z;
}

void XYZfromRgb709(double& X, double& Y, double& Z, double rl, double gl, double bl) {
  X = 0.412390799 * rl + 0.357584339 * gl + 0.180480788 * bl;
  Y = 0.212639006 * rl + 0.715168679 * gl + 0.072192315 * bl;
  Z = 0.019330819 * rl + 0.119194780 * gl + 0.950532152 * bl;
}
void rgb709fromXYZ(double& rl, double& gl, double& bl, double X, double Y, double Z) {
  rl = 3.240969944 * X - 1.537383176 * Y - 0.4986107602 * Z;
  gl = -0.9692436371 * X + 1.875967501 * Y + 0.04155505794 * Z;
  bl = 0.05563007901 * X - 0.2039769588 * Y + 1.056971515 * Z;
}

/*void LsAsBsFromXYZ(double& Ls, double& As, double& Bs, double X, double Y, double Z) {
  // conversion to CIELAB without all the constants
  constexpr double Xw = 0.950455927;
  constexpr double Yw = 1;
  constexpr double Zw = 1.089057751;
  Ls = Y / Yw;
  As = X / Xw - Ls;
  Bs = Ls - Z / Zw;
}
void XYZfromLsAsBs(double& X, double& Y, double& Z, double Ls, double As, double Bs) {
  constexpr double Xw = 0.950455927;
  constexpr double Yw = 1;
  constexpr double Zw = 1.089057751;
  Y = Ls * Yw;
  X = (As + Ls) * Xw;
  Z = (Ls - Bs) * Zw;
}
void LsAsBsFromRgb709(double& Ls, double& As, double& Bs, double rl, double gl, double bl) {
  double X, Y, Z;
  XYZfromRgb709(X, Y, Z, rl, gl, bl);
  LsAsBsFromXYZ(Ls, As, Bs, X, Y, Z);
}
void LsAsBsFromRgb2020(double& Ls, double& As, double& Bs, double rl, double gl, double bl) {
  double X, Y, Z;
  XYZfromRgb2020(X, Y, Z, rl, gl, bl);
  LsAsBsFromXYZ(Ls, As, Bs, X, Y, Z);
}
void rgb709fromLsAsBs(double& rl, double& gl, double& bl, double Ls, double As, double Bs) {
  double X, Y, Z;
  XYZfromLsAsBs(X, Y, Z, Ls, As, Bs);
  rgb709fromXYZ(rl, gl, bl, X, Y, Z);
}
void rgb2020fromLsAsBs(double& rl, double& gl, double& bl, double Ls, double As, double Bs) {
  double X, Y, Z;
  XYZfromLsAsBs(X, Y, Z, Ls, As, Bs);
  rgb2020fromXYZ(rl, gl, bl, X, Y, Z);
}*/

void lmsFromXYZ(double& l, double& m, double& s, double X, double Y, double Z) {
  l = 0.8189330101 * X + 0.3618667424 * Y - 0.1288597137 * Z;
  m = 0.0329845436 * X + 0.9293118715 * Y + 0.0361456387 * Z;
  s = 0.0482003018 * X + 0.2643662691 * Y + 0.6338517070 * Z;
}
void XYZfromLms(double& X, double& Y, double& Z, double l, double m, double s) {
  X = 1.227013851 * l - 0.5577999807 * m + 0.281256149 * s;
  Y = -0.04058017842 * l + 1.11225687 * m - 0.07167667867 * s;
  Z = -0.07638128451 * l - 0.4214819784 * m + 1.58616322 * s;
}

void labFromLms(double& L, double& a, double& b, double l, double m, double s) {
  // the original exponent of 3.0 produces a blue color line which curls into itself
  constexpr double exponent = 2.3605929376;
  l = std::pow(l, 1/exponent);
  m = std::pow(m, 1/exponent);
  s = std::pow(s, 1/exponent);
  L = 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s;
  a = 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s;
  b = 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s;
}
void lmsFromLab(double& l, double& m, double& s, double L, double a, double b) {
  l = 0.9999999985 * L + 0.3963377922 * a + 0.2158037581 * b;
  m = 1.000000009 * L - 0.1055613423 * a - 0.06385417477 * b;
  s = 1.000000055 * L - 0.08948418209 * a - 1.291485538 * b;
  // the original exponent of 3.0 produces a blue color line which curls into itself
  constexpr double exponent = 2.3605929376;
  l = std::pow(l, exponent);
  m = std::pow(m, exponent);
  s = std::pow(s, exponent);
}
void XYZfromLab(double& X, double& Y, double& Z, double L, double a, double b) {
  double l, m, s;
  lmsFromLab(l, m, s, L, a, b);
  XYZfromLms(X, Y, Z, l, m, s);
}
void labFromXYZ(double& L, double& a, double& b, double X, double Y, double Z) {
  double l, m, s;
  lmsFromXYZ(l, m, s, X, Y, Z);
  labFromLms(L, a, b, l, m, s);
}

void labFromRgb2020(double& L, double& a, double& b, double rl, double gl, double bl) {
  double X, Y, Z;
  XYZfromRgb2020(X, Y, Z, rl, gl, bl);
  labFromXYZ(L, a, b, X, Y, Z);
}
void labFromRgb709(double& L, double& a, double& b, double rl, double gl, double bl) {
  double X, Y, Z;
  XYZfromRgb709(X, Y, Z, rl, gl, bl);
  labFromXYZ(L, a, b, X, Y, Z);
}
void rgb2020fromLab(double& rl, double& gl, double& bl, double L, double a, double b) {
  double X, Y, Z;
  XYZfromLab(X, Y, Z, L, a, b);
  rgb2020fromXYZ(rl, gl, bl, X, Y, Z);
}
void rgb709fromLab(double& rl, double& gl, double& bl, double L, double a, double b) {
  double X, Y, Z;
  XYZfromLab(X, Y, Z, L, a, b);
  rgb709fromXYZ(rl, gl, bl, X, Y, Z);
}

void convert2020To709(double& rl, double& gl, double& bl) {
  double x, y, z;
  XYZfromRgb2020(x, y, z, rl, gl, bl);
  rgb709fromXYZ(rl, gl, bl, x, y, z);
}

void reduceChroma(double& rl, double& gl, double& bl, double f = 0.9418) {
  double L, a, b;
  labFromRgb2020(L, a, b, rl, gl, bl);
  rgb2020fromLab(rl, gl, bl, L, a * f, b * f);
}

typedef void (*rgbFromLab)(double& rl, double& gl, double& bl, double L, double a, double b);
typedef void (*labFromRgb)(double& L, double& a, double& b, double rl, double gl, double bl);

bool outOfLumaRange(double r, double g, double b) {
  return r > 1. || g > 1. || b > 1.;
}
bool outOfChromaRange(double r, double g, double b) {
  return r < 0 || g < 0 || b < 0;
}
bool outOfRange(double r, double g, double b) {
  return outOfChromaRange(r, g, b) || outOfLumaRange(r, g, b);
}
bool outOfChromaRangeApprox(double r, double g, double b) {
  return r < -0.000001 || g < -0.000001 || b < -0.000001;
}
bool outOfRangeApprox(double r, double g, double b) {
  return outOfChromaRangeApprox(r, g, b) || outOfLumaRange(r, g, b);
}

// This finction finds the border of the color space in direction along the saturation of the given color
// The value returned is the position of the given color as a percentage of the distance to the border
double findBorderInside(double Li, double& Ai, double& Bi, rgbFromLab rgbFromLabFunc) {
  double f = 1;
  double r, g, b;
  for (int_fast32_t d = 2; d < (1 << 29) + 1; d *= 2) {
    rgbFromLabFunc(r, g, b, Li, Ai * f, Bi * f);
    if (outOfRangeApprox(r, g, b)) {
      f -= 1.0 / d;
    }
    else {
      f += 1.0 / d;
    }
  }
  Ai *= f;
  Bi *= f;

  //double x, y, X, Y, Z;
  //XYZfromLab(X, Y, Z, Li, Ai, Bi);
  //xyFromXYZ(x, y, X, Y, Z);

  return 1/f;
}
double findBorderOutside(double Li, double& Ai, double& Bi, rgbFromLab rgbFromLabFunc) {
  double f = 1;
  double r, g, b;
  for (int_fast32_t d = 2; d < (1 << 29) + 1; d *= 2) {
    rgbFromLabFunc(r, g, b, Li, Ai / f, Bi / f);
    if (outOfRangeApprox(r, g, b)) {
      f += 1.0 / d;
    }
    else {
      f -= 1.0 / d;
    }
  }
  Ai /= f;
  Bi /= f;

  //double x, y, X, Y, Z;
  //XYZfromLab(X, Y, Z, Li, Ai, Bi);
  //xyFromXYZ(x, y, X, Y, Z);

  return f;
}

// Returns the bt709 color for a given bt2020 color preserving the look as much as possible while preventing any clipping.
double hybridColorMapping(double& ri, double& gi, double& bi) {
  double rc = ri, gc = gi, bc = bi;
  convert2020To709(rc, gc, bc);

  double L, A, B;
  labFromRgb2020(L, A, B, ri, gi, bi); // same as labFromRgb709(L, A, B, rc, gc, bc);
  //double C = sqrt(A * A + B * B);
  //double h = atan2(B, A);

  double fracSrc = 0;
  double aHullSrc = A;
  double bHullSrc = B;
  if (outOfRange(ri, gi, bi))
    fracSrc = findBorderInside(L, aHullSrc, bHullSrc, rgb2020fromLab);
  else
    fracSrc = findBorderOutside(L, aHullSrc, bHullSrc, rgb2020fromLab);
 
  double fracDst = 0;
  double aHullDst = A;
  double bHullDst = B;
  if (outOfRange(rc, gc, bc))
    fracDst = findBorderInside(L, aHullDst, bHullDst, rgb709fromLab);
  else
    fracDst = findBorderOutside(L, aHullDst, bHullDst, rgb709fromLab);
  
  // the amount of reduction along the given color saturaion line
  double reductionFactor = fracDst / fracSrc;
  if (fracSrc > 1) {
    // this should only be numerical flukes
    fracSrc = 1;
  }

  double aLinScaled = aHullDst * fracSrc; //same as A/reduction;
  double bLinScaled = bHullDst * fracSrc; //same as B/reduction;

  double weight1 = 0;
  double weight2 = 1;
  if (reductionFactor > 1) {
    weight1 = fracSrc;
    weight2 = (1 - fracSrc);
    if (reductionFactor < 2) {
      // the weight can be increased so that we remain nearer to the original color
      weight2 = ((fracSrc - 1) * ((reductionFactor - 2) * (reductionFactor - 1) * fracSrc - reductionFactor)) / reductionFactor;
    }
  }
  if (reductionFactor > 2) {
    // the weight must be decreased such that we don't overshoot
    double weight2sup = (1 + (reductionFactor / 2 - 1) * fracSrc);
    weight2 /= weight2sup;
  }
  double aWeighted = aLinScaled * weight1 + A * weight2;
  double bWeighted = bLinScaled * weight1 + B * weight2;
  //double Co = sqrt(aWeighted * aWeighted + bWeighted * bWeighted);
  //double ho = atan2(bWeighted, aWeighted);

  //double x, y, X, Y, Z;
  //XYZfromLab(X, Y, Z, L, aHullDst, bHullDst);
  //xyFromXYZ(x, y, X, Y, Z);

  rgb709fromLab(ri, gi, bi, L, aWeighted, bWeighted);

  return reductionFactor;
}

void selfTestHybridConversion() {
  double r, g, b, f;

  for (int i = 0; i < 101; i++) {
    r = 0 + (i * 0.01); g = 0 + (i * 0.01); b = 1 - (i * 0.01);
    f = hybridColorMapping(r, g, b);
  }
  r = 0; g = 0; b = 0;
  f = hybridColorMapping(r, g, b);
  r = 1; g = 1; b = 1;
  f = hybridColorMapping(r, g, b);
  r = 1; g = 0; b = 0;
  f = hybridColorMapping(r, g, b);
  r = 0; g = 1; b = 0;
  f = hybridColorMapping(r, g, b);
  r = 0; g = 0; b = 1;
  f = hybridColorMapping(r, g, b);
  r = 0; g = 1; b = 1;
  f = hybridColorMapping(r, g, b);
  r = 1; g = 0; b = 1;
  f = hybridColorMapping(r, g, b);
  r = 1; g = 1; b = 0;
  f = hybridColorMapping(r, g, b);
  r = 0.2; g = 0.4; b = 0.5;
  f = hybridColorMapping(r, g, b);
}

void showBestLutSizes() {
  printf("For non-renormalized PQ inputs only the following LUT sizes should be used:\n");
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

void showUsage(std::string& execname) {
  printf("USAGE:\n");
  printf("%s <string> -s|--size <string> [-n|--normalized] [-d|--sdr <int>] [-g|--gain <float>] [-c|--compression <float>] [-r|--reduction <float>]\n", execname.c_str());
  printf("\n");
  printf("ARGUMENTS:\n");
  printf("<string>                   name of the to-be-generated LUT file\n");
  printf("-s, --size <int>           size along one dimension of the LUT\n");
  printf("-i, --input <int>          expected input format to convert from:\n");
  printf("                             0: BT.2100 PQ (default)\n");
  printf("                             1: BT.2100 PQ which was re-normalized to 1000nits\n");
  printf("                             2: BT.2100 HLG\n");
  printf("                             3: BT.2020 SDR\n");
  printf("-o, --output <int>         output format to convert to:\n");
  printf("                             0: BT.2100 HLG (default)\n");
  printf("                             1: BT.2020 SDR\n");
  printf("                             2: BT.709 SDR with hard color clipping\n");
  printf("                             3: BT.709 SDR with advanced color mapping\n");
  printf("-g, --gain <float>         gain of light midtones in conversions to SDR (0.0)\n");
  printf("-c, --compression <float>  highlight compression in conversions to SDR (0.0)\n");
  printf("-r, --reduction <float>    chroma reduction factor in SDR conversions (1.0)\n");
  printf("\n");
  printf("examples for PQ to HLG conversions:\n");
  printf("%s pq2hlg.cube -s 65 -i 0 -o 0\n", execname.c_str());
  printf("%s pq2hlg_normalizedInput.cube -s 50 -i 1 -o 0\n", execname.c_str());
  printf("\n");
  printf("examples for PQ to BT.709 SDR conversions:\n");
  printf("%s pq2sdr709.cube -s 65 -i 0 -o 3\n", execname.c_str());
  printf("%s pq2sdr709_normalizedInput.cube -s 50 -i 1 -o 3\n", execname.c_str());
  printf("\n");
  printf("example for HLG to BT.709 SDR conversions:\n");
  printf("%s hlg2sdr709.cube -s 50 -i 2 -o 3\n", execname.c_str());
  printf("\n");
  printf("example for BT.2020 SDR to BT.709 SDR conversions:\n");
  printf("%s bt2020to709.cube -s 50 -i 3 -o 3\n", execname.c_str());
  printf("\n");
  showBestLutSizes();
}

bool hasArgument(
  std::vector<std::string>& args,
  const std::string& option_name,
  const std::string& option_alt_name) {
  for (auto it = args.begin(), end = args.end(); it != end; ++it) {
    if (*it == option_name || *it == option_alt_name) {
      return true;
    }
  }
  return false;
}
bool getFlag(
  std::vector<std::string>& args,
  const std::string& option_name,
  const std::string& option_alt_name) {
  for (auto it = args.begin(), end = args.end(); it != end; ++it) {
    if (*it == option_name || *it == option_alt_name) {
      args.erase(it);
      return true;
    }
  }
  return false;
}
std::string getArgumentValue(
  std::vector<std::string>& args,
  const std::string& option_name,
  const std::string& option_alt_name) {
  for (auto it = args.begin(), end = args.end(); it != end; ++it) {
    if (*it == option_name || *it == option_alt_name) {
      std::string value = *(it + 1);
      args.erase(it + 1);
      args.erase(it);
      return value;
    }
  }
  return "";
}
std::string getPositionalArgument(std::vector<std::string>& args) {
  std::string value = *(args.begin());
  args.erase(args.begin());
  return value;
}

#ifdef DOVI_LUTGEN
int main(int argc, char* argv[])
{
  //selfTestHybridConversion();

  std::string execname = std::filesystem::path(*argv).filename().string();
  if (argc < 5) {
    showUsage(execname);
    return 1;
  }

  std::vector<std::string> args(argv + 1, argv + argc);
  std::string lutFileName;
  int lutSize;
  bool reNormalizedInput = false;
  int inputFormat = 0;
  int outputFormat = 0;
  double sdrGain = 0.0;
  double sdrCompression = 0.0;
  double chromaReduction = 1.0;

  if (!hasArgument(args, "-s", "--size")) {
    printf("LUT size not given!\n\n");
    showUsage(execname);
    return 1;
  }
  lutSize = std::stoi(getArgumentValue(args, "-s", "--size"));
  if (hasArgument(args, "-i", "--input")) {
    inputFormat = std::stoi(getArgumentValue(args, "-i", "--input"));
  }
  if (hasArgument(args, "-o", "--output")) {
    outputFormat = std::stoi(getArgumentValue(args, "-o", "--output"));
  }
  if (hasArgument(args, "-g", "--gain")) {
    sdrGain = std::stof(getArgumentValue(args, "-g", "--gain"));
  }
  if (hasArgument(args, "-c", "--compression")) {
    sdrCompression = std::stof(getArgumentValue(args, "-c", "--compression"));
  }
  if (hasArgument(args, "-r", "--reduction")) {
    chromaReduction = std::stof(getArgumentValue(args, "-r", "--reduction"));
  }
  if (args.begin() == args.end()) {
    printf("No output filename given!\n\n");
    showUsage(execname);
    return 2;
  }
  lutFileName = getPositionalArgument(args);
  if (args.size() > 0) {
    printf("Unknown arguments given!\n\n");
    showUsage(execname);
    return 2;
  }

  if (lutSize < 1 || lutSize > 200) {
    printf("LUT size not supported");
    return 2;
  }
  std::string inFormatName = "BT.2100 PQ";
  switch (inputFormat)
  {
  case 0: inFormatName = "BT.2100 PQ"; break;
  case 1: inFormatName = "BT.2100 PQ"; reNormalizedInput = true; break;
  case 2: inFormatName = "BT.2100 HLG"; break;
  case 3: inFormatName = "BT.2020 SDR"; break;
  default:
    break;
  }
  std::string outFormatName = "BT.2100 HLG";
  switch (outputFormat)
  {
  case 0: outFormatName = "BT.2100 HLG"; break;
  case 1: outFormatName = "BT.2020 SDR"; break;
  case 2: outFormatName = "BT.709 SDR"; break;
  case 3: outFormatName = "BT.709 SDR"; break;
  default:
    break;
  }
  if (!inputFormat) {
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
  if (inputFormat > 2 && outputFormat < 2) {
    printf("Unsupported conversion");
    return 2;
  }
  if (inputFormat > 1 && outputFormat < 1) {
    printf("Unsupported conversion");
    return 2;
  }
  if (sdrGain < 0) {
    printf("Gain cannot be negative");
    return 2;
  }
  if (sdrCompression > 1) {
    printf("Compression cannot be above 1");
    return 2;
  }
  if (chromaReduction < 0) {
    printf("Reduction cannot be negative");
    return 2;
  }

  std::ofstream fsCube;
  fsCube.open(lutFileName, std::ifstream::out);
  if (!fsCube.is_open()) {
    fsCube.close();
    printf("Unable to open output file\n");
    return 2;
  }

  double sdrKneeStart = std::sqrt(sdrGain) * 0.2103131223 + 0.5;

  fsCube << "# Generated by DoViLutGen" << std::endl;
  fsCube << "# in format: " << inFormatName << std::endl;
  if (inputFormat < 2) {
    fsCube << "# re-normalized input: " << reNormalizedInput << std::endl;
  }
  fsCube << "# out format: " << outFormatName << std::endl;
  if (outputFormat) {
    if (inputFormat < 3) {
      fsCube << "# midtone gain: " << sdrGain << " (kS=" << sdrKneeStart << ")" << std::endl;
      fsCube << "# highlight compression: " << sdrCompression << " (=m1Factor)" << std::endl;
    }
    if (outputFormat > 1) {
      fsCube << "# chroma reduction factor: " << chromaReduction << std::endl;
      fsCube << "# color conversion: ";
      if (outputFormat > 2) {
        fsCube << "advanced mapping" << std::endl;
      }
      else {
        fsCube << "hard clipping" << std::endl;
      }
    }
    printf("Producing a LUT for %s -> %s conversions of size %i\n", inFormatName.c_str(), outFormatName.c_str(), lutSize);
    fsCube << "# LUT for conversions from " << inFormatName << " to " << outFormatName << std::endl;
  }
  else {
    printf("Producing LUT for PQ -> HLG conversions of size %i\n", lutSize);
    fsCube << "# LUT for conversions from BT.2100 HDR PQ to BT.2100 HDR HLG following the BBC specs with 75% reference white" << std::endl;
  }

  double inputScale = 1;
  if (reNormalizedInput) {
    printf("The LUT will expect the PQ input to be re-normalized to 1000 nits\n");
    inputScale = 0.7518271;
    fsCube << "# ATTENTION: This special LUT expects the PQ input to be re-normalized to 1000nits max brightness!" << std::endl;
    fsCube << "TITLE \"PQ renorm 1000nits to " << outFormatName << "\"" << std::endl;
  }
  else {
    fsCube << "TITLE \"" << inFormatName << " to " << outFormatName << "\"" << std::endl;
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
        //if (gi != ri)continue;
        //if (bi != lutSize - 1) continue;
        //if (gi != ri)continue;
        //if (gi != 25-(bi-25)) continue;
        //if (ri != 23 || gi != 29 || bi != 36) continue;
        double erp = double(ri) / (lutSize - 1) * inputScale;
        double rd = EOTFpq(erp) * 10000 / 1000;
        rd = (rd > 1) ? 1 : rd;

        double yd = 0.2627 * rd + 0.6780 * gd + 0.0593 * bd;
        double rs = OOTFhlgInv(rd, yd);
        double gs = OOTFhlgInv(gd, yd);
        double bs = OOTFhlgInv(bd, yd);

        if (rs < 0)rs = 0;
        if (gs < 0)gs = 0;
        if (bs < 0)bs = 0;

        double rg = OETFhlg(rs);
        double gg = OETFhlg(gs);
        double bg = OETFhlg(bs);

        if (rg > 1)rg = 1;
        if (gg > 1)gg = 1;
        if (bg > 1)bg = 1;

        if (inputFormat == 2) {
          rg = double(ri) / (lutSize - 1);
          gg = double(gi) / (lutSize - 1);
          bg = double(bi) / (lutSize - 1);
        }

        if (outputFormat) {
          rg = hlg2sdr(rg, sdrKneeStart, sdrCompression);
          gg = hlg2sdr(gg, sdrKneeStart, sdrCompression);
          bg = hlg2sdr(bg, sdrKneeStart, sdrCompression);

          if (inputFormat == 3) {
            rg = double(ri) / (lutSize - 1);
            gg = double(gi) / (lutSize - 1);
            bg = double(bi) / (lutSize - 1);
          }

          if (outputFormat > 1) {
            double rl = EOTFsdr(rg);
            double gl = EOTFsdr(gg);
            double bl = EOTFsdr(bg);
            if (chromaReduction != 1.0)
              reduceChroma(rl, gl, bl, chromaReduction);
            if(outputFormat > 2)
              hybridColorMapping(rl, gl, bl);
            else
              convert2020To709(rl, gl, bl);
            
            //fsCube << ri << " " << gi << " " << bi << " " << rl << " " << gl << " " << bl << std::endl;

            if (rl < 0)rl = 0;
            if (gl < 0)gl = 0;
            if (bl < 0)bl = 0;

            rg = EOTFsdrInv(rl);
            gg = EOTFsdrInv(gl);
            bg = EOTFsdrInv(bl);
          }

          if (rg > 1)rg = 1;
          if (gg > 1)gg = 1;
          if (bg > 1)bg = 1;
        }

        //fsCube << std::format("{:.7f} {:.7f} {:.7f}",rg, gg, bg) << std::endl;
        fsCube << rg << " " << gg << " " << bg << std::endl;
      }
    }
  }
  fsCube.close();
}
#endif // DOVI_LUTGEN