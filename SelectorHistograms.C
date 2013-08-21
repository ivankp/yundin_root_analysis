
#include "SelectorHistograms.h"

// --------------------------------------------------------------------------- //
// Histograms
// --------------------------------------------------------------------------- //

Histogram::Histogram(const TString& filename_, const TString& name_,
                     int nbin_, double x1_, double x2_)
  : filename(filename_), name(name_),
    x1(x1_), x2(x2_), x12(x2-x1),
    nbin(nbin_), prevevt(-1), lastidx(0),
    curidx(nbin), events(nbin), curwgt(nbin),
    wgt(nbin), wgt2(nbin), bwidth(nbin), edge(nbin+1)
{ }

TString Histogram::getFile() const
{
  return filename;
}

void Histogram::print(std::ostream& stream, const TString& runname, double count, bool unweight)
{
  flush();
  stream << "# BEGIN HISTOGRAM /" << runname << "/" << name << std::endl;
  stream << "AidaPath=/" << runname << "/" << name << std::endl;
  stream << "Title=" << std::endl;

  double area = 0.;
  for (int i=0; i<nbin; i++) {
    area += wgt[i];
  }
  area /= count;
  stream << "## Area: " << area << std::endl;
  stream << "## Num bins: " << nbin << std::endl;
  stream << "## xlow  \txhigh   \tval    \terrminus\terrplus" << std::endl;

  for (int i=0; i<nbin; i++) {
    const double value = wgt[i]/count;
    double error = 0;
    if (unweight) {
      double variance = wgt2[i]; //events[i] > 1 ? (wgt2[i]-wgt[i]*wgt[i]/events[i])/(events[i]-1) : wgt[i]*wgt[i];
      error = variance > 0. ? sqrt(variance)/count : 0.;
    } else {
      double variance = (wgt2[i]-wgt[i]*wgt[i]/count)/(count-1);
      error = variance > 0. ? sqrt(variance/count) : 0.;
    }

    stream << edge[i] << "\t"
          << edge[i+1] << "\t"
          << value/bwidth[i] << "\t"
          << error/bwidth[i] << "\t"
          << error/bwidth[i] << std::endl;
  }
  stream << "# END HISTOGRAM" << std::endl << std::endl << std::endl;
}

void Histogram::flush()
{
  for (int i=0; i<lastidx; i++) {
    int idx = curidx[i];
    wgt[idx] += curwgt[i];
    wgt2[idx] += curwgt[i]*curwgt[i];
    events[idx]++;
    curwgt[i] = 0.;
  }
  lastidx = 0;
}

void Histogram::setedges() {
  edge[0] = x1;
  for (int i=0; i<nbin; i++) {
    edge[i+1] = edge[i] + bwidth[i];
  }
}

inline
void Histogram::fill(int evt, int n, double w)
{
  assert(0 <= n && n < nbin);

  if (evt != prevevt) {
    flush();
    prevevt = evt;
  }

  int i;
  for (i=0; i<lastidx; i++) {
    if (curidx[i] == n) {
      break;
    }
  }
  curidx[i] = n;
  curwgt[i] += w;
  if (i == lastidx) {
    lastidx++;
  }
}

LinearHistogram::LinearHistogram(const TString& filename_, const TString& name_,
                                 int nbin_, double x1_, double x2_,
                                 double /*param1*/, double /*param2*/, double /*param3*/)
  : Histogram(filename_, name_, nbin_, x1_, x2_), step(x12/nbin)
{
  for (int i=0; i<nbin; i++) {
    bwidth[i] = step;
  }
  setedges();
}

void LinearHistogram::bin(int nextevt, double x, double w)
{
//       std::cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << std::endl; std::cout.flush();
  if (x < x1 or x > x2) return;
  int n = static_cast<int>(nbin*(x-x1)/x12);
  assert(0 <= n && n < nbin);
  fill(nextevt, n, w);
}

SmearedLinearHistogram::SmearedLinearHistogram(const TString& filename_, const TString& name_,
                                               int nbin_, double x1_, double x2_,
                                               double smear_, double /*param2*/, double /*param3*/)
  : LinearHistogram(filename_, name_, nbin_, x1_, x2_), smear(smear_)
{
}

void SmearedLinearHistogram::bin(int nextevt, double x, double w)
{
//       std::cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << std::endl; std::cout.flush();
  if (x < x1 or x > x2) return;
  const double dn = (x-x1)/step;
  const int n = static_cast<int>(dn);
  assert(0 <= n && n < nbin);

  double dist = dn - n; // [0, 1)
  if (dist < 0.5*smear && n > 0) {
    fill(nextevt, n, 0.5*w);
    fill(nextevt, n-1, 0.5*w);
  } else if ((1. - dist) < 0.5*smear && n < nbin-1) {
    fill(nextevt, n, 0.5*w);
    fill(nextevt, n+1, 0.5*w);
  } else {
    fill(nextevt, n, w);
  }
}

QuadraticHistogram::QuadraticHistogram(const TString& filename_, const TString& name_,
                                       int nbin_, double x1_, double x2_,
                                       double f, double /*param2*/, double /*param3*/)
  : Histogram(filename_, name_, nbin_, x1_, x2_),
    step(2*x12/((1 + f)*nbin)),
    slope((f - 1)/(nbin - 1))
{
  for (int i=0; i<nbin; i++) {
    bwidth[i] = step*(1 + i*slope);
  }
  setedges();
}

void QuadraticHistogram::bin(int nextevt, double x, double w)
{
//       std::cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << std::endl; std::cout.flush();
  if (x < x1 or x > x2) return;
  const double dn = (slope-2 + sqrt((slope-2)*(slope-2) + (8*slope*(x - x1))/step))/(2*slope);
  const int n = static_cast<int>(dn);
  assert(0 <= n && n < nbin);
  fill(nextevt, n, w);
}

SmearedQuadraticHistogram::SmearedQuadraticHistogram(const TString& filename_, const TString& name_,
                                       int nbin_, double x1_, double x2_,
                                       double f, double smear_, double /*param3*/)
  : QuadraticHistogram(filename_, name_, nbin_, x1_, x2_, f), smear(smear_)
{
}

void SmearedQuadraticHistogram::bin(int nextevt, double x, double w)
{
//       std::cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << std::endl; std::cout.flush();
  if (x < x1 or x > x2) return;
  const double dn = (slope-2 + sqrt((slope-2)*(slope-2) + (8*slope*(x - x1))/step))/(2*slope);
  const int n = static_cast<int>(dn);
  assert(0 <= n && n < nbin);

  double dist = dn - n; // [0, 1)
  if (dist < 0.5*smear && n > 0) {
    fill(nextevt, n, 0.5*w);
    fill(nextevt, n-1, 0.5*w);
  } else if ((1. - dist) < 0.5*smear && n < nbin-1) {
    fill(nextevt, n, 0.5*w);
    fill(nextevt, n+1, 0.5*w);
  } else {
    fill(nextevt, n, w);
  }
}
