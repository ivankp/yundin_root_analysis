
#include "SelectorHistograms.h"

// --------------------------------------------------------------------------- //
// Histograms
// --------------------------------------------------------------------------- //

Histogram::Histogram(const TString& filename_, const TString& name_,
                     int nbin_, double x1_, double x2_)
  : filename(filename_), name(name_),
    x1(x1_), x2(x2_), x12(x2-x1),
    nbin(nbin_), prevevt(-1), lastidx(0),
    curidxwgt(nbin), events(nbin),
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
    int idx = curidxwgt[i].first;
    double w = curidxwgt[i].second;
    curidxwgt[i].second = 0.;
    wgt[idx] += w;
    wgt2[idx] += w*w;
    events[idx]++;
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
    if (curidxwgt[i].first == n) {
      break;
    }
  }
  curidxwgt[i].first = n;
  curidxwgt[i].second += w;
  if (i == lastidx) {
    lastidx++;
  }
}

SmearedHistogram::SmearedHistogram(const TString& filename_, const TString& name_,
                                   int nbin_, double x1_, double x2_,
                                   double smear_, double /*param2*/, double /*param3*/)
  : Histogram(filename_, name_, nbin_, x1_, x2_),
    smear(smear_), wgtvec(nbin), xwgtvec(nbin)
{
}

void SmearedHistogram::flush()
{
  if (lastidx > 1) {
    std::sort(curidxwgt.begin(), curidxwgt.begin()+lastidx);
    for (int i=1; i<lastidx; i++) {
      if (curidxwgt[i-1].first + 1 == curidxwgt[i].first) {
        int idx0 = curidxwgt[i-1].first;
        int idx1 = curidxwgt[i].first;
        double s0 = (edge[idx1] - xwgtvec[idx0]/wgtvec[idx0])/bwidth[idx0];
        double s1 = (xwgtvec[idx1]/wgtvec[idx1] - edge[idx1])/bwidth[idx1];
        if (s1 < s0 and s1 < smear) { // swap to merge i to i-1
          std::swap(curidxwgt[i-1].first, curidxwgt[i].first);
          std::swap(curidxwgt[i-1].second, curidxwgt[i].second);
          std::swap(idx0, idx1);
          std::swap(s0, s1);
        }
        if (s0 < s1 and s0 < smear) {
          curidxwgt[i].second += curidxwgt[i-1].second;
          wgtvec[idx1] += wgtvec[idx0];
          xwgtvec[idx1] += xwgtvec[idx0];
          curidxwgt[i-1].second = 0;
        }
      }
    }
  }
  wgtvec.assign(wgtvec.size(), 0.);
  xwgtvec.assign(wgtvec.size(), 0.);
  Histogram::flush();
}

inline
void SmearedHistogram::fill(int evt, int n, double w, double x)
{
  Histogram::fill(evt, n, w);
  wgtvec[n] += abs(w);
  xwgtvec[n] += x*abs(w);
}

// specific histograms

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
//  assert(0 <= n && n < nbin);
  if (n<0 or n>=nbin) {
    std::cout << "WARN: " << name << " 0 <= " << n << " " << nbin << " x=" << x << " w=" << w << std::endl;
    return;
  }
  fill(nextevt, n, w);
}

SmearedLinearHistogram::SmearedLinearHistogram(const TString& filename_, const TString& name_,
                                       int nbin_, double x1_, double x2_,
                                       double smear_, double /*param2*/, double /*param3*/)
  : SmearedHistogram(filename_, name_, nbin_, x1_, x2_, smear_),
    step(x12/nbin)
{
  for (int i=0; i<nbin; i++) {
    bwidth[i] = step;
  }
  setedges();
}

void SmearedLinearHistogram::bin(int nextevt, double x, double w)
{
//       std::cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << std::endl; std::cout.flush();
  if (x < x1 or x > x2) return;
  int n = static_cast<int>(nbin*(x-x1)/x12);
  assert(0 <= n && n < nbin);
  fill(nextevt, n, w, x);
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
  : SmearedHistogram(filename_, name_, nbin_, x1_, x2_, smear_),
    step(2*x12/((1 + f)*nbin)),
    slope((f - 1)/(nbin - 1))
{
  for (int i=0; i<nbin; i++) {
    bwidth[i] = step*(1 + i*slope);
  }
  setedges();
}

void SmearedQuadraticHistogram::bin(int nextevt, double x, double w)
{
//       std::cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << std::endl; std::cout.flush();
  if (x < x1 or x > x2) return;
  const double dn = (slope-2 + sqrt((slope-2)*(slope-2) + (8*slope*(x - x1))/step))/(2*slope);
  const int n = static_cast<int>(dn);
  assert(0 <= n && n < nbin);
  fill(nextevt, n, w, x);
}

// APPLgrid histograms

double Grid::aparam = 5.;
std::string Grid::pdf_function = "ntuplejets";
ntuple_pdf* Grid::pdf_object = 0;
bool Grid::pdfWeight = false;
int Grid::born_alphapower = 0;
int Grid::nloops = 0;

GridOpts Grid::def_opts = GridOpts(50, 0., 16e6, 5,
                                   50, 1e-6, 1., 5);

bool Grid::valid = false;

void Grid::static_init()
{
  if (valid) {
    return;
  }
  valid = true;

  appl::grid::transformvar(aparam);
  if (pdf_function == "ntuplejets") {
    pdf_object = new ntuplejets_pdf();
  } else if (pdf_function == "ntuplephjets") {
    pdf_object = new ntuplephjets_pdf();
  } else if (pdf_function == "ntupleall") {
    pdf_object = new ntupleall_pdf();
  } else {
    std::cout << "Unknown pdf_function " << pdf_function << std::endl;
    exit(0);
  }
}

Grid::Grid(const std::string& name)
  : warmup(false), filename(name)
{
  init();

  m_grid = new appl::grid(filename);

  if (m_grid->isOptimised()) {
    delete m_grid;
    std::cout << "Grid is aready optimised. Quitting ..." << std::endl;
    exit(0);
  }

  // reseting reference histgram
  TH1D* htemp = m_grid->getReference();
  for (int i = 0; i <= htemp->GetNbinsX()+1; i++) {
    htemp->SetBinContent(i, 0);
  }

  m_grid->optimise();
}

Grid::Grid(const std::string& name,
           const std::vector<double>& edges, const GridOpts& opts)
  : warmup(true), filename(name)
{
  init();

  m_grid = new appl::grid(edges,
                          opts.Q2bins, opts.Q2low, opts.Q2high, opts.Q2order,
                          opts.Xbins, opts.Xlow, opts.Xhigh, opts.Xorder,
                          pdf_function, born_alphapower, nloops);
  m_grid->reweight(pdfWeight);
}

Grid::~Grid()
{
  delete m_grid;
}

void Grid::init()
{
  static_init();
  weights.resize(pdf_object->Nproc());
  if (warmup) {
    weights.assign(weights.size(), 1.);
  }
}

void Grid::fill(int /*id*/, int id1, int id2,
                double x1, double x2, double Q,
                const double* fA, const double* fB,
                double obs, double g_w, double h_w, int order)
{
  if (x1 == 1. or x2 == 1.) {
    return;
  }
  if (not warmup) {
    const int ch = pdf_object->channel(id1, id2);
    weights.assign(weights.size(), 0.);
    weights[ch] = pdf_object->reweight(g_w, ch, id1, id2, fA, fB);
  }
  m_grid->fill_grid(x1, x2, Q*Q, obs, weights.data(), order);

  double binwidth = m_grid->deltaobs(m_grid->obsbin(obs));
  m_grid->getReference()->Fill(obs, h_w/binwidth);
}

void Grid::write(double count)
{
  if (not warmup) {
    (*m_grid) *= 1./count;
  }
  m_grid->Write(filename);
}
