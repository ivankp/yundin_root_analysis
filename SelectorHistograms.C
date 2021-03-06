
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

void Histogram::setedges()
{
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

// ------------------------- SmearedHistogram base ----------------------------

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

// --------------------------------------------------------------------------- //
// Specific histograms
// --------------------------------------------------------------------------- //

// ----------------------- ListHistogram -------------------------

ListHistogram::ListHistogram(const TString& filename_, const TString& name_,
                             const std::vector<double>& edge_)
  : Histogram(filename_, name_, edge_.size()-1, edge_.front(), edge_.back())
{
  for (int i=0; i<nbin; i++) {
    bwidth[i] = edge_[i+1] - edge_[i];
    assert(bwidth[i] > 0.);
  }
  edge.assign(edge_.begin(), edge_.end());
}

int ListHistogram::getbinid(double x)
{
  int lo = 0;
  int hi = nbin - 1;
  int n = (hi - lo)/2;
  // do a binary search of the bin index
  while (true) {
    if (edge[n] > x) {
      hi = n;
      n -= (n - lo + 1)/2;
      continue;
    } else {
      lo = n;
    }
    if (edge[n+1] <= x) {
      lo = n;
      n += (hi - n + 1)/2;
      continue;
    } else {
      hi = n;
    }
    return n;
  }
}

void ListHistogram::bin(int nextevt, double x, double w)
{
//       std::cout << name << ": E(" << evt << ") LE (" << lastevt << ") LI(" << lastidx << ")" << std::endl; std::cout.flush();
  if (x < x1 or x >= x2) return;
  int n = getbinid(x);
  assert(0 <= n && n < nbin);
  fill(nextevt, n, w);}

// ---------------------- LinearHistogram ------------------------

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
  if (x < x1 or x >= x2) return;
  int n = static_cast<int>(nbin*(x-x1)/x12);
  assert(0 <= n && n < nbin);
  fill(nextevt, n, w);
}

// ---------------------- SmearedLinearHistogram --------------------------

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
  if (x < x1 or x >= x2) return;
  int n = static_cast<int>(nbin*(x-x1)/x12);
  assert(0 <= n && n < nbin);
  fill(nextevt, n, w, x);
}

// -------------------------- QuadraticHistogram --------------------------

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
  if (x < x1 or x >= x2) return;
  const double dn = (slope-2 + sqrt((slope-2)*(slope-2) + (8*slope*(x - x1))/step))/(2*slope);
  const int n = static_cast<int>(dn);
  if (n < 0 or n >= nbin) {
    std::cout << "WARNING " << 0 << " <= " << n << " < " << nbin
              << " while " << x1 << " <= " << x << " < " << x2 << std::endl;
    return;
  }
  fill(nextevt, n, w);
}

// ---------------------- SmearedQuadraticHistogram -----------------------

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
  if (x < x1 or x >= x2) return;
  const double dn = (slope-2 + sqrt((slope-2)*(slope-2) + (8*slope*(x - x1))/step))/(2*slope);
  const int n = static_cast<int>(dn);
  if (n < 0 or n >= nbin) {
    std::cout << "WARNING " << 0 << " <= " << n << " < " << nbin
              << " while " << x1 << " <= " << x << " < " << x2 << std::endl;
    return;
  }
  fill(nextevt, n, w, x);
}

// --------------------------------------------------------------------------- //
// 2D Histograms
// --------------------------------------------------------------------------- //

// ---------------------- Histogram2D ------------------------

template <typename Hist1D>
Histogram2D<Hist1D>::Histogram2D(const TString& filename_, const TString& name_,
                                 int nbinx_, double x1_, double x2_,
                                 int nbiny_, double y1_, double y2_)
  : name(name_), y1(y1_), y2(y2_), y12(y2-y1),
    nbin(nbiny_), bwidth(nbin), edge(nbin+1)
{
  hists1d.reserve(nbin);
  for (int i=0; i<nbin; i++) {
    hists1d.push_back(Hist1D(filename_, name_ + TString::Format("_%d", i), nbinx_, x1_, x2_));
  }
  assert(nbin > 0);
}

template <typename Hist1D>
TString Histogram2D<Hist1D>::getFile() const
{
  return hists1d[0].getFile();
}

template <typename Hist1D>
void Histogram2D<Hist1D>::print(std::ostream& stream, const TString& runname,
                                double count, bool unweight)
{
  stream << "# BEGIN HISTOGRAM2D /" << runname << "/" << name << std::endl;

  stream << "#";
  for (unsigned i=0; i<edge.size(); i++) {
    stream << " " << edge[i];
  }
  stream << std::endl;

  for (unsigned i=0; i<hists1d.size(); i++) {
    hists1d[i].print(stream, runname, count, unweight);
  }
  stream << "# END HISTOGRAM2D" << std::endl << std::endl << std::endl;
}

template <typename Hist1D>
void Histogram2D<Hist1D>::setedges()
{
  edge[0] = y1;
  for (int i=0; i<nbin; i++) {
    edge[i+1] = edge[i] + bwidth[i];
  }
}

template <typename Hist1D>
inline
void Histogram2D<Hist1D>::fill(double x, int evt, int n, double w)
{
  assert(0 <= n && n < nbin);
  hists1d[n].bin(evt, x, w);
}

// ---------------------- LinearHistogram2D ------------------------

LinearHistogram2D::LinearHistogram2D(const TString& filename_, const TString& name_,
                                     int nbinx_, double x1_, double x2_,
                                     int nbiny_, double y1_, double y2_,
                                     double /*param1*/, double /*param2*/, double /*param3*/)
  : BaseClass(filename_, name_, nbinx_, x1_, x2_, nbiny_, y1_, y2_), step(y12/nbin)
{
  for (int i=0; i<nbin; i++) {
    bwidth[i] = step;
  }
  setedges();
}

void LinearHistogram2D::bin2d(int nextevt, double x, double y, double w)
{
  if (y < y1 or y >= y2) return;
  int n = static_cast<int>(nbin*(y-y1)/y12);
  assert(0 <= n && n < nbin);
  fill(x, nextevt, n, w);
}


// --------------------------------------------------------------------------- //
// APPLgrid histograms
// --------------------------------------------------------------------------- //

double Grid::aparam = 5.;
std::string Grid::pdf_function = "ntuplejets";
bool Grid::pdfWeight = false;
int Grid::born_alphaspower = -1;
int Grid::nloops = 0;

GridOpts Grid::def_opts = GridOpts(50, 0., 16e6, 5,
                                   50, 1e-6, 1., 5);
bool Grid::valid = false;

#ifndef DISABLE_APPLGRID

ntuple_pdf* Grid::pdf_object = 0;


void Grid::static_init()
{
  if (valid) {
    return;
  }
  valid = true;

  if (born_alphaspower < 0) {
    std::cout << "Invalid born_alphaspower = " << born_alphaspower
              << " Check your file for typos" << std::endl;
    exit(0);
  }

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
                          pdf_function, born_alphaspower, nloops);
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

#endif // DISABLE_APPLGRID
// --------------------------------------------------------------------------- //
