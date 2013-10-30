
#ifndef SELECTOR_HISTOGRAMS_H
#define SELECTOR_HISTOGRAMS_H

#include <utility>

class Histogram
{
  public:
    Histogram(const TString& filename_, const TString& name_,
              int nbin_, double x1_, double x2_);
    virtual ~Histogram() {};

    virtual void bin(int nextevt, double x, double w) = 0;
    virtual void print(std::ostream& stream, const TString& runname, double count, bool unweight=true);

    TString getFile() const;

  protected:
    typedef std::pair<int, double> TIdxWgt;

    virtual void flush();
    void setedges();
    void fill(int evt, int n, double w);

    TString filename;
    TString name;

    double x1, x2, x12;
    const int nbin;
    int prevevt, lastidx;
    std::vector<TIdxWgt> curidxwgt;
    std::vector<int> events;
    std::vector<double> wgt;
    std::vector<double> wgt2;
    std::vector<double> bwidth;
    std::vector<double> edge;
};


class SmearedHistogram : public Histogram
{
  public:
    SmearedHistogram(const TString& filename_, const TString& name_,
                     int nbin_, double x1_, double x2_,
                     double smear_=1., double param2=0, double param3=0);

    virtual void bin(int nextevt, double x, double w) = 0;

  protected:
    void flush();
    void fill(int evt, int n, double w, double x);

    const double smear;

    std::vector<double> wgtvec;
    std::vector<double> xwgtvec;
};


class ListHistogram : public Histogram
{
  public:
    ListHistogram(const TString& filename_, const TString& name_,
                  const std::vector<double>& edge_);
    void bin(int nextevt, double x, double w);

  protected:
    int getbinid(double x);
};


class LinearHistogram : public Histogram
{
  public:
    LinearHistogram(const TString& filename_, const TString& name_,
                    int nbin_, double x1_, double x2_,
                    double param1=0, double param2=0, double param3=0);
    void bin(int nextevt, double x, double w);

  protected:
    const double step;
};


class QuadraticHistogram : public Histogram
{
  public:
    QuadraticHistogram(const TString& filename_, const TString& name_,
                       int nbin_, double x1_, double x2_,
                       double f, double param2=0, double param3=0);
    void bin(int nextevt, double x, double w);

  protected:
    const double step;
    const double slope;
};

class SmearedLinearHistogram : public SmearedHistogram
{
  public:
    SmearedLinearHistogram(const TString& filename_, const TString& name_,
                       int nbin_, double x1_, double x2_,
                       double smear_=1., double param2=0, double param3=0);
    void bin(int nextevt, double x, double w);

  protected:
    const double step;
};

class SmearedQuadraticHistogram : public SmearedHistogram
{
  public:
    SmearedQuadraticHistogram(const TString& filename_, const TString& name_,
                       int nbin_, double x1_, double x2_,
                       double f, double smear_=1., double param3=0);
    void bin(int nextevt, double x, double w);

  protected:
    const double step;
    const double slope;
};

#include <appl_grid/appl_grid.h>
#include "ntuple_pdf.h"

struct GridOpts
{
  GridOpts(int qbin, double qlo, double qhi, int qord,
           int xbin, double xlo, double xhi, int xord)
    : Q2bins(qbin), Q2low(qlo), Q2high(qhi), Q2order(qord),
      Xbins(xbin), Xlow(xlo), Xhigh(xhi), Xorder(xord)
  { }

  int Q2bins;
  double Q2low, Q2high;
  int Q2order;

  int Xbins;
  double  Xlow, Xhigh;
  int Xorder;
};

class Grid
{
  public:
    static double aparam;
    static std::string pdf_function;
    static ntuple_pdf* pdf_object;
    static bool pdfWeight;

    static int born_alphapower;
    static int nloops;

    static GridOpts def_opts;

    static bool valid;
    static void static_init();

    static Grid* capture(Grid* p) { return p; }

    Grid(const std::string& name);
    Grid(const std::string& name,
         const std::vector<double>& edges, const GridOpts& opts = Grid::def_opts);
    void init();
    ~Grid();

    void fill(int id, int id1, int id2,
              double x1, double x2, double Q2,
              const double* fA, const double* fB,
              double obs, double g_w, double h_w, int order);

    bool isWarmup() { return warmup; }
    void setFilename(const std::string& filename_) { filename = filename_; }
    void write(double count);

    bool warmup;
    appl::grid* m_grid;
    std::string filename;
    std::vector<double> weights;
};

#if defined(__MAKECINT__)
#pragma link C++ class Histogram;
#pragma link C++ class ListHistogram;
#pragma link C++ class LinearHistogram;
#pragma link C++ class QuadraticHistogram;
#pragma link C++ class SmearedHistogram;
#pragma link C++ class SmearedLinearHistogram;
#pragma link C++ class SmearedQuadraticHistogram;
#pragma link C++ class GridOpts;
#pragma link C++ class Grid;
#endif

#endif // SELECTOR_HISTOGRAMS_H
