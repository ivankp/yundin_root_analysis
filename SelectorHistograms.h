
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

#if defined(__MAKECINT__)
#pragma link C++ class Histogram;
#pragma link C++ class LinearHistogram;
#pragma link C++ class QuadraticHistogram;
#pragma link C++ class SmearedLinearHistogram;
#pragma link C++ class SmearedQuadraticHistogram;
#endif

#endif // SELECTOR_HISTOGRAMS_H
