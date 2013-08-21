
#ifndef SELECTOR_HISTOGRAMS_H
#define SELECTOR_HISTOGRAMS_H

class Histogram
{
  public:
    Histogram(const TString& filename_, const TString& name_,
              int nbin_, double x1_, double x2_);
    virtual ~Histogram() {};

    virtual void bin(int nextevt, double x, double w) = 0;
    void print(std::ostream& stream, const TString& runname, double count, bool unweight=true);

    TString getFile() const;

  protected:
    void flush();
    void setedges();
    void fill(int evt, int n, double w);

    TString filename;
    TString name;

    double x1, x2, x12;
    const int nbin;
    int prevevt, lastidx;
    std::vector<int> curidx;
    std::vector<int> events;
    std::vector<double> curwgt;
    std::vector<double> wgt;
    std::vector<double> wgt2;
    std::vector<double> bwidth;
    std::vector<double> edge;
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

class SmearedLinearHistogram : public LinearHistogram
{
  public:
    SmearedLinearHistogram(const TString& filename_, const TString& name_,
                    int nbin_, double x1_, double x2_,
                    double smear_=1., double param2=0, double param3=0);
    void bin(int nextevt, double x, double w);

  protected:
    const double smear;
};

class SmearedQuadraticHistogram : public QuadraticHistogram
{
  public:
    SmearedQuadraticHistogram(const TString& filename_, const TString& name_,
                       int nbin_, double x1_, double x2_,
                       double f, double smear_=1., double param3=0);
    void bin(int nextevt, double x, double w);

  protected:
    const double smear;
};

#if defined(__MAKECINT__)
#pragma link C++ class Histogram;
#pragma link C++ class LinearHistogram;
#pragma link C++ class QuadraticHistogram;
#pragma link C++ class SmearedLinearHistogram;
#pragma link C++ class SmearedQuadraticHistogram;
#endif

#endif // SELECTOR_HISTOGRAMS_H
