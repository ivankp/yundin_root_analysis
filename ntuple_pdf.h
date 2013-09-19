
#ifndef NTUPLEJETS_PDF_H
#define NTUPLEJETS_PDF_H

#include "appl_grid/appl_pdf.h"
#include <cstdlib>
#include <cassert>

class ntuple_pdf : public appl::appl_pdf
{
  public:
    ntuple_pdf(const std::string& name)
      : appl::appl_pdf(name)
    { }

    virtual int channel(int id1, int id2) = 0;
    virtual double reweight(double w, int ch, int id1, int id2,
                            const double* fA, const double* fB) = 0;
};

class ntuplejets_pdf : public ntuple_pdf
{
  public:
    ntuplejets_pdf()
      : ntuple_pdf("ntuplejets")
    {
      m_Nproc = TOT_CHANNELS;
    }

    enum CHANNELS {
      chGG,
      chQG,
      chGQ,
      chQR,
      chQQ,
      chQQx,
      chQRx,
      TOT_CHANNELS
    };

    void evaluate(const double* fA, const double* fB, double* H)
    {
      fA += 6;
      fB += 6;

      double pdfG1 = fA[0];
      double pdfQ1 = 0;
      double pdfQ1x = 0;

      double pdfG2 = fB[0];
      double pdfQ2 = 0;
      double pdfQ2x = 0;

      double pdfD = 0;
      double pdfDx = 0;

      for (int i=1; i<=5; i++) {
        pdfQ1 += fA[i];
        pdfQ1x += fA[-i];
        pdfQ2 += fB[i];
        pdfQ2x += fB[-i];
        pdfD += fA[i]*fB[i] + fA[-i]*fB[-i];
        pdfDx += fA[-i]*fB[i] + fA[i]*fB[-i];
      }

      H[chGG] = pdfG1*pdfG2;
      H[chQG] = (pdfQ1 + pdfQ1x)*pdfG2;
      H[chGQ] = pdfG1*(pdfQ2 + pdfQ2x);
      H[chQR] = (pdfQ1*pdfQ2 + pdfQ1x*pdfQ2x - pdfD);
      H[chQQ] = pdfD;
      H[chQQx] = pdfDx;
      H[chQRx] = (pdfQ1*pdfQ2x + pdfQ1x*pdfQ2 - pdfDx);
    }

    int channel(int id1, int id2)
    {
      id1 = id1 == 21 ? 0. : id1;
      id2 = id2 == 21 ? 0. : id2;

      if (id1 == 0 and id2 == 0) {
        return chGG;
      }
      if (id2 == 0) {
        return chQG;
      }
      if (id1 == 0) {
        return chGQ;
      }
      if (id1 != id2 and id1*id2 > 0) {
        return chQR;
      }
      if (id1 == id2) {
        return chQQ;
      }
      if (abs(id1) == abs(id2)) {
        return chQQx;
      }
      if (abs(id1) != abs(id2)) {
        return chQRx;
      }
      std::cout << "Unknown channel " << id1 << " " << id2 << std::endl;
      exit(0);
    }

    double reweight(double w, int ch, int id1, int id2,
                    const double* fA, const double* fB)
    {
      return w;
    }
};

class ntuplephjets_pdf : public ntuple_pdf
{
  public:
    ntuplephjets_pdf()
      : ntuple_pdf("ntuplephjets")
    {
      m_Nproc = TOT_CHANNELS;
    }

    enum CHANNELS {
      chGG,
      chQdG,
      chQuG,
      chGQd,
      chGQu,
      chQdRd,
      chQuRu,
      chQdRu,
      chQuRd,
      chQdQd,
      chQuQu,
      chQdQdx,
      chQuQux,
      chQdRdx,
      chQuRux,
      chQdRux,
      chQuRdx,
      TOT_CHANNELS
    };

    void evaluate(const double* fA, const double* fB, double* H)
    {
      fA += 6;
      fB += 6;

      double pdfG1 = fA[0];
      double pdfQ1d = 0;
      double pdfQ1u = 0;
      double pdfQ1dx = 0;
      double pdfQ1ux = 0;

      double pdfG2 = fB[0];
      double pdfQ2d = 0;
      double pdfQ2u = 0;
      double pdfQ2dx = 0;
      double pdfQ2ux = 0;

      double pdfDdd = 0;
      double pdfDuu = 0;
      double pdfDddx = 0;
      double pdfDuux = 0;

      for (int i=1; i<=6; i+=2) {
        const int id = i;
        const int iu = i+1;
        pdfQ1d += fA[id];
        pdfQ1u += fA[iu];
        pdfQ1dx += fA[-id];
        pdfQ1ux += fA[-iu];
        pdfQ2d += fB[id];
        pdfQ2u += fB[iu];
        pdfQ2dx += fB[-id];
        pdfQ2ux += fB[-iu];
        pdfDdd += fA[id]*fB[id] + fA[-id]*fB[-id];
        pdfDuu += fA[iu]*fB[iu] + fA[-iu]*fB[-iu];
        pdfDddx += fA[id]*fB[-id] + fA[-id]*fB[id];
        pdfDuux += fA[iu]*fB[-iu] + fA[-iu]*fB[iu];
      }

      H[chGG] = pdfG1*pdfG2;
      H[chQdG] = (pdfQ1d + pdfQ1dx)*pdfG2;
      H[chQuG] = (pdfQ1u + pdfQ1ux)*pdfG2;
      H[chGQd] = pdfG1*(pdfQ2d + pdfQ2dx);
      H[chGQu] = pdfG1*(pdfQ2u + pdfQ2ux);
      H[chQdRd] = pdfQ1d*pdfQ2d + pdfQ1dx*pdfQ2dx - pdfDdd;
      H[chQuRu] = pdfQ1u*pdfQ2u + pdfQ1ux*pdfQ2ux - pdfDuu;
      H[chQdRu] = pdfQ1d*pdfQ2u + pdfQ1dx*pdfQ2ux;
      H[chQuRd] = pdfQ1u*pdfQ2d + pdfQ1ux*pdfQ2dx;
      H[chQdQd] = pdfDdd;
      H[chQuQu] = pdfDuu;
      H[chQdQdx] = pdfDddx;
      H[chQuQux] = pdfDuux;
      H[chQdRdx] = pdfQ1d*pdfQ2dx + pdfQ1dx*pdfQ2d - pdfDddx;
      H[chQuRux] = pdfQ1u*pdfQ2ux + pdfQ1ux*pdfQ2u - pdfDuux;
      H[chQdRux] = pdfQ1d*pdfQ2ux + pdfQ1dx*pdfQ2u;
      H[chQuRdx] = pdfQ1u*pdfQ2dx + pdfQ1ux*pdfQ2d;
    }

    int channel(int id1, int id2)
    {
      id1 = id1 == 21 ? 0. : id1;
      id2 = id2 == 21 ? 0. : id2;

      bool down1 = abs(id1) % 2 == 1;
      bool down2 = abs(id2) % 2 == 1;

      if (id1 == 0 and id2 == 0) {
        return chGG;
      }
      if (id2 == 0) {
        return down1 ? chQdG : chQuG;
      }
      if (id1 == 0) {
        return down2 ? chGQd : chGQu;
      }
      if (id1 != id2 and id1*id2 > 0) {
        if (down1 == down2) {
          return down1 ? chQdRd : chQuRu;
        } else {
          return down1 ? chQdRu : chQuRd;
        }
      }
      if (id1 == id2) {
        return down1 ? chQdQd : chQuQu;
      }
      if (abs(id1) == abs(id2)) {
        return down1 ? chQdQdx : chQuQux;
      }
      if (abs(id1) != abs(id2)) {
        if (down1 == down2) {
          return down1 ? chQdRdx : chQuRux;
        } else {
          return down1 ? chQdRux : chQuRdx;
        }
      }
      std::cout << "Unknown channel " << id1 << " " << id2 << std::endl;
      exit(0);
    }

    double reweight(double w, int ch, int id1, int id2,
                    const double* fA, const double* fB)
    {
      fA += 6;
      fB += 6;
      double factor = fA[id1]*fB[id2];
      switch (ch) {
        case chGG:
          factor = 1.;
          break;
        case chQdG:
          factor /= (fA[5] + fA[3] + fA[1] + fA[-5] + fA[-3] + fA[-1])*fB[0];
          break;
        case chQuG:
          factor /= (fA[4] + fA[2] + fA[-4] + fA[-2])*fB[0];
          break;
        case chGQd:
          factor /= fA[0]*(fB[5] + fB[3] + fB[1] + fB[-5] + fB[-3] + fB[-1]);
          break;
        case chGQu:
          factor /= fA[0]*(fB[4] + fB[2] + fB[-4] + fB[-2]);
          break;
        case chQdRd:
          factor /= (fA[5]*fB[3] + fA[5]*fB[1] +
                     fA[3]*fB[5] + fA[3]*fB[1] +
                     fA[1]*fB[5] + fA[1]*fB[3] +
                     fA[-5]*fB[-3] + fA[-5]*fB[-1] +
                     fA[-3]*fB[-5] + fA[-3]*fB[-1] +
                     fA[-1]*fB[-5] + fA[-1]*fB[-3]);
          break;
        case chQuRu:
          factor /= (fA[4]*fB[2] + fA[2]*fB[4] +
                     fA[-4]*fB[-2] + fA[-2]*fB[-4]);
          break;
        case chQdRu:
          factor /= (fA[5]*fB[4] + fA[5]*fB[2] +
                     fA[3]*fB[4] + fA[3]*fB[2] +
                     fA[1]*fB[4] + fA[1]*fB[2] +
                     fA[-5]*fB[-4] + fA[-5]*fB[-2] +
                     fA[-3]*fB[-4] + fA[-3]*fB[-2] +
                     fA[-1]*fB[-4] + fA[-1]*fB[-2]);
          break;
        case chQuRd:
          factor /= (
                     fA[4]*fB[5] + fA[2]*fB[5] +
                     fA[4]*fB[3] + fA[2]*fB[3] +
                     fA[4]*fB[1] + fA[2]*fB[1] +
                     fA[-4]*fB[-5] + fA[-2]*fB[-5] +
                     fA[-4]*fB[-3] + fA[-2]*fB[-3] +
                     fA[-4]*fB[-1] + fA[-2]*fB[-1]);
          break;
        case chQdQd:
          factor /= (fA[5]*fB[5] + fA[3]*fB[3] + fA[1]*fB[1] +
                     fA[-5]*fB[-5] + fA[-3]*fB[-3] + fA[-1]*fB[-1]);
          break;
        case chQuQu:
          factor /= (fA[4]*fB[4] + fA[2]*fB[2] +
                     fA[-4]*fB[-4] + fA[-2]*fB[-2]);
          break;
        case chQdQdx:
          factor /= (fA[5]*fB[-5] + fA[3]*fB[-3] + fA[1]*fB[-1] +
                     fA[-5]*fB[5] + fA[-3]*fB[3] + fA[-1]*fB[1]);
          break;
        case chQuQux:
          factor /= (fA[4]*fB[-4] + fA[2]*fB[-2] +
                     fA[-4]*fB[4] + fA[-2]*fB[2]);
          break;
        case chQdRdx:
          factor /= (fA[5]*fB[-3] + fA[5]*fB[-1] +
                     fA[3]*fB[-5] + fA[3]*fB[-1] +
                     fA[1]*fB[-5] + fA[1]*fB[-3] +
                     fA[-5]*fB[3] + fA[-5]*fB[1] +
                     fA[-3]*fB[5] + fA[-3]*fB[1] +
                     fA[-1]*fB[5] + fA[-1]*fB[3]);
          break;
        case chQuRux:
          factor /= (fA[4]*fB[-2] + fA[2]*fB[-4] +
                     fA[-4]*fB[2] + fA[-2]*fB[4]);
          break;
        case chQdRux:
          factor /= (fA[5]*fB[-4] + fA[5]*fB[-2] +
                     fA[3]*fB[-4] + fA[3]*fB[-2] +
                     fA[1]*fB[-4] + fA[1]*fB[-2] +
                     fA[-5]*fB[4] + fA[-5]*fB[2] +
                     fA[-3]*fB[4] + fA[-3]*fB[2] +
                     fA[-1]*fB[4] + fA[-1]*fB[2]);
          break;
        case chQuRdx:
          factor /= (fA[4]*fB[-5] + fA[2]*fB[-5] +
                     fA[4]*fB[-3] + fA[2]*fB[-3] +
                     fA[4]*fB[-1] + fA[2]*fB[-1] +
                     fA[-4]*fB[5] + fA[-2]*fB[5] +
                     fA[-4]*fB[3] + fA[-2]*fB[3] +
                     fA[-4]*fB[1] + fA[-2]*fB[1]);
          break;
      }
      return w*factor;
    }
};

class ntupleall_pdf : public ntuple_pdf
{
  public:
    ntupleall_pdf()
      : ntuple_pdf("ntupleall")
    {
      m_Nproc = 121;
    }

    void evaluate(const double* fA, const double* fB, double* H)
    {
      fA += 6;
      fB += 6;

      int h = 0;
      for (int ia=-5; ia<=5; ia++) {
        for (int ib=-5; ib<=5; ib++) {
          H[h++] = fA[ia]*fB[ib];
        }
      }
    }

    int channel(int id1, int id2)
    {
      id1 = id1 == 21 ? 0. : id1;
      id2 = id2 == 21 ? 0. : id2;

      return 11*(id1 + 5) + (id2 + 5);
    }

    double reweight(double w, int ch, int id1, int id2,
                    const double* fA, const double* fB)
    {
      return w;
    }
};

#endif // NTUPLEJETS_PDF_H
