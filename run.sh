#!/bin/bash

#CPLUS_INCLUDE_PATH=/home/ivanp/local/LHAPDF5/include:$CPLUS_INCLUDE_PATH
#LIBRARY_PATH=/home/ivanp/local/LHAPDF5/lib:$LIBRARY_PATH
#LD_LIBRARY_PATH=/home/ivanp/local/LHAPDF5/lib:$LD_LIBRARY_PATH

#source /msu/opt/cern/rootSL6/v5.34.18_64/bin/thisroot.sh

#__LHADIR=/home/mondrag3/software/lhapdf
__LHADIR=/home/ivanp/local/LHAPDF5
PATH=$__LHADIR/bin:$PATH
CPLUS_INCLUDE_PATH=$__LHADIR/include:$CPLUS_INCLUDE_PATH
LIBRARY_PATH=$__LHADIR/lib:$LIBRARY_PATH
LD_LIBRARY_PATH=$__LHADIR/lib:$LD_LIBRARY_PATH
#LHAPATH=$__LHADIR/share/lhapdf/PDFsets/

JETS=3
HAM=./hammer.py     # full path to hammer.py
ANA=./ham-higgs.py  # full path to ham-higgs.py
SEED=100_100
FPDF=CT10nlo      # PDF set used in NTuples
TPDF=CT10nlo      # PDF to use in analysis
for RESC in mult_higgs fixed_higgs; do
  for PART in B V I RS; do
    for SCL in 1 0.5 2; do
      FILE=/msu/data/t3work4/luisonig/Ntuples/H${JETS}.0j_GGF_${PART}_4000_pt25.0_eta4.5_r${SEED}.root
      OUT="ATLAS-H$JETS-r$SEED-$TPDF-$RESC-$PART-$SCL"
      $HAM --noapplgrid --noloopsim --analysis $ANA -o "$OUT-%s" \
           -n $JETS -s $SCL -f $FPDF -t $TPDF --rescaler=$RESC $FILE > ${OUT}.log
    done
    wait
  done
done
