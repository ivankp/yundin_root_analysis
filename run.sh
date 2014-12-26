#!/bin/bash

JETS=3
HAM=hammer.py     # full path to hammer.py
ANA=ham-higgs.py  # full path to ham-higgs.py
SEED=100_100
FPDF=CT10nlo      # PDF set used in NTuples
TPDF=CT10nlo      # PDF to use in analysis
for RESC in mult_higgs fixed_higgs; do
  for PART in B V I RS; do
    for SCL in 1 0.5 2; do
      FILE=H${JETS}.0j_GGF_${PART}_4000_pt25.0_eta4.5_r${SEED}.root
      OUT="ATLAS-H$JETS-r$SEED-$TPDF-$RESC-$PART-$SCL"
      $HAM --noapplgrid --noloopsim --analysis $ANA -o "$OUT-%s" \
           -n $JETS -s $SCL -f $FPDF -t $TPDF --rescaler=$RESC $FILE > ${OUT}.log
    done
    wait
  done
done
