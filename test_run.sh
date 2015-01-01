#!/bin/bash

__LHADIR=/home/ivanp/local/LHAPDF5
PATH=$__LHADIR/bin:$PATH
CPLUS_INCLUDE_PATH=$__LHADIR/include:$CPLUS_INCLUDE_PATH
LIBRARY_PATH=$__LHADIR/lib:$LIBRARY_PATH
LD_LIBRARY_PATH=$__LHADIR/lib:$LD_LIBRARY_PATH

./hammer.py --noapplgrid --noloopsim --analysis "./ham-higgs.py" -o "born-%s" \
            -n 3 -s 1 -f CT10nlo -t CT10nlo --rescaler="fixed_higgs" \
            /home/ivanp/work/bh_analysis/data/ten/born_bh.root
