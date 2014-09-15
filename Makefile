
all:
	./hammer.py --noapplgrid -a ham-jets --warmup -o 'dummy-%s.tmp' -n 2 -s 1. -p 2 -f 'CT10.LHgrid' -t 'CT10.LHgrid' dummy.root
	rm -f dummy-*.tmp.hist dummy-*.tmp.root

clean:
	rm -f *.d *.o *.so *_ACLiC_*

root-analysis: root-analysis.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ `python-config --includes` `python-config --libs`

libroot-analysis.so: root-analysis.cpp
	$(CXX) -shared -fPIC $(CXXFLAGS) -o $@ $^ `python-config --includes` `python-config --libs` -lCore
