
all:
	./hammer.py -a ham-jets --warmup -o 'dummy-%s.tmp' -n 2 -s 1. -p 2 -f 'CT10.LHgrid' -t 'CT10.LHgrid' dummy.root
	rm -f dummy-*.tmp.hist dummy-*.tmp.root

clean:
	rm -f *.d *.o *.so
