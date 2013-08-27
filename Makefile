
all:
	./hammer.py -a ham-jets -o 'dummy-%s.tmp' -n 2 -s 1. -p 2 -d -f 'CT10.LHgrid' -t 'CT10.LHgrid' dummy.root
	rm -f dummy-*.tmp

clean:
	rm -f *.d *.o *.so
