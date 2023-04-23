LIBS=-lMinuit2Math -larmadillo -lMinuit2
CFLAGS=-O0
INCLUDE=-I/usr/local/include/Minuit2

test: test.o toyjets/gaus.o toyjets/gen.o
	g++ $^ -o $@ $(LIBS) $(CFLAGS)

%.o: %.cc
	g++ -c -o $@ $^ $(INCLUDE) $(CFLAGS)

clean: 
	rm -f *.o
	find . -maxdepth 1 -type f -executable -exec rm {} +
