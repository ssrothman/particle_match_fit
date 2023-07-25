LIBS=-lMinuit2Math -larmadillo -lMinuit2
CFLAGS=-O0 -fPIC
INCLUDE=-I/usr/local/include/Minuit2 -I//home/simon/EECdev/usercode -I/home/simon/miniforge3/envs/EEC/include
OBJECTS=chisqLossFCN.o matcher.o MatchingFilter.o matchingUtil.o ParticleUncertainty.o prefit.o refinePrefit.o

matching.so:$(OBJECTS)
	g++ --shared -o $@ $(OBJECTS)

%.o: %.cc *.h
	g++ -c -o $@ $< $(INCLUDE) $(CFLAGS)
clean: 
	rm -f *.o
	rm -f *.so
	find . -maxdepth 1 -type f -executable -exec rm {} +
