LIBS="-larmadillo"
CFLAGS="-O1"
INCLUDE="-I/usr/local/include/Minuit2"

test: test.o
	g++ $^ -o $@ $(LIBS) $(CFLAGS)

%.o: %.cc
	g++ -c -o $@ $^ $(INCLUDE) $(CFLAGS)

clean: 
	rm -f *.o
	find . -maxdepth 1 -type f -executable -exec rm {} +
