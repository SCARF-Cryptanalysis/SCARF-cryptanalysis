CFLAGS   = -Wall -Wextra -O3 -march=native -g -fopenmp
CXXFLAGS = -Wall -Wextra -O3 -march=native -g -fopenmp -std=c++20
LDFLAGS  = -fopenmp
LDLIBS   = -lm


all: 6r_distinguisher 8r_keyrecovery

6r_distinguisher: aes-rng.o 6r_distinguisher.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	rm -f 6r_distinguisher 8r_keyrecovery *.o

.PHONY: all clean
