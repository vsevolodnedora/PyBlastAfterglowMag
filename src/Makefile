# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
# CFLAGS  = `pkg-config --cflags hdf5-serial` src/main.cpp `pkg-config --libs hdf5-serial` -lhdf5_cpp -O2 -std=c++17 -Wall -g -lstdc++fs
# CFLAGS = -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_hl -lhdf5 -lhdf5_cpp -v

#main: pba.out
#	g++ `pkg-config --cflags hdf5-serial` main.cpp `pkg-config --libs hdf5-serial` -lhdf5_cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp
	
	
main: pba.out
	g++ `pkg-config --cflags hdf5-serial` main.cpp `pkg-config --libs hdf5-serial` -lhdf5_cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp

alternative: pba.out
	g++ main.cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_hl -lhdf5 -lhdf5_cpp -v
