CXX = clang++

run: bin/main.o bin/random_normal.o bin/european_option.o bin/barrier_option.o
	g++ bin/main.o bin/random_normal.o bin/european_option.o bin/barrier_option.o -o run

bin/main.o: main.cpp
	g++ -c main.cpp -o bin/main.o

bin/random_normal.o: random_normal.hpp random_normal.cpp
	g++ -c random_normal.cpp -o bin/random_normal.o

bin/european_option.o: european_option.hpp european_option.cpp
	g++ -c european_option.cpp -o bin/european_option.o

bin/barrier_option.o: barrier_option.hpp barrier_option.cpp
	g++ -c barrier_option.cpp -o bin/barrier_option.o
