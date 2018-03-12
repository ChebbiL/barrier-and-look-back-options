CXX = clang++

run: bin/main.o bin/random_normal.o bin/european_option.o bin/barrier_option.o bin/lookback_option.o
	g++ bin/main.o bin/random_normal.o bin/european_option.o bin/barrier_option.o bin/lookback_option.o -o run -std=c++11
	#g++ -dynamiclib -fPIC -o ../libs/BarrierLookBackOptions.dylib bin/random_normal.o bin/european_option.o bin/barrier_option.o bin/lookback_option.o -std=c++11
	#g++ -shared -o ../libs/libBarrierLookBackOptions.so bin/random_normal.o bin/european_option.o bin/barrier_option.o bin/lookback_option.o
	#ar rvs ../libs/BarrierLookBackOptions.a  bin/random_normal.o bin/european_option.o bin/barrier_option.o bin/lookback_option.o -std=c++11

bin/main.o: main.cpp
	g++ -c main.cpp -o bin/main.o -std=c++11

bin/random_normal.o: random_normal.hpp random_normal.cpp
	g++ -c random_normal.cpp -o bin/random_normal.o -std=c++11

bin/european_option.o: european_option.hpp european_option.cpp
	g++ -c european_option.cpp -o bin/european_option.o -std=c++11

bin/barrier_option.o: barrier_option.hpp barrier_option.cpp
	g++ -c barrier_option.cpp -o bin/barrier_option.o -std=c++11

bin/lookback_option.o: lookback_option.hpp lookback_option.cpp
	g++ -c lookback_option.cpp -o bin/lookback_option.o -std=c++11

clean:
	rm -f bin/*.o
	rm -f ../libs/*.so
	rm -f ../libs/*.dylib
	rm -f ../libs/*.dll
	rm -f run
