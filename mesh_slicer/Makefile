
all: build/CMakeCache.txt
	(cd build; make -j4)

build/CMakeCache.txt: build
	(cd build; cmake ..)

build:
	mkdir -p build

clean:
	rm -rf build bin




