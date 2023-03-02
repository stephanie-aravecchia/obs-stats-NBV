

all: build/CMakeCache.txt
	(cd build; $(MAKE) all)

build: 
	mkdir -p build

build/CMakeCache.txt: build
	(cd build; cmake ..)

clean:
	-(cd build; $(MAKE) clean)
	rm -rf build


