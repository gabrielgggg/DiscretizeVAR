#
# "gnu" or "intel"
#
tools = gnu

ifeq ($(tools),gnu)
  compiler = g++
  flags = -Ofast -fopenmp -pipe -flto -std=c++17
  fWarn = -Wall -pedantic -g
  libs = -larmadillo 
else ifeq ($(tools),intel)
  compiler = icpc
  flags = -O3 -ipo -parallel -align -march=native -mtune=native -mcpu=native -qopenmp
  fWarn = -warn all
  libs = -larmadillo -mkl
else
  $(error Unknown toolchain?)
endif

.PHONY: all clean

all:	main.cpp TimeSeries.o
	$(compiler) $(flags) $(fWarn) -o main *.o main.cpp ${libs}

TimeSeries.o:	TimeSeries.cpp TimeSeries.hpp
	$(compiler) $(flags) $(fWarn) -c TimeSeries.cpp

clean:
	rm -f ./main ./*.o

cleanresults:
	rm -f results/*.tab results/*.bin

cleanmlab:
	rm -f results/*.mat results/*.fig

ccc:	clean cleanresults

format:	
	# clang-format --style=chromium -i ./*.cpp ./*.hpp
	clang-format -i ./*.cpp ./*.hpp
