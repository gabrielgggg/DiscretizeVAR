#
# "gnu" or "intel"
#
tools = intel

ifeq ($(tools),gnu)
  compiler = g++
  flags = -Ofast -fopenmp -pipe -flto -std=c++2a -march=native -mcpu=native
  fWarn = -Wall -pedantic # -g
  libs = -larmadillo 
else ifeq ($(tools),intel)
  compiler = icpc
  flags = -Ofast -parallel -align -march=native -mcpu=native -qopenmp -std=c++20 -ipo # -xHost 
  fWarn = -W -pedantic # -g
  libs = -larmadillo -mkl
else
  $(error Unknown toolchain?)
endif

.PHONY: all clean cleanresults ccc format

all:	main.cpp TimeSeries.o
	$(compiler) $(flags) $(fWarn) -o main *.o main.cpp ${libs}

TimeSeries.o:	TimeSeries.cpp TimeSeries.hpp
	$(compiler) $(flags) $(fWarn) -c TimeSeries.cpp

clean:
	rm -f ./main ./*.o

cleanresults:
	rm -f ./*.h5

ccc:	clean cleanresults

format:	
	clang-format -i ./*.cpp ./*.hpp
