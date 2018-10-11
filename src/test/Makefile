CXX = /usr/bin/c++

CXXFLAGS = -std=c++11 -O3 -g -fmessage-length=0 

OALignUtility = ../../bin/utility/AlignUtility.o 

TAlignUtility = ../../bin/AlignUtility

$(OAlignUtility): ../tr/GlobAlignE.h
	$(CXX) $(CXXFlAGS) -c ../tr/AlignUtility.cpp -o $(OAlignUtility)

util: $(OAlignUtility)