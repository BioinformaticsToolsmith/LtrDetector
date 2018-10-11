# CXX = /usr/bin/c++
CXX = /usr/bin/g++

#CXXFLAGS = -std=c++11 -O3 -g -fmessage-length=0 #-Wall
CXXFLAGS = -std=c++11 -g -fmessage-length=0 #-Wall
#
# Objects
#

ORed  = ../bin/Red.o
JProj = ../bin/Joseph.o
OTestTr = ../bin/Tester.o

# Exception
OInvalidInputException  = ../bin/exception/InvalidInputException.o
OInvalidStateException  = ../bin/exception/InvalidStateException.o
OFileDoesNotExistException  = ../bin/exception/FileDoesNotExistException.o
OInvalidOrderOfOperationsException  = ../bin/exception/InvalidOrderOfOperationsException.o
OInvalidScoreException  = ../bin/exception/InvalidScoreException.o
OInvalidOperationException = ../bin/exception/InvalidOperationException.o

# Utility
OUtil = ../bin/utility/Util.o
OLocation = ../bin/utility/Location.o
OEmptyLocation = ../bin/utility/EmptyLocation.o
OLCS = ../bin/utility/LCS.o
OLCSLen = ../bin/utility/LCSLen.o
OLCSubStr = ../bin/utility/LCSubStr.o
OTSD = ../bin/utility/TSD.o
OEmptyTSD = ../bin/utility/EmptyTSD.o
OTailFinder = ../bin/utility/TailFinder.o
OTail = ../bin/utility/Tail.o
OITail = ../bin/utility/ITail.o
OEmptyTail = ../bin/utility/EmptyTail.o
OITSD = ../bin/utility/ITSD.o



# Non TR
OChromosome  = ../bin/nonltr/Chromosome.o
OChromosomeOneDigit = ../bin/nonltr/ChromosomeOneDigit.o
OChromosomeRandom = ../bin/nonltr/ChromosomeRandom.o
OChromListMaker = ../bin/nonltr/ChromListMaker.o
OTableBuilder = ../bin/nonltr/TableBuilder.o
OScorer = ../bin/nonltr/Scorer.o
ODetectorMaxima = ../bin/nonltr/DetectorMaxima.o
OChromDetectorMaxima =  ../bin/nonltr/ChromDetectorMaxima.o
OHMM = ../bin/nonltr/HMM.o
OScanner = ../bin/nonltr/Scanner.o
OTrainer = ../bin/nonltr/Trainer.o
OLocationList = ../bin/nonltr/LocationList.o
OLocationListCollection = ../bin/nonltr/LocationListTailCollection.o
OJosephClient = ../bin/ltr/JosephClient.o
OBJS = $(ORed) $(OInvalidInputException) $(OInvalidStateException) $(OFileDoesNotExistException) $(OInvalidOrderOfOperationsException) $(OInvalidOperationException) $(OInvalidScoreException) $(OUtil) $(OLocation) $(OEmptyLocation) $(OChromosome) $(OChromosomeOneDigit) $(OChromosomeRandom) $(OChromListMaker) $(OTableBuilder) $(OScorer) $(ODetectorMaxima) $(OChromDetector) $(OChromDetectorMaxima) $(OHMM) $(OScanner) $(OTrainer) $(OLocationList) $(OLocationListCollection)
#$(OUtil) $(OLocation) $(OChromosome) $(OChromosomeOneDigit) $(OChromosomeRandom) $(OScorerTr) $(ODetectorTr) $(OTr) $(OForwardTr) $(OBackwardTr) $(OTrCollector) $(OLtrTe) $(OTrKVisitor) $(OTrCsVisitor) $(OTrSineVisitor)

#LTR
OScorerTr = ../bin/ltr/ScorerTr.o
ODetectorTr = ../bin/ltr/DetectorTr.o
OFilterTr = ../bin/ltr/FilterTr.o
OTr = ../bin/ltr/Tr.o
OForwardTr = ../bin/ltr/ForwardTr.o
OBackwardTr = ../bin/ltr/BackwardTr.o

OTrCollector = ../bin/ltr/TrCollector.o
OLtrTe = ../bin/ltr/LtrTe.o
OTrKVisitor = ../bin/ltr/TrKVisitor.o
OTrCsVisitor = ../bin/ltr/TrCsVisitor.o
OTrSineVisitor = ../bin/ltr/TrSineVisitor.o


TestObs = $(OTestTr) $(OChromosome) $(OChromosomeOneDigit) $(OChromosomeRandom) $(OFileDoesNotExistException) $(OInvalidInputException) $(OEmptyTail) $(OInvalidStateException) $(OInvalidOrderOfOperationsException) $(OInvalidOperationException) $(OInvalidScoreException) $(OLocation) $(OUtil) $(OEmptyLocation) $(OLCS) $(OLCSLen) $(OLCSubStr) $(OTSD) $(OEmptyTSD) $(OTailFinder) $(OTail) $(OITail) $(OScorerTr) $(ODetectorTr) $(OFilterTr) $(OTr) $(OForwardTr) $(OBackwardTr) $(OLtrTe) $(OTrCollector) $(OTrKVisitor) $(OTrCsVisitor) $(OTrSineVisitor)
#


#
# TargetTail
#The main function must be declared as a non-member function in the global namespace. This means that it cannot be a static or non-static member function of a class, nor can it be placed in a namespace (even the unnamed namespace). 
TRed = ../bin/Red

JStuff = ../bin/Joseph

LongTerm = ../bin/LongTerm

#
# Make RepeatsDetector
#

$(TRed): $(OBJS)
	$(CXX) -o $(TRed) $(OBJS)

#
# LTR COMPILE
#

$(OJosephClient): ltr/JosephClient.cpp nonltr/Chromosome.h nonltr/ChromListMaker.h
	$(CXX) $(CXXFLAGS) -c ltr/JosephClient.cpp -o $(OJosephClient)

$(JStuff): $(OJosephClient) $(OChromosome) $(OChromListMaker) $(OChromosomeOneDigit) $(OInvalidInputException) $(OInvalidOperationException) 
	$(CXX) -o $(JStuff) $(OJosephClient) $(OChromosome) $(OChromListMaker) $(OChromosomeOneDigit) $(OInvalidInputException) $(OInvalidOperationException) 

$(OTestTr): test/TestTr.cpp nonltr/ChromosomeOneDigit.h nonltr/Chromosome.h nonltr/ChromosomeRandom.h ltr/ScorerTr.h ltr/DetectorTr.h ltr/FilterTr.h ltr/Tr.h ltr/ForwardTr.h ltr/BackwardTr.h ltr/TrCollector.h utility/Util.h utility/LCS.h utility/LCSLen.h utility/Location.h utility/ILocation.h utility/LCSubStr.h utility/TSD.h
	$(CXX) $(CXXFLAGS) -c test/TestTr.cpp -o $(OTestTr)

$(LongTerm): $(TestObs)
	$(CXX) $(CXXFLAGS) -o $(LongTerm) $(TestObs)
#
# RepeatsDetector
#

$(ORed): RepeatsDetector.cpp nonltr/KmerHashTable.h nonltr/KmerHashTable.cpp nonltr/TableBuilder.h nonltr/HMM.h nonltr/Scanner.h nonltr/Trainer.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c RepeatsDetector.cpp -o $(ORed) 

#
# ExceptionTTest
#
$(OInvalidInputException): exception/InvalidInputException.cpp exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidInputException.cpp -o $(OInvalidInputException)
	
$(OInvalidStateException): exception/InvalidStateException.cpp exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidStateException.cpp -o $(OInvalidStateException)	
	
$(OFileDoesNotExistException): exception/FileDoesNotExistException.cpp exception/FileDoesNotExistException.h
	$(CXX) $(CXXFLAGS) -c exception/FileDoesNotExistException.cpp -o $(OFileDoesNotExistException)

$(OInvalidOrderOfOperationsException): exception/InvalidOrderOfOperationsException.cpp exception/InvalidOrderOfOperationsException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidOrderOfOperationsException.cpp -o $(OInvalidOrderOfOperationsException)

$(OInvalidScoreException): exception/InvalidScoreException.cpp exception/InvalidScoreException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidScoreException.cpp -o $(OInvalidScoreException)

$(OInvalidOperationException): exception/InvalidOperationException.cpp exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c exception/InvalidOperationException.cpp -o $(OInvalidOperationException)

#utility/Util.h exception/InvalidInputException
# Utility
#

$(OUtil): utility/Util.cpp utility/Util.h utility/Location.h exception/FileDoesNotExistException.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c utility/Util.cpp -o $(OUtil)

	
$(OLocation): utility/Location.cpp utility/Location.h utility/ILocation.h exception/InvalidInputException.h utility/Util.h 
	$(CXX) $(CXXFLAGS) -c utility/Location.cpp -o $(OLocation)
	
$(OEmptyLocation): utility/EmptyLocation.cpp utility/EmptyLocation.h utility/ILocation.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c utility/EmptyLocation.cpp -o $(OEmptyLocation)	

$(OLCS): utility/LCS.cpp utility/LCS.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c utility/LCS.cpp -o $(OLCS)

$(OLCSLen): utility/LCSLen.cpp utility/LCSLen.h utility/Util.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c utility/LCSLen.cpp -o $(OLCSLen)


$(OLCSubStr): utility/LCSubStr.cpp utility/LCSubStr.h utility/Location.h utility/ILocation.h utility/Util.h 
	$(CXX) $(CXXFLAGS) -c utility/LCSubStr.cpp -o $(OLCSubStr)

$(OTSD): utility/TSD.cpp utility/TSD.h utility/ITSD.h utility/Location.h utility/EmptyLocation.h utility/Util.h utility/LCSubStr.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c utility/TSD.cpp -o $(OTSD)

$(OEmptyTSD): utility/EmptyTSD.cpp utility/EmptyTSD.h utility/Util.h exception/InvalidOperationException.h utility/ITSD.h utility/ILocation.h
	$(CXX) $(CXXFLAGS) -c utility/EmptyTSD.cpp -o $(OEmptyTSD)

$(OTailFinder): utility/TailFinder.cpp utility/TailFinder.h	utility/ILocation.h utility/Util.h exception/InvalidInputException.h exception/InvalidStateException.h 
	$(CXX) $(CXXFLAGS) -c utility/TailFinder.cpp -o $(OTailFinder)

$(OTail): utility/Tail.cpp utility/Tail.h utility/ITail.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c utility/Tail.cpp -o $(OTail)	

$(OITail): utility/ITail.cpp utility/ITail.h utility/ILocation.h
	$(CXX) $(CXXFLAGS) -c utility/ITail.cpp -o $(OITail)


$(OEmptyTail): utility/EmptyTail.h utility/EmptyTail.cpp utility/ITail.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c utility/EmptyTail.cpp -o $(OEmptyTail)





#
# Non LTR
#

$(OChromosome): nonltr/Chromosome.cpp nonltr/Chromosome.h nonltr/IChromosome.h utility/Util.h exception/InvalidInputException.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c nonltr/Chromosome.cpp -o $(OChromosome)
	
$(OChromosomeOneDigit): nonltr/ChromosomeOneDigit.cpp nonltr/ChromosomeOneDigit.h nonltr/Chromosome.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromosomeOneDigit.cpp	-o $(OChromosomeOneDigit)
	
$(OChromosomeRandom): nonltr/ChromosomeRandom.cpp nonltr/ChromosomeRandom.h nonltr/IChromosome.h exception/InvalidInputException.h exception/InvalidStateException.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromosomeRandom.cpp -o $(OChromosomeRandom)
	
$(OTableBuilder): nonltr/TableBuilder.cpp nonltr/TableBuilder.h utility/Util.h nonltr/ChromosomeOneDigit.h nonltr/ITableView.h nonltr/KmerHashTable.h nonltr/KmerHashTable.cpp nonltr/EnrichmentMarkovView.h nonltr/EnrichmentMarkovView.cpp exception/InvalidStateException.h nonltr/ChromListMaker.h nonltr/IChromosome.h
	$(CXX) $(CXXFLAGS) -c nonltr/TableBuilder.cpp -o $(OTableBuilder)
	
$(OScorer): nonltr/Scorer.cpp nonltr/Scorer.h nonltr/ChromosomeOneDigit.h utility/Util.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/Scorer.cpp -o $(OScorer)
		
$(ODetectorMaxima): nonltr/DetectorMaxima.cpp nonltr/DetectorMaxima.h utility/ILocation.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/DetectorMaxima.cpp -o $(ODetectorMaxima)

$(OChromDetectorMaxima): nonltr/ChromDetectorMaxima.cpp nonltr/ChromDetectorMaxima.h nonltr/DetectorMaxima.h nonltr/ChromosomeOneDigit.h utility/Util.h utility/ILocation.h utility/Location.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromDetectorMaxima.cpp -o $(OChromDetectorMaxima)	
	
$(OHMM): nonltr/HMM.cpp nonltr/HMM.h utility/ILocation.h exception/InvalidStateException.h exception/InvalidInputException.h exception/FileDoesNotExistException.h exception/InvalidOperationException.h
	$(CXX) $(CXXFLAGS) -c nonltr/HMM.cpp -o $(OHMM)
	
$(OScanner): nonltr/Scanner.cpp nonltr/Scanner.h nonltr/Chromosome.h nonltr/ChromosomeOneDigit.h nonltr/HMM.h nonltr/ITableView.h nonltr/Scorer.h utility/Util.h utility/ILocation.h exception/InvalidInputException.h exception/InvalidStateException.h exception/FileDoesNotExistException.h exception/InvalidOperationException.h	
	$(CXX) $(CXXFLAGS) -c nonltr/Scanner.cpp -o $(OScanner)
	
$(OTrainer): nonltr/Trainer.cpp nonltr/Trainer.h nonltr/TableBuilder.h nonltr/KmerHashTable.h nonltr/KmerHashTable.cpp nonltr/HMM.h nonltr/ChromDetectorMaxima.h nonltr/Scorer.h nonltr/ChromListMaker.h utility/Util.h nonltr/LocationListCollection.h
	$(CXX) $(CXXFLAGS) -c nonltr/Trainer.cpp -o $(OTrainer)
	
$(OChromListMaker): nonltr/ChromListMaker.cpp nonltr/ChromListMaker.h nonltr/Chromosome.h nonltr/ChromosomeOneDigit.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromListMaker.cpp -o $(OChromListMaker) 
	
$(OCluster): nonltr/Cluster.cpp nonltr/Cluster.h utility/Util.h exception/InvalidStateException.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c nonltr/Cluster.cpp -o $(OCluster)


$(OLocationList): nonltr/LocationList.cpp nonltr/LocationList.h utility/ILocation.h utility/Location.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/LocationList.cpp -o $(OLocationList)

$(OLocationListCollection): nonltr/LocationListCollection.cpp nonltr/LocationListCollection.h utility/Location.h exception/InvalidStateException.h
	$(CXX) $(CXXFLAGS) -c nonltr/LocationListCollection.cpp -o $(OLocationListCollection)
	ltr/ChromListMaker.cpp nonltr/ChromListMaker.h nonltr/Chromosome.h nonltr/ChromosomeOneDigit.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c nonltr/ChromListMaker.cpp -o $(OChromListMaker) 
	

#LTR

$(OScorerTr): ltr/ScorerTr.cpp ltr/ScorerTr.h nonltr/KmerHashTable.h nonltr/ChromosomeOneDigit.h utility/Util.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c ltr/ScorerTr.cpp -o $(OScorerTr)

$(ODetectorTr): ltr/DetectorTr.cpp  ltr/DetectorTr.cpp ltr/BackwardTr.h ltr/ForwardTr.h utility/Util.h exception/InvalidStateException.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c ltr/DetectorTr.cpp -o $(ODetectorTr)

$(OFilterTr): ltr/FilterTr.cpp ltr/FilterTr.h ltr/BackwardTr.h ltr/LtrTe.h ltr/TrKVisitor.h ltr/TrCsVisitor.h ltr/TrSineVisitor.h utility/Util.h utility/TSD.h utility/EmptyTSD.h utility/TailFinder.h utility/Tail.h utility/EmptyTail.h exception/InvalidStateException.h exception/InvalidInputException.h
	$(CXX) $(CXXFLAGS) -c ltr/FilterTr.cpp -o $(OFilterTr)

$(OTr): ltr/Tr.cpp ltr/Tr.h ltr/ITrVisitor.h exception/InvalidStateException.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c ltr/Tr.cpp -o $(OTr)

$(OForwardTr): ltr/ForwardTr.cpp ltr/ForwardTr.h ltr/BackwardTr.h ltr/Tr.h exception/InvalidOperationException.h exception/InvalidStateException.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c ltr/ForwardTr.cpp -o $(OForwardTr)
 
$(OBackwardTr): ltr/BackwardTr.cpp ltr/BackwardTr.h ltr/Tr.h ltr/ForwardTr.h exception/InvalidStateException.h exception/InvalidOperationException.h utility/Util.h
	$(CXX) $(CXXFLAGS) -c ltr/BackwardTr.cpp -o $(OBackwardTr)

$(OTrKVisitor): ltr/TrKVisitor.cpp ltr/TrKVisitor.h ltr/ITrVisitor.h
	$(CXX) $(CXXFLAGS) -c ltr/TrKVisitor.cpp -o $(OTrKVisitor)

$(OTrCsVisitor): ltr/TrCsVisitor.cpp ltr/TrCsVisitor.h ltr/ITrVisitor.h
	$(CXX) $(CXXFLAGS) -c ltr/TrCsVisitor.cpp -o $(OTrCsVisitor)

$(OTrSineVisitor): ltr/TrSineVisitor.cpp ltr/TrSineVisitor.h ltr/ITrVisitor.h
	$(CXX) $(CXXFLAGS) -c ltr/TrSineVisitor.cpp -o $(OTrSineVisitor)

$(OTrCollector): ltr/TrCollector.cpp ltr/TrCollector.h ltr/ScorerTr.h ltr/DetectorTr.h ltr/FilterTr.h nonltr/ChromosomeOneDigit.h
	$(CXX) $(CXXFLAGS) -c ltr/TrCollector.cpp -o $(OTrCollector)

$(OLtrTe): ltr/LtrTe.cpp ltr/LtrTe.h exception/InvalidStateException.h exception/InvalidOperationException.h utility/Util.h utility/ITail.h utility/Tail.h utility/EmptyTail.h utility/EmptyTSD.h utility/ITSD.h
	$(CXX) $(CXXFLAGS) -c ltr/LtrTe.cpp -o $(OLtrTe)
	


#
# Make binary directories
#

red: $(TRed)

#
# Make Red
#




# Make Long Term Repeats
#
#
ltr: $(LongTerm)

this: $(OTr)



bin: 
	mkdir ../bin
	mkdir ../bin/exception			
	mkdir ../bin/utility
	mkdir ../bin/nonltr	
	mkdir ../bin/ltr	

#
# Make clean
#

clean:
	rm -f ../bin/*.o ../bin/exception/*.o ../bin/ms/*.o ../bin/nonltr/*.o ../bin/test/*.o ../bin/utility/*.o ../bin/ltr/*.o*.o $(LongTerm)

