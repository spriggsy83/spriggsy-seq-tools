# Builds programs:
# getSeqSizeStats
# getSeqQCStats
# getSeqCGstats
# getSeqSizeList
# getSeqSizeChart
# getSeqCountTable
# filterSeqSize
# excludeEveryXSeq
# extractEveryXSeq
# excludeSeqsBySAM
# tallySNPs2
# reverseComplement
# splitSeqsIntoXFiles
# tallyGeneCoverageSamGZ

#Requires Boost C++ Libraries and OpenMPI
#module load boost
#module load openmpi

cd CppLibrary
g++ -lboost_iostreams -o ../getSeqSizeStats getSeqSizeStats.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../getSeqQCStats getSeqQCStats.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../getSeqCGstats getSeqCGstats.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../getSeqSizeList getSeqSizeList.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../getSeqSizeChart getSeqSizeChart.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../filterSeqSize filterSeqSize.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../getSeqCountTable getSeqCountTable.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../excludeEveryXSeq excludeEveryXSeq.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../extractEveryXSeq extractEveryXSeq.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../excludeSeqsBySAM excludeSeqsBySAM.cpp SeqReader.cpp
g++ -fopenmp -lboost_iostreams -o ../tallySNPs2 tallySNPs2.cpp SNPTallyer2.cpp SeqReader.cpp AlignedRead.cpp
g++ -lboost_iostreams -o ../reverseComplement reverseComplement.cpp SeqReader.cpp
g++ -lboost_iostreams -o ../splitSeqsIntoXFiles splitSeqsIntoXFiles.cpp SeqReader.cpp
g++ -fopenmp -lboost_iostreams -o ../tallyGeneCoverageSamGZ tallyGeneCoverageSamGZ.cpp GeneCoverageTallyerSamGZ.cpp
