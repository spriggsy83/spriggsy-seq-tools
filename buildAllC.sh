# Builds programs:
# getSeqSizeStats
# getSeqSizeStatsT
# getSeqQCStats
# getSeqCGstats
# getSeqSizeList
# getSeqSizeChart
# getSeqCountTable
# filterSeqSize
# getSubSeqs
# extractSeqSubsets
# excludeSeqsBySAM
# tallySNPs2
# reverseComplement
# splitSeqsIntoXFiles
# tallyGeneCoverageSamGZ
# mergeKmerCounts

#Requires Boost C++ Libraries and OpenMPI
#module load boost
#module load openmpi

cd CppLibrary
g++ -o ../getSeqSizeStats getSeqSizeStats.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../getSeqSizeStatsT getSeqSizeStatsT.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../getSeqQCStats getSeqQCStats.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../getSeqCGstats getSeqCGstats.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../getSeqSizeList getSeqSizeList.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../getSeqSizeChart getSeqSizeChart.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../filterSeqSize filterSeqSize.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../getSubSeqs getSubSeqs.cpp SeqReader.cpp -lboost_iostreams -lz -lboost_regex
g++ -o ../getSeqCountTable getSeqCountTable.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../extractSeqSubsets extractSeqSubsets.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../excludeSeqsBySAM excludeSeqsBySAM.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../tallySNPs2 tallySNPs2.cpp SNPTallyer2.cpp SeqReader.cpp AlignedRead.cpp -fopenmp -lboost_iostreams -lz
g++ -o ../reverseComplement reverseComplement.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../splitSeqsIntoXFiles splitSeqsIntoXFiles.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -o ../tallyGeneCoverageSamGZ tallyGeneCoverageSamGZ.cpp GeneCoverageTallyerSamGZ.cpp -fopenmp -lboost_iostreams -lz
g++ -o ../mergeKmerCounts mergeKmerCounts.cpp KmerCountMerger.cpp SeqReader.cpp -lboost_iostreams -lz
g++ -std=c++0x -o ../buildKmerChaos buildKmerChaos.cpp KmerChaosGen.cpp SeqReader.cpp -fopenmp -lboost_iostreams -lz
