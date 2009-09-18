c==============================================================================
c Configuration of harmonic folds to use 
c==============================================================================
c
c 'make clean; make' after modifying this file
c
      integer nfolds
      parameter(nfolds=5)
      integer foldvals(nfolds)
      data foldvals /1,2,4,8,16/
