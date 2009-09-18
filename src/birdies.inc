C @(#)birdies.inc	3.1 12/17/92

      integer nbrd, nba, nbb, irange_a, irange_b, nhbrd, zapped
      real bfa, bfb, tot, bffa, bffb, perzapped
      logical beambirdie

      common/birdie/nbrd,bfa(20),bfb(20),nba(20),nbb(20),tot,
     +    irange_a(20),irange_b(20),nhbrd(20),bffa(20), bffb(20)
     +    ,beambirdie(20),zapped,perzapped
