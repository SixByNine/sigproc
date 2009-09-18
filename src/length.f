c     DECK LENGTH
c     
c     
c     
c     
c     RETURNS THE LENGTH OF 'STRING' EXCLUDING ANY TRAILING SPACES.
c     
      integer function length(string)
      implicit none
      character string*(*)
c     
c     OBTAIN THE LOCATION OF THE LAST NON-SPACE CHARACTER.
c     
      integer ilen,i
c     search for the first null
      ilen = len(string)
c     use the position of the first null
      do 1 i = ilen, 1, -1
c     
c     LENGTH FOUND.
c     
         if (string(i:i) .ne. char(32) .and.
     &       string(i:i).ne.char(0)) then
            length = i
            return 
         end if
c     
c     STRING IS ALL SPACES OR ZERO LENGTH.
c     
    1 continue
      length = 0
c     
c     END OF INTEGER FUNCTION LENGTH.
c     
      return 
      end










