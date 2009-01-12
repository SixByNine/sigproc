C @(#)find_summary.inc	3.1 12/17/92
c=================================================================
c     VLSA : Include file for the  " .s. "  output file
c=================================================================

C n.b. remember to update the find_s_read.f file when new parameters
C are added

C split off the common block with the characters - pah - The HP optimizer
C will produce garbage if characters in commons.

      real*8 s_tsmp
      
      real s_dm,s_rms,s_avs,s_ave,s_sprms,s_alfa,
     &    s_spectrum(2,1024,2),
     &    s_pamp(100,5),s_prds(100,5),s_ixhst(15),s_iyhst(15),
     &    s_bfa(20),s_bfb(20),snbrd(20)

      integer s_lbig(100,5),s_nbrd,s_nba(20),s_nbb(20),s_ndat,
     &    s_nskip,s_nclip,s_nrm,s_nfreq,s_ndatfft,s_nhbrd(20),
     &    s_nkill,s_chk(16),s_loop,s_nssum,s_dm_step,s_rebin_step

      character s_lb1*130,s_lb2*130,s_lb3*130
	
      byte bufc(32768), bufh(1024)

      common/summary/s_tsmp,s_dm,s_rms,s_avs,s_ave,s_sprms,s_alfa,
     &    s_ixhst,s_iyhst,s_nbrd,s_nba,s_nbb,s_nssum,
     &    s_bfa,s_bfb,s_nhbrd,s_ndatfft,s_nkill,s_chk,
     &    s_ndat,s_nskip,s_nclip,s_nrm,s_nfreq,
     &    s_loop,s_dm_step,s_rebin_step,
     &    s_spectrum,s_pamp,s_prds,s_lbig,snbrd

      common/summaryc/s_lb1,s_lb2,s_lb3

      equivalence(s_tsmp,bufc(1)),(s_tsmp,bufh(1))
c==================================================================




