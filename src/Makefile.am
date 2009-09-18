core_src=add_channels.c add_samples.c aliases.c angle_split.c aoftm.c bandfactors.c baseline.f bit.c bpp2fb.c cel2gal.c chan_freqs.c clock.f close_log.c cpower.f cprofc.c decimate_data.c decimate_header.c dedisperse_data.c dedisperse_data.c.bak dedisperse_header.c deg2dms.c dialog.C dialog.h dmdelay.c dmshift.c dosearch.f epnpack.f eraseDP.c error_message.c fast_mean.f fetch_hdrval.c ffft.f ffreq.c fftdata.f file_exists.c filterbank_header.c find_norm.f fold_data.c fold_header.c formspec.f freq.f fshift_prof.c getddis.f getfast_rmea.f getmeanrms.f getmjd.f getrmea.f getrrms.f glun.f gmrt2fb.c help.csh help.c histmax.f head_parse.c ignored_channels.c indexxf77.f indexx.c indexxf77.f inv_cerf.c inputclip.c length.f lookup.c machine2prf.c meanrms.f minmax.f mjd.c mmzap.f normal.f norm_prof.c np2.c nrselect.c nrsort.c nrutil.c nrutil.h nsamples.c ooty2fb.c open_file.c open_files.c open_files.c open_log.c pack_unpack.c phcalc.c print_version.c process.c prof_adds.c prof_ddis.c prof_sbas.c prof_sumc.c prof_sumifs.c pspm2fb.c pspm_decode.c pspm_prof.c pspm_tstart.c pulsar2k2fb.c pulse.h pulse.c pulse.inc put.c quikgray.f ralphs_mean.f random.c rdfbtab.f rdhead.f read_aoscan.c read_block.c read_blocksubzero.c subzero.c readdat.f read_header.c read_polyco.c readspec.f readsus.f readtim.c rebin.f recipes.c recon_prof.f resample.f resample_old.f reverse_baseline_normalise.c rotate_time.f rwepn.f scaledata.c scale_prof.c scamp2fb.c seekin.f select.c send_stuff.c shift_prof.c short.f single_ch.c singlepulse.f sizeof_file.c slalib.c slalib.h slfit.f smooth.c spcsnr.f sprof.f ssm.c strings_equal.c submean.c submedian.c sumhrm.f swap_bytes.c thresh_1d.c timer.f typeof_inputdata.c update_log.c uttime.c vanvleck.c vmax.c vmin.c wappcorrect.c wapp_prof.c whiten.c write_epn.c writeepn.f write_profiles.c write_pulses.c writespec.f y.tab.c y.tab.h zap_birdies.f zapit.f dialog.C dialog.h libplotfil.h libplotfil.C gtools.C gtools.h newnewsumhrm.f newoldsumhrm.f 

fortran_src=baseline.f best.f chaninfo.f clock.f cpower.f dosearch.f epnpack.f fast_mean.f ffa.f ffft.f fftdata.f fftwdata.f fftwdummy.f find_norm.f formspec.f freq.f getddis.f getfast_rmea.f getmeanrms.f getmjd.f getrmea.f getrrms.f glun.f grey.f histmax.f indexxf77.f length.f mask.f meanrms.f minmax.f mmzap.f newnewsumhrm.f newoldsumhrm.f normal.f peak.f plotpulses.f plotpulses.old.f postproc.f ralphs_mean.f rdfbtab.f rdhead.f readdat.f readspec.f readsus.f rebin.f recon_prof.f resample.f resample_old.f rotate_time.f rwepn.f seekin.f short.f singlepulse.f slfit.f sort3.f spcsnr.f sprof.f step.f sumhrm.f tbary.f timer.f tssb.f writeepn.f writespec.f zap_birdies.f zapit.f


pgpl_src=quickplot.f quikgray.f



bin_PROGRAMS=seek spec

seek_SOURCES=seek.f dosearch.f ffft.f fftdata.f  find_norm.f formspec.f freq.f glun.f sumhrm.f seekin.f baseline.f readspec.f writespec.f singlepulse.f getrmea.f getmeanrms.f zapit.f timer.f indexxf77.f slfit.f meanrms.f getrrms.f readdat.f resample.f rotate_time.f zap_birdies.f mmzap.f recon_prof.f length.f rdhead.f rdfbtab.f getddis.f getmjd.f  sigproc.h readtim.c pulse.c read_header.c angle_split.c clock.f open_file.c strings_equal.c aliases.c open_files.c single_ch.c write_pulses.c thresh_1d.c smooth.c

if HAVE_FFTW
seek_SOURCES+=fftwdata.f
seek_LDFLAGS=@FFTW_LIBS@
INCLUDES=@FFTW_CFLAGS@
else
seek_SOURCES+=fftwdummy.f
endif

spec_SOURCES=spec.f $(fortran_src)

sigproc.h: include.csh
	./include.csh
