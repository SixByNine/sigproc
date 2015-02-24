bin_PROGRAMS=seek header dedisperse extract fake flatten clip downsample fold profile flux header bandpass reader readchunk splitter snrdm decimate splice dice depolyco blanker barycentre tree polyco2period makePsrXml ffa grey chaninfo step postproc filedit chop_fil dedisperse_all shredder filterbank fast_fake zerodm 
bin_SCRIPTS=accn  ahunt  csearch  foldsignals  hunt makedummy  monitor  quicklook spsearch dmsampdelay
EXTRA_DIST=accn.csh ahunt.csh csearch.csh differ.csh exporter.csh foldsignals.csh help.csh hunt.csh include.csh makedummy.csh mkmonitor.csh quicklook.csh subint2prf.csh spsearch.pl dmsampdelay.pl
MOSTLYCLEANFILES=sigproc.h accn  ahunt  csearch  foldsignals  hunt makedummy  monitor  quicklook spsearch dmsampdelay
noinst_LIBRARIES=libsigproc.a

EXTRA_PROGRAMS=fddtw

if HAVE_T2PRED
bin_PROGRAMS += inject_pulsar
if HAVE_PGPLOT
bin_PROGRAMS += rfi_analyse
endif
endif

BUILT_SOURCES=sigproc.h

fast_fake_SOURCES=fast_fake.c
fast_fake_CPPFLAGS=-std=gnu99

AM_FFLAGS=-fno-second-underscore
AM_CPPFLAGS=.

f_headers=birdies.inc csamp.inc epnhdr.inc find_summary.inc folds.inc head.inc npts.inc peak.inc pulse.inc seek.inc time.inc vers.inc 
c_headers=bacspin.h bpp.h bpphdr.h cpgplot.h decimate.h dedisperse.h dialog.h epn.h filterbank.h find_baseline.h find.h fold.h gtools.h header.h key.h libplotfil.h longnam.h nr.h nrutil.h older_wapp_header.h old_wapp_header.h ooty.h polyco.h prec_proto.h profile.h pspmhdr.h pulse.h scamp_header.h slalib.h slamac.h version.h wapp_header.h y.tab.h mjklog.h mjk_random.h mjk_cmd.h

libsigproc_a_SOURCES=sigproc.h $(c_headers) add_channels.c add_samples.c aliases.c angle_split.c aoftm.c bandfactors.c bit.c cel2gal.c chan_freqs.c close_log.c cprofc.c decimate_data.c decimate_header.c deg2dms.c dmdelay.c dmshift.c eraseDP.c error_message.c fetch_hdrval.c ffreq.c file_exists.c filterbank_header.c fold_data.c fold_header.c fshift_prof.c gmrt2fb.c head_parse.c help.c ignored_channels.c indexx.c inputclip.c inv_cerf.c lookup.c mjd.c mkheader.c norm_prof.c np2.c nrselect.c nrsort.c nrutil.c nsamples.c ooty2fb.c open_file.c open_files.c open_log.c pack_unpack.c phcalc.c print_version.c process.c prof_adds.c prof_cent.c prof_ddis.c prof_sbas.c prof_sumc.c prof_sumifs.c pspm2fb.c pspm_decode.c pspm_prof.c pspm_tstart.c pulsar2k2fb.c pulse.c put.c random.c read_aoscan.c read_block.c read_blocksubzero.c read_header.c read_polyco.c readtim.c recipes.c reverse_baseline_normalise.c scaledata.c scale_prof.c scamp2fb.c select.c send_stuff.c shift_prof.c single_ch.c sizeof_file.c slalib.c smooth.c ssm.c strings_equal.c submean.c submedian.c subzero.c swap_bytes.c thresh_1d.c typeof_inputdata.c update_log.c uttime.c vanvleck.c vmax.c vmin.c wapp2fb.c wappcorrect.c wapp_prof.c whiten.c write_epn.c write_profiles.c write_pulses.c ffft.f dedisperse_data.c dedisperse_header.c machine2prf.c mkheaderlex.c alfa_position.c bpp2fb.c mjk_random.c mjk_cmd.c newgmrt2fb.c

if HAVE_T2PRED
inject_pulsar_LDADD=-ltempo2pred libsigproc.a @FFTW_LIBS@
inject_pulsar_CFLAGS=-std=gnu99
inject_pulsar_SOURCES=inject_pulsar.c 
if HAVE_PGPLOT
rfi_analyse_LDADD=-lt2toolkit libsigproc.a @FFTW_LIBS@ @PGPLOT_LIBS@
rfi_analyse_CFLAGS=-std=gnu99
rfi_analyse_SOURCES=rfi_analyse.c
endif
endif

LDADD=libsigproc.a -lm

if HAVE_CFITSIO
libsigproc_a_SOURCES+=psrfits2fb.c readpsrfits_hd.c 
LDADD+=@CFITSIO_LIBS@
AM_CPPFLAGS+=@CFITSIO_CFLAGS@ -DPSRFITS
endif

if HAVE_PGPLOT
bin_PROGRAMS+=spec tune giant spexplorer pgplotter getpulse peak best quickplot plotpulses quickgplot extragplot
# PGPLOT BASED CODE HERE...
spec_SOURCES=spec.f readspec.f freq.f minmax.f histmax.f glun.f
spec_LDADD=@PGPLOT_LIBS@

tune_SOURCES=tune.c
tune_LDADD=@PGPLOT_LIBS@ $(LDADD)

giant_SOURCES=giant.C Creadspec.C  dialog.C  find_baseline.C  find_fft.C  find_formspec.C  gtools.C  libplotfil.C fftdata.f newoldsumhrm.f toolkit.h toolkit.C TKfit.h TKfit.C
giant_LDADD=@PGPLOT_LIBS@ $(LDADD)

spexplorer_SOURCES=spexplorer.C Creadspec.C  dialog.C  find_baseline.C  find_fft.C  find_formspec.C  gtools.C  libplotfil.C fftdata.f newoldsumhrm.f toolkit.h toolkit.C TKfit.h TKfit.C
spexplorer_LDADD=@PGPLOT_LIBS@ $(LDADD)

extragplot_SOURCES=extragplot.C
extragplot_LDADD=@PGPLOT_LIBS@ $(LDADD)

quickgplot_SOURCES=quickgplot.C find_baseline.C libplotfil.C toolkit.h toolkit.C
quickgplot_LDADD=@PGPLOT_LIBS@ $(LDADD)

pgplotter_LDADD=@PGPLOT_LIBS@ $(LDADD)

getpulse_LDADD=@PGPLOT_LIBS@ $(LDADD)

peak_SOURCES=sigproc.h peak.f glun.f writeepn.f timer.f rwepn.f readdat.f getmjd.f readtim.c read_header.c epnpack.f strings_equal.c angle_split.c  rdhead.f clock.f baseline.f aliases.c open_files.c open_file.c length.f getddis.f rdfbtab.f slfit.f sprof.f
peak_LDADD=@PGPLOT_LIBS@ 

best_SOURCES=best.f length.f glun.f indexxf77.f readsus.f
best_LDADD=@PGPLOT_LIBS@

quickplot_SOURCES=quickplot.f
quickplot_LDADD=@PGPLOT_LIBS@

plotpulses_SOURCES=plotpulses.f
plotpulses_LDADD=@PGPLOT_LIBS@

endif

filedit_SOURCES=filedit.c header.h read_header.c strings_equal.c
filedit_LDADD=-lm

# seek
seek_SOURCES=sigproc.h $(f_headers) seek.f dosearch.f ffft.f fftdata.f  find_norm.f formspec.f freq.f glun.f sumhrm.f seekin.f baseline.f readspec.f writespec.f singlepulse.f getrmea.f getmeanrms.f zapit.f timer.f indexxf77.f slfit.f meanrms.f getrrms.f readdat.f resample.f rotate_time.f zap_birdies.f mmzap.f recon_prof.f length.f rdhead.f rdfbtab.f getddis.f getmjd.f readtim.c pulse.c read_header.c angle_split.c clock.f open_file.c strings_equal.c aliases.c open_files.c single_ch.c write_pulses.c thresh_1d.c smooth.c

ffa_SOURCES=ffa.f glun.f timer.f clock.f getmjd.f

grey_SOURCES=grey.f glun.f writeepn.f timer.f rwepn.f readdat.f getmjd.f  clock.f epnpack.f sprof.f rdfbtab.f getddis.f rdhead.f length.f

chaninfo_SOURCES=chaninfo.f

step_SOURCES=step.f
postproc_SOURCES=postproc.f

dedisperse_all_SOURCES=dedisperse_all.C dedisperse_all.h getDMtable.C getDMtable.h gtools.C

shredder_SOURCES=shredder.C dedisperse_all.h getDMtable.C getDMtable.h gtools.C


fddtw_SOURCES=fddtw.C

fddtw_CPPFLAGS=-std=gnu++11 -Wno-write-strings

if HAVE_FFTW
seek_SOURCES+=fftwdata.f
seek_LDFLAGS=@FFTW_LIBS@
AM_CPPFLAGS+=@FFTW_CFLAGS@
else
seek_SOURCES+=fftwdummy.f
endif

src=@top_srcdir@/src

sigproc.h: include.csh
	$(src)/include.csh

monitor : mkmonitor.csh
	$(src)/mkmonitor.csh >monitor
	chmod +x monitor

polyco : polyco.tcl
	echo "#!`which tclsh`" >  polyco
	cat $(src)/polyco.tcl         >> polyco
	chmod +x polyco

makedmlist : makedmlist.tcl
	echo "#!`which tclsh`" >             makedmlist
	cat $(src)/makedmlist.tcl     >>            /makedmlist
	chmod +x makedmlist

quicklook : quicklook.csh
	echo "#!`which csh`" >               quicklook
	cat $(src)/quicklook.csh    >>              quicklook
	chmod +x quicklook

foldsignals : foldsignals.csh
	echo "#!`which csh`" >            foldsignals
	echo 'set bin        = `dirname $$0`' >> foldsignals
	cat $(src)/foldsignals.csh            >> foldsignals
	chmod +x foldsignals

hunt : hunt.csh
	echo "#!`which csh`" >             hunt
	echo 'set bin        = `dirname $$0`'  >> hunt
	cat $(src)/hunt.csh                    >> hunt
	chmod +x hunt

ahunt : ahunt.csh
	echo "#!`which csh`" >             ahunt
	echo 'set bin        = `dirname $$0`'  >> ahunt
	cat $(src)/ahunt.csh                    >> ahunt
	chmod +x ahunt

spsearch : spsearch.pl
	cp $(src)/spsearch.pl spsearch
	chmod +x spsearch

dmsampdelay : dmsampdelay.pl
	cp $(src)/dmsampdelay.pl dmsampdelay
	chmod +x dmsampdelay

csearch : csearch.csh
	echo "#!`which csh`" >             csearch
	echo 'set bin        = `dirname $$0`'  >> csearch
	cat $(src)/csearch.csh                 >> csearch
	chmod +x csearch

accn : accn.csh
	echo "#!`which csh`" >             accn
	echo 'set bin        = `dirname $$0`'  >> accn
	cat $(src)/accn.csh                    >> accn
	chmod +x accn

makedummy : makedummy.csh
	cp $(src)/makedummy.csh makedummy
	chmod +x makedummy






