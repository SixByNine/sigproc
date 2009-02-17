###############################################################################
# "makefile" for SIGPROC - the generic pulsar signal processing software 
# Report any problems compiling to drl@jb.man.ac.uk
# YOU SHOULD NOT NORMALLY HAVE TO CHANGE ANYTHING IN THIS FILE !!!!!!!!!!!!!
###############################################################################
include makefile.$(OSTYPE)
CC = $(CCC) $(DFITS) $(DFFTW) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
###############################################################################
LIB = libsigproc_$(OSTYPE).a
LIBOBJECTS = \
$(LIB)(add_channels.o)\
$(LIB)(add_samples.o)\
$(LIB)(aliases.o)\
$(LIB)(angle_split.o)\
$(LIB)(aoftm.o)\
$(LIB)(bandfactors.o)\
$(LIB)(baseline.o)\
$(LIB)(bit.o)\
$(LIB)(bpp2fb.o)\
$(LIB)(cel2gal.o)\
$(LIB)(chan_freqs.o)\
$(LIB)(clock.o)\
$(LIB)(close_log.o)\
$(LIB)(cpower.o)\
$(LIB)(cprofc.o)\
$(LIB)(decimate_data.o)\
$(LIB)(decimate_header.o)\
$(LIB)(dedisperse_data.o)\
$(LIB)(dedisperse_header.o)\
$(LIB)(deg2dms.o)\
$(LIB)(dialog.o)\
$(LIB)(dmdelay.o)\
$(LIB)(dmshift.o)\
$(LIB)(dosearch.o)\
$(LIB)(epnpack.o)\
$(LIB)(eraseDP.o)\
$(LIB)(error_message.o)\
$(LIB)(fast_mean.o)\
$(LIB)(fetch_hdrval.o)\
$(LIB)(ffft.o)\
$(LIB)(ffreq.o)\
$(LIB)(fftdata.o)\
$(LIB)($(FFTWF))\
$(LIB)(file_exists.o)\
$(LIB)(filterbank_header.o)\
$(LIB)(find_norm.o)\
$(LIB)(fold_data.o)\
$(LIB)(fold_header.o)\
$(LIB)(formspec.o)\
$(LIB)(freq.o)\
$(LIB)(fshift_prof.o)\
$(LIB)(getddis.o)\
$(LIB)(getfast_rmea.o)\
$(LIB)(getmeanrms.o)\
$(LIB)(getmjd.o)\
$(LIB)(getrmea.o)\
$(LIB)(getrrms.o)\
$(LIB)(glun.o)\
$(LIB)(gmrt2fb.o)\
$(LIB)(help.o)\
$(LIB)(histmax.o)\
$(LIB)(ignored_channels.o)\
$(LIB)(indexxf77.o)\
$(LIB)(indexx.o)\
$(LIB)(inv_cerf.o)\
$(LIB)(length.o)\
$(LIB)(lookup.o)\
$(LIB)(machine2prf.o)\
$(LIB)(meanrms.o)\
$(LIB)(minmax.o)\
$(LIB)(mjd.o)\
$(LIB)(mmzap.o)\
$(LIB)(normal.o)\
$(LIB)(norm_prof.o)\
$(LIB)(np2.o)\
$(LIB)(nrselect.o)\
$(LIB)(nrsort.o)\
$(LIB)(nrutil.o)\
$(LIB)(nsamples.o)\
$(LIB)(ooty2fb.o)\
$(LIB)(open_file.o)\
$(LIB)(open_files.o)\
$(LIB)(open_log.o)\
$(LIB)(pack_unpack.o)\
$(LIB)(phcalc.o)\
$(LIB)(print_version.o)\
$(LIB)(process.o)\
$(LIB)(prof_adds.o)\
$(LIB)(prof_ddis.o)\
$(LIB)(prof_sbas.o)\
$(LIB)(prof_sumc.o)\
$(LIB)(prof_sumifs.o)\
$(LIB)(pspm2fb.o)\
$(LIB)(pspm_decode.o)\
$(LIB)(pspm_prof.o)\
$(LIB)(pspm_tstart.o)\
$(LIB)(pulsar2k2fb.o)\
$(LIB)(pulse.o)\
$(LIB)(put.o)\
$(LIB)(quikgray.o)\
$(LIB)(ralphs_mean.o)\
$(LIB)(random.o)\
$(LIB)(rdfbtab.o)\
$(LIB)(rdhead.o)\
$(LIB)(read_aoscan.o)\
$(LIB)(read_block.o)\
$(LIB)(readdat.o)\
$(LIB)(read_header.o)\
$(LIB)(read_polyco.o)\
$(LIB)(readspec.o)\
$(LIB)(readsus.o)\
$(LIB)(readtim.o)\
$(LIB)(rebin.o)\
$(LIB)(recipes.o)\
$(LIB)(recon_prof.o)\
$(LIB)(resample.o)\
$(LIB)(rotate_time.o)\
$(LIB)(rwepn.o)\
$(LIB)(scaledata.o)\
$(LIB)(scale_prof.o)\
$(LIB)(scamp2fb.o)\
$(LIB)(seekin.o)\
$(LIB)(select.o)\
$(LIB)(send_stuff.o)\
$(LIB)(shift_prof.o)\
$(LIB)(short.o)\
$(LIB)(single_ch.o)\
$(LIB)(singlepulse.o)\
$(LIB)(sizeof_file.o)\
$(LIB)(slalib.o)\
$(LIB)(slfit.o)\
$(LIB)(smooth.o)\
$(LIB)(spcsnr.o)\
$(LIB)(sprof.o)\
$(LIB)(ssm.o)\
$(LIB)(strings_equal.o)\
$(LIB)(submean.o)\
$(LIB)(submedian.o)\
$(LIB)(sumhrm.o)\
$(LIB)(swap_bytes.o)\
$(LIB)(thresh_1d.o)\
$(LIB)(timer.o)\
$(LIB)(typeof_inputdata.o)\
$(LIB)(update_log.o)\
$(LIB)(uttime.o)\
$(LIB)(vanvleck.o)\
$(LIB)(vmax.o)\
$(LIB)(vmin.o)\
$(LIB)(wappcorrect.o)\
$(LIB)(wapp_prof.o)\
$(LIB)(whiten.o)\
$(LIB)(write_epn.o)\
$(LIB)(writeepn.o)\
$(LIB)(write_profiles.o)\
$(LIB)(write_pulses.o)\
$(LIB)(writespec.o)\
$(LIB)(y.tab.o)\
$(LIB)(zap_birdies.o)\
$(LIB)(zapit.o)
###############################################################################
all: library programs scripts 

programs: filterbank fake header bandpass reader decimate \
	dedisperse fold polyco2period profile flux splice barycentre \
	seek best step postproc grey mask csearch polyco depolyco tune \
	chaninfo plotpulses getpulse

scripts: monitor quicklook polyco makedmlist foldsignals hunt accn ahunt

sigproc.h : *.c
	./include.csh

library : sigproc.h $(LIBOBJECTS)

filterbank  : filterbank.o head_parse.o wapp2fb.o library 
	$(CC) -o $(BIN)/filterbank filterbank.o  $(LIB) $(LFFTW) -lm
	rm -f filterbank.o

extract  : extract.o library 
	$(CC) -o $(BIN)/extract extract.o  $(LIB) -lm
	rm -f extract.o

fake  : fake.o library 
	$(CC) -o $(BIN)/fake fake.o  $(LIB) -lm
	rm -f fake.o

gmrt2fil  : gmrt2fil.o library 
	$(CC) -o $(BIN)/gmrt2fil gmrt2fil.o  $(LIB) -lm
	rm -f gmrt2fil.o

tune  : tune.o library
	$(CC)  -c tune.c 
	$(F77) -o $(BIN)/tune tune.o $(LIB) -lm $(LFITS) $(LPGPLOT)
	rm -f tune.o
flatten  : flatten.o library 
	$(CC) -o $(BIN)/flatten flatten.o  $(LIB) -lm $(LFITS)
	rm -f flatten.o
 
clip  : clip.o library 
	$(CC) -o $(BIN)/clip clip.o  $(LIB) -lm $(LFITS)
	rm -f clip.o
downsample  : downsample.o library 
	$(CC) -o $(BIN)/downsample downsample.o  $(LIB) -lm $(LFITS)
	rm -f downsample.o



fold  : fold.o library 
	$(CC) -o $(BIN)/fold fold.o  $(LIB) -lm $(LFITS)
	rm -f fold.o

profile  : profile.o library 
	$(CC) -o $(BIN)/profile profile.o  $(LIB) -lm
	rm -f profile.o

flux  : flux.o library 
	$(CC) -o $(BIN)/flux flux.o  $(LIB) -lm
	rm -f flux.o

header  : header.o library
	$(CC) -o $(BIN)/header header.o $(LIB) -lm
	rm -f header.o

bandpass  : bandpass.o library
	$(CC) -o $(BIN)/bandpass bandpass.o  $(LIB) -lm
	rm -f bandpass.o

reader  : reader.o library
	$(CC) -o $(BIN)/reader reader.o $(LIB) -lm
	rm -f reader.o

readchunk  : readchunk.o library
	$(CC) -o $(BIN)/readchunk readchunk.o $(LIB) -lm
	rm -f readchunk.o

splitter  : splitter.o library
	$(CC) -o $(BIN)/splitter splitter.o $(LIB) -lm
	rm -f splitter.o

snrdm  : snrdm.o library
	$(CC) -o $(BIN)/snrdm snrdm.o $(LIB) -lm
	rm -f snrdm.o

decimate  : decimate.o library
	$(CC) -o $(BIN)/decimate decimate.o $(LIB) -lm
	rm -f decimate.o

splice  : splice.o library
	$(CC) -o $(BIN)/splice splice.o $(LIB) -lm
	rm -f splice.o

dice  : dice.o library
	$(CC) -o $(BIN)/dice dice.o $(LIB) -lm
	rm -f dice.o

depolyco  : depolyco.o library 
	$(CC) -o $(BIN)/depolyco depolyco.o $(LIB) -lm  $(LFITS) $(LFFTW)
	rm -f depolyco.o

blanker  : blanker.o library 
	$(CC) -o $(BIN)/blanker blanker.o $(LIB) -lm 
	rm -f blanker.o

clipper  : clipper.o library 
	$(CC) -o $(BIN)/clipper clipper.o $(LIB) -lm 
	rm -f clipper.o

brute  : brute.o library 
	$(CC) -o $(BIN)/brute brute.o $(LIB) -lm 
	rm -f brute.o

giant  : giant.o library 
	$(CXX) -o $(BIN)/giant giant.o $(LIB) $(LPGPLOT) -lm 
	rm -f giant.o

barycentre  : barycentre.o library 
	$(CC) -o $(BIN)/barycentre barycentre.o $(LIB) -lm  $(LFITS) $(LFFTW)
	rm -f barycentre.o

dedisperse  : dedisperse.o library
	$(CC) -o $(BIN)/dedisperse dedisperse.o $(LIB) -lm $(SUNLM)
	rm -f dedisperse.o

tree  : tree.o library
	$(CC) -o $(BIN)/tree tree.o $(LIB) -lm $(SUNLM)
	rm -f tree.o

polyco2period  : polyco2period.o library
	$(CC) -o $(BIN)/polyco2period polyco2period.o $(LIB) -lm $(LFITS)
	rm -f polyco2period.o

monitor :
	./mkmonitor.csh > $(BIN)/monitor
	chmod +x $(BIN)/monitor

pgplotter: pgplotter.o library
	$(CC) -c pgplotter.c
	$(F77) -o $(BIN)/pgplotter pgplotter.o $(LIB) $(LPGPLOT) -lm

polyco :
	echo "#!`which tclsh`" >  $(BIN)/polyco
	cat polyco.tcl         >> $(BIN)/polyco
	chmod +x $(BIN)/polyco

makedmlist :
	echo "#!`which tclsh`" >             $(BIN)/makedmlist
	cat makedmlist.tcl     >>            $(BIN)/makedmlist
	chmod +x $(BIN)/makedmlist

quicklook :
	echo "#!`which csh`" >               $(BIN)/quicklook
	cat quicklook.csh    >>              $(BIN)/quicklook
	chmod +x $(BIN)/quicklook

foldsignals :
	echo "#!`which csh`" >            $(BIN)/foldsignals
	echo "set bin        = $(BIN)" >> $(BIN)/foldsignals
	cat foldsignals.csh            >> $(BIN)/foldsignals
	chmod +x $(BIN)/foldsignals

seek : seek.o library
	$(LINK.f) -o $(BIN)/seek seek.o $(LIB) $(LFFTW)
	rm -f seek.o

sumfft : sumfft.o library
	$(LINK.f) -o $(BIN)/sumfft sumfft.o $(LIB)
	rm -f sumfft.o

peak : peak.o library
	$(LINK.f) -o $(BIN)/peak peak.o $(LIB) $(LPGPLOT)
	rm -f peak.o

ffa : ffa.o library
	$(LINK.f) -o $(BIN)/ffa ffa.o $(LIB) 
	rm -f ffa.o

best : best.o library
	$(LINK.f) -o $(BIN)/best best.o $(LIB) $(LPGPLOT)
	rm -f best.o

mask : mask.o library
	$(LINK.f) -o $(BIN)/mask mask.o $(LIB) $(LPGPLOT)
	rm -f mask.o

spec : spec.o library
	$(LINK.f) -o $(BIN)/spec spec.o $(LIB) $(LPGPLOT)
	rm -f spec.o

grey : grey.o library
	$(LINK.f) -o $(BIN)/grey grey.o $(LIB) $(LPGPLOT)
	rm -f grey.o

chaninfo : chaninfo.o 
	$(LINK.f) -o $(BIN)/chaninfo chaninfo.o
	rm -f chaninfo.o

step : step.o 
	$(LINK.f) -o $(BIN)/step step.o
	rm -f step.o

postproc : postproc.o 
	$(LINK.f) -o $(BIN)/postproc postproc.o
	rm -f postproc.o

hunt :
	echo "#!`which csh`" >             $(BIN)/hunt
	echo "set bin        = $(BIN)"  >> $(BIN)/hunt
	cat hunt.csh                    >> $(BIN)/hunt
	chmod +x $(BIN)/hunt

ahunt :
	echo "#!`which csh`" >             $(BIN)/ahunt
	echo "set bin        = $(BIN)"  >> $(BIN)/ahunt
	cat ahunt.csh                    >> $(BIN)/ahunt
	chmod +x $(BIN)/ahunt

csearch :
	echo "#!`which csh`" >             $(BIN)/csearch
	echo "set bin        = $(BIN)"  >> $(BIN)/csearch
	cat csearch.csh                 >> $(BIN)/csearch
	chmod +x $(BIN)/csearch

accn :
	echo "#!`which csh`" >             $(BIN)/accn
	echo "set bin        = $(BIN)"  >> $(BIN)/accn
	cat accn.csh                    >> $(BIN)/accn
	chmod +x $(BIN)/accn

makedummy :
	cp makedummy.csh $(BIN)/makedummy
	chmod +x $(BIN)/makedummy

quickplot : quickplot.o 
	$(F77) -o $(BIN)/quickplot quickplot.o $(LPGPLOT)

plotpulses : plotpulses.o 
	$(F77) -o $(BIN)/plotpulses plotpulses.o $(LPGPLOT)

getpulse : getpulse.o library
	$(CC) -c getpulse.c
	$(F77) -o $(BIN)/getpulse getpulse.o $(LIB) $(LPGPLOT) -lm

head_parse.o: head_parse.c mkheaderlex.c
	$(CC) -I. -I/opt/local/include -D$(OSTYPE) -c head_parse.c
	ar rv $(LIB) head_parse.o
	rm -f head_parse.o

wapp2fb.o: wapp2fb.c 
	$(CC) -I. -D$(OSTYPE) -c wapp2fb.c
	ar rv $(LIB) wapp2fb.o
	rm -f wapp2fb.o

alfa2fb.o: alfa2fb.c 
	$(CC) -I. -D$(OSTYPE) -c alfa2fb.c
	ar rv $(LIB) alfa2fb.o
	rm -f alfa2fb.o

y.tab.c: mkheader.y
	yacc -vdt mkheader.y

mkheaderlex.c: mkheader.l y.tab.c
	lex -t mkheader.l >mkheaderlex.c

clean :
	rm -f $(LIB) *~ *.o sigproc.h \
	sigproc.aux sigproc.dvi sigproc.log sigproc.ps sigproc.toc \
	sigproc.idx sigproc.ilg sigproc.ind

help :
	./help.csh $(BIN)/filterbank
	./help.csh $(BIN)/fake
	./help.csh $(BIN)/header
	./help.csh $(BIN)/bandpass
	./help.csh $(BIN)/reader
	./help.csh $(BIN)/decimate
	./help.csh $(BIN)/splice
	./help.csh $(BIN)/dedisperse
	./help.csh $(BIN)/fold
	./help.csh $(BIN)/profile
	./help.csh $(BIN)/flux
	./help.csh $(BIN)/polyco
	./help.csh $(BIN)/barycentre
	./help.csh $(BIN)/quicklook

documentation: help
	latex sigproc
	makeindex sigproc
	latex sigproc
	latex sigproc
	dvipdf sigproc

export :
	./exporter.csh

###############################################################################
