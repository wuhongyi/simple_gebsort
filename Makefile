

PLTFMOPT=Makefile.$(shell uname)
include $(PLTFMOPT)

#                  +------- link to where EPICS is on this machine
#                  |        see below
EPICS_ROOT="curEPICS"

# at ANL we use
#rm curEPICS; ln -s /global/devel6_sl62/base/base-3.14.12.1 curEPICS
#or maybe
#rm curEPICS; ln -s /global/devel6_sl62/base/R3.14.10 curEPICS

# at gretina @ MSU we use
#rm curEPICS; ln -s /global/devel/base/R3.14.11 curEPICS

EPICS_LIBSUBDIR=linux-x86
EPICS_INCSUBDIR=Linux

#------------------------------------------ end modifications ---

EPICSLIBDIR=$(EPICS_ROOT)/lib
EPICSARCHLIBDIR=$(EPICSLIBDIR)/$(EPICS_LIBSUBDIR)
EPICSINCDIR=$(EPICS_ROOT)/include
EPICSOSINCDIR=$(EPICS_ROOT)/include/os/$(EPICS_INCSUBDIR)

CFLAGS=-I$(EPICSINCDIR) -I$(EPICSOSINCDIR) -I.

EPICSFLAGS=-L$(EPICSARCHLIBDIR) -lCom -Wl,"-rpath=$(EPICSARCHLIBDIR)"

#-------------------

all:		
		make GEBSort_nogeb
#		make G4
#		make mkMap
#		make GEBMerge
#		go lr_60Co_ascii
offline:	GEBSort_nogeb
online:		GEBSort
G4:		
		make G4toMode2_EU
		make G4toMode2_USa_1
		make G4toMode2_USa_2
		make G4toMode2_USa_3
		make G4toMode2_USb

#-------------------


GEBHeader.c:	GTMerge.h gdecomp.h
GTPrint.c:	GTMerge.h gdecomp.h GEBSort.h
GTUtilsCC.c:	GTMerge.h gdecomp.h

GEBSort:	GTPrintCC.o GEBSort.o time_stampCC.o GEBHeaderCC.o get_a_seedCC.o jta.o \
		bin_mode2.o bin_mode1.o bin_mode3.o bin_template.o bin_final.o bin_dgs.o bin_XA.o bin_dfma.o spe_funCC.o \
		GEBClient.o gretTapClient.o bin_gtcal.o str_decompCC.o validate.o bin_linpol.o \
		bin_angcor_GT.o bin_angcor_DGS.o bin_angdis.o \
                utilsCC.o
		$(CC) $^ -DHAVE_VXWORKS=0 $(ROOTLIB) $(EPICSFLAGS) -o $@  -lNew
#		go sim_map
#		go sim_root

#-------------------
# Channel Access (CA)
# must be compile and used on 32 bit Linux
# __ because that is where EPICS is compiled

test_CA:	test_CA.c
		$(cc)  -IcurEPICS/include -IcurEPICS/include/os/Linux/ \
		-LcurEPICS/lib/linux-x86  -lca -lCom -lm -o $@ $^
		$@ tank_display


temp_ge:	temp_ge.c
		$(cc)  -IcurEPICS/include -IcurEPICS/include/os/Linux/ \
		-LcurEPICS/lib/linux-x86  -lca -lCom -lm -o $@ $^
#		$@ > x.html; scp x.html tlauritsen@www.phy.anl.gov:webgroups/gammasphere/ln_snapshot_new.html
#		$@

#-------------------
		

GEBSort_test:	GEBSort_test.cxx
		$(CC) $^ $(ROOTLIB) -o $@ -lz -lNew
		./$@


GEBSort.o:	GEBSort.cxx GTMerge.h GTmode3.h GEBSort.h GEBClient.h GEBLink.h gretTapClient.h\
		UserInclude.h UserChat.h UserStat.h UserFunctions.h \
		UserDeclare.h UserInit.h UserRawEv.h UserGoodEv.h \
		UserPreCond.h UserEv.h UserExit.h gdecomp.h veto_pos.h \
		permtable1.h permtable2.h permtable3.h permtable4.h \
		permtable5.h permtable6.h permtable7.h permtable8.h 
		$(CC) -DHAVE_VXWORKS=0  $(ROOTINC) $(EPICSFLAGS) $(CFLAGS) -c GEBSort.cxx

GEBSort_nogeb:	GTPrintCC.o GEBSort_nogeb.o time_stampCC.o GEBHeaderCC.o get_a_seedCC.o jta.o \
		bin_dgs.o bin_XA.o  bin_dfma.o\
		spe_funCC.o bin_gtcal.o str_decompCC.o validate.o bin_angcor_GT.o bin_angcor_DGS.o bin_angdis.c\
                utilsCC.o trig_funCC.o LinkDict.o
		$(CC) $^ -DHAVE_VXWORKS=0 $(ROOTLIB)  -o $@ 
#		./GEBSort_nogeb -input disk GTDATA/g4_mode1.out -rootfile GTDATA/test.root RECREATE -chat GEBSort.chat 
#		go mode3

get_a_seedCC.o:	get_a_seed.c
		$(CC)  -c  -o $@ $^

LinkDict.o : LinkDict.cxx
	    $(CC) -c $(ROOTINC) LinkDict.cxx -o LinkDict.o  

LinkDict.cxx : GTMerge.h Linkdef.h
	    #root 6
		#rootcling -f LinkDict.cxx -c $^
		# root 5
		rootcint -f LinkDict.cxx -c $^

str_decompCC.o:	str_decomp.c
		$(CC)  -c  -o $@ $^

utilsCC.o:	utils.c
		$(CC)  -c  -o $@ $^

util.o:		utils.c
		$(cc)  -c  -o $@ $^

str_decomp.c:	ctk.h

GEBSort.cxx:	GTMerge.h GEBSort.h \
		UserInclude.h UserChat.h UserStat.h UserFunctions.h \
		UserDeclare.h UserInit.h UserRawEv.h UserGoodEv.h \
		UserPreCond.h UserEv.h UserExit.h gdecomp.h  veto_pos.h \
		permtable1.h permtable2.h permtable3.h permtable4.h \
		permtable5.h permtable6.h permtable7.h permtable8.h

GEBSort_nogeb.o:	GEBSort.cxx 
		$(CC) -DHAVE_VXWORKS=0  $(ROOTINC)  $(CFLAGS) -c GEBSort.cxx -o GEBSort_nogeb.o

GEBSort.cxx:	gdecomp.h veto_pos.h GEBSort.h GTMerge.h 

GSUtil:		GSUtil.cc spe_funCC.o 2d_funCC.o
		$(cc) $(ROOTINC) $(CFLAGS) -o $@ $^
		
time_stampCC.o:	time_stamp.c
		$(CC)  -c  -o time_stampCC.o time_stamp.c  

GEBHeaderCC.o:	GEBHeader.c
		$(CC) -c $(ROOTINC) -o  $@ $^

2d_funCC.o:	2d_fun.c
		$(cc) -c -o $@ $^

GTPrintCC.o:	GTPrint.c
		$(CC) -c -o $@  $^

jta.o:	jta.c
		$(CC) -c -o $@  $^

bin_mode3.o:	bin_mode3.c
		$(CC) -c -o $@  $^

bin_gtcal.o:	bin_gtcal.c
		$(CC) -c -o $@  $^

bin_mode2.o:	bin_mode2.c
		$(CC) -c -o $@  $^

bin_mode2.c:	veto_pos.h

bin_mode1.o:	bin_mode1.c
		$(CC) -c -o $@  $^

bin_template.o:	bin_template.c
		$(CC) -c -o $@  $^

bin_final.o:	bin_final.c
		$(CC) -c -o $@  $^

bin_angcor_GT.o:	bin_angcor_GT.c
		$(CC) -c -o $@  $^

bin_angcor_DGS.o:	bin_angcor_DGS.c
		$(CC) -c -o $@  $^

bin_angdis.o:	bin_angdis.c
		$(CC) -c -o $@  $^

trig_funCC.o:	trig_fun.c
		$(CC) $(CFLAGS) -c -o $@  $^

bin_linpol.o:	bin_linpol.c 
		$(CC) -c -o $@  $^ 

validate.o:	validate.c
		$(CC) -c -o $@  $^

bin_dgs.o:	bin_dgs.c
		$(CC) -c $(ROOTINC)  -o $@ $^

bin_XA.o:	bin_XA.c
		$(CC) -c -o $@ $^

bin_dfma.o:	bin_dfma.c
		$(CC) -c -o $@  $^

GTPrint.o:	GTPrint.c GTMerge.h
		$(CC) -c $^

time_stamp.o:	time_stamp.c
		$(cc)  -c  $^

GTUtilsCC.o:	GTUtils.c
		$(cc) -c $^

GEBClient.o:	GEBClient.c 
		$(CC) $(EPICSFLAGS) $(CFLAGS) -c -o $@  $^

GEBClient.c:	GEBClient.h

gretTapClient.o:	gretTapClient.c 
			$(CC) $(EPICSFLAGS) $(CFLAGS) -c -o $@  $^

gretTapClient.c:	GEBClient.h gretTapClient.h

GEBSumGate:	GEBSumGate.c getEvent.o writeTrack_repeat.o printEvent.o
		$(cc) -o $@ $^ -lm
#		$@ /media/20140424_1131/user/gsfma315_GT/Run0239/Global.dat /media/r20150911_1334/user/166Ho/s1787/Global_1787.dat 1700 1800

GEBbunch:	GEBbunch.c getEvent.o writeTrack_repeat.o printEvent.o
		$(cc) -o $@ $^ -lm
#		$@ /media/r20150911_1334/user/20140424_1131/gsfma315_GT/Run0240/merged.gtd_000 MM 30 1831 1840 20000

GEBmkhkrsp:	GEBmkhkrsp.c spe_fun.o
		$(cc) -o $@ $^ -lm
		$@

GEBcombine:	GEBcombine.c getEvent.o  printEvent.o
		$(cc) -o $@ $^ -lm
#		$@ GTDATA/MM5.dat GTDATA/MM5_c.dat 30 5
#		ls -l GTDATA/MM5.dat GTDATA/MM5_c.dat

#-------------------

str_decomp.o:	str_decomp.c 
		$(cc) -c  $^

printEvent.c:	ctk.h
printEvent.o:	printEvent.c 
		$(cc) -c $^

#-------------------


GEBMerge.c:	GTMerge.h GTMerge_readnew.h pbuf.h ctk.h gdecomp.h
GEBMerge:	GEBMerge.c spe_fun.o \
		get_a_seed.o time_stamp.o\
		spe_fun.o printEvent.o str_decomp.o
		$(cc) -o $@ -lm -lz $^
#		rm core.*; GEBMerge gtmerge.chat DATA/merged.gtd DATA/run_39.chico_000		

#-------------------

GEBFilter:	GEBFilter.c GF_veto_cube.o get_a_seed.o
		$(cc) -o $@ -lm -lz $^
#		$@  GEBFilter.chat GTDATA/merged.gtd_000 filtertest.gtd

GEBFilter.c:	ctk.h gdecomp.h GEBFilter.h

GF_veto_cube.c:	ctk.h gdecomp.h veto_pos.h GEBFilter.h

GF_veto_cube.o:	GF_veto_cube.c
		$(cc) -c $^

#-------------------
# /media/r20150911_1334/user/AG_mode3/event_mezzdata.cdat.0003
# ^^^^ has compressed traces
# /media/d20150930_1152/user/AG2GT/1/event_mezzdata.bdat.0000
# ^^^^ has uncompressed traces

AG_mode3_to_GT:	AG_mode3_to_GT.c util.o
		$(cc) -o $@ -lm  $^
#		$@ /media/r20150911_1334/user/AG_mode3/event_mezzdata.cdat.0003 
#		$@ /media/d20150930_1152/user/AG2GT/1/event_mezzdata.bdat.0000
#		$@ /media/20171129_1623/user/AG2GT/1/event_mezzdata.bdat.0000  > x.log
#		go AG_mode3_to_GT

#-------------------

decimate:	decimate.c
		$(cc) -o $@ $^
#		$@ read

listTS:		listTS.c gdecomp.h
		$(cc) -o $@ $^
#		$@ /media/120621a/user/gtdata/anl1/Global.dat 1000 300

# EU/AGATA G4 ascii data
# put -g on G4toMode2_EU to overcome unknown problem on mac

G4toMode2_EU:		G4toMode2.c spe_fun.o printEvent.o get_a_seed.o str_decomp.o
			$(cc) -g -D USGEANTFORMAT=0 -D ASCIIDATA=1 -D NEWFORMAT=1 -o $@ $^ -lm

# LR UCGretina G4 ascii data, format 1 is old and 'new format 2 as of july 2014'

G4toMode2_USa_1:	G4toMode2.c spe_fun.o printEvent.o get_a_seed.o str_decomp.o
			$(cc) -D USGEANTFORMAT=1 -D ASCIIDATA=1 -D NEWFORMAT=1 -o $@ $^ -lm

G4toMode2_USa_2:	G4toMode2.c spe_fun.o printEvent.o get_a_seed.o str_decomp.o
			$(cc) -D USGEANTFORMAT=1 -D ASCIIDATA=1 -D NEWFORMAT=2 -o $@ $^ -lm
#			go lr_60Co_ascii.
#			$@  GTDATA/g4.dat  GTDATA/g4_mode2.out  G4toMode2.chat   1000000
			
G4toMode2_USa_3:	G4toMode2.c spe_fun.o printEvent.o get_a_seed.o str_decomp.o
			$(cc) -D USGEANTFORMAT=1 -D ASCIIDATA=1 -D NEWFORMAT=3 -o $@ $^ -lm
#			$@ /home/tl/tmp/g4.dat /home/tl/tmp/4_mode2.out  G4toMode2.chat   1000000

# LR UCGretina G4, writing mode 2 format

G4toMode2_USb:	G4toMode2.c spe_fun.o printEvent.o get_a_seed.o str_decomp.o
		$(cc) -D USGEANTFORMAT=1 -D ASCIIDATA=0 -D NEWFORMAT=1 -o $@ $^ -lm

spe_fun.o:	spe_fun.c
		$(cc) -c $^

spe_funCC.o:	spe_fun.c
		$(CC) -c -o spe_funCC.o  $^

get_a_seed.o:	get_a_seed.c 
		$(cc) -c  $^

G4toMode2.c:	gdecomp.h

#-------------------

zzip:		zzip.c
		$(cc) -o $@ $^ -lz

zunzip:		zunzip.c
		$(cc) -o $@ $^ -lz

mkMap:		mkMap.c
		$(cc) -o $@ $^ -lz

inv3x3:		inv3x3.c
		$(cc) -o $@ $^ -lm
		$@

#-------------------

format_2_3:	format_2_3.c
		$(cc) -o $@ $^

format_1_3:	format_1_3.c
		$(cc) -o $@ $^
#		$@ GTDATA/Global.dat Global_mod.dat

#-------------------

dgs_pz:		dgs_pz.c spe_fun.o
		$(cc) -o $@ $^ -lm
#		$@ 350 141 dgs_pz.cal 1.0 | grep pz096

dgs_pz.c:	GTMerge.h

dgs_ecal:	dgs_ecal.c spe_fun.o
		$(cc) -o $@ $^ -lm
#		$@ dgs_ehi.cal 207Bi
dgs_ecal2:	dgs_ecal2.c spe_fun.o pt_1sp.c
		$(cc) -o $@ $^ -lm
#		$@ dgs_ehi.cal 88Y 600 1.0

dgs_ecal.c:	GTMerge.h

#-------------------

clean:		

		rm -f *.o
		rm -f core.*
		rm -f *%
		rm -f a.out
		rm -f *.gch
		rm -f G4toMode2_EU G4toMode2_US GEBMerge GEBSort 
		rm -f GEBSort_nogeb listTS mkMap

debug:		
#		gdb -d ./ G4toMode2_EU `ls core.*`
#		gdb -d ./ GEBSort_nogeb `ls core.*`
#		gdb -d ./ GEBMerge `ls core.*`
#		gdb -d ./ GEBbunch `ls core.*`
#		gdb -d ./ AG_mode3_to_GT `ls core.*`
		gdb -d ./ dgs_ecal2 `ls core.*`


