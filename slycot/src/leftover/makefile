####################################################################
#  SLICOT routines makefile                                        #
#  Makefile for creating/updating the SLICOT Library object file   # 
#  on Unix machines.                                               #
#  SLICOT, Release 5.0                      ./slicot/src/makefile  #
#  Vasile Sima, KU Leuven                                          #
#  October 31, 1996.                                               #
#  Revised June 25, 1998.                                          # 
#  Revised December 7, 1999, September 5, 2003, Jan. 9, 2009.      #
####################################################################
#
#  This is the makefile to create/update the library for SLICOT.
#  The SLICOT Library routines are written for double precision only.
#
#  The command
#       make
#  without any arguments creates or updates a library called
#       slicot.a
#  in the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make double FRC=FRC
#
#######################################################################
 
include ../make.inc

DSLSRC = \
    AB01MD.o AB01ND.o AB01OD.o AB04MD.o AB05MD.o AB05ND.o AB05OD.o \
    AB05PD.o AB05QD.o AB05RD.o AB05SD.o AB07MD.o AB07ND.o AB08MD.o \
    AB08ND.o AB08NX.o AB09AD.o AB09AX.o AB09BD.o AB09BX.o AB09CD.o \
    AB09CX.o AB09DD.o AB09ED.o AB09FD.o AB09GD.o AB09HD.o AB09HX.o \
    AB09HY.o AB09ID.o AB09IX.o AB09IY.o AB09JD.o AB09JV.o AB09JW.o \
    AB09JX.o AB09KD.o AB09KX.o AB09MD.o AB09ND.o AB13AD.o AB13AX.o \
    AB13BD.o AB13CD.o AB13DD.o AB13DX.o AB13ED.o AB13FD.o AB13MD.o \
    AG07BD.o AG08BD.o AG08BY.o \
    BB01AD.o BB02AD.o BB03AD.o BB04AD.o BD01AD.o BD02AD.o \
    DE01OD.o DE01PD.o DF01MD.o DG01MD.o DG01ND.o DG01NY.o DG01OD.o \
    DK01MD.o \
    FB01QD.o FB01RD.o FB01SD.o FB01TD.o FB01VD.o FD01AD.o \
    IB01AD.o IB01BD.o IB01CD.o IB01MD.o IB01MY.o IB01ND.o IB01OD.o \
    IB01OY.o IB01PD.o IB01PX.o IB01PY.o IB01QD.o IB01RD.o \
    IB03AD.o IB03BD.o \
    MA01AD.o MA02AD.o MA02BD.o MA02CD.o MA02DD.o MA02ED.o MA02FD.o \
    MA02GD.o MA02HD.o MB01PD.o MB01QD.o MB01RD.o MB01RU.o MB01RW.o \
    MB01RX.o MB01RY.o MB01SD.o MB01TD.o MB01UD.o MB01UW.o MB01VD.o \
    MB01WD.o MB01XD.o MB01XY.o MB01YD.o MB01ZD.o MB02CD.o MB02CU.o \
    MB02CV.o MB02CX.o MB02CY.o MB02DD.o MB02ED.o MB02FD.o MB02GD.o \
    MB02HD.o MB02ID.o MB02JD.o MB02JX.o MB02KD.o MB02MD.o MB02ND.o \
    MB02NY.o MB02OD.o MB02PD.o MB02QD.o MB02QY.o MB02RD.o MB02RZ.o \
    MB02SD.o MB02SZ.o MB02TD.o MB02TZ.o MB02UD.o MB02UU.o MB02UV.o \
    MB02VD.o MB02WD.o MB02XD.o MB02YD.o \
    MB03MD.o MB03MY.o MB03ND.o MB03NY.o MB03OD.o MB03OY.o \
    MB03PD.o MB03PY.o MB03QD.o MB03QX.o MB03QY.o MB03RD.o MB03RX.o \
    MB03RY.o MB03SD.o MB03UD.o MB03VD.o MB03VY.o MB03WD.o MB03WX.o \
    MB04DY.o MB04GD.o MB04ID.o MB04IY.o MB04JD.o MB04KD.o MB04LD.o \
    MB04MD.o MB04ND.o MB04NY.o MB04OD.o MB04OW.o MB04OX.o MB04OY.o \
    MB04PY.o MB04TT.o MB04TU.o MB04TV.o MB04TW.o MB04TX.o MB04TY.o \
    MB04UD.o MB04VD.o MB04VX.o MB04XD.o MB04XY.o MB04YD.o MB04YW.o \
    MB04ZD.o MB05MD.o MB05MY.o MB05ND.o MB05OD.o MB05OY.o MC01MD.o \
    MC01ND.o MC01OD.o MC01PD.o MC01PY.o MC01QD.o MC01RD.o MC01SD.o \
    MC01SW.o MC01SX.o MC01SY.o MC01TD.o MC01VD.o MC01WD.o MC03MD.o \
    MC03ND.o MC03NX.o MC03NY.o \
    MD03AD.o MD03BA.o MD03BB.o MD03BD.o MD03BF.o MD03BX.o MD03BY.o \
    NF01AD.o NF01AY.o NF01BA.o NF01BB.o NF01BD.o NF01BE.o NF01BF.o \
    NF01BP.o NF01BQ.o NF01BR.o NF01BS.o NF01BU.o NF01BV.o NF01BW.o \
    NF01BX.o NF01BY.o \
    SB01BD.o SB01BX.o SB01BY.o SB01DD.o SB01FY.o SB01MD.o SB02CX.o \
    SB02MD.o SB02MR.o SB02MS.o SB02MT.o SB02MU.o SB02MV.o SB02MW.o \
    SB02ND.o SB02OD.o SB02OU.o SB02OV.o SB02OW.o SB02OX.o SB02OY.o \
    SB02PD.o SB02QD.o SB02RD.o SB02RU.o SB02SD.o SB03MD.o SB03MU.o \
    SB03MV.o SB03MW.o SB03MX.o SB03MY.o SB03OD.o SB03OR.o SB03OT.o \
    SB03OU.o SB03OV.o SB03OY.o SB03PD.o SB03QD.o SB03QX.o SB03QY.o \
    SB03RD.o SB03SD.o SB03SX.o SB03SY.o SB03TD.o SB03UD.o SB04MD.o \
    SB04MR.o SB04MU.o SB04MW.o SB04MY.o SB04ND.o SB04NV.o SB04NW.o \
    SB04NX.o SB04NY.o SB04OD.o SB04PD.o SB04PX.o SB04PY.o SB04QD.o \
    SB04QR.o SB04QU.o SB04QY.o SB04RD.o SB04RV.o SB04RW.o SB04RX.o \
    SB04RY.o SB06ND.o SB08CD.o SB08DD.o SB08ED.o SB08FD.o SB08GD.o \
    SB08HD.o SB08MD.o SB08MY.o SB08ND.o SB08NY.o SB09MD.o SB10AD.o \
    SB10DD.o SB10ED.o SB10FD.o SB10HD.o SB10ID.o SB10JD.o SB10KD.o \
    SB10MD.o SB10LD.o SB10PD.o SB10QD.o SB10RD.o SB10SD.o SB10TD.o \
    SB10UD.o SB10VD.o SB10WD.o SB10YD.o SB10ZD.o SB10ZP.o SB16AD.o \
    SB16AY.o SB16BD.o SB16CD.o SB16CY.o SG02AD.o SG03AD.o SG03AX.o \
    SG03AY.o SG03BD.o SG03BU.o SG03BV.o SG03BW.o SG03BX.o SG03BY.o \
    TB01ID.o TB01KD.o TB01LD.o TB01MD.o TB01ND.o TB01PD.o TB01TD.o \
    TB01TY.o TB01UD.o TB01VD.o TB01VY.o TB01WD.o TB01XD.o TB01YD.o \
    TB01ZD.o TB03AD.o TB03AY.o TB04AD.o TB04AY.o TB04BD.o TB04BV.o \
    TB04BW.o TB04BX.o TB04CD.o TB05AD.o TC01OD.o TC04AD.o TC05AD.o \
    TD03AD.o TD03AY.o TD04AD.o TD05AD.o TF01MD.o TF01MX.o TF01MY.o \
    TF01ND.o TF01OD.o TF01PD.o TF01QD.o TF01RD.o TG01AD.o TG01BD.o \
    TG01CD.o TG01DD.o TG01ED.o TG01FD.o TG01HD.o TG01HX.o TG01ID.o \
    TG01JD.o TG01WD.o \
    UD01BD.o UD01CD.o UD01DD.o UD01MD.o UD01ND.o select.o delctg.o \
    MA02ID.o MA02JD.o MB01MD.o MB01ND.o MB01UX.o MB03TD.o MB03TS.o \
    MB03WA.o MB03XD.o MB03XP.o MB03XU.o MB03YA.o MB03YD.o MB03YT.o \
    MB03ZA.o MB03ZD.o MB04DD.o MB04DI.o MB04DS.o MB04PA.o MB04PB.o \
    MB04PU.o MB04QB.o MB04QC.o MB04QF.o MB04QU.o MB04TB.o MB04TS.o \
    MB04WD.o MB04WP.o MB04WR.o MB04WU.o SB04OW.o UE01MD.o UD01MZ.o \
    AB08MZ.o AB08NZ.o AB8NXZ.o AG08BZ.o AG8BYZ.o MA02BZ.o MA02CZ.o \
    MB04IZ.o MB3OYZ.o MB3PYZ.o TB01IZ.o TB01XZ.o TG01AZ.o TG01FZ.o

all: double

double: $(DSLSRC)
	$(ARCH) $(ARCHFLAGS) $(SLICOTLIB) $(DSLSRC)

$(DSLSRC): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

.f.o: 
	$(FORTRAN) $(OPTS) -c $<
