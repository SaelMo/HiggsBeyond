      SUBROUTINE ML5_0_0_1_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=5)
      INTEGER    NCOMB
      PARAMETER (NCOMB=4)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=276, NLOOPGROUPS=276, NCTAMPS=78)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=354)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=44,NLOOPWAVEFUNCS=416)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*16 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)

      LOGICAL DUMMYFALSE
      DATA DUMMYFALSE/.FALSE./
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
      INCLUDE 'mp_coupl.inc'

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/ML5_0_0_1_FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/ML5_0_0_1_INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/ML5_0_0_1_SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_0_0_1_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/ML5_0_0_1_I_SO/I_SO
      INTEGER I_LIB
      COMMON/ML5_0_0_1_I_LIB/I_LIB

      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/ML5_0_0_1_W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/ML5_0_0_1_WL/WL,PL

      COMPLEX*16 AMPL(3,NLOOPAMPS)
      COMMON/ML5_0_0_1_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL SXXXXX(P(0,4),+1*IC(4),W(1,4))
      CALL SXXXXX(P(0,5),+1*IC(5),W(1,5))
      CALL SSS1_1(W(1,3),W(1,4),GC_HHETA0,MDL_META,MDL_WETA,W(1,6))
      CALL SSS1_2(W(1,6),W(1,5),GC_HHETA0,MDL_MH,MDL_WH,W(1,7))
C     Counter-term amplitude(s) for loop diagram number 7
      CALL VVS1_0(W(1,1),W(1,2),W(1,7),R2_GGHB,AMPL(1,1))
      CALL SSS1_1(W(1,6),W(1,5),GC_HETA0ETA0,MDL_META,MDL_WETA,W(1,8))
C     Counter-term amplitude(s) for loop diagram number 8
      CALL VVS1_0(W(1,1),W(1,2),W(1,8),R2_GGETA0B,AMPL(1,2))
      CALL SSS1_2(W(1,6),W(1,5),GC_HIOTA0ETA0,MDL_MIOTA,MDL_WIOTA,W(1
     $ ,9))
C     Counter-term amplitude(s) for loop diagram number 9
      CALL VVS1_0(W(1,1),W(1,2),W(1,9),R2_GGIOTA0B,AMPL(1,3))
C     Counter-term amplitude(s) for loop diagram number 10
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,6),R2_GGHETA0B,AMPL(1,4))
      CALL SSS1_1(W(1,3),W(1,4),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1,10))
      CALL SSS1_2(W(1,10),W(1,5),GC_HHIOTA0,MDL_MH,MDL_WH,W(1,11))
C     Counter-term amplitude(s) for loop diagram number 12
      CALL VVS1_0(W(1,1),W(1,2),W(1,11),R2_GGHB,AMPL(1,5))
      CALL SSS1_1(W(1,10),W(1,5),GC_HIOTA0ETA0,MDL_META,MDL_WETA,W(1
     $ ,12))
C     Counter-term amplitude(s) for loop diagram number 13
      CALL VVS1_0(W(1,1),W(1,2),W(1,12),R2_GGETA0B,AMPL(1,6))
      CALL SSS1_1(W(1,10),W(1,5),GC_HIOTA0IOTA0,MDL_MIOTA,MDL_WIOTA
     $ ,W(1,13))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL VVS1_0(W(1,1),W(1,2),W(1,13),R2_GGIOTA0B,AMPL(1,7))
C     Counter-term amplitude(s) for loop diagram number 15
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,10),R2_GGHIOTA0B,AMPL(1
     $ ,8))
      CALL SSS1_1(W(1,3),W(1,4),GC_30,MDL_MH,MDL_WH,W(1,14))
      CALL SSS1_1(W(1,5),W(1,14),GC_30,MDL_MH,MDL_WH,W(1,15))
C     Counter-term amplitude(s) for loop diagram number 17
      CALL VVS1_0(W(1,1),W(1,2),W(1,15),R2_GGHB,AMPL(1,9))
      CALL SSS1_1(W(1,5),W(1,14),GC_HHETA0,MDL_META,MDL_WETA,W(1,16))
C     Counter-term amplitude(s) for loop diagram number 18
      CALL VVS1_0(W(1,1),W(1,2),W(1,16),R2_GGETA0B,AMPL(1,10))
      CALL SSS1_1(W(1,5),W(1,14),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1,17)
     $ )
C     Counter-term amplitude(s) for loop diagram number 19
      CALL VVS1_0(W(1,1),W(1,2),W(1,17),R2_GGIOTA0B,AMPL(1,11))
C     Counter-term amplitude(s) for loop diagram number 20
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,14),R2_GGHHB,AMPL(1,12))
      CALL SSS1_1(W(1,3),W(1,5),GC_HHETA0,MDL_META,MDL_WETA,W(1,18))
      CALL SSS1_2(W(1,18),W(1,4),GC_HHETA0,MDL_MH,MDL_WH,W(1,19))
C     Counter-term amplitude(s) for loop diagram number 22
      CALL VVS1_0(W(1,1),W(1,2),W(1,19),R2_GGHB,AMPL(1,13))
      CALL SSS1_1(W(1,18),W(1,4),GC_HETA0ETA0,MDL_META,MDL_WETA,W(1,20)
     $ )
C     Counter-term amplitude(s) for loop diagram number 23
      CALL VVS1_0(W(1,1),W(1,2),W(1,20),R2_GGETA0B,AMPL(1,14))
      CALL SSS1_2(W(1,18),W(1,4),GC_HIOTA0ETA0,MDL_MIOTA,MDL_WIOTA,W(1
     $ ,21))
C     Counter-term amplitude(s) for loop diagram number 24
      CALL VVS1_0(W(1,1),W(1,2),W(1,21),R2_GGIOTA0B,AMPL(1,15))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,18),R2_GGHETA0B,AMPL(1
     $ ,16))
      CALL SSS1_1(W(1,3),W(1,5),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1,22))
      CALL SSS1_2(W(1,22),W(1,4),GC_HHIOTA0,MDL_MH,MDL_WH,W(1,23))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL VVS1_0(W(1,1),W(1,2),W(1,23),R2_GGHB,AMPL(1,17))
      CALL SSS1_1(W(1,22),W(1,4),GC_HIOTA0ETA0,MDL_META,MDL_WETA,W(1
     $ ,24))
C     Counter-term amplitude(s) for loop diagram number 28
      CALL VVS1_0(W(1,1),W(1,2),W(1,24),R2_GGETA0B,AMPL(1,18))
      CALL SSS1_1(W(1,22),W(1,4),GC_HIOTA0IOTA0,MDL_MIOTA,MDL_WIOTA
     $ ,W(1,25))
C     Counter-term amplitude(s) for loop diagram number 29
      CALL VVS1_0(W(1,1),W(1,2),W(1,25),R2_GGIOTA0B,AMPL(1,19))
C     Counter-term amplitude(s) for loop diagram number 30
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,22),R2_GGHIOTA0B,AMPL(1
     $ ,20))
      CALL SSS1_1(W(1,3),W(1,5),GC_30,MDL_MH,MDL_WH,W(1,26))
      CALL SSS1_1(W(1,4),W(1,26),GC_30,MDL_MH,MDL_WH,W(1,27))
C     Counter-term amplitude(s) for loop diagram number 32
      CALL VVS1_0(W(1,1),W(1,2),W(1,27),R2_GGHB,AMPL(1,21))
      CALL SSS1_1(W(1,4),W(1,26),GC_HHETA0,MDL_META,MDL_WETA,W(1,28))
C     Counter-term amplitude(s) for loop diagram number 33
      CALL VVS1_0(W(1,1),W(1,2),W(1,28),R2_GGETA0B,AMPL(1,22))
      CALL SSS1_1(W(1,4),W(1,26),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1,29)
     $ )
C     Counter-term amplitude(s) for loop diagram number 34
      CALL VVS1_0(W(1,1),W(1,2),W(1,29),R2_GGIOTA0B,AMPL(1,23))
C     Counter-term amplitude(s) for loop diagram number 35
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,26),R2_GGHHB,AMPL(1,24))
      CALL SSS1_1(W(1,4),W(1,5),GC_HHETA0,MDL_META,MDL_WETA,W(1,30))
      CALL SSS1_2(W(1,30),W(1,3),GC_HHETA0,MDL_MH,MDL_WH,W(1,31))
C     Counter-term amplitude(s) for loop diagram number 37
      CALL VVS1_0(W(1,1),W(1,2),W(1,31),R2_GGHB,AMPL(1,25))
      CALL SSS1_1(W(1,30),W(1,3),GC_HETA0ETA0,MDL_META,MDL_WETA,W(1,32)
     $ )
C     Counter-term amplitude(s) for loop diagram number 38
      CALL VVS1_0(W(1,1),W(1,2),W(1,32),R2_GGETA0B,AMPL(1,26))
      CALL SSS1_2(W(1,30),W(1,3),GC_HIOTA0ETA0,MDL_MIOTA,MDL_WIOTA,W(1
     $ ,33))
C     Counter-term amplitude(s) for loop diagram number 39
      CALL VVS1_0(W(1,1),W(1,2),W(1,33),R2_GGIOTA0B,AMPL(1,27))
C     Counter-term amplitude(s) for loop diagram number 40
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,3),W(1,30),R2_GGHETA0B,AMPL(1
     $ ,28))
      CALL SSS1_1(W(1,4),W(1,5),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1,34))
      CALL SSS1_2(W(1,34),W(1,3),GC_HHIOTA0,MDL_MH,MDL_WH,W(1,35))
C     Counter-term amplitude(s) for loop diagram number 42
      CALL VVS1_0(W(1,1),W(1,2),W(1,35),R2_GGHB,AMPL(1,29))
      CALL SSS1_1(W(1,34),W(1,3),GC_HIOTA0ETA0,MDL_META,MDL_WETA,W(1
     $ ,36))
C     Counter-term amplitude(s) for loop diagram number 43
      CALL VVS1_0(W(1,1),W(1,2),W(1,36),R2_GGETA0B,AMPL(1,30))
      CALL SSS1_1(W(1,34),W(1,3),GC_HIOTA0IOTA0,MDL_MIOTA,MDL_WIOTA
     $ ,W(1,37))
C     Counter-term amplitude(s) for loop diagram number 44
      CALL VVS1_0(W(1,1),W(1,2),W(1,37),R2_GGIOTA0B,AMPL(1,31))
C     Counter-term amplitude(s) for loop diagram number 45
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,3),W(1,34),R2_GGHIOTA0B,AMPL(1
     $ ,32))
      CALL SSS1_1(W(1,4),W(1,5),GC_30,MDL_MH,MDL_WH,W(1,38))
      CALL SSS1_1(W(1,3),W(1,38),GC_30,MDL_MH,MDL_WH,W(1,39))
C     Counter-term amplitude(s) for loop diagram number 47
      CALL VVS1_0(W(1,1),W(1,2),W(1,39),R2_GGHB,AMPL(1,33))
      CALL SSS1_1(W(1,3),W(1,38),GC_HHETA0,MDL_META,MDL_WETA,W(1,40))
C     Counter-term amplitude(s) for loop diagram number 48
      CALL VVS1_0(W(1,1),W(1,2),W(1,40),R2_GGETA0B,AMPL(1,34))
      CALL SSS1_1(W(1,3),W(1,38),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1,41)
     $ )
C     Counter-term amplitude(s) for loop diagram number 49
      CALL VVS1_0(W(1,1),W(1,2),W(1,41),R2_GGIOTA0B,AMPL(1,35))
C     Counter-term amplitude(s) for loop diagram number 50
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,38),W(1,3),R2_GGHHB,AMPL(1,36))
      CALL SSSS1_1(W(1,3),W(1,4),W(1,5),GC_HHHH,MDL_MH,MDL_WH,W(1,42))
C     Counter-term amplitude(s) for loop diagram number 52
      CALL VVS1_0(W(1,1),W(1,2),W(1,42),R2_GGHB,AMPL(1,37))
      CALL SSSS1_4(W(1,3),W(1,4),W(1,5),GC_HHHETA0,MDL_META,MDL_WETA
     $ ,W(1,43))
C     Counter-term amplitude(s) for loop diagram number 53
      CALL VVS1_0(W(1,1),W(1,2),W(1,43),R2_GGETA0B,AMPL(1,38))
      CALL SSSS1_4(W(1,3),W(1,4),W(1,5),GC_HHHIOTA0,MDL_MIOTA
     $ ,MDL_WIOTA,W(1,44))
C     Counter-term amplitude(s) for loop diagram number 54
      CALL VVS1_0(W(1,1),W(1,2),W(1,44),R2_GGIOTA0B,AMPL(1,39))
C     Counter-term amplitude(s) for loop diagram number 145
      CALL VVS1_0(W(1,1),W(1,2),W(1,7),R2_GGHT,AMPL(1,40))
C     Counter-term amplitude(s) for loop diagram number 146
      CALL VVS1_0(W(1,1),W(1,2),W(1,8),R2_GGETA0T,AMPL(1,41))
C     Counter-term amplitude(s) for loop diagram number 147
      CALL VVS1_0(W(1,1),W(1,2),W(1,9),R2_GGIOTA0T,AMPL(1,42))
C     Counter-term amplitude(s) for loop diagram number 148
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,6),R2_GGHETA0T,AMPL(1,43)
     $ )
C     Counter-term amplitude(s) for loop diagram number 150
      CALL VVS1_0(W(1,1),W(1,2),W(1,11),R2_GGHT,AMPL(1,44))
C     Counter-term amplitude(s) for loop diagram number 151
      CALL VVS1_0(W(1,1),W(1,2),W(1,12),R2_GGETA0T,AMPL(1,45))
C     Counter-term amplitude(s) for loop diagram number 152
      CALL VVS1_0(W(1,1),W(1,2),W(1,13),R2_GGIOTA0T,AMPL(1,46))
C     Counter-term amplitude(s) for loop diagram number 153
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,10),R2_GGHIOTA0T,AMPL(1
     $ ,47))
C     Counter-term amplitude(s) for loop diagram number 155
      CALL VVS1_0(W(1,1),W(1,2),W(1,15),R2_GGHT,AMPL(1,48))
C     Counter-term amplitude(s) for loop diagram number 156
      CALL VVS1_0(W(1,1),W(1,2),W(1,16),R2_GGETA0T,AMPL(1,49))
C     Counter-term amplitude(s) for loop diagram number 157
      CALL VVS1_0(W(1,1),W(1,2),W(1,17),R2_GGIOTA0T,AMPL(1,50))
C     Counter-term amplitude(s) for loop diagram number 158
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,14),R2_GGHHT,AMPL(1,51))
C     Counter-term amplitude(s) for loop diagram number 160
      CALL VVS1_0(W(1,1),W(1,2),W(1,19),R2_GGHT,AMPL(1,52))
C     Counter-term amplitude(s) for loop diagram number 161
      CALL VVS1_0(W(1,1),W(1,2),W(1,20),R2_GGETA0T,AMPL(1,53))
C     Counter-term amplitude(s) for loop diagram number 162
      CALL VVS1_0(W(1,1),W(1,2),W(1,21),R2_GGIOTA0T,AMPL(1,54))
C     Counter-term amplitude(s) for loop diagram number 163
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,18),R2_GGHETA0T,AMPL(1
     $ ,55))
C     Counter-term amplitude(s) for loop diagram number 165
      CALL VVS1_0(W(1,1),W(1,2),W(1,23),R2_GGHT,AMPL(1,56))
C     Counter-term amplitude(s) for loop diagram number 166
      CALL VVS1_0(W(1,1),W(1,2),W(1,24),R2_GGETA0T,AMPL(1,57))
C     Counter-term amplitude(s) for loop diagram number 167
      CALL VVS1_0(W(1,1),W(1,2),W(1,25),R2_GGIOTA0T,AMPL(1,58))
C     Counter-term amplitude(s) for loop diagram number 168
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,22),R2_GGHIOTA0T,AMPL(1
     $ ,59))
C     Counter-term amplitude(s) for loop diagram number 170
      CALL VVS1_0(W(1,1),W(1,2),W(1,27),R2_GGHT,AMPL(1,60))
C     Counter-term amplitude(s) for loop diagram number 171
      CALL VVS1_0(W(1,1),W(1,2),W(1,28),R2_GGETA0T,AMPL(1,61))
C     Counter-term amplitude(s) for loop diagram number 172
      CALL VVS1_0(W(1,1),W(1,2),W(1,29),R2_GGIOTA0T,AMPL(1,62))
C     Counter-term amplitude(s) for loop diagram number 173
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,26),R2_GGHHT,AMPL(1,63))
C     Counter-term amplitude(s) for loop diagram number 175
      CALL VVS1_0(W(1,1),W(1,2),W(1,31),R2_GGHT,AMPL(1,64))
C     Counter-term amplitude(s) for loop diagram number 176
      CALL VVS1_0(W(1,1),W(1,2),W(1,32),R2_GGETA0T,AMPL(1,65))
C     Counter-term amplitude(s) for loop diagram number 177
      CALL VVS1_0(W(1,1),W(1,2),W(1,33),R2_GGIOTA0T,AMPL(1,66))
C     Counter-term amplitude(s) for loop diagram number 178
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,3),W(1,30),R2_GGHETA0T,AMPL(1
     $ ,67))
C     Counter-term amplitude(s) for loop diagram number 180
      CALL VVS1_0(W(1,1),W(1,2),W(1,35),R2_GGHT,AMPL(1,68))
C     Counter-term amplitude(s) for loop diagram number 181
      CALL VVS1_0(W(1,1),W(1,2),W(1,36),R2_GGETA0T,AMPL(1,69))
C     Counter-term amplitude(s) for loop diagram number 182
      CALL VVS1_0(W(1,1),W(1,2),W(1,37),R2_GGIOTA0T,AMPL(1,70))
C     Counter-term amplitude(s) for loop diagram number 183
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,3),W(1,34),R2_GGHIOTA0T,AMPL(1
     $ ,71))
C     Counter-term amplitude(s) for loop diagram number 185
      CALL VVS1_0(W(1,1),W(1,2),W(1,39),R2_GGHT,AMPL(1,72))
C     Counter-term amplitude(s) for loop diagram number 186
      CALL VVS1_0(W(1,1),W(1,2),W(1,40),R2_GGETA0T,AMPL(1,73))
C     Counter-term amplitude(s) for loop diagram number 187
      CALL VVS1_0(W(1,1),W(1,2),W(1,41),R2_GGIOTA0T,AMPL(1,74))
C     Counter-term amplitude(s) for loop diagram number 188
      CALL R2_GGHH_0(W(1,1),W(1,2),W(1,38),W(1,3),R2_GGHHT,AMPL(1,75))
C     Counter-term amplitude(s) for loop diagram number 190
      CALL VVS1_0(W(1,1),W(1,2),W(1,42),R2_GGHT,AMPL(1,76))
C     Counter-term amplitude(s) for loop diagram number 191
      CALL VVS1_0(W(1,1),W(1,2),W(1,43),R2_GGETA0T,AMPL(1,77))
C     Counter-term amplitude(s) for loop diagram number 192
      CALL VVS1_0(W(1,1),W(1,2),W(1,44),R2_GGIOTA0T,AMPL(1,78))
C     At this point, all CT amps needed for (QCD=4), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

