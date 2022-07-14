      SUBROUTINE ML5_0_0_1_MP_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
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
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/ML5_0_0_1_FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/ML5_0_0_1_SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_0_0_1_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/ML5_0_0_1_MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/ML5_0_0_1_MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NLOOPAMPS)
      COMMON/ML5_0_0_1_MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL MP_VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL MP_VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL MP_SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL MP_SXXXXX(P(0,4),+1*IC(4),W(1,4))
      CALL MP_SXXXXX(P(0,5),+1*IC(5),W(1,5))
      CALL MP_SSS1_1(W(1,3),W(1,4),GC_HHETA0,MDL_META,MDL_WETA,W(1,6))
      CALL MP_SSS1_2(W(1,6),W(1,5),GC_HHETA0,MDL_MH,MDL_WH,W(1,7))
C     Counter-term amplitude(s) for loop diagram number 7
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,7),R2_GGHB,AMPL(1,1))
      CALL MP_SSS1_1(W(1,6),W(1,5),GC_HETA0ETA0,MDL_META,MDL_WETA,W(1
     $ ,8))
C     Counter-term amplitude(s) for loop diagram number 8
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,8),R2_GGETA0B,AMPL(1,2))
      CALL MP_SSS1_2(W(1,6),W(1,5),GC_HIOTA0ETA0,MDL_MIOTA,MDL_WIOTA
     $ ,W(1,9))
C     Counter-term amplitude(s) for loop diagram number 9
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,9),R2_GGIOTA0B,AMPL(1,3))
C     Counter-term amplitude(s) for loop diagram number 10
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,6),R2_GGHETA0B,AMPL(1
     $ ,4))
      CALL MP_SSS1_1(W(1,3),W(1,4),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1
     $ ,10))
      CALL MP_SSS1_2(W(1,10),W(1,5),GC_HHIOTA0,MDL_MH,MDL_WH,W(1,11))
C     Counter-term amplitude(s) for loop diagram number 12
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,11),R2_GGHB,AMPL(1,5))
      CALL MP_SSS1_1(W(1,10),W(1,5),GC_HIOTA0ETA0,MDL_META,MDL_WETA
     $ ,W(1,12))
C     Counter-term amplitude(s) for loop diagram number 13
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,12),R2_GGETA0B,AMPL(1,6))
      CALL MP_SSS1_1(W(1,10),W(1,5),GC_HIOTA0IOTA0,MDL_MIOTA,MDL_WIOTA
     $ ,W(1,13))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,13),R2_GGIOTA0B,AMPL(1,7))
C     Counter-term amplitude(s) for loop diagram number 15
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,10),R2_GGHIOTA0B
     $ ,AMPL(1,8))
      CALL MP_SSS1_1(W(1,3),W(1,4),GC_30,MDL_MH,MDL_WH,W(1,14))
      CALL MP_SSS1_1(W(1,5),W(1,14),GC_30,MDL_MH,MDL_WH,W(1,15))
C     Counter-term amplitude(s) for loop diagram number 17
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,15),R2_GGHB,AMPL(1,9))
      CALL MP_SSS1_1(W(1,5),W(1,14),GC_HHETA0,MDL_META,MDL_WETA,W(1,16)
     $ )
C     Counter-term amplitude(s) for loop diagram number 18
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,16),R2_GGETA0B,AMPL(1,10))
      CALL MP_SSS1_1(W(1,5),W(1,14),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1
     $ ,17))
C     Counter-term amplitude(s) for loop diagram number 19
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,17),R2_GGIOTA0B,AMPL(1,11))
C     Counter-term amplitude(s) for loop diagram number 20
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,14),R2_GGHHB,AMPL(1
     $ ,12))
      CALL MP_SSS1_1(W(1,3),W(1,5),GC_HHETA0,MDL_META,MDL_WETA,W(1,18))
      CALL MP_SSS1_2(W(1,18),W(1,4),GC_HHETA0,MDL_MH,MDL_WH,W(1,19))
C     Counter-term amplitude(s) for loop diagram number 22
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,19),R2_GGHB,AMPL(1,13))
      CALL MP_SSS1_1(W(1,18),W(1,4),GC_HETA0ETA0,MDL_META,MDL_WETA,W(1
     $ ,20))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,20),R2_GGETA0B,AMPL(1,14))
      CALL MP_SSS1_2(W(1,18),W(1,4),GC_HIOTA0ETA0,MDL_MIOTA,MDL_WIOTA
     $ ,W(1,21))
C     Counter-term amplitude(s) for loop diagram number 24
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,21),R2_GGIOTA0B,AMPL(1,15))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,18),R2_GGHETA0B
     $ ,AMPL(1,16))
      CALL MP_SSS1_1(W(1,3),W(1,5),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1
     $ ,22))
      CALL MP_SSS1_2(W(1,22),W(1,4),GC_HHIOTA0,MDL_MH,MDL_WH,W(1,23))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,23),R2_GGHB,AMPL(1,17))
      CALL MP_SSS1_1(W(1,22),W(1,4),GC_HIOTA0ETA0,MDL_META,MDL_WETA
     $ ,W(1,24))
C     Counter-term amplitude(s) for loop diagram number 28
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,24),R2_GGETA0B,AMPL(1,18))
      CALL MP_SSS1_1(W(1,22),W(1,4),GC_HIOTA0IOTA0,MDL_MIOTA,MDL_WIOTA
     $ ,W(1,25))
C     Counter-term amplitude(s) for loop diagram number 29
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,25),R2_GGIOTA0B,AMPL(1,19))
C     Counter-term amplitude(s) for loop diagram number 30
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,22),R2_GGHIOTA0B
     $ ,AMPL(1,20))
      CALL MP_SSS1_1(W(1,3),W(1,5),GC_30,MDL_MH,MDL_WH,W(1,26))
      CALL MP_SSS1_1(W(1,4),W(1,26),GC_30,MDL_MH,MDL_WH,W(1,27))
C     Counter-term amplitude(s) for loop diagram number 32
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,27),R2_GGHB,AMPL(1,21))
      CALL MP_SSS1_1(W(1,4),W(1,26),GC_HHETA0,MDL_META,MDL_WETA,W(1,28)
     $ )
C     Counter-term amplitude(s) for loop diagram number 33
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,28),R2_GGETA0B,AMPL(1,22))
      CALL MP_SSS1_1(W(1,4),W(1,26),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1
     $ ,29))
C     Counter-term amplitude(s) for loop diagram number 34
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,29),R2_GGIOTA0B,AMPL(1,23))
C     Counter-term amplitude(s) for loop diagram number 35
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,26),R2_GGHHB,AMPL(1
     $ ,24))
      CALL MP_SSS1_1(W(1,4),W(1,5),GC_HHETA0,MDL_META,MDL_WETA,W(1,30))
      CALL MP_SSS1_2(W(1,30),W(1,3),GC_HHETA0,MDL_MH,MDL_WH,W(1,31))
C     Counter-term amplitude(s) for loop diagram number 37
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,31),R2_GGHB,AMPL(1,25))
      CALL MP_SSS1_1(W(1,30),W(1,3),GC_HETA0ETA0,MDL_META,MDL_WETA,W(1
     $ ,32))
C     Counter-term amplitude(s) for loop diagram number 38
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,32),R2_GGETA0B,AMPL(1,26))
      CALL MP_SSS1_2(W(1,30),W(1,3),GC_HIOTA0ETA0,MDL_MIOTA,MDL_WIOTA
     $ ,W(1,33))
C     Counter-term amplitude(s) for loop diagram number 39
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,33),R2_GGIOTA0B,AMPL(1,27))
C     Counter-term amplitude(s) for loop diagram number 40
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,3),W(1,30),R2_GGHETA0B
     $ ,AMPL(1,28))
      CALL MP_SSS1_1(W(1,4),W(1,5),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1
     $ ,34))
      CALL MP_SSS1_2(W(1,34),W(1,3),GC_HHIOTA0,MDL_MH,MDL_WH,W(1,35))
C     Counter-term amplitude(s) for loop diagram number 42
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,35),R2_GGHB,AMPL(1,29))
      CALL MP_SSS1_1(W(1,34),W(1,3),GC_HIOTA0ETA0,MDL_META,MDL_WETA
     $ ,W(1,36))
C     Counter-term amplitude(s) for loop diagram number 43
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,36),R2_GGETA0B,AMPL(1,30))
      CALL MP_SSS1_1(W(1,34),W(1,3),GC_HIOTA0IOTA0,MDL_MIOTA,MDL_WIOTA
     $ ,W(1,37))
C     Counter-term amplitude(s) for loop diagram number 44
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,37),R2_GGIOTA0B,AMPL(1,31))
C     Counter-term amplitude(s) for loop diagram number 45
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,3),W(1,34),R2_GGHIOTA0B
     $ ,AMPL(1,32))
      CALL MP_SSS1_1(W(1,4),W(1,5),GC_30,MDL_MH,MDL_WH,W(1,38))
      CALL MP_SSS1_1(W(1,3),W(1,38),GC_30,MDL_MH,MDL_WH,W(1,39))
C     Counter-term amplitude(s) for loop diagram number 47
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,39),R2_GGHB,AMPL(1,33))
      CALL MP_SSS1_1(W(1,3),W(1,38),GC_HHETA0,MDL_META,MDL_WETA,W(1,40)
     $ )
C     Counter-term amplitude(s) for loop diagram number 48
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,40),R2_GGETA0B,AMPL(1,34))
      CALL MP_SSS1_1(W(1,3),W(1,38),GC_HHIOTA0,MDL_MIOTA,MDL_WIOTA,W(1
     $ ,41))
C     Counter-term amplitude(s) for loop diagram number 49
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,41),R2_GGIOTA0B,AMPL(1,35))
C     Counter-term amplitude(s) for loop diagram number 50
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,38),W(1,3),R2_GGHHB,AMPL(1
     $ ,36))
      CALL MP_SSSS1_1(W(1,3),W(1,4),W(1,5),GC_HHHH,MDL_MH,MDL_WH,W(1
     $ ,42))
C     Counter-term amplitude(s) for loop diagram number 52
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,42),R2_GGHB,AMPL(1,37))
      CALL MP_SSSS1_4(W(1,3),W(1,4),W(1,5),GC_HHHETA0,MDL_META
     $ ,MDL_WETA,W(1,43))
C     Counter-term amplitude(s) for loop diagram number 53
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,43),R2_GGETA0B,AMPL(1,38))
      CALL MP_SSSS1_4(W(1,3),W(1,4),W(1,5),GC_HHHIOTA0,MDL_MIOTA
     $ ,MDL_WIOTA,W(1,44))
C     Counter-term amplitude(s) for loop diagram number 54
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,44),R2_GGIOTA0B,AMPL(1,39))
C     Counter-term amplitude(s) for loop diagram number 145
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,7),R2_GGHT,AMPL(1,40))
C     Counter-term amplitude(s) for loop diagram number 146
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,8),R2_GGETA0T,AMPL(1,41))
C     Counter-term amplitude(s) for loop diagram number 147
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,9),R2_GGIOTA0T,AMPL(1,42))
C     Counter-term amplitude(s) for loop diagram number 148
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,6),R2_GGHETA0T,AMPL(1
     $ ,43))
C     Counter-term amplitude(s) for loop diagram number 150
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,11),R2_GGHT,AMPL(1,44))
C     Counter-term amplitude(s) for loop diagram number 151
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,12),R2_GGETA0T,AMPL(1,45))
C     Counter-term amplitude(s) for loop diagram number 152
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,13),R2_GGIOTA0T,AMPL(1,46))
C     Counter-term amplitude(s) for loop diagram number 153
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,10),R2_GGHIOTA0T
     $ ,AMPL(1,47))
C     Counter-term amplitude(s) for loop diagram number 155
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,15),R2_GGHT,AMPL(1,48))
C     Counter-term amplitude(s) for loop diagram number 156
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,16),R2_GGETA0T,AMPL(1,49))
C     Counter-term amplitude(s) for loop diagram number 157
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,17),R2_GGIOTA0T,AMPL(1,50))
C     Counter-term amplitude(s) for loop diagram number 158
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,5),W(1,14),R2_GGHHT,AMPL(1
     $ ,51))
C     Counter-term amplitude(s) for loop diagram number 160
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,19),R2_GGHT,AMPL(1,52))
C     Counter-term amplitude(s) for loop diagram number 161
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,20),R2_GGETA0T,AMPL(1,53))
C     Counter-term amplitude(s) for loop diagram number 162
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,21),R2_GGIOTA0T,AMPL(1,54))
C     Counter-term amplitude(s) for loop diagram number 163
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,18),R2_GGHETA0T
     $ ,AMPL(1,55))
C     Counter-term amplitude(s) for loop diagram number 165
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,23),R2_GGHT,AMPL(1,56))
C     Counter-term amplitude(s) for loop diagram number 166
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,24),R2_GGETA0T,AMPL(1,57))
C     Counter-term amplitude(s) for loop diagram number 167
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,25),R2_GGIOTA0T,AMPL(1,58))
C     Counter-term amplitude(s) for loop diagram number 168
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,22),R2_GGHIOTA0T
     $ ,AMPL(1,59))
C     Counter-term amplitude(s) for loop diagram number 170
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,27),R2_GGHT,AMPL(1,60))
C     Counter-term amplitude(s) for loop diagram number 171
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,28),R2_GGETA0T,AMPL(1,61))
C     Counter-term amplitude(s) for loop diagram number 172
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,29),R2_GGIOTA0T,AMPL(1,62))
C     Counter-term amplitude(s) for loop diagram number 173
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,4),W(1,26),R2_GGHHT,AMPL(1
     $ ,63))
C     Counter-term amplitude(s) for loop diagram number 175
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,31),R2_GGHT,AMPL(1,64))
C     Counter-term amplitude(s) for loop diagram number 176
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,32),R2_GGETA0T,AMPL(1,65))
C     Counter-term amplitude(s) for loop diagram number 177
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,33),R2_GGIOTA0T,AMPL(1,66))
C     Counter-term amplitude(s) for loop diagram number 178
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,3),W(1,30),R2_GGHETA0T
     $ ,AMPL(1,67))
C     Counter-term amplitude(s) for loop diagram number 180
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,35),R2_GGHT,AMPL(1,68))
C     Counter-term amplitude(s) for loop diagram number 181
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,36),R2_GGETA0T,AMPL(1,69))
C     Counter-term amplitude(s) for loop diagram number 182
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,37),R2_GGIOTA0T,AMPL(1,70))
C     Counter-term amplitude(s) for loop diagram number 183
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,3),W(1,34),R2_GGHIOTA0T
     $ ,AMPL(1,71))
C     Counter-term amplitude(s) for loop diagram number 185
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,39),R2_GGHT,AMPL(1,72))
C     Counter-term amplitude(s) for loop diagram number 186
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,40),R2_GGETA0T,AMPL(1,73))
C     Counter-term amplitude(s) for loop diagram number 187
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,41),R2_GGIOTA0T,AMPL(1,74))
C     Counter-term amplitude(s) for loop diagram number 188
      CALL MP_R2_GGHH_0(W(1,1),W(1,2),W(1,38),W(1,3),R2_GGHHT,AMPL(1
     $ ,75))
C     Counter-term amplitude(s) for loop diagram number 190
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,42),R2_GGHT,AMPL(1,76))
C     Counter-term amplitude(s) for loop diagram number 191
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,43),R2_GGETA0T,AMPL(1,77))
C     Counter-term amplitude(s) for loop diagram number 192
      CALL MP_VVS1_0(W(1,1),W(1,2),W(1,44),R2_GGIOTA0T,AMPL(1,78))
C     At this point, all CT amps needed for (QCD=4), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      MP_CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

