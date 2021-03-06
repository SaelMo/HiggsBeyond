      SUBROUTINE ML5_0_0_1_COEF_CONSTRUCTION_1(P,NHEL,H,IC)
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
      IF (FILTER_SO.AND.LOOP_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Coefficient construction for loop diagram with ID 1
      CALL FFV1L1_2(PL(0,0),W(1,1),GC_5,MDL_MB,ZERO,PL(0,1),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,1))
      CALL FFV1L1_2(PL(0,1),W(1,2),GC_5,MDL_MB,ZERO,PL(0,2),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,2))
      CALL FFS1L1_2(PL(0,2),W(1,4),GC_33,MDL_MB,ZERO,PL(0,3),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,3))
      CALL FFS1L1_2(PL(0,3),W(1,5),GC_33,MDL_MB,ZERO,PL(0,4),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,3),4,COEFS,4,4,WL(1,0,1,4))
      CALL FFS1L1_2(PL(0,4),W(1,3),GC_33,MDL_MB,ZERO,PL(0,5),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,5))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,5),5,4,1,1,1,79)
C     Coefficient construction for loop diagram with ID 2
      CALL FFS1L1_2(PL(0,2),W(1,5),GC_33,MDL_MB,ZERO,PL(0,6),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,6))
      CALL FFS1L1_2(PL(0,6),W(1,4),GC_33,MDL_MB,ZERO,PL(0,7),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,6),4,COEFS,4,4,WL(1,0,1,7))
      CALL FFS1L1_2(PL(0,7),W(1,3),GC_33,MDL_MB,ZERO,PL(0,8),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,8))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,8),5,4,2,1,1,80)
C     Coefficient construction for loop diagram with ID 3
      CALL FFS1L1_2(PL(0,2),W(1,3),GC_33,MDL_MB,ZERO,PL(0,9),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,9))
      CALL FFS1L1_2(PL(0,9),W(1,5),GC_33,MDL_MB,ZERO,PL(0,10),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,10)
     $ )
      CALL FFS1L1_2(PL(0,10),W(1,4),GC_33,MDL_MB,ZERO,PL(0,11),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,10),4,COEFS,4,4,WL(1,0,1
     $ ,11))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,11),5,4,3,1,1,81)
C     Coefficient construction for loop diagram with ID 4
      CALL FFS1L1_2(PL(0,6),W(1,3),GC_33,MDL_MB,ZERO,PL(0,12),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,6),4,COEFS,4,4,WL(1,0,1,12)
     $ )
      CALL FFS1L1_2(PL(0,12),W(1,4),GC_33,MDL_MB,ZERO,PL(0,13),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,12),4,COEFS,4,4,WL(1,0,1
     $ ,13))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,13),5,4,4,1,1,82)
C     Coefficient construction for loop diagram with ID 5
      CALL FFS1L1_2(PL(0,9),W(1,4),GC_33,MDL_MB,ZERO,PL(0,14),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,14)
     $ )
      CALL FFS1L1_2(PL(0,14),W(1,5),GC_33,MDL_MB,ZERO,PL(0,15),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1
     $ ,15))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,15),5,4,5,1,1,83)
C     Coefficient construction for loop diagram with ID 6
      CALL FFS1L1_2(PL(0,3),W(1,3),GC_33,MDL_MB,ZERO,PL(0,16),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,3),4,COEFS,4,4,WL(1,0,1,16)
     $ )
      CALL FFS1L1_2(PL(0,16),W(1,5),GC_33,MDL_MB,ZERO,PL(0,17),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,16),4,COEFS,4,4,WL(1,0,1
     $ ,17))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,17),5,4,6,1,1,84)
C     Coefficient construction for loop diagram with ID 7
      CALL FFS1L1_2(PL(0,2),W(1,7),GC_33,MDL_MB,ZERO,PL(0,18),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,18)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,18),3,4,7,1,1,85)
C     Coefficient construction for loop diagram with ID 8
      CALL FFS1L1_2(PL(0,2),W(1,8),GC_933,MDL_MB,ZERO,PL(0,19),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,19)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,19),3,4,8,1,1,86)
C     Coefficient construction for loop diagram with ID 9
      CALL FFS1L1_2(PL(0,2),W(1,9),GC_733,MDL_MB,ZERO,PL(0,20),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,20)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,20),3,4,9,1,1,87)
C     Coefficient construction for loop diagram with ID 10
      CALL FFS1L1_2(PL(0,6),W(1,6),GC_933,MDL_MB,ZERO,PL(0,21),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,6),4,COEFS,4,4,WL(1,0,1,21)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,21),4,4,10,1,1,88)
C     Coefficient construction for loop diagram with ID 11
      CALL FFS1L1_2(PL(0,2),W(1,6),GC_933,MDL_MB,ZERO,PL(0,22),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,22)
     $ )
      CALL FFS1L1_2(PL(0,22),W(1,5),GC_33,MDL_MB,ZERO,PL(0,23),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,22),4,COEFS,4,4,WL(1,0,1
     $ ,23))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,23),4,4,11,1,1,89)
C     Coefficient construction for loop diagram with ID 12
      CALL FFS1L1_2(PL(0,2),W(1,11),GC_33,MDL_MB,ZERO,PL(0,24),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,24)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,24),3,4,12,1,1,90)
C     Coefficient construction for loop diagram with ID 13
      CALL FFS1L1_2(PL(0,2),W(1,12),GC_933,MDL_MB,ZERO,PL(0,25),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,25)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,25),3,4,13,1,1,91)
C     Coefficient construction for loop diagram with ID 14
      CALL FFS1L1_2(PL(0,2),W(1,13),GC_733,MDL_MB,ZERO,PL(0,26),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,26)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,26),3,4,14,1,1,92)
C     Coefficient construction for loop diagram with ID 15
      CALL FFS1L1_2(PL(0,6),W(1,10),GC_733,MDL_MB,ZERO,PL(0,27),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,6),4,COEFS,4,4,WL(1,0,1,27)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,27),4,4,15,1,1,93)
C     Coefficient construction for loop diagram with ID 16
      CALL FFS1L1_2(PL(0,2),W(1,10),GC_733,MDL_MB,ZERO,PL(0,28),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,28)
     $ )
      CALL FFS1L1_2(PL(0,28),W(1,5),GC_33,MDL_MB,ZERO,PL(0,29),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,28),4,COEFS,4,4,WL(1,0,1
     $ ,29))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,29),4,4,16,1,1,94)
C     Coefficient construction for loop diagram with ID 17
      CALL FFS1L1_2(PL(0,2),W(1,15),GC_33,MDL_MB,ZERO,PL(0,30),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,30)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,30),3,4,17,1,1,95)
C     Coefficient construction for loop diagram with ID 18
      CALL FFS1L1_2(PL(0,2),W(1,16),GC_933,MDL_MB,ZERO,PL(0,31),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,31)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,31),3,4,18,1,1,96)
C     Coefficient construction for loop diagram with ID 19
      CALL FFS1L1_2(PL(0,2),W(1,17),GC_733,MDL_MB,ZERO,PL(0,32),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,32)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,32),3,4,19,1,1,97)
C     Coefficient construction for loop diagram with ID 20
      CALL FFS1L1_2(PL(0,6),W(1,14),GC_33,MDL_MB,ZERO,PL(0,33),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,6),4,COEFS,4,4,WL(1,0,1,33)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,33),4,4,20,1,1,98)
C     Coefficient construction for loop diagram with ID 21
      CALL FFS1L1_2(PL(0,2),W(1,14),GC_33,MDL_MB,ZERO,PL(0,34),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,34)
     $ )
      CALL FFS1L1_2(PL(0,34),W(1,5),GC_33,MDL_MB,ZERO,PL(0,35),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1
     $ ,35))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,35),4,4,21,1,1,99)
C     Coefficient construction for loop diagram with ID 22
      CALL FFS1L1_2(PL(0,2),W(1,19),GC_33,MDL_MB,ZERO,PL(0,36),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,36)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,36),3,4,22,1,1,100)
C     Coefficient construction for loop diagram with ID 23
      CALL FFS1L1_2(PL(0,2),W(1,20),GC_933,MDL_MB,ZERO,PL(0,37),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,37)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,37),3,4,23,1,1,101)
C     Coefficient construction for loop diagram with ID 24
      CALL FFS1L1_2(PL(0,2),W(1,21),GC_733,MDL_MB,ZERO,PL(0,38),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,38)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,38),3,4,24,1,1,102)
C     Coefficient construction for loop diagram with ID 25
      CALL FFS1L1_2(PL(0,3),W(1,18),GC_933,MDL_MB,ZERO,PL(0,39),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,3),4,COEFS,4,4,WL(1,0,1,39)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,39),4,4,25,1,1,103)
C     Coefficient construction for loop diagram with ID 26
      CALL FFS1L1_2(PL(0,2),W(1,18),GC_933,MDL_MB,ZERO,PL(0,40),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,40)
     $ )
      CALL FFS1L1_2(PL(0,40),W(1,4),GC_33,MDL_MB,ZERO,PL(0,41),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,40),4,COEFS,4,4,WL(1,0,1
     $ ,41))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,41),4,4,26,1,1,104)
C     Coefficient construction for loop diagram with ID 27
      CALL FFS1L1_2(PL(0,2),W(1,23),GC_33,MDL_MB,ZERO,PL(0,42),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,42)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,42),3,4,27,1,1,105)
C     Coefficient construction for loop diagram with ID 28
      CALL FFS1L1_2(PL(0,2),W(1,24),GC_933,MDL_MB,ZERO,PL(0,43),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,43)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,43),3,4,28,1,1,106)
C     Coefficient construction for loop diagram with ID 29
      CALL FFS1L1_2(PL(0,2),W(1,25),GC_733,MDL_MB,ZERO,PL(0,44),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,44)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,44),3,4,29,1,1,107)
C     Coefficient construction for loop diagram with ID 30
      CALL FFS1L1_2(PL(0,3),W(1,22),GC_733,MDL_MB,ZERO,PL(0,45),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,3),4,COEFS,4,4,WL(1,0,1,45)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,45),4,4,30,1,1,108)
C     Coefficient construction for loop diagram with ID 31
      CALL FFS1L1_2(PL(0,2),W(1,22),GC_733,MDL_MB,ZERO,PL(0,46),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,46)
     $ )
      CALL FFS1L1_2(PL(0,46),W(1,4),GC_33,MDL_MB,ZERO,PL(0,47),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,46),4,COEFS,4,4,WL(1,0,1
     $ ,47))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,47),4,4,31,1,1,109)
C     Coefficient construction for loop diagram with ID 32
      CALL FFS1L1_2(PL(0,2),W(1,27),GC_33,MDL_MB,ZERO,PL(0,48),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,48)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,48),3,4,32,1,1,110)
C     Coefficient construction for loop diagram with ID 33
      CALL FFS1L1_2(PL(0,2),W(1,28),GC_933,MDL_MB,ZERO,PL(0,49),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,49)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,49),3,4,33,1,1,111)
C     Coefficient construction for loop diagram with ID 34
      CALL FFS1L1_2(PL(0,2),W(1,29),GC_733,MDL_MB,ZERO,PL(0,50),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,50)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,50),3,4,34,1,1,112)
C     Coefficient construction for loop diagram with ID 35
      CALL FFS1L1_2(PL(0,3),W(1,26),GC_33,MDL_MB,ZERO,PL(0,51),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,3),4,COEFS,4,4,WL(1,0,1,51)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,51),4,4,35,1,1,113)
C     Coefficient construction for loop diagram with ID 36
      CALL FFS1L1_2(PL(0,2),W(1,26),GC_33,MDL_MB,ZERO,PL(0,52),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,52)
     $ )
      CALL FFS1L1_2(PL(0,52),W(1,4),GC_33,MDL_MB,ZERO,PL(0,53),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,52),4,COEFS,4,4,WL(1,0,1
     $ ,53))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,53),4,4,36,1,1,114)
C     Coefficient construction for loop diagram with ID 37
      CALL FFS1L1_2(PL(0,2),W(1,31),GC_33,MDL_MB,ZERO,PL(0,54),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,54)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,54),3,4,37,1,1,115)
C     Coefficient construction for loop diagram with ID 38
      CALL FFS1L1_2(PL(0,2),W(1,32),GC_933,MDL_MB,ZERO,PL(0,55),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,55)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,55),3,4,38,1,1,116)
C     Coefficient construction for loop diagram with ID 39
      CALL FFS1L1_2(PL(0,2),W(1,33),GC_733,MDL_MB,ZERO,PL(0,56),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,56)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,56),3,4,39,1,1,117)
C     Coefficient construction for loop diagram with ID 40
      CALL FFS1L1_2(PL(0,2),W(1,30),GC_933,MDL_MB,ZERO,PL(0,57),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,57)
     $ )
      CALL FFS1L1_2(PL(0,57),W(1,3),GC_33,MDL_MB,ZERO,PL(0,58),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,57),4,COEFS,4,4,WL(1,0,1
     $ ,58))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,58),4,4,40,1,1,118)
C     Coefficient construction for loop diagram with ID 41
      CALL FFS1L1_2(PL(0,9),W(1,30),GC_933,MDL_MB,ZERO,PL(0,59),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,59)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,59),4,4,41,1,1,119)
C     Coefficient construction for loop diagram with ID 42
      CALL FFS1L1_2(PL(0,2),W(1,35),GC_33,MDL_MB,ZERO,PL(0,60),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,60)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,60),3,4,42,1,1,120)
C     Coefficient construction for loop diagram with ID 43
      CALL FFS1L1_2(PL(0,2),W(1,36),GC_933,MDL_MB,ZERO,PL(0,61),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,61)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,61),3,4,43,1,1,121)
C     Coefficient construction for loop diagram with ID 44
      CALL FFS1L1_2(PL(0,2),W(1,37),GC_733,MDL_MB,ZERO,PL(0,62),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,62)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,62),3,4,44,1,1,122)
C     Coefficient construction for loop diagram with ID 45
      CALL FFS1L1_2(PL(0,2),W(1,34),GC_733,MDL_MB,ZERO,PL(0,63),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,63)
     $ )
      CALL FFS1L1_2(PL(0,63),W(1,3),GC_33,MDL_MB,ZERO,PL(0,64),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,63),4,COEFS,4,4,WL(1,0,1
     $ ,64))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,64),4,4,45,1,1,123)
C     Coefficient construction for loop diagram with ID 46
      CALL FFS1L1_2(PL(0,9),W(1,34),GC_733,MDL_MB,ZERO,PL(0,65),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,65)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,65),4,4,46,1,1,124)
C     Coefficient construction for loop diagram with ID 47
      CALL FFS1L1_2(PL(0,2),W(1,39),GC_33,MDL_MB,ZERO,PL(0,66),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,66)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,66),3,4,47,1,1,125)
C     Coefficient construction for loop diagram with ID 48
      CALL FFS1L1_2(PL(0,2),W(1,40),GC_933,MDL_MB,ZERO,PL(0,67),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,67)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,67),3,4,48,1,1,126)
C     Coefficient construction for loop diagram with ID 49
      CALL FFS1L1_2(PL(0,2),W(1,41),GC_733,MDL_MB,ZERO,PL(0,68),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,68)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,68),3,4,49,1,1,127)
C     Coefficient construction for loop diagram with ID 50
      CALL FFS1L1_2(PL(0,2),W(1,38),GC_33,MDL_MB,ZERO,PL(0,69),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,69)
     $ )
      CALL FFS1L1_2(PL(0,69),W(1,3),GC_33,MDL_MB,ZERO,PL(0,70),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,69),4,COEFS,4,4,WL(1,0,1
     $ ,70))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,70),4,4,50,1,1,128)
C     Coefficient construction for loop diagram with ID 51
      CALL FFS1L1_2(PL(0,9),W(1,38),GC_33,MDL_MB,ZERO,PL(0,71),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,71)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,71),4,4,51,1,1,129)
C     Coefficient construction for loop diagram with ID 52
      CALL FFS1L1_2(PL(0,2),W(1,42),GC_33,MDL_MB,ZERO,PL(0,72),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,72)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,72),3,4,52,1,1,130)
C     Coefficient construction for loop diagram with ID 53
      CALL FFS1L1_2(PL(0,2),W(1,43),GC_933,MDL_MB,ZERO,PL(0,73),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,73)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,73),3,4,53,1,1,131)
C     Coefficient construction for loop diagram with ID 54
      CALL FFS1L1_2(PL(0,2),W(1,44),GC_733,MDL_MB,ZERO,PL(0,74),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,74)
     $ )
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,74),3,4,54,1,1,132)
C     Coefficient construction for loop diagram with ID 55
      CALL FFV1L2_1(PL(0,0),W(1,1),GC_5,MDL_MB,ZERO,PL(0,75),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,75)
     $ )
      CALL FFV1L2_1(PL(0,75),W(1,2),GC_5,MDL_MB,ZERO,PL(0,76),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,75),4,COEFS,4,4,WL(1,0,1
     $ ,76))
      CALL FFS1L2_1(PL(0,76),W(1,7),GC_33,MDL_MB,ZERO,PL(0,77),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,77))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,77),3,4,55,1,1,133)
C     Coefficient construction for loop diagram with ID 56
      CALL FFS1L2_1(PL(0,76),W(1,8),GC_933,MDL_MB,ZERO,PL(0,78),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,78))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,78),3,4,56,1,1,134)
C     Coefficient construction for loop diagram with ID 57
      CALL FFS1L2_1(PL(0,76),W(1,9),GC_733,MDL_MB,ZERO,PL(0,79),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,79))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,79),3,4,57,1,1,135)
C     Coefficient construction for loop diagram with ID 58
      CALL FFS1L2_1(PL(0,76),W(1,5),GC_33,MDL_MB,ZERO,PL(0,80),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,80))
      CALL FFS1L2_1(PL(0,80),W(1,6),GC_933,MDL_MB,ZERO,PL(0,81),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,80),4,COEFS,4,4,WL(1,0,1
     $ ,81))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,81),4,4,58,1,1,136)
C     Coefficient construction for loop diagram with ID 59
      CALL FFS1L2_1(PL(0,75),W(1,5),GC_33,MDL_MB,ZERO,PL(0,82),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,75),4,COEFS,4,4,WL(1,0,1
     $ ,82))
      CALL FFV1L2_1(PL(0,82),W(1,2),GC_5,MDL_MB,ZERO,PL(0,83),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,82),4,COEFS,4,4,WL(1,0,1
     $ ,83))
      CALL FFS1L2_1(PL(0,83),W(1,6),GC_933,MDL_MB,ZERO,PL(0,84),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,83),4,COEFS,4,4,WL(1,0,1
     $ ,84))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,84),4,4,59,1,1,137)
C     Coefficient construction for loop diagram with ID 60
      CALL FFS1L2_1(PL(0,76),W(1,11),GC_33,MDL_MB,ZERO,PL(0,85),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,85))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,85),3,4,60,1,1,138)
C     Coefficient construction for loop diagram with ID 61
      CALL FFS1L2_1(PL(0,76),W(1,13),GC_733,MDL_MB,ZERO,PL(0,86),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,86))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,86),3,4,61,1,1,139)
C     Coefficient construction for loop diagram with ID 62
      CALL FFS1L2_1(PL(0,76),W(1,12),GC_933,MDL_MB,ZERO,PL(0,87),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,87))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,87),3,4,62,1,1,140)
C     Coefficient construction for loop diagram with ID 63
      CALL FFS1L2_1(PL(0,80),W(1,10),GC_733,MDL_MB,ZERO,PL(0,88),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,80),4,COEFS,4,4,WL(1,0,1
     $ ,88))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,88),4,4,63,1,1,141)
C     Coefficient construction for loop diagram with ID 64
      CALL FFS1L2_1(PL(0,83),W(1,10),GC_733,MDL_MB,ZERO,PL(0,89),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,83),4,COEFS,4,4,WL(1,0,1
     $ ,89))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,89),4,4,64,1,1,142)
C     Coefficient construction for loop diagram with ID 65
      CALL FFS1L2_1(PL(0,76),W(1,16),GC_933,MDL_MB,ZERO,PL(0,90),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,90))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,90),3,4,65,1,1,143)
C     Coefficient construction for loop diagram with ID 66
      CALL FFS1L2_1(PL(0,76),W(1,17),GC_733,MDL_MB,ZERO,PL(0,91),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,91))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,91),3,4,66,1,1,144)
C     Coefficient construction for loop diagram with ID 67
      CALL FFS1L2_1(PL(0,76),W(1,15),GC_33,MDL_MB,ZERO,PL(0,92),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,92))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,92),3,4,67,1,1,145)
C     Coefficient construction for loop diagram with ID 68
      CALL FFS1L2_1(PL(0,80),W(1,14),GC_33,MDL_MB,ZERO,PL(0,93),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,80),4,COEFS,4,4,WL(1,0,1
     $ ,93))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,93),4,4,68,1,1,146)
C     Coefficient construction for loop diagram with ID 69
      CALL FFS1L2_1(PL(0,83),W(1,14),GC_33,MDL_MB,ZERO,PL(0,94),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,83),4,COEFS,4,4,WL(1,0,1
     $ ,94))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,94),4,4,69,1,1,147)
C     Coefficient construction for loop diagram with ID 70
      CALL FFS1L2_1(PL(0,76),W(1,6),GC_933,MDL_MB,ZERO,PL(0,95),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,95))
      CALL FFS1L2_1(PL(0,95),W(1,5),GC_33,MDL_MB,ZERO,PL(0,96),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,95),4,COEFS,4,4,WL(1,0,1
     $ ,96))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,96),4,4,70,1,1,148)
C     Coefficient construction for loop diagram with ID 71
      CALL FFS1L1_2(PL(0,1),W(1,5),GC_33,MDL_MB,ZERO,PL(0,97),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,97)
     $ )
      CALL FFV1L1_2(PL(0,97),W(1,2),GC_5,MDL_MB,ZERO,PL(0,98),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,97),4,COEFS,4,4,WL(1,0,1
     $ ,98))
      CALL FFS1L1_2(PL(0,98),W(1,6),GC_933,MDL_MB,ZERO,PL(0,99),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,98),4,COEFS,4,4,WL(1,0,1
     $ ,99))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,99),4,4,71,1,1,149)
C     Coefficient construction for loop diagram with ID 72
      CALL FFS1L2_1(PL(0,76),W(1,10),GC_733,MDL_MB,ZERO,PL(0,100)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,100))
      CALL FFS1L2_1(PL(0,100),W(1,5),GC_33,MDL_MB,ZERO,PL(0,101),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,100),4,COEFS,4,4,WL(1,0,1
     $ ,101))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,101),4,4,72,1,1,150)
C     Coefficient construction for loop diagram with ID 73
      CALL FFS1L1_2(PL(0,98),W(1,10),GC_733,MDL_MB,ZERO,PL(0,102)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,98),4,COEFS,4,4,WL(1,0,1
     $ ,102))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,102),4,4,73,1,1,151)
C     Coefficient construction for loop diagram with ID 74
      CALL FFS1L2_1(PL(0,76),W(1,14),GC_33,MDL_MB,ZERO,PL(0,103),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,103))
      CALL FFS1L2_1(PL(0,103),W(1,5),GC_33,MDL_MB,ZERO,PL(0,104),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,103),4,COEFS,4,4,WL(1,0,1
     $ ,104))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,104),4,4,74,1,1,152)
C     Coefficient construction for loop diagram with ID 75
      CALL FFS1L1_2(PL(0,98),W(1,14),GC_33,MDL_MB,ZERO,PL(0,105),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,98),4,COEFS,4,4,WL(1,0,1
     $ ,105))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,105),4,4,75,1,1,153)
C     Coefficient construction for loop diagram with ID 76
      CALL FFS1L2_1(PL(0,76),W(1,19),GC_33,MDL_MB,ZERO,PL(0,106),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,106))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,106),3,4,76,1,1,154)
C     Coefficient construction for loop diagram with ID 77
      CALL FFS1L2_1(PL(0,76),W(1,20),GC_933,MDL_MB,ZERO,PL(0,107)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,107))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,107),3,4,77,1,1,155)
C     Coefficient construction for loop diagram with ID 78
      CALL FFS1L2_1(PL(0,76),W(1,21),GC_733,MDL_MB,ZERO,PL(0,108)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,108))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,108),3,4,78,1,1,156)
C     Coefficient construction for loop diagram with ID 79
      CALL FFS1L2_1(PL(0,76),W(1,4),GC_33,MDL_MB,ZERO,PL(0,109),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,109))
      CALL FFS1L2_1(PL(0,109),W(1,18),GC_933,MDL_MB,ZERO,PL(0,110)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,109),4,COEFS,4,4,WL(1,0,1
     $ ,110))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,110),4,4,79,1,1,157)
C     Coefficient construction for loop diagram with ID 80
      CALL FFS1L2_1(PL(0,75),W(1,4),GC_33,MDL_MB,ZERO,PL(0,111),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,75),4,COEFS,4,4,WL(1,0,1
     $ ,111))
      CALL FFV1L2_1(PL(0,111),W(1,2),GC_5,MDL_MB,ZERO,PL(0,112),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,111),4,COEFS,4,4,WL(1,0,1
     $ ,112))
      CALL FFS1L2_1(PL(0,112),W(1,18),GC_933,MDL_MB,ZERO,PL(0,113)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,112),4,COEFS,4,4,WL(1,0,1
     $ ,113))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,113),4,4,80,1,1,158)
C     Coefficient construction for loop diagram with ID 81
      CALL FFS1L2_1(PL(0,76),W(1,23),GC_33,MDL_MB,ZERO,PL(0,114),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,114))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,114),3,4,81,1,1,159)
C     Coefficient construction for loop diagram with ID 82
      CALL FFS1L2_1(PL(0,76),W(1,25),GC_733,MDL_MB,ZERO,PL(0,115)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,115))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,115),3,4,82,1,1,160)
C     Coefficient construction for loop diagram with ID 83
      CALL FFS1L2_1(PL(0,76),W(1,24),GC_933,MDL_MB,ZERO,PL(0,116)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,116))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,116),3,4,83,1,1,161)
C     Coefficient construction for loop diagram with ID 84
      CALL FFS1L2_1(PL(0,109),W(1,22),GC_733,MDL_MB,ZERO,PL(0,117)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,109),4,COEFS,4,4,WL(1,0,1
     $ ,117))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,117),4,4,84,1,1,162)
C     Coefficient construction for loop diagram with ID 85
      CALL FFS1L2_1(PL(0,112),W(1,22),GC_733,MDL_MB,ZERO,PL(0,118)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,112),4,COEFS,4,4,WL(1,0,1
     $ ,118))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,118),4,4,85,1,1,163)
C     Coefficient construction for loop diagram with ID 86
      CALL FFS1L2_1(PL(0,76),W(1,28),GC_933,MDL_MB,ZERO,PL(0,119)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,119))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,119),3,4,86,1,1,164)
C     Coefficient construction for loop diagram with ID 87
      CALL FFS1L2_1(PL(0,76),W(1,29),GC_733,MDL_MB,ZERO,PL(0,120)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,120))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,120),3,4,87,1,1,165)
C     Coefficient construction for loop diagram with ID 88
      CALL FFS1L2_1(PL(0,76),W(1,27),GC_33,MDL_MB,ZERO,PL(0,121),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,121))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,121),3,4,88,1,1,166)
C     Coefficient construction for loop diagram with ID 89
      CALL FFS1L2_1(PL(0,109),W(1,26),GC_33,MDL_MB,ZERO,PL(0,122)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,109),4,COEFS,4,4,WL(1,0,1
     $ ,122))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,122),4,4,89,1,1,167)
C     Coefficient construction for loop diagram with ID 90
      CALL FFS1L2_1(PL(0,112),W(1,26),GC_33,MDL_MB,ZERO,PL(0,123)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,112),4,COEFS,4,4,WL(1,0,1
     $ ,123))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,123),4,4,90,1,1,168)
C     Coefficient construction for loop diagram with ID 91
      CALL FFS1L2_1(PL(0,76),W(1,18),GC_933,MDL_MB,ZERO,PL(0,124)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,124))
      CALL FFS1L2_1(PL(0,124),W(1,4),GC_33,MDL_MB,ZERO,PL(0,125),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,124),4,COEFS,4,4,WL(1,0,1
     $ ,125))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,125),4,4,91,1,1,169)
C     Coefficient construction for loop diagram with ID 92
      CALL FFS1L1_2(PL(0,1),W(1,4),GC_33,MDL_MB,ZERO,PL(0,126),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1
     $ ,126))
      CALL FFV1L1_2(PL(0,126),W(1,2),GC_5,MDL_MB,ZERO,PL(0,127),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,126),4,COEFS,4,4,WL(1,0,1
     $ ,127))
      CALL FFS1L1_2(PL(0,127),W(1,18),GC_933,MDL_MB,ZERO,PL(0,128)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,127),4,COEFS,4,4,WL(1,0,1
     $ ,128))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,128),4,4,92,1,1,170)
C     Coefficient construction for loop diagram with ID 93
      CALL FFS1L2_1(PL(0,76),W(1,22),GC_733,MDL_MB,ZERO,PL(0,129)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,129))
      CALL FFS1L2_1(PL(0,129),W(1,4),GC_33,MDL_MB,ZERO,PL(0,130),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,129),4,COEFS,4,4,WL(1,0,1
     $ ,130))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,130),4,4,93,1,1,171)
C     Coefficient construction for loop diagram with ID 94
      CALL FFS1L1_2(PL(0,127),W(1,22),GC_733,MDL_MB,ZERO,PL(0,131)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,127),4,COEFS,4,4,WL(1,0,1
     $ ,131))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,131),4,4,94,1,1,172)
C     Coefficient construction for loop diagram with ID 95
      CALL FFS1L2_1(PL(0,76),W(1,26),GC_33,MDL_MB,ZERO,PL(0,132),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,132))
      CALL FFS1L2_1(PL(0,132),W(1,4),GC_33,MDL_MB,ZERO,PL(0,133),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,132),4,COEFS,4,4,WL(1,0,1
     $ ,133))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,133),4,4,95,1,1,173)
C     Coefficient construction for loop diagram with ID 96
      CALL FFS1L1_2(PL(0,127),W(1,26),GC_33,MDL_MB,ZERO,PL(0,134)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,127),4,COEFS,4,4,WL(1,0,1
     $ ,134))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,134),4,4,96,1,1,174)
C     Coefficient construction for loop diagram with ID 97
      CALL FFS1L2_1(PL(0,80),W(1,4),GC_33,MDL_MB,ZERO,PL(0,135),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,80),4,COEFS,4,4,WL(1,0,1
     $ ,135))
      CALL FFS1L2_1(PL(0,135),W(1,3),GC_33,MDL_MB,ZERO,PL(0,136),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,135),4,COEFS,4,4,WL(1,0,1
     $ ,136))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,136),5,4,97,1,1,175)
C     Coefficient construction for loop diagram with ID 98
      CALL FFS1L2_1(PL(0,109),W(1,5),GC_33,MDL_MB,ZERO,PL(0,137),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,109),4,COEFS,4,4,WL(1,0,1
     $ ,137))
      CALL FFS1L2_1(PL(0,137),W(1,3),GC_33,MDL_MB,ZERO,PL(0,138),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,137),4,COEFS,4,4,WL(1,0,1
     $ ,138))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,138),5,4,98,1,1,176)
C     Coefficient construction for loop diagram with ID 99
      CALL FFS1L1_2(PL(0,1),W(1,3),GC_33,MDL_MB,ZERO,PL(0,139),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1
     $ ,139))
      CALL FFV1L1_2(PL(0,139),W(1,2),GC_5,MDL_MB,ZERO,PL(0,140),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,139),4,COEFS,4,4,WL(1,0,1
     $ ,140))
      CALL FFS1L1_2(PL(0,140),W(1,5),GC_33,MDL_MB,ZERO,PL(0,141),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,140),4,COEFS,4,4,WL(1,0,1
     $ ,141))
      CALL FFS1L1_2(PL(0,141),W(1,4),GC_33,MDL_MB,ZERO,PL(0,142),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,141),4,COEFS,4,4,WL(1,0,1
     $ ,142))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,142),5,4,99,1,1,177)
C     Coefficient construction for loop diagram with ID 100
      CALL FFS1L1_2(PL(0,139),W(1,5),GC_33,MDL_MB,ZERO,PL(0,143),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,139),4,COEFS,4,4,WL(1,0,1
     $ ,143))
      CALL FFV1L1_2(PL(0,143),W(1,2),GC_5,MDL_MB,ZERO,PL(0,144),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,143),4,COEFS,4,4,WL(1,0,1
     $ ,144))
      CALL FFS1L1_2(PL(0,144),W(1,4),GC_33,MDL_MB,ZERO,PL(0,145),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,144),4,COEFS,4,4,WL(1,0,1
     $ ,145))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,145),5,4,100,1,1,178)
C     Coefficient construction for loop diagram with ID 101
      CALL FFS1L1_2(PL(0,140),W(1,4),GC_33,MDL_MB,ZERO,PL(0,146),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,140),4,COEFS,4,4,WL(1,0,1
     $ ,146))
      CALL FFS1L1_2(PL(0,146),W(1,5),GC_33,MDL_MB,ZERO,PL(0,147),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,146),4,COEFS,4,4,WL(1,0,1
     $ ,147))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,147),5,4,101,1,1,179)
C     Coefficient construction for loop diagram with ID 102
      CALL FFS1L1_2(PL(0,139),W(1,4),GC_33,MDL_MB,ZERO,PL(0,148),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,139),4,COEFS,4,4,WL(1,0,1
     $ ,148))
      CALL FFV1L1_2(PL(0,148),W(1,2),GC_5,MDL_MB,ZERO,PL(0,149),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,148),4,COEFS,4,4,WL(1,0,1
     $ ,149))
      CALL FFS1L1_2(PL(0,149),W(1,5),GC_33,MDL_MB,ZERO,PL(0,150),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,149),4,COEFS,4,4,WL(1,0,1
     $ ,150))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,150),5,4,102,1,1,180)
C     Coefficient construction for loop diagram with ID 103
      CALL FFS1L2_1(PL(0,76),W(1,30),GC_933,MDL_MB,ZERO,PL(0,151)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,151))
      CALL FFS1L2_1(PL(0,151),W(1,3),GC_33,MDL_MB,ZERO,PL(0,152),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,151),4,COEFS,4,4,WL(1,0,1
     $ ,152))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,152),4,4,103,1,1,181)
C     Coefficient construction for loop diagram with ID 104
      CALL FFS1L1_2(PL(0,140),W(1,30),GC_933,MDL_MB,ZERO,PL(0,153)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,140),4,COEFS,4,4,WL(1,0,1
     $ ,153))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,153),4,4,104,1,1,182)
C     Coefficient construction for loop diagram with ID 105
      CALL FFS1L2_1(PL(0,76),W(1,34),GC_733,MDL_MB,ZERO,PL(0,154)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,154))
      CALL FFS1L2_1(PL(0,154),W(1,3),GC_33,MDL_MB,ZERO,PL(0,155),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,154),4,COEFS,4,4,WL(1,0,1
     $ ,155))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,155),4,4,105,1,1,183)
C     Coefficient construction for loop diagram with ID 106
      CALL FFS1L1_2(PL(0,140),W(1,34),GC_733,MDL_MB,ZERO,PL(0,156)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,140),4,COEFS,4,4,WL(1,0,1
     $ ,156))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,156),4,4,106,1,1,184)
C     Coefficient construction for loop diagram with ID 107
      CALL FFS1L2_1(PL(0,76),W(1,38),GC_33,MDL_MB,ZERO,PL(0,157),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,157))
      CALL FFS1L2_1(PL(0,157),W(1,3),GC_33,MDL_MB,ZERO,PL(0,158),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,157),4,COEFS,4,4,WL(1,0,1
     $ ,158))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,158),4,4,107,1,1,185)
C     Coefficient construction for loop diagram with ID 108
      CALL FFS1L1_2(PL(0,140),W(1,38),GC_33,MDL_MB,ZERO,PL(0,159)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,140),4,COEFS,4,4,WL(1,0,1
     $ ,159))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,159),4,4,108,1,1,186)
C     Coefficient construction for loop diagram with ID 109
      CALL FFS1L2_1(PL(0,76),W(1,31),GC_33,MDL_MB,ZERO,PL(0,160),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,160))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,160),3,4,109,1,1,187)
C     Coefficient construction for loop diagram with ID 110
      CALL FFS1L2_1(PL(0,76),W(1,32),GC_933,MDL_MB,ZERO,PL(0,161)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,161))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,161),3,4,110,1,1,188)
C     Coefficient construction for loop diagram with ID 111
      CALL FFS1L2_1(PL(0,76),W(1,33),GC_733,MDL_MB,ZERO,PL(0,162)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,162))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,162),3,4,111,1,1,189)
C     Coefficient construction for loop diagram with ID 112
      CALL FFS1L2_1(PL(0,76),W(1,3),GC_33,MDL_MB,ZERO,PL(0,163),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,163))
      CALL FFS1L2_1(PL(0,163),W(1,30),GC_933,MDL_MB,ZERO,PL(0,164)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,163),4,COEFS,4,4,WL(1,0,1
     $ ,164))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,164),4,4,112,1,1,190)
C     Coefficient construction for loop diagram with ID 113
      CALL FFS1L2_1(PL(0,75),W(1,3),GC_33,MDL_MB,ZERO,PL(0,165),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,75),4,COEFS,4,4,WL(1,0,1
     $ ,165))
      CALL FFV1L2_1(PL(0,165),W(1,2),GC_5,MDL_MB,ZERO,PL(0,166),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,165),4,COEFS,4,4,WL(1,0,1
     $ ,166))
      CALL FFS1L2_1(PL(0,166),W(1,30),GC_933,MDL_MB,ZERO,PL(0,167)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,166),4,COEFS,4,4,WL(1,0,1
     $ ,167))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,167),4,4,113,1,1,191)
C     Coefficient construction for loop diagram with ID 114
      CALL FFS1L2_1(PL(0,76),W(1,35),GC_33,MDL_MB,ZERO,PL(0,168),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,168))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,168),3,4,114,1,1,192)
C     Coefficient construction for loop diagram with ID 115
      CALL FFS1L2_1(PL(0,76),W(1,37),GC_733,MDL_MB,ZERO,PL(0,169)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,169))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,169),3,4,115,1,1,193)
C     Coefficient construction for loop diagram with ID 116
      CALL FFS1L2_1(PL(0,76),W(1,36),GC_933,MDL_MB,ZERO,PL(0,170)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,170))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,170),3,4,116,1,1,194)
C     Coefficient construction for loop diagram with ID 117
      CALL FFS1L2_1(PL(0,163),W(1,34),GC_733,MDL_MB,ZERO,PL(0,171)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,163),4,COEFS,4,4,WL(1,0,1
     $ ,171))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,171),4,4,117,1,1,195)
C     Coefficient construction for loop diagram with ID 118
      CALL FFS1L2_1(PL(0,166),W(1,34),GC_733,MDL_MB,ZERO,PL(0,172)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,166),4,COEFS,4,4,WL(1,0,1
     $ ,172))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,172),4,4,118,1,1,196)
C     Coefficient construction for loop diagram with ID 119
      CALL FFS1L2_1(PL(0,76),W(1,40),GC_933,MDL_MB,ZERO,PL(0,173)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,173))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,173),3,4,119,1,1,197)
C     Coefficient construction for loop diagram with ID 120
      CALL FFS1L2_1(PL(0,76),W(1,41),GC_733,MDL_MB,ZERO,PL(0,174)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,174))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,174),3,4,120,1,1,198)
C     Coefficient construction for loop diagram with ID 121
      CALL FFS1L2_1(PL(0,76),W(1,39),GC_33,MDL_MB,ZERO,PL(0,175),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,175))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,175),3,4,121,1,1,199)
C     Coefficient construction for loop diagram with ID 122
      CALL FFS1L2_1(PL(0,163),W(1,38),GC_33,MDL_MB,ZERO,PL(0,176)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,163),4,COEFS,4,4,WL(1,0,1
     $ ,176))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,176),4,4,122,1,1,200)
C     Coefficient construction for loop diagram with ID 123
      CALL FFS1L2_1(PL(0,166),W(1,38),GC_33,MDL_MB,ZERO,PL(0,177)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,166),4,COEFS,4,4,WL(1,0,1
     $ ,177))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,177),4,4,123,1,1,201)
C     Coefficient construction for loop diagram with ID 124
      CALL FFS1L2_1(PL(0,80),W(1,3),GC_33,MDL_MB,ZERO,PL(0,178),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,80),4,COEFS,4,4,WL(1,0,1
     $ ,178))
      CALL FFS1L2_1(PL(0,178),W(1,4),GC_33,MDL_MB,ZERO,PL(0,179),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,178),4,COEFS,4,4,WL(1,0,1
     $ ,179))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,179),5,4,124,1,1,202)
C     Coefficient construction for loop diagram with ID 125
      CALL FFS1L2_1(PL(0,163),W(1,5),GC_33,MDL_MB,ZERO,PL(0,180),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,163),4,COEFS,4,4,WL(1,0,1
     $ ,180))
      CALL FFS1L2_1(PL(0,180),W(1,4),GC_33,MDL_MB,ZERO,PL(0,181),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,180),4,COEFS,4,4,WL(1,0,1
     $ ,181))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,181),5,4,125,1,1,203)
C     Coefficient construction for loop diagram with ID 126
      CALL FFS1L2_1(PL(0,165),W(1,5),GC_33,MDL_MB,ZERO,PL(0,182),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,165),4,COEFS,4,4,WL(1,0,1
     $ ,182))
      CALL FFV1L2_1(PL(0,182),W(1,2),GC_5,MDL_MB,ZERO,PL(0,183),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,182),4,COEFS,4,4,WL(1,0,1
     $ ,183))
      CALL FFS1L2_1(PL(0,183),W(1,4),GC_33,MDL_MB,ZERO,PL(0,184),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,183),4,COEFS,4,4,WL(1,0,1
     $ ,184))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,184),5,4,126,1,1,204)
C     Coefficient construction for loop diagram with ID 127
      CALL FFS1L2_1(PL(0,166),W(1,5),GC_33,MDL_MB,ZERO,PL(0,185),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,166),4,COEFS,4,4,WL(1,0,1
     $ ,185))
      CALL FFS1L2_1(PL(0,185),W(1,4),GC_33,MDL_MB,ZERO,PL(0,186),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,185),4,COEFS,4,4,WL(1,0,1
     $ ,186))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,186),5,4,127,1,1,205)
C     Coefficient construction for loop diagram with ID 128
      CALL FFS1L1_2(PL(0,127),W(1,3),GC_33,MDL_MB,ZERO,PL(0,187),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,127),4,COEFS,4,4,WL(1,0,1
     $ ,187))
      CALL FFS1L1_2(PL(0,187),W(1,5),GC_33,MDL_MB,ZERO,PL(0,188),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,187),4,COEFS,4,4,WL(1,0,1
     $ ,188))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,188),5,4,128,1,1,206)
C     Coefficient construction for loop diagram with ID 129
      CALL FFS1L1_2(PL(0,126),W(1,3),GC_33,MDL_MB,ZERO,PL(0,189),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,126),4,COEFS,4,4,WL(1,0,1
     $ ,189))
      CALL FFV1L1_2(PL(0,189),W(1,2),GC_5,MDL_MB,ZERO,PL(0,190),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,189),4,COEFS,4,4,WL(1,0,1
     $ ,190))
      CALL FFS1L1_2(PL(0,190),W(1,5),GC_33,MDL_MB,ZERO,PL(0,191),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,190),4,COEFS,4,4,WL(1,0,1
     $ ,191))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,191),5,4,129,1,1,207)
C     Coefficient construction for loop diagram with ID 130
      CALL FFS1L2_1(PL(0,109),W(1,3),GC_33,MDL_MB,ZERO,PL(0,192),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,109),4,COEFS,4,4,WL(1,0,1
     $ ,192))
      CALL FFS1L2_1(PL(0,192),W(1,5),GC_33,MDL_MB,ZERO,PL(0,193),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,192),4,COEFS,4,4,WL(1,0,1
     $ ,193))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,193),5,4,130,1,1,208)
C     Coefficient construction for loop diagram with ID 131
      CALL FFS1L2_1(PL(0,163),W(1,4),GC_33,MDL_MB,ZERO,PL(0,194),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,163),4,COEFS,4,4,WL(1,0,1
     $ ,194))
      CALL FFS1L2_1(PL(0,194),W(1,5),GC_33,MDL_MB,ZERO,PL(0,195),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,194),4,COEFS,4,4,WL(1,0,1
     $ ,195))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,195),5,4,131,1,1,209)
C     Coefficient construction for loop diagram with ID 132
      CALL FFS1L2_1(PL(0,165),W(1,4),GC_33,MDL_MB,ZERO,PL(0,196),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,165),4,COEFS,4,4,WL(1,0,1
     $ ,196))
      CALL FFV1L2_1(PL(0,196),W(1,2),GC_5,MDL_MB,ZERO,PL(0,197),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,196),4,COEFS,4,4,WL(1,0,1
     $ ,197))
      CALL FFS1L2_1(PL(0,197),W(1,5),GC_33,MDL_MB,ZERO,PL(0,198),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,197),4,COEFS,4,4,WL(1,0,1
     $ ,198))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,198),5,4,132,1,1,210)
C     Coefficient construction for loop diagram with ID 133
      CALL FFS1L2_1(PL(0,166),W(1,4),GC_33,MDL_MB,ZERO,PL(0,199),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,166),4,COEFS,4,4,WL(1,0,1
     $ ,199))
      CALL FFS1L2_1(PL(0,199),W(1,5),GC_33,MDL_MB,ZERO,PL(0,200),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,199),4,COEFS,4,4,WL(1,0,1
     $ ,200))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,200),5,4,133,1,1,211)
C     Coefficient construction for loop diagram with ID 134
      CALL FFS1L2_1(PL(0,111),W(1,3),GC_33,MDL_MB,ZERO,PL(0,201),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,111),4,COEFS,4,4,WL(1,0,1
     $ ,201))
      CALL FFV1L2_1(PL(0,201),W(1,2),GC_5,MDL_MB,ZERO,PL(0,202),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,201),4,COEFS,4,4,WL(1,0,1
     $ ,202))
      CALL FFS1L2_1(PL(0,202),W(1,5),GC_33,MDL_MB,ZERO,PL(0,203),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,202),4,COEFS,4,4,WL(1,0,1
     $ ,203))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,203),5,4,134,1,1,212)
C     Coefficient construction for loop diagram with ID 135
      CALL FFS1L2_1(PL(0,112),W(1,3),GC_33,MDL_MB,ZERO,PL(0,204),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,112),4,COEFS,4,4,WL(1,0,1
     $ ,204))
      CALL FFS1L2_1(PL(0,204),W(1,5),GC_33,MDL_MB,ZERO,PL(0,205),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,204),4,COEFS,4,4,WL(1,0,1
     $ ,205))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,205),5,4,135,1,1,213)
C     Coefficient construction for loop diagram with ID 136
      CALL FFS1L2_1(PL(0,76),W(1,42),GC_33,MDL_MB,ZERO,PL(0,206),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,206))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,206),3,4,136,1,1,214)
C     Coefficient construction for loop diagram with ID 137
      CALL FFS1L2_1(PL(0,76),W(1,43),GC_933,MDL_MB,ZERO,PL(0,207)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,207))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,207),3,4,137,1,1,215)
C     Coefficient construction for loop diagram with ID 138
      CALL FFS1L2_1(PL(0,76),W(1,44),GC_733,MDL_MB,ZERO,PL(0,208)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,76),4,COEFS,4,4,WL(1,0,1
     $ ,208))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,208),3,4,138,1,1,216)
C     Coefficient construction for loop diagram with ID 139
      CALL FFV1L1_2(PL(0,0),W(1,1),GC_5,MDL_MT,MDL_WT,PL(0,209),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1
     $ ,209))
      CALL FFV1L1_2(PL(0,209),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,210)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,209),4,COEFS,4,4,WL(1,0,1
     $ ,210))
      CALL FFS1L1_2(PL(0,210),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,211)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,211))
      CALL FFS1L1_2(PL(0,211),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,212)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,211),4,COEFS,4,4,WL(1,0,1
     $ ,212))
      CALL FFS1L1_2(PL(0,212),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,213)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,212),4,COEFS,4,4,WL(1,0,1
     $ ,213))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,213),5,4,139,1,1,217)
C     Coefficient construction for loop diagram with ID 140
      CALL FFS1L1_2(PL(0,210),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,214)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,214))
      CALL FFS1L1_2(PL(0,214),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,215)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,214),4,COEFS,4,4,WL(1,0,1
     $ ,215))
      CALL FFS1L1_2(PL(0,215),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,216)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,215),4,COEFS,4,4,WL(1,0,1
     $ ,216))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,216),5,4,140,1,1,218)
C     Coefficient construction for loop diagram with ID 141
      CALL FFS1L1_2(PL(0,210),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,217)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,217))
      CALL FFS1L1_2(PL(0,217),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,218)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,217),4,COEFS,4,4,WL(1,0,1
     $ ,218))
      CALL FFS1L1_2(PL(0,218),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,219)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,218),4,COEFS,4,4,WL(1,0,1
     $ ,219))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,219),5,4,141,1,1,219)
C     Coefficient construction for loop diagram with ID 142
      CALL FFS1L1_2(PL(0,214),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,220)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,214),4,COEFS,4,4,WL(1,0,1
     $ ,220))
      CALL FFS1L1_2(PL(0,220),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,221)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,220),4,COEFS,4,4,WL(1,0,1
     $ ,221))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,221),5,4,142,1,1,220)
C     Coefficient construction for loop diagram with ID 143
      CALL FFS1L1_2(PL(0,217),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,222)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,217),4,COEFS,4,4,WL(1,0,1
     $ ,222))
      CALL FFS1L1_2(PL(0,222),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,223)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,222),4,COEFS,4,4,WL(1,0,1
     $ ,223))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,223),5,4,143,1,1,221)
C     Coefficient construction for loop diagram with ID 144
      CALL FFS1L1_2(PL(0,211),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,224)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,211),4,COEFS,4,4,WL(1,0,1
     $ ,224))
      CALL FFS1L1_2(PL(0,224),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,225)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,224),4,COEFS,4,4,WL(1,0,1
     $ ,225))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,225),5,4,144,1,1,222)
C     Coefficient construction for loop diagram with ID 145
      CALL FFS1L1_2(PL(0,210),W(1,7),GC_37,MDL_MT,MDL_WT,PL(0,226)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,226))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,226),3,4,145,1,1,223)
C     Coefficient construction for loop diagram with ID 146
      CALL FFS1L1_2(PL(0,210),W(1,8),GC_937,MDL_MT,MDL_WT,PL(0,227)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,227))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,227),3,4,146,1,1,224)
C     Coefficient construction for loop diagram with ID 147
      CALL FFS1L1_2(PL(0,210),W(1,9),GC_737,MDL_MT,MDL_WT,PL(0,228)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,228))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,228),3,4,147,1,1,225)
C     Coefficient construction for loop diagram with ID 148
      CALL FFS1L1_2(PL(0,214),W(1,6),GC_937,MDL_MT,MDL_WT,PL(0,229)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,214),4,COEFS,4,4,WL(1,0,1
     $ ,229))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,229),4,4,148,1,1,226)
C     Coefficient construction for loop diagram with ID 149
      CALL FFS1L1_2(PL(0,210),W(1,6),GC_937,MDL_MT,MDL_WT,PL(0,230)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,230))
      CALL FFS1L1_2(PL(0,230),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,231)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,230),4,COEFS,4,4,WL(1,0,1
     $ ,231))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,231),4,4,149,1,1,227)
C     Coefficient construction for loop diagram with ID 150
      CALL FFS1L1_2(PL(0,210),W(1,11),GC_37,MDL_MT,MDL_WT,PL(0,232)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,232))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,232),3,4,150,1,1,228)
C     Coefficient construction for loop diagram with ID 151
      CALL FFS1L1_2(PL(0,210),W(1,12),GC_937,MDL_MT,MDL_WT,PL(0,233)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,233))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,233),3,4,151,1,1,229)
C     Coefficient construction for loop diagram with ID 152
      CALL FFS1L1_2(PL(0,210),W(1,13),GC_737,MDL_MT,MDL_WT,PL(0,234)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,234))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,234),3,4,152,1,1,230)
C     Coefficient construction for loop diagram with ID 153
      CALL FFS1L1_2(PL(0,214),W(1,10),GC_737,MDL_MT,MDL_WT,PL(0,235)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,214),4,COEFS,4,4,WL(1,0,1
     $ ,235))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,235),4,4,153,1,1,231)
C     Coefficient construction for loop diagram with ID 154
      CALL FFS1L1_2(PL(0,210),W(1,10),GC_737,MDL_MT,MDL_WT,PL(0,236)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,236))
      CALL FFS1L1_2(PL(0,236),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,237)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,236),4,COEFS,4,4,WL(1,0,1
     $ ,237))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,237),4,4,154,1,1,232)
C     Coefficient construction for loop diagram with ID 155
      CALL FFS1L1_2(PL(0,210),W(1,15),GC_37,MDL_MT,MDL_WT,PL(0,238)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,238))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,238),3,4,155,1,1,233)
C     Coefficient construction for loop diagram with ID 156
      CALL FFS1L1_2(PL(0,210),W(1,16),GC_937,MDL_MT,MDL_WT,PL(0,239)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,239))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,239),3,4,156,1,1,234)
C     Coefficient construction for loop diagram with ID 157
      CALL FFS1L1_2(PL(0,210),W(1,17),GC_737,MDL_MT,MDL_WT,PL(0,240)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,240))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,240),3,4,157,1,1,235)
C     Coefficient construction for loop diagram with ID 158
      CALL FFS1L1_2(PL(0,214),W(1,14),GC_37,MDL_MT,MDL_WT,PL(0,241)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,214),4,COEFS,4,4,WL(1,0,1
     $ ,241))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,241),4,4,158,1,1,236)
C     Coefficient construction for loop diagram with ID 159
      CALL FFS1L1_2(PL(0,210),W(1,14),GC_37,MDL_MT,MDL_WT,PL(0,242)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,242))
      CALL FFS1L1_2(PL(0,242),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,243)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,242),4,COEFS,4,4,WL(1,0,1
     $ ,243))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,243),4,4,159,1,1,237)
C     Coefficient construction for loop diagram with ID 160
      CALL FFS1L1_2(PL(0,210),W(1,19),GC_37,MDL_MT,MDL_WT,PL(0,244)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,244))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,244),3,4,160,1,1,238)
C     Coefficient construction for loop diagram with ID 161
      CALL FFS1L1_2(PL(0,210),W(1,20),GC_937,MDL_MT,MDL_WT,PL(0,245)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,245))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,245),3,4,161,1,1,239)
C     Coefficient construction for loop diagram with ID 162
      CALL FFS1L1_2(PL(0,210),W(1,21),GC_737,MDL_MT,MDL_WT,PL(0,246)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,246))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,246),3,4,162,1,1,240)
C     Coefficient construction for loop diagram with ID 163
      CALL FFS1L1_2(PL(0,211),W(1,18),GC_937,MDL_MT,MDL_WT,PL(0,247)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,211),4,COEFS,4,4,WL(1,0,1
     $ ,247))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,247),4,4,163,1,1,241)
C     Coefficient construction for loop diagram with ID 164
      CALL FFS1L1_2(PL(0,210),W(1,18),GC_937,MDL_MT,MDL_WT,PL(0,248)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,248))
      CALL FFS1L1_2(PL(0,248),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,249)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,248),4,COEFS,4,4,WL(1,0,1
     $ ,249))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,249),4,4,164,1,1,242)
C     Coefficient construction for loop diagram with ID 165
      CALL FFS1L1_2(PL(0,210),W(1,23),GC_37,MDL_MT,MDL_WT,PL(0,250)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,250))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,250),3,4,165,1,1,243)
C     Coefficient construction for loop diagram with ID 166
      CALL FFS1L1_2(PL(0,210),W(1,24),GC_937,MDL_MT,MDL_WT,PL(0,251)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,251))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,251),3,4,166,1,1,244)
C     Coefficient construction for loop diagram with ID 167
      CALL FFS1L1_2(PL(0,210),W(1,25),GC_737,MDL_MT,MDL_WT,PL(0,252)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,252))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,252),3,4,167,1,1,245)
C     Coefficient construction for loop diagram with ID 168
      CALL FFS1L1_2(PL(0,211),W(1,22),GC_737,MDL_MT,MDL_WT,PL(0,253)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,211),4,COEFS,4,4,WL(1,0,1
     $ ,253))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,253),4,4,168,1,1,246)
C     Coefficient construction for loop diagram with ID 169
      CALL FFS1L1_2(PL(0,210),W(1,22),GC_737,MDL_MT,MDL_WT,PL(0,254)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,254))
      CALL FFS1L1_2(PL(0,254),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,255)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,254),4,COEFS,4,4,WL(1,0,1
     $ ,255))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,255),4,4,169,1,1,247)
C     Coefficient construction for loop diagram with ID 170
      CALL FFS1L1_2(PL(0,210),W(1,27),GC_37,MDL_MT,MDL_WT,PL(0,256)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,256))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,256),3,4,170,1,1,248)
C     Coefficient construction for loop diagram with ID 171
      CALL FFS1L1_2(PL(0,210),W(1,28),GC_937,MDL_MT,MDL_WT,PL(0,257)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,257))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,257),3,4,171,1,1,249)
C     Coefficient construction for loop diagram with ID 172
      CALL FFS1L1_2(PL(0,210),W(1,29),GC_737,MDL_MT,MDL_WT,PL(0,258)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,258))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,258),3,4,172,1,1,250)
C     Coefficient construction for loop diagram with ID 173
      CALL FFS1L1_2(PL(0,211),W(1,26),GC_37,MDL_MT,MDL_WT,PL(0,259)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,211),4,COEFS,4,4,WL(1,0,1
     $ ,259))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,259),4,4,173,1,1,251)
C     Coefficient construction for loop diagram with ID 174
      CALL FFS1L1_2(PL(0,210),W(1,26),GC_37,MDL_MT,MDL_WT,PL(0,260)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,260))
      CALL FFS1L1_2(PL(0,260),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,261)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,260),4,COEFS,4,4,WL(1,0,1
     $ ,261))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,261),4,4,174,1,1,252)
C     Coefficient construction for loop diagram with ID 175
      CALL FFS1L1_2(PL(0,210),W(1,31),GC_37,MDL_MT,MDL_WT,PL(0,262)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,262))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,262),3,4,175,1,1,253)
C     Coefficient construction for loop diagram with ID 176
      CALL FFS1L1_2(PL(0,210),W(1,32),GC_937,MDL_MT,MDL_WT,PL(0,263)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,263))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,263),3,4,176,1,1,254)
C     Coefficient construction for loop diagram with ID 177
      CALL FFS1L1_2(PL(0,210),W(1,33),GC_737,MDL_MT,MDL_WT,PL(0,264)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,264))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,264),3,4,177,1,1,255)
C     Coefficient construction for loop diagram with ID 178
      CALL FFS1L1_2(PL(0,210),W(1,30),GC_937,MDL_MT,MDL_WT,PL(0,265)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,265))
      CALL FFS1L1_2(PL(0,265),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,266)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,265),4,COEFS,4,4,WL(1,0,1
     $ ,266))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,266),4,4,178,1,1,256)
C     Coefficient construction for loop diagram with ID 179
      CALL FFS1L1_2(PL(0,217),W(1,30),GC_937,MDL_MT,MDL_WT,PL(0,267)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,217),4,COEFS,4,4,WL(1,0,1
     $ ,267))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,267),4,4,179,1,1,257)
C     Coefficient construction for loop diagram with ID 180
      CALL FFS1L1_2(PL(0,210),W(1,35),GC_37,MDL_MT,MDL_WT,PL(0,268)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,268))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,268),3,4,180,1,1,258)
C     Coefficient construction for loop diagram with ID 181
      CALL FFS1L1_2(PL(0,210),W(1,36),GC_937,MDL_MT,MDL_WT,PL(0,269)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,269))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,269),3,4,181,1,1,259)
C     Coefficient construction for loop diagram with ID 182
      CALL FFS1L1_2(PL(0,210),W(1,37),GC_737,MDL_MT,MDL_WT,PL(0,270)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,270))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,270),3,4,182,1,1,260)
C     Coefficient construction for loop diagram with ID 183
      CALL FFS1L1_2(PL(0,210),W(1,34),GC_737,MDL_MT,MDL_WT,PL(0,271)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,271))
      CALL FFS1L1_2(PL(0,271),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,272)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,271),4,COEFS,4,4,WL(1,0,1
     $ ,272))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,272),4,4,183,1,1,261)
C     Coefficient construction for loop diagram with ID 184
      CALL FFS1L1_2(PL(0,217),W(1,34),GC_737,MDL_MT,MDL_WT,PL(0,273)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,217),4,COEFS,4,4,WL(1,0,1
     $ ,273))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,273),4,4,184,1,1,262)
C     Coefficient construction for loop diagram with ID 185
      CALL FFS1L1_2(PL(0,210),W(1,39),GC_37,MDL_MT,MDL_WT,PL(0,274)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,274))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,274),3,4,185,1,1,263)
C     Coefficient construction for loop diagram with ID 186
      CALL FFS1L1_2(PL(0,210),W(1,40),GC_937,MDL_MT,MDL_WT,PL(0,275)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,275))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,275),3,4,186,1,1,264)
C     Coefficient construction for loop diagram with ID 187
      CALL FFS1L1_2(PL(0,210),W(1,41),GC_737,MDL_MT,MDL_WT,PL(0,276)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,276))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,276),3,4,187,1,1,265)
C     Coefficient construction for loop diagram with ID 188
      CALL FFS1L1_2(PL(0,210),W(1,38),GC_37,MDL_MT,MDL_WT,PL(0,277)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,277))
      CALL FFS1L1_2(PL(0,277),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,278)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,277),4,COEFS,4,4,WL(1,0,1
     $ ,278))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,278),4,4,188,1,1,266)
C     Coefficient construction for loop diagram with ID 189
      CALL FFS1L1_2(PL(0,217),W(1,38),GC_37,MDL_MT,MDL_WT,PL(0,279)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,217),4,COEFS,4,4,WL(1,0,1
     $ ,279))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,279),4,4,189,1,1,267)
C     Coefficient construction for loop diagram with ID 190
      CALL FFS1L1_2(PL(0,210),W(1,42),GC_37,MDL_MT,MDL_WT,PL(0,280)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,280))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,280),3,4,190,1,1,268)
C     Coefficient construction for loop diagram with ID 191
      CALL FFS1L1_2(PL(0,210),W(1,43),GC_937,MDL_MT,MDL_WT,PL(0,281)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,281))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,281),3,4,191,1,1,269)
C     Coefficient construction for loop diagram with ID 192
      CALL FFS1L1_2(PL(0,210),W(1,44),GC_737,MDL_MT,MDL_WT,PL(0,282)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,210),4,COEFS,4,4,WL(1,0,1
     $ ,282))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,282),3,4,192,1,1,270)
C     Coefficient construction for loop diagram with ID 193
      CALL FFV1L2_1(PL(0,0),W(1,1),GC_5,MDL_MT,MDL_WT,PL(0,283),COEFS)
      CALL ML5_0_0_1_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1
     $ ,283))
      CALL FFV1L2_1(PL(0,283),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,284)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,283),4,COEFS,4,4,WL(1,0,1
     $ ,284))
      CALL FFS1L2_1(PL(0,284),W(1,7),GC_37,MDL_MT,MDL_WT,PL(0,285)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,285))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,285),3,4,193,1,1,271)
C     Coefficient construction for loop diagram with ID 194
      CALL FFS1L2_1(PL(0,284),W(1,8),GC_937,MDL_MT,MDL_WT,PL(0,286)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,286))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,286),3,4,194,1,1,272)
C     Coefficient construction for loop diagram with ID 195
      CALL FFS1L2_1(PL(0,284),W(1,9),GC_737,MDL_MT,MDL_WT,PL(0,287)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,287))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,287),3,4,195,1,1,273)
C     Coefficient construction for loop diagram with ID 196
      CALL FFS1L2_1(PL(0,284),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,288)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,288))
      CALL FFS1L2_1(PL(0,288),W(1,6),GC_937,MDL_MT,MDL_WT,PL(0,289)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,288),4,COEFS,4,4,WL(1,0,1
     $ ,289))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,289),4,4,196,1,1,274)
C     Coefficient construction for loop diagram with ID 197
      CALL FFS1L2_1(PL(0,283),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,290)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,283),4,COEFS,4,4,WL(1,0,1
     $ ,290))
      CALL FFV1L2_1(PL(0,290),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,291)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,290),4,COEFS,4,4,WL(1,0,1
     $ ,291))
      CALL FFS1L2_1(PL(0,291),W(1,6),GC_937,MDL_MT,MDL_WT,PL(0,292)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,291),4,COEFS,4,4,WL(1,0,1
     $ ,292))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,292),4,4,197,1,1,275)
C     Coefficient construction for loop diagram with ID 198
      CALL FFS1L2_1(PL(0,284),W(1,11),GC_37,MDL_MT,MDL_WT,PL(0,293)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,293))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,293),3,4,198,1,1,276)
C     Coefficient construction for loop diagram with ID 199
      CALL FFS1L2_1(PL(0,284),W(1,13),GC_737,MDL_MT,MDL_WT,PL(0,294)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,294))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,294),3,4,199,1,1,277)
C     Coefficient construction for loop diagram with ID 200
      CALL FFS1L2_1(PL(0,284),W(1,12),GC_937,MDL_MT,MDL_WT,PL(0,295)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,295))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,295),3,4,200,1,1,278)
C     Coefficient construction for loop diagram with ID 201
      CALL FFS1L2_1(PL(0,288),W(1,10),GC_737,MDL_MT,MDL_WT,PL(0,296)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,288),4,COEFS,4,4,WL(1,0,1
     $ ,296))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,296),4,4,201,1,1,279)
C     Coefficient construction for loop diagram with ID 202
      CALL FFS1L2_1(PL(0,291),W(1,10),GC_737,MDL_MT,MDL_WT,PL(0,297)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,291),4,COEFS,4,4,WL(1,0,1
     $ ,297))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,297),4,4,202,1,1,280)
C     Coefficient construction for loop diagram with ID 203
      CALL FFS1L2_1(PL(0,284),W(1,16),GC_937,MDL_MT,MDL_WT,PL(0,298)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,298))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,298),3,4,203,1,1,281)
C     Coefficient construction for loop diagram with ID 204
      CALL FFS1L2_1(PL(0,284),W(1,17),GC_737,MDL_MT,MDL_WT,PL(0,299)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,299))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,299),3,4,204,1,1,282)
C     Coefficient construction for loop diagram with ID 205
      CALL FFS1L2_1(PL(0,284),W(1,15),GC_37,MDL_MT,MDL_WT,PL(0,300)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,300))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,300),3,4,205,1,1,283)
C     Coefficient construction for loop diagram with ID 206
      CALL FFS1L2_1(PL(0,288),W(1,14),GC_37,MDL_MT,MDL_WT,PL(0,301)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,288),4,COEFS,4,4,WL(1,0,1
     $ ,301))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,301),4,4,206,1,1,284)
C     Coefficient construction for loop diagram with ID 207
      CALL FFS1L2_1(PL(0,291),W(1,14),GC_37,MDL_MT,MDL_WT,PL(0,302)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,291),4,COEFS,4,4,WL(1,0,1
     $ ,302))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,302),4,4,207,1,1,285)
C     Coefficient construction for loop diagram with ID 208
      CALL FFS1L2_1(PL(0,284),W(1,6),GC_937,MDL_MT,MDL_WT,PL(0,303)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,303))
      CALL FFS1L2_1(PL(0,303),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,304)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,303),4,COEFS,4,4,WL(1,0,1
     $ ,304))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,304),4,4,208,1,1,286)
C     Coefficient construction for loop diagram with ID 209
      CALL FFS1L1_2(PL(0,209),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,305)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,209),4,COEFS,4,4,WL(1,0,1
     $ ,305))
      CALL FFV1L1_2(PL(0,305),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,306)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,305),4,COEFS,4,4,WL(1,0,1
     $ ,306))
      CALL FFS1L1_2(PL(0,306),W(1,6),GC_937,MDL_MT,MDL_WT,PL(0,307)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,306),4,COEFS,4,4,WL(1,0,1
     $ ,307))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,307),4,4,209,1,1,287)
C     Coefficient construction for loop diagram with ID 210
      CALL FFS1L2_1(PL(0,284),W(1,10),GC_737,MDL_MT,MDL_WT,PL(0,308)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,308))
      CALL FFS1L2_1(PL(0,308),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,309)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,308),4,COEFS,4,4,WL(1,0,1
     $ ,309))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,309),4,4,210,1,1,288)
C     Coefficient construction for loop diagram with ID 211
      CALL FFS1L1_2(PL(0,306),W(1,10),GC_737,MDL_MT,MDL_WT,PL(0,310)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,306),4,COEFS,4,4,WL(1,0,1
     $ ,310))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,310),4,4,211,1,1,289)
C     Coefficient construction for loop diagram with ID 212
      CALL FFS1L2_1(PL(0,284),W(1,14),GC_37,MDL_MT,MDL_WT,PL(0,311)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,311))
      CALL FFS1L2_1(PL(0,311),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,312)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,311),4,COEFS,4,4,WL(1,0,1
     $ ,312))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,312),4,4,212,1,1,290)
C     Coefficient construction for loop diagram with ID 213
      CALL FFS1L1_2(PL(0,306),W(1,14),GC_37,MDL_MT,MDL_WT,PL(0,313)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,306),4,COEFS,4,4,WL(1,0,1
     $ ,313))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,313),4,4,213,1,1,291)
C     Coefficient construction for loop diagram with ID 214
      CALL FFS1L2_1(PL(0,284),W(1,19),GC_37,MDL_MT,MDL_WT,PL(0,314)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,314))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,314),3,4,214,1,1,292)
C     Coefficient construction for loop diagram with ID 215
      CALL FFS1L2_1(PL(0,284),W(1,20),GC_937,MDL_MT,MDL_WT,PL(0,315)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,315))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,315),3,4,215,1,1,293)
C     Coefficient construction for loop diagram with ID 216
      CALL FFS1L2_1(PL(0,284),W(1,21),GC_737,MDL_MT,MDL_WT,PL(0,316)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,316))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,316),3,4,216,1,1,294)
C     Coefficient construction for loop diagram with ID 217
      CALL FFS1L2_1(PL(0,284),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,317)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,317))
      CALL FFS1L2_1(PL(0,317),W(1,18),GC_937,MDL_MT,MDL_WT,PL(0,318)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,317),4,COEFS,4,4,WL(1,0,1
     $ ,318))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,318),4,4,217,1,1,295)
C     Coefficient construction for loop diagram with ID 218
      CALL FFS1L2_1(PL(0,283),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,319)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,283),4,COEFS,4,4,WL(1,0,1
     $ ,319))
      CALL FFV1L2_1(PL(0,319),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,320)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,319),4,COEFS,4,4,WL(1,0,1
     $ ,320))
      CALL FFS1L2_1(PL(0,320),W(1,18),GC_937,MDL_MT,MDL_WT,PL(0,321)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,320),4,COEFS,4,4,WL(1,0,1
     $ ,321))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,321),4,4,218,1,1,296)
C     Coefficient construction for loop diagram with ID 219
      CALL FFS1L2_1(PL(0,284),W(1,23),GC_37,MDL_MT,MDL_WT,PL(0,322)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,322))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,322),3,4,219,1,1,297)
C     Coefficient construction for loop diagram with ID 220
      CALL FFS1L2_1(PL(0,284),W(1,25),GC_737,MDL_MT,MDL_WT,PL(0,323)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,323))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,323),3,4,220,1,1,298)
C     Coefficient construction for loop diagram with ID 221
      CALL FFS1L2_1(PL(0,284),W(1,24),GC_937,MDL_MT,MDL_WT,PL(0,324)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,324))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,324),3,4,221,1,1,299)
C     Coefficient construction for loop diagram with ID 222
      CALL FFS1L2_1(PL(0,317),W(1,22),GC_737,MDL_MT,MDL_WT,PL(0,325)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,317),4,COEFS,4,4,WL(1,0,1
     $ ,325))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,325),4,4,222,1,1,300)
C     Coefficient construction for loop diagram with ID 223
      CALL FFS1L2_1(PL(0,320),W(1,22),GC_737,MDL_MT,MDL_WT,PL(0,326)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,320),4,COEFS,4,4,WL(1,0,1
     $ ,326))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,326),4,4,223,1,1,301)
C     Coefficient construction for loop diagram with ID 224
      CALL FFS1L2_1(PL(0,284),W(1,28),GC_937,MDL_MT,MDL_WT,PL(0,327)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,327))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,327),3,4,224,1,1,302)
C     Coefficient construction for loop diagram with ID 225
      CALL FFS1L2_1(PL(0,284),W(1,29),GC_737,MDL_MT,MDL_WT,PL(0,328)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,328))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,328),3,4,225,1,1,303)
C     Coefficient construction for loop diagram with ID 226
      CALL FFS1L2_1(PL(0,284),W(1,27),GC_37,MDL_MT,MDL_WT,PL(0,329)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,329))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,329),3,4,226,1,1,304)
C     Coefficient construction for loop diagram with ID 227
      CALL FFS1L2_1(PL(0,317),W(1,26),GC_37,MDL_MT,MDL_WT,PL(0,330)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,317),4,COEFS,4,4,WL(1,0,1
     $ ,330))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,330),4,4,227,1,1,305)
C     Coefficient construction for loop diagram with ID 228
      CALL FFS1L2_1(PL(0,320),W(1,26),GC_37,MDL_MT,MDL_WT,PL(0,331)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,320),4,COEFS,4,4,WL(1,0,1
     $ ,331))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,331),4,4,228,1,1,306)
C     Coefficient construction for loop diagram with ID 229
      CALL FFS1L2_1(PL(0,284),W(1,18),GC_937,MDL_MT,MDL_WT,PL(0,332)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,332))
      CALL FFS1L2_1(PL(0,332),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,333)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,332),4,COEFS,4,4,WL(1,0,1
     $ ,333))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,333),4,4,229,1,1,307)
C     Coefficient construction for loop diagram with ID 230
      CALL FFS1L1_2(PL(0,209),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,334)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,209),4,COEFS,4,4,WL(1,0,1
     $ ,334))
      CALL FFV1L1_2(PL(0,334),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,335)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,334),4,COEFS,4,4,WL(1,0,1
     $ ,335))
      CALL FFS1L1_2(PL(0,335),W(1,18),GC_937,MDL_MT,MDL_WT,PL(0,336)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,335),4,COEFS,4,4,WL(1,0,1
     $ ,336))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,336),4,4,230,1,1,308)
C     Coefficient construction for loop diagram with ID 231
      CALL FFS1L2_1(PL(0,284),W(1,22),GC_737,MDL_MT,MDL_WT,PL(0,337)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,337))
      CALL FFS1L2_1(PL(0,337),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,338)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,337),4,COEFS,4,4,WL(1,0,1
     $ ,338))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,338),4,4,231,1,1,309)
C     Coefficient construction for loop diagram with ID 232
      CALL FFS1L1_2(PL(0,335),W(1,22),GC_737,MDL_MT,MDL_WT,PL(0,339)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,335),4,COEFS,4,4,WL(1,0,1
     $ ,339))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,339),4,4,232,1,1,310)
C     Coefficient construction for loop diagram with ID 233
      CALL FFS1L2_1(PL(0,284),W(1,26),GC_37,MDL_MT,MDL_WT,PL(0,340)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,340))
      CALL FFS1L2_1(PL(0,340),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,341)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,340),4,COEFS,4,4,WL(1,0,1
     $ ,341))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,341),4,4,233,1,1,311)
C     Coefficient construction for loop diagram with ID 234
      CALL FFS1L1_2(PL(0,335),W(1,26),GC_37,MDL_MT,MDL_WT,PL(0,342)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,335),4,COEFS,4,4,WL(1,0,1
     $ ,342))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,342),4,4,234,1,1,312)
C     Coefficient construction for loop diagram with ID 235
      CALL FFS1L2_1(PL(0,288),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,343)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,288),4,COEFS,4,4,WL(1,0,1
     $ ,343))
      CALL FFS1L2_1(PL(0,343),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,344)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,343),4,COEFS,4,4,WL(1,0,1
     $ ,344))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,344),5,4,235,1,1,313)
C     Coefficient construction for loop diagram with ID 236
      CALL FFS1L2_1(PL(0,317),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,345)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,317),4,COEFS,4,4,WL(1,0,1
     $ ,345))
      CALL FFS1L2_1(PL(0,345),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,346)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,345),4,COEFS,4,4,WL(1,0,1
     $ ,346))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,346),5,4,236,1,1,314)
C     Coefficient construction for loop diagram with ID 237
      CALL FFS1L1_2(PL(0,209),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,347)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,209),4,COEFS,4,4,WL(1,0,1
     $ ,347))
      CALL FFV1L1_2(PL(0,347),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,348)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,347),4,COEFS,4,4,WL(1,0,1
     $ ,348))
      CALL FFS1L1_2(PL(0,348),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,349)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,348),4,COEFS,4,4,WL(1,0,1
     $ ,349))
      CALL FFS1L1_2(PL(0,349),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,350)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,349),4,COEFS,4,4,WL(1,0,1
     $ ,350))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,350),5,4,237,1,1,315)
C     Coefficient construction for loop diagram with ID 238
      CALL FFS1L1_2(PL(0,347),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,351)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,347),4,COEFS,4,4,WL(1,0,1
     $ ,351))
      CALL FFV1L1_2(PL(0,351),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,352)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,351),4,COEFS,4,4,WL(1,0,1
     $ ,352))
      CALL FFS1L1_2(PL(0,352),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,353)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,352),4,COEFS,4,4,WL(1,0,1
     $ ,353))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,353),5,4,238,1,1,316)
C     Coefficient construction for loop diagram with ID 239
      CALL FFS1L1_2(PL(0,348),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,354)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,348),4,COEFS,4,4,WL(1,0,1
     $ ,354))
      CALL FFS1L1_2(PL(0,354),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,355)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,354),4,COEFS,4,4,WL(1,0,1
     $ ,355))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,355),5,4,239,1,1,317)
C     Coefficient construction for loop diagram with ID 240
      CALL FFS1L1_2(PL(0,347),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,356)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,347),4,COEFS,4,4,WL(1,0,1
     $ ,356))
      CALL FFV1L1_2(PL(0,356),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,357)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,356),4,COEFS,4,4,WL(1,0,1
     $ ,357))
      CALL FFS1L1_2(PL(0,357),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,358)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,357),4,COEFS,4,4,WL(1,0,1
     $ ,358))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,358),5,4,240,1,1,318)
C     Coefficient construction for loop diagram with ID 241
      CALL FFS1L2_1(PL(0,284),W(1,30),GC_937,MDL_MT,MDL_WT,PL(0,359)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,359))
      CALL FFS1L2_1(PL(0,359),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,360)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,359),4,COEFS,4,4,WL(1,0,1
     $ ,360))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,360),4,4,241,1,1,319)
C     Coefficient construction for loop diagram with ID 242
      CALL FFS1L1_2(PL(0,348),W(1,30),GC_937,MDL_MT,MDL_WT,PL(0,361)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,348),4,COEFS,4,4,WL(1,0,1
     $ ,361))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,361),4,4,242,1,1,320)
C     Coefficient construction for loop diagram with ID 243
      CALL FFS1L2_1(PL(0,284),W(1,34),GC_737,MDL_MT,MDL_WT,PL(0,362)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,362))
      CALL FFS1L2_1(PL(0,362),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,363)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,362),4,COEFS,4,4,WL(1,0,1
     $ ,363))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,363),4,4,243,1,1,321)
C     Coefficient construction for loop diagram with ID 244
      CALL FFS1L1_2(PL(0,348),W(1,34),GC_737,MDL_MT,MDL_WT,PL(0,364)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,348),4,COEFS,4,4,WL(1,0,1
     $ ,364))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,364),4,4,244,1,1,322)
C     Coefficient construction for loop diagram with ID 245
      CALL FFS1L2_1(PL(0,284),W(1,38),GC_37,MDL_MT,MDL_WT,PL(0,365)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,365))
      CALL FFS1L2_1(PL(0,365),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,366)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,365),4,COEFS,4,4,WL(1,0,1
     $ ,366))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,366),4,4,245,1,1,323)
C     Coefficient construction for loop diagram with ID 246
      CALL FFS1L1_2(PL(0,348),W(1,38),GC_37,MDL_MT,MDL_WT,PL(0,367)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,348),4,COEFS,4,4,WL(1,0,1
     $ ,367))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,367),4,4,246,1,1,324)
C     Coefficient construction for loop diagram with ID 247
      CALL FFS1L2_1(PL(0,284),W(1,31),GC_37,MDL_MT,MDL_WT,PL(0,368)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,368))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,368),3,4,247,1,1,325)
C     Coefficient construction for loop diagram with ID 248
      CALL FFS1L2_1(PL(0,284),W(1,32),GC_937,MDL_MT,MDL_WT,PL(0,369)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,369))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,369),3,4,248,1,1,326)
C     Coefficient construction for loop diagram with ID 249
      CALL FFS1L2_1(PL(0,284),W(1,33),GC_737,MDL_MT,MDL_WT,PL(0,370)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,370))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,370),3,4,249,1,1,327)
C     Coefficient construction for loop diagram with ID 250
      CALL FFS1L2_1(PL(0,284),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,371)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,371))
      CALL FFS1L2_1(PL(0,371),W(1,30),GC_937,MDL_MT,MDL_WT,PL(0,372)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,371),4,COEFS,4,4,WL(1,0,1
     $ ,372))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,372),4,4,250,1,1,328)
C     Coefficient construction for loop diagram with ID 251
      CALL FFS1L2_1(PL(0,283),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,373)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_1_1(WL(1,0,1,283),4,COEFS,4,4,WL(1,0,1
     $ ,373))
      CALL FFV1L2_1(PL(0,373),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,374)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,373),4,COEFS,4,4,WL(1,0,1
     $ ,374))
      CALL FFS1L2_1(PL(0,374),W(1,30),GC_937,MDL_MT,MDL_WT,PL(0,375)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,374),4,COEFS,4,4,WL(1,0,1
     $ ,375))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,375),4,4,251,1,1,329)
C     Coefficient construction for loop diagram with ID 252
      CALL FFS1L2_1(PL(0,284),W(1,35),GC_37,MDL_MT,MDL_WT,PL(0,376)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,376))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,376),3,4,252,1,1,330)
C     Coefficient construction for loop diagram with ID 253
      CALL FFS1L2_1(PL(0,284),W(1,37),GC_737,MDL_MT,MDL_WT,PL(0,377)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,377))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,377),3,4,253,1,1,331)
C     Coefficient construction for loop diagram with ID 254
      CALL FFS1L2_1(PL(0,284),W(1,36),GC_937,MDL_MT,MDL_WT,PL(0,378)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,378))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,378),3,4,254,1,1,332)
C     Coefficient construction for loop diagram with ID 255
      CALL FFS1L2_1(PL(0,371),W(1,34),GC_737,MDL_MT,MDL_WT,PL(0,379)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,371),4,COEFS,4,4,WL(1,0,1
     $ ,379))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,379),4,4,255,1,1,333)
C     Coefficient construction for loop diagram with ID 256
      CALL FFS1L2_1(PL(0,374),W(1,34),GC_737,MDL_MT,MDL_WT,PL(0,380)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,374),4,COEFS,4,4,WL(1,0,1
     $ ,380))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,380),4,4,256,1,1,334)
C     Coefficient construction for loop diagram with ID 257
      CALL FFS1L2_1(PL(0,284),W(1,40),GC_937,MDL_MT,MDL_WT,PL(0,381)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,381))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,381),3,4,257,1,1,335)
C     Coefficient construction for loop diagram with ID 258
      CALL FFS1L2_1(PL(0,284),W(1,41),GC_737,MDL_MT,MDL_WT,PL(0,382)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,382))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,382),3,4,258,1,1,336)
C     Coefficient construction for loop diagram with ID 259
      CALL FFS1L2_1(PL(0,284),W(1,39),GC_37,MDL_MT,MDL_WT,PL(0,383)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,383))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,383),3,4,259,1,1,337)
C     Coefficient construction for loop diagram with ID 260
      CALL FFS1L2_1(PL(0,371),W(1,38),GC_37,MDL_MT,MDL_WT,PL(0,384)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,371),4,COEFS,4,4,WL(1,0,1
     $ ,384))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,384),4,4,260,1,1,338)
C     Coefficient construction for loop diagram with ID 261
      CALL FFS1L2_1(PL(0,374),W(1,38),GC_37,MDL_MT,MDL_WT,PL(0,385)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,374),4,COEFS,4,4,WL(1,0,1
     $ ,385))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,385),4,4,261,1,1,339)
C     Coefficient construction for loop diagram with ID 262
      CALL FFS1L2_1(PL(0,288),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,386)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,288),4,COEFS,4,4,WL(1,0,1
     $ ,386))
      CALL FFS1L2_1(PL(0,386),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,387)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,386),4,COEFS,4,4,WL(1,0,1
     $ ,387))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,387),5,4,262,1,1,340)
C     Coefficient construction for loop diagram with ID 263
      CALL FFS1L2_1(PL(0,371),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,388)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,371),4,COEFS,4,4,WL(1,0,1
     $ ,388))
      CALL FFS1L2_1(PL(0,388),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,389)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,388),4,COEFS,4,4,WL(1,0,1
     $ ,389))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,389),5,4,263,1,1,341)
C     Coefficient construction for loop diagram with ID 264
      CALL FFS1L2_1(PL(0,373),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,390)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,373),4,COEFS,4,4,WL(1,0,1
     $ ,390))
      CALL FFV1L2_1(PL(0,390),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,391)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,390),4,COEFS,4,4,WL(1,0,1
     $ ,391))
      CALL FFS1L2_1(PL(0,391),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,392)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,391),4,COEFS,4,4,WL(1,0,1
     $ ,392))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,392),5,4,264,1,1,342)
C     Coefficient construction for loop diagram with ID 265
      CALL FFS1L2_1(PL(0,374),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,393)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,374),4,COEFS,4,4,WL(1,0,1
     $ ,393))
      CALL FFS1L2_1(PL(0,393),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,394)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,393),4,COEFS,4,4,WL(1,0,1
     $ ,394))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,394),5,4,265,1,1,343)
C     Coefficient construction for loop diagram with ID 266
      CALL FFS1L1_2(PL(0,335),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,395)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,335),4,COEFS,4,4,WL(1,0,1
     $ ,395))
      CALL FFS1L1_2(PL(0,395),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,396)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,395),4,COEFS,4,4,WL(1,0,1
     $ ,396))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,396),5,4,266,1,1,344)
C     Coefficient construction for loop diagram with ID 267
      CALL FFS1L1_2(PL(0,334),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,397)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,334),4,COEFS,4,4,WL(1,0,1
     $ ,397))
      CALL FFV1L1_2(PL(0,397),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,398)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,397),4,COEFS,4,4,WL(1,0,1
     $ ,398))
      CALL FFS1L1_2(PL(0,398),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,399)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,398),4,COEFS,4,4,WL(1,0,1
     $ ,399))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,399),5,4,267,1,1,345)
C     Coefficient construction for loop diagram with ID 268
      CALL FFS1L2_1(PL(0,317),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,400)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,317),4,COEFS,4,4,WL(1,0,1
     $ ,400))
      CALL FFS1L2_1(PL(0,400),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,401)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,400),4,COEFS,4,4,WL(1,0,1
     $ ,401))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,401),5,4,268,1,1,346)
C     Coefficient construction for loop diagram with ID 269
      CALL FFS1L2_1(PL(0,371),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,402)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,371),4,COEFS,4,4,WL(1,0,1
     $ ,402))
      CALL FFS1L2_1(PL(0,402),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,403)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,402),4,COEFS,4,4,WL(1,0,1
     $ ,403))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,403),5,4,269,1,1,347)
C     Coefficient construction for loop diagram with ID 270
      CALL FFS1L2_1(PL(0,373),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,404)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,373),4,COEFS,4,4,WL(1,0,1
     $ ,404))
      CALL FFV1L2_1(PL(0,404),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,405)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,404),4,COEFS,4,4,WL(1,0,1
     $ ,405))
      CALL FFS1L2_1(PL(0,405),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,406)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,405),4,COEFS,4,4,WL(1,0,1
     $ ,406))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,406),5,4,270,1,1,348)
C     Coefficient construction for loop diagram with ID 271
      CALL FFS1L2_1(PL(0,374),W(1,4),GC_37,MDL_MT,MDL_WT,PL(0,407)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,374),4,COEFS,4,4,WL(1,0,1
     $ ,407))
      CALL FFS1L2_1(PL(0,407),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,408)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,407),4,COEFS,4,4,WL(1,0,1
     $ ,408))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,408),5,4,271,1,1,349)
C     Coefficient construction for loop diagram with ID 272
      CALL FFS1L2_1(PL(0,319),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,409)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,319),4,COEFS,4,4,WL(1,0,1
     $ ,409))
      CALL FFV1L2_1(PL(0,409),W(1,2),GC_5,MDL_MT,MDL_WT,PL(0,410)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,409),4,COEFS,4,4,WL(1,0,1
     $ ,410))
      CALL FFS1L2_1(PL(0,410),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,411)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,410),4,COEFS,4,4,WL(1,0,1
     $ ,411))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,411),5,4,272,1,1,350)
C     Coefficient construction for loop diagram with ID 273
      CALL FFS1L2_1(PL(0,320),W(1,3),GC_37,MDL_MT,MDL_WT,PL(0,412)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_3_1(WL(1,0,1,320),4,COEFS,4,4,WL(1,0,1
     $ ,412))
      CALL FFS1L2_1(PL(0,412),W(1,5),GC_37,MDL_MT,MDL_WT,PL(0,413)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_4_1(WL(1,0,1,412),4,COEFS,4,4,WL(1,0,1
     $ ,413))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,413),5,4,273,1,1,351)
C     Coefficient construction for loop diagram with ID 274
      CALL FFS1L2_1(PL(0,284),W(1,42),GC_37,MDL_MT,MDL_WT,PL(0,414)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,414))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,414),3,4,274,1,1,352)
C     Coefficient construction for loop diagram with ID 275
      CALL FFS1L2_1(PL(0,284),W(1,43),GC_937,MDL_MT,MDL_WT,PL(0,415)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,415))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,415),3,4,275,1,1,353)
C     Coefficient construction for loop diagram with ID 276
      CALL FFS1L2_1(PL(0,284),W(1,44),GC_737,MDL_MT,MDL_WT,PL(0,416)
     $ ,COEFS)
      CALL ML5_0_0_1_UPDATE_WL_2_1(WL(1,0,1,284),4,COEFS,4,4,WL(1,0,1
     $ ,416))
      CALL ML5_0_0_1_CREATE_LOOP_COEFS(WL(1,0,1,416),3,4,276,1,1,354)
C     At this point, all loop coefficients needed for (QCD=4), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 4000

      GOTO 1001
 4000 CONTINUE
      LOOP_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

