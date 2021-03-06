      MODULE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INCLUDE 'coef_specs.inc'
      INCLUDE 'loop_max_coefs.inc'

C     Map associating a rank to each coefficient position
      INTEGER COEFTORANK_MAP(0:LOOPMAXCOEFS-1)
      DATA COEFTORANK_MAP(0:0)/1*0/
      DATA COEFTORANK_MAP(1:4)/4*1/
      DATA COEFTORANK_MAP(5:14)/10*2/
      DATA COEFTORANK_MAP(15:34)/20*3/
      DATA COEFTORANK_MAP(35:69)/35*4/
      DATA COEFTORANK_MAP(70:125)/56*5/

C     Map defining the number of coefficients for a symmetric tensor
C      of a given rank
      INTEGER NCOEF_R(0:5)
      DATA NCOEF_R/1,5,15,35,70,126/

C     Map defining the coef position resulting from the multiplication
C      of two lower rank coefs.
      INTEGER COMB_COEF_POS(0:LOOPMAXCOEFS-1,0:4)
      DATA COMB_COEF_POS(  0,  0:  4) /  0,  1,  2,  3,  4/
      DATA COMB_COEF_POS(  1,  0:  4) /  1,  5,  6,  8, 11/
      DATA COMB_COEF_POS(  2,  0:  4) /  2,  6,  7,  9, 12/
      DATA COMB_COEF_POS(  3,  0:  4) /  3,  8,  9, 10, 13/
      DATA COMB_COEF_POS(  4,  0:  4) /  4, 11, 12, 13, 14/
      DATA COMB_COEF_POS(  5,  0:  4) /  5, 15, 16, 19, 25/
      DATA COMB_COEF_POS(  6,  0:  4) /  6, 16, 17, 20, 26/
      DATA COMB_COEF_POS(  7,  0:  4) /  7, 17, 18, 21, 27/
      DATA COMB_COEF_POS(  8,  0:  4) /  8, 19, 20, 22, 28/
      DATA COMB_COEF_POS(  9,  0:  4) /  9, 20, 21, 23, 29/
      DATA COMB_COEF_POS( 10,  0:  4) / 10, 22, 23, 24, 30/
      DATA COMB_COEF_POS( 11,  0:  4) / 11, 25, 26, 28, 31/
      DATA COMB_COEF_POS( 12,  0:  4) / 12, 26, 27, 29, 32/
      DATA COMB_COEF_POS( 13,  0:  4) / 13, 28, 29, 30, 33/
      DATA COMB_COEF_POS( 14,  0:  4) / 14, 31, 32, 33, 34/
      DATA COMB_COEF_POS( 15,  0:  4) / 15, 35, 36, 40, 50/
      DATA COMB_COEF_POS( 16,  0:  4) / 16, 36, 37, 41, 51/
      DATA COMB_COEF_POS( 17,  0:  4) / 17, 37, 38, 42, 52/
      DATA COMB_COEF_POS( 18,  0:  4) / 18, 38, 39, 43, 53/
      DATA COMB_COEF_POS( 19,  0:  4) / 19, 40, 41, 44, 54/
      DATA COMB_COEF_POS( 20,  0:  4) / 20, 41, 42, 45, 55/
      DATA COMB_COEF_POS( 21,  0:  4) / 21, 42, 43, 46, 56/
      DATA COMB_COEF_POS( 22,  0:  4) / 22, 44, 45, 47, 57/
      DATA COMB_COEF_POS( 23,  0:  4) / 23, 45, 46, 48, 58/
      DATA COMB_COEF_POS( 24,  0:  4) / 24, 47, 48, 49, 59/
      DATA COMB_COEF_POS( 25,  0:  4) / 25, 50, 51, 54, 60/
      DATA COMB_COEF_POS( 26,  0:  4) / 26, 51, 52, 55, 61/
      DATA COMB_COEF_POS( 27,  0:  4) / 27, 52, 53, 56, 62/
      DATA COMB_COEF_POS( 28,  0:  4) / 28, 54, 55, 57, 63/
      DATA COMB_COEF_POS( 29,  0:  4) / 29, 55, 56, 58, 64/
      DATA COMB_COEF_POS( 30,  0:  4) / 30, 57, 58, 59, 65/
      DATA COMB_COEF_POS( 31,  0:  4) / 31, 60, 61, 63, 66/
      DATA COMB_COEF_POS( 32,  0:  4) / 32, 61, 62, 64, 67/
      DATA COMB_COEF_POS( 33,  0:  4) / 33, 63, 64, 65, 68/
      DATA COMB_COEF_POS( 34,  0:  4) / 34, 66, 67, 68, 69/
      DATA COMB_COEF_POS( 35,  0:  4) / 35, 70, 71, 76, 91/
      DATA COMB_COEF_POS( 36,  0:  4) / 36, 71, 72, 77, 92/
      DATA COMB_COEF_POS( 37,  0:  4) / 37, 72, 73, 78, 93/
      DATA COMB_COEF_POS( 38,  0:  4) / 38, 73, 74, 79, 94/
      DATA COMB_COEF_POS( 39,  0:  4) / 39, 74, 75, 80, 95/
      DATA COMB_COEF_POS( 40,  0:  4) / 40, 76, 77, 81, 96/
      DATA COMB_COEF_POS( 41,  0:  4) / 41, 77, 78, 82, 97/
      DATA COMB_COEF_POS( 42,  0:  4) / 42, 78, 79, 83, 98/
      DATA COMB_COEF_POS( 43,  0:  4) / 43, 79, 80, 84, 99/
      DATA COMB_COEF_POS( 44,  0:  4) / 44, 81, 82, 85,100/
      DATA COMB_COEF_POS( 45,  0:  4) / 45, 82, 83, 86,101/
      DATA COMB_COEF_POS( 46,  0:  4) / 46, 83, 84, 87,102/
      DATA COMB_COEF_POS( 47,  0:  4) / 47, 85, 86, 88,103/
      DATA COMB_COEF_POS( 48,  0:  4) / 48, 86, 87, 89,104/
      DATA COMB_COEF_POS( 49,  0:  4) / 49, 88, 89, 90,105/
      DATA COMB_COEF_POS( 50,  0:  4) / 50, 91, 92, 96,106/
      DATA COMB_COEF_POS( 51,  0:  4) / 51, 92, 93, 97,107/
      DATA COMB_COEF_POS( 52,  0:  4) / 52, 93, 94, 98,108/
      DATA COMB_COEF_POS( 53,  0:  4) / 53, 94, 95, 99,109/
      DATA COMB_COEF_POS( 54,  0:  4) / 54, 96, 97,100,110/
      DATA COMB_COEF_POS( 55,  0:  4) / 55, 97, 98,101,111/
      DATA COMB_COEF_POS( 56,  0:  4) / 56, 98, 99,102,112/
      DATA COMB_COEF_POS( 57,  0:  4) / 57,100,101,103,113/
      DATA COMB_COEF_POS( 58,  0:  4) / 58,101,102,104,114/
      DATA COMB_COEF_POS( 59,  0:  4) / 59,103,104,105,115/
      DATA COMB_COEF_POS( 60,  0:  4) / 60,106,107,110,116/
      DATA COMB_COEF_POS( 61,  0:  4) / 61,107,108,111,117/
      DATA COMB_COEF_POS( 62,  0:  4) / 62,108,109,112,118/
      DATA COMB_COEF_POS( 63,  0:  4) / 63,110,111,113,119/
      DATA COMB_COEF_POS( 64,  0:  4) / 64,111,112,114,120/
      DATA COMB_COEF_POS( 65,  0:  4) / 65,113,114,115,121/
      DATA COMB_COEF_POS( 66,  0:  4) / 66,116,117,119,122/
      DATA COMB_COEF_POS( 67,  0:  4) / 67,117,118,120,123/
      DATA COMB_COEF_POS( 68,  0:  4) / 68,119,120,121,124/
      DATA COMB_COEF_POS( 69,  0:  4) / 69,122,123,124,125/
      DATA COMB_COEF_POS( 70,  0:  4) / 70,126,127,133,154/
      DATA COMB_COEF_POS( 71,  0:  4) / 71,127,128,134,155/
      DATA COMB_COEF_POS( 72,  0:  4) / 72,128,129,135,156/
      DATA COMB_COEF_POS( 73,  0:  4) / 73,129,130,136,157/
      DATA COMB_COEF_POS( 74,  0:  4) / 74,130,131,137,158/
      DATA COMB_COEF_POS( 75,  0:  4) / 75,131,132,138,159/
      DATA COMB_COEF_POS( 76,  0:  4) / 76,133,134,139,160/
      DATA COMB_COEF_POS( 77,  0:  4) / 77,134,135,140,161/
      DATA COMB_COEF_POS( 78,  0:  4) / 78,135,136,141,162/
      DATA COMB_COEF_POS( 79,  0:  4) / 79,136,137,142,163/
      DATA COMB_COEF_POS( 80,  0:  4) / 80,137,138,143,164/
      DATA COMB_COEF_POS( 81,  0:  4) / 81,139,140,144,165/
      DATA COMB_COEF_POS( 82,  0:  4) / 82,140,141,145,166/
      DATA COMB_COEF_POS( 83,  0:  4) / 83,141,142,146,167/
      DATA COMB_COEF_POS( 84,  0:  4) / 84,142,143,147,168/
      DATA COMB_COEF_POS( 85,  0:  4) / 85,144,145,148,169/
      DATA COMB_COEF_POS( 86,  0:  4) / 86,145,146,149,170/
      DATA COMB_COEF_POS( 87,  0:  4) / 87,146,147,150,171/
      DATA COMB_COEF_POS( 88,  0:  4) / 88,148,149,151,172/
      DATA COMB_COEF_POS( 89,  0:  4) / 89,149,150,152,173/
      DATA COMB_COEF_POS( 90,  0:  4) / 90,151,152,153,174/
      DATA COMB_COEF_POS( 91,  0:  4) / 91,154,155,160,175/
      DATA COMB_COEF_POS( 92,  0:  4) / 92,155,156,161,176/
      DATA COMB_COEF_POS( 93,  0:  4) / 93,156,157,162,177/
      DATA COMB_COEF_POS( 94,  0:  4) / 94,157,158,163,178/
      DATA COMB_COEF_POS( 95,  0:  4) / 95,158,159,164,179/
      DATA COMB_COEF_POS( 96,  0:  4) / 96,160,161,165,180/
      DATA COMB_COEF_POS( 97,  0:  4) / 97,161,162,166,181/
      DATA COMB_COEF_POS( 98,  0:  4) / 98,162,163,167,182/
      DATA COMB_COEF_POS( 99,  0:  4) / 99,163,164,168,183/
      DATA COMB_COEF_POS(100,  0:  4) /100,165,166,169,184/
      DATA COMB_COEF_POS(101,  0:  4) /101,166,167,170,185/
      DATA COMB_COEF_POS(102,  0:  4) /102,167,168,171,186/
      DATA COMB_COEF_POS(103,  0:  4) /103,169,170,172,187/
      DATA COMB_COEF_POS(104,  0:  4) /104,170,171,173,188/
      DATA COMB_COEF_POS(105,  0:  4) /105,172,173,174,189/
      DATA COMB_COEF_POS(106,  0:  4) /106,175,176,180,190/
      DATA COMB_COEF_POS(107,  0:  4) /107,176,177,181,191/
      DATA COMB_COEF_POS(108,  0:  4) /108,177,178,182,192/
      DATA COMB_COEF_POS(109,  0:  4) /109,178,179,183,193/
      DATA COMB_COEF_POS(110,  0:  4) /110,180,181,184,194/
      DATA COMB_COEF_POS(111,  0:  4) /111,181,182,185,195/
      DATA COMB_COEF_POS(112,  0:  4) /112,182,183,186,196/
      DATA COMB_COEF_POS(113,  0:  4) /113,184,185,187,197/
      DATA COMB_COEF_POS(114,  0:  4) /114,185,186,188,198/
      DATA COMB_COEF_POS(115,  0:  4) /115,187,188,189,199/
      DATA COMB_COEF_POS(116,  0:  4) /116,190,191,194,200/
      DATA COMB_COEF_POS(117,  0:  4) /117,191,192,195,201/
      DATA COMB_COEF_POS(118,  0:  4) /118,192,193,196,202/
      DATA COMB_COEF_POS(119,  0:  4) /119,194,195,197,203/
      DATA COMB_COEF_POS(120,  0:  4) /120,195,196,198,204/
      DATA COMB_COEF_POS(121,  0:  4) /121,197,198,199,205/
      DATA COMB_COEF_POS(122,  0:  4) /122,200,201,203,206/
      DATA COMB_COEF_POS(123,  0:  4) /123,201,202,204,207/
      DATA COMB_COEF_POS(124,  0:  4) /124,203,204,205,208/
      DATA COMB_COEF_POS(125,  0:  4) /125,206,207,208,209/

      END MODULE ML5_0_0_1_POLYNOMIAL_CONSTANTS


C     The subroutine to create the loop coefficients form the last
C      loop wf.
C     In this case of loop-induced process, the reduction is performed
C      at the loop
C     amplitude level so that no multiplication is performed.

      SUBROUTINE ML5_0_0_1_CREATE_LOOP_COEFS(LOOP_WF,RANK,LCUT_SIZE
     $ ,LOOP_GROUP_NUMBER,SYMFACT,MULTIPLIER)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      REAL*8 ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(ZERO,ONE))
      COMPLEX*16 CMPLX_ZERO
      PARAMETER (CMPLX_ZERO=(ZERO,ZERO))
      INTEGER    NLOOPGROUPS
      PARAMETER (NLOOPGROUPS=276)
      INTEGER    NCOMB
      PARAMETER (NCOMB=4)
C     
C     ARGUMENTS 
C     
      COMPLEX*16 LOOP_WF(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER RANK, SYMFACT, LCUT_SIZE, LOOP_GROUP_NUMBER, MULTIPLIER
C     
C     LOCAL VARIABLES 
C     
      COMPLEX*16 CONST
      INTEGER I,J
C     
C     GLOBAL VARIABLES
C     
      COMPLEX*16 LOOPCOEFS(0:LOOPMAXCOEFS-1,NLOOPGROUPS)
      COMMON/ML5_0_0_1_LCOEFS/LOOPCOEFS
C     
C     BEGIN CODE
C     
      CONST=CMPLX((ONE*MULTIPLIER)/SYMFACT,ZERO,KIND=16)

      CALL ML5_0_0_1_MERGE_WL(LOOP_WF,RANK,LCUT_SIZE,CONST,LOOPCOEFS(0
     $ ,LOOP_GROUP_NUMBER))

      END

      SUBROUTINE ML5_0_0_1_INVERT_MOMENTA_IN_POLYNOMIAL(NCOEFS
     $ ,POLYNOMIAL)
C     Just a handy subroutine to modify the coefficients for the
C     tranformation q_loop -> -q_loop
C     It is only used for the NINJA interface
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE

      INTEGER I, NCOEFS

      COMPLEX*16 POLYNOMIAL(0:NCOEFS-1)

      DO I=0,NCOEFS-1
        IF (MOD(COEFTORANK_MAP(I),2).EQ.1) THEN
          POLYNOMIAL(I)=-POLYNOMIAL(I)
        ENDIF
      ENDDO

      END

C     Now the routines to update the wavefunctions



C     The subroutine to create the loop coefficients form the last
C      loop wf.
C     In this case of loop-induced process, the reduction is performed
C      at the loop
C     amplitude level so that no multiplication is performed.

      SUBROUTINE MP_ML5_0_0_1_CREATE_LOOP_COEFS(LOOP_WF,RANK,LCUT_SIZE
     $ ,LOOP_GROUP_NUMBER,SYMFACT,MULTIPLIER)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      REAL*16 ZERO,ONE
      PARAMETER (ZERO=0.0E0_16,ONE=1.0E0_16)
      COMPLEX*32 IMAG1
      PARAMETER (IMAG1=(ZERO,ONE))
      COMPLEX*32 CMPLX_ZERO
      PARAMETER (CMPLX_ZERO=(ZERO,ZERO))
      INTEGER    NLOOPGROUPS
      PARAMETER (NLOOPGROUPS=276)
      INTEGER    NCOMB
      PARAMETER (NCOMB=4)
C     
C     ARGUMENTS 
C     
      COMPLEX*32 LOOP_WF(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER RANK, SYMFACT, LCUT_SIZE, LOOP_GROUP_NUMBER, MULTIPLIER
C     
C     LOCAL VARIABLES 
C     
      COMPLEX*32 CONST
      INTEGER I,J
C     
C     GLOBAL VARIABLES
C     
      COMPLEX*32 LOOPCOEFS(0:LOOPMAXCOEFS-1,NLOOPGROUPS)
      COMMON/ML5_0_0_1_MP_LCOEFS/LOOPCOEFS
C     
C     BEGIN CODE
C     
      CONST=CMPLX((ONE*MULTIPLIER)/SYMFACT,ZERO,KIND=16)

      CALL MP_ML5_0_0_1_MERGE_WL(LOOP_WF,RANK,LCUT_SIZE,CONST
     $ ,LOOPCOEFS(0,LOOP_GROUP_NUMBER))

      END

      SUBROUTINE MP_ML5_0_0_1_INVERT_MOMENTA_IN_POLYNOMIAL(NCOEFS
     $ ,POLYNOMIAL)
C     Just a handy subroutine to modify the coefficients for the
C     tranformation q_loop -> -q_loop
C     It is only used for the NINJA interface
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE

      INTEGER I, NCOEFS

      COMPLEX*32 POLYNOMIAL(0:NCOEFS-1)

      DO I=0,NCOEFS-1
        IF (MOD(COEFTORANK_MAP(I),2).EQ.1) THEN
          POLYNOMIAL(I)=-POLYNOMIAL(I)
        ENDIF
      ENDDO

      END

C     Now the routines to update the wavefunctions



      SUBROUTINE ML5_0_0_1_EVAL_POLY(C,R,Q,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      COMPLEX*16 C(0:LOOPMAXCOEFS-1)
      INTEGER R
      COMPLEX*16 Q(0:3)
      COMPLEX*16 OUT

      OUT=C(0)
      IF (R.GE.1) THEN
        OUT=OUT+C(1)*Q(0)+C(2)*Q(1)+C(3)*Q(2)+C(4)*Q(3)
      ENDIF
      IF (R.GE.2) THEN
        OUT=OUT+C(5)*Q(0)*Q(0)+C(6)*Q(0)*Q(1)+C(7)*Q(1)*Q(1)+C(8)*Q(0)
     $   *Q(2)+C(9)*Q(1)*Q(2)+C(10)*Q(2)*Q(2)+C(11)*Q(0)*Q(3)+C(12)
     $   *Q(1)*Q(3)+C(13)*Q(2)*Q(3)+C(14)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.3) THEN
        OUT=OUT+C(15)*Q(0)*Q(0)*Q(0)+C(16)*Q(0)*Q(0)*Q(1)+C(17)*Q(0)
     $   *Q(1)*Q(1)+C(18)*Q(1)*Q(1)*Q(1)+C(19)*Q(0)*Q(0)*Q(2)+C(20)
     $   *Q(0)*Q(1)*Q(2)+C(21)*Q(1)*Q(1)*Q(2)+C(22)*Q(0)*Q(2)*Q(2)
     $   +C(23)*Q(1)*Q(2)*Q(2)+C(24)*Q(2)*Q(2)*Q(2)+C(25)*Q(0)*Q(0)
     $   *Q(3)+C(26)*Q(0)*Q(1)*Q(3)+C(27)*Q(1)*Q(1)*Q(3)+C(28)*Q(0)
     $   *Q(2)*Q(3)+C(29)*Q(1)*Q(2)*Q(3)+C(30)*Q(2)*Q(2)*Q(3)+C(31)
     $   *Q(0)*Q(3)*Q(3)+C(32)*Q(1)*Q(3)*Q(3)+C(33)*Q(2)*Q(3)*Q(3)
     $   +C(34)*Q(3)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.4) THEN
        OUT=OUT+C(35)*Q(0)*Q(0)*Q(0)*Q(0)+C(36)*Q(0)*Q(0)*Q(0)*Q(1)
     $   +C(37)*Q(0)*Q(0)*Q(1)*Q(1)+C(38)*Q(0)*Q(1)*Q(1)*Q(1)+C(39)
     $   *Q(1)*Q(1)*Q(1)*Q(1)+C(40)*Q(0)*Q(0)*Q(0)*Q(2)+C(41)*Q(0)*Q(0)
     $   *Q(1)*Q(2)+C(42)*Q(0)*Q(1)*Q(1)*Q(2)+C(43)*Q(1)*Q(1)*Q(1)*Q(2)
     $   +C(44)*Q(0)*Q(0)*Q(2)*Q(2)+C(45)*Q(0)*Q(1)*Q(2)*Q(2)+C(46)
     $   *Q(1)*Q(1)*Q(2)*Q(2)+C(47)*Q(0)*Q(2)*Q(2)*Q(2)+C(48)*Q(1)*Q(2)
     $   *Q(2)*Q(2)+C(49)*Q(2)*Q(2)*Q(2)*Q(2)+C(50)*Q(0)*Q(0)*Q(0)*Q(3)
     $   +C(51)*Q(0)*Q(0)*Q(1)*Q(3)+C(52)*Q(0)*Q(1)*Q(1)*Q(3)+C(53)
     $   *Q(1)*Q(1)*Q(1)*Q(3)+C(54)*Q(0)*Q(0)*Q(2)*Q(3)+C(55)*Q(0)*Q(1)
     $   *Q(2)*Q(3)+C(56)*Q(1)*Q(1)*Q(2)*Q(3)+C(57)*Q(0)*Q(2)*Q(2)*Q(3)
     $   +C(58)*Q(1)*Q(2)*Q(2)*Q(3)+C(59)*Q(2)*Q(2)*Q(2)*Q(3)+C(60)
     $   *Q(0)*Q(0)*Q(3)*Q(3)+C(61)*Q(0)*Q(1)*Q(3)*Q(3)+C(62)*Q(1)*Q(1)
     $   *Q(3)*Q(3)+C(63)*Q(0)*Q(2)*Q(3)*Q(3)+C(64)*Q(1)*Q(2)*Q(3)*Q(3)
        OUT=OUT+C(65)*Q(2)*Q(2)*Q(3)*Q(3)+C(66)*Q(0)*Q(3)*Q(3)*Q(3)
     $   +C(67)*Q(1)*Q(3)*Q(3)*Q(3)+C(68)*Q(2)*Q(3)*Q(3)*Q(3)+C(69)
     $   *Q(3)*Q(3)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.5) THEN
        OUT=OUT+C(70)*Q(0)*Q(0)*Q(0)*Q(0)*Q(0)+C(71)*Q(0)*Q(0)*Q(0)
     $   *Q(0)*Q(1)+C(72)*Q(0)*Q(0)*Q(0)*Q(1)*Q(1)+C(73)*Q(0)*Q(0)*Q(1)
     $   *Q(1)*Q(1)+C(74)*Q(0)*Q(1)*Q(1)*Q(1)*Q(1)+C(75)*Q(1)*Q(1)*Q(1)
     $   *Q(1)*Q(1)+C(76)*Q(0)*Q(0)*Q(0)*Q(0)*Q(2)+C(77)*Q(0)*Q(0)*Q(0)
     $   *Q(1)*Q(2)+C(78)*Q(0)*Q(0)*Q(1)*Q(1)*Q(2)+C(79)*Q(0)*Q(1)*Q(1)
     $   *Q(1)*Q(2)+C(80)*Q(1)*Q(1)*Q(1)*Q(1)*Q(2)+C(81)*Q(0)*Q(0)*Q(0)
     $   *Q(2)*Q(2)+C(82)*Q(0)*Q(0)*Q(1)*Q(2)*Q(2)+C(83)*Q(0)*Q(1)*Q(1)
     $   *Q(2)*Q(2)+C(84)*Q(1)*Q(1)*Q(1)*Q(2)*Q(2)+C(85)*Q(0)*Q(0)*Q(2)
     $   *Q(2)*Q(2)+C(86)*Q(0)*Q(1)*Q(2)*Q(2)*Q(2)+C(87)*Q(1)*Q(1)*Q(2)
     $   *Q(2)*Q(2)+C(88)*Q(0)*Q(2)*Q(2)*Q(2)*Q(2)+C(89)*Q(1)*Q(2)*Q(2)
     $   *Q(2)*Q(2)+C(90)*Q(2)*Q(2)*Q(2)*Q(2)*Q(2)+C(91)*Q(0)*Q(0)*Q(0)
     $   *Q(0)*Q(3)+C(92)*Q(0)*Q(0)*Q(0)*Q(1)*Q(3)+C(93)*Q(0)*Q(0)*Q(1)
     $   *Q(1)*Q(3)+C(94)*Q(0)*Q(1)*Q(1)*Q(1)*Q(3)+C(95)*Q(1)*Q(1)*Q(1)
     $   *Q(1)*Q(3)+C(96)*Q(0)*Q(0)*Q(0)*Q(2)*Q(3)+C(97)*Q(0)*Q(0)*Q(1)
     $   *Q(2)*Q(3)+C(98)*Q(0)*Q(1)*Q(1)*Q(2)*Q(3)+C(99)*Q(1)*Q(1)*Q(1)
     $   *Q(2)*Q(3)
        OUT=OUT+C(100)*Q(0)*Q(0)*Q(2)*Q(2)*Q(3)+C(101)*Q(0)*Q(1)*Q(2)
     $   *Q(2)*Q(3)+C(102)*Q(1)*Q(1)*Q(2)*Q(2)*Q(3)+C(103)*Q(0)*Q(2)
     $   *Q(2)*Q(2)*Q(3)+C(104)*Q(1)*Q(2)*Q(2)*Q(2)*Q(3)+C(105)*Q(2)
     $   *Q(2)*Q(2)*Q(2)*Q(3)+C(106)*Q(0)*Q(0)*Q(0)*Q(3)*Q(3)+C(107)
     $   *Q(0)*Q(0)*Q(1)*Q(3)*Q(3)+C(108)*Q(0)*Q(1)*Q(1)*Q(3)*Q(3)
     $   +C(109)*Q(1)*Q(1)*Q(1)*Q(3)*Q(3)+C(110)*Q(0)*Q(0)*Q(2)*Q(3)
     $   *Q(3)+C(111)*Q(0)*Q(1)*Q(2)*Q(3)*Q(3)+C(112)*Q(1)*Q(1)*Q(2)
     $   *Q(3)*Q(3)+C(113)*Q(0)*Q(2)*Q(2)*Q(3)*Q(3)+C(114)*Q(1)*Q(2)
     $   *Q(2)*Q(3)*Q(3)+C(115)*Q(2)*Q(2)*Q(2)*Q(3)*Q(3)+C(116)*Q(0)
     $   *Q(0)*Q(3)*Q(3)*Q(3)+C(117)*Q(0)*Q(1)*Q(3)*Q(3)*Q(3)+C(118)
     $   *Q(1)*Q(1)*Q(3)*Q(3)*Q(3)+C(119)*Q(0)*Q(2)*Q(3)*Q(3)*Q(3)
     $   +C(120)*Q(1)*Q(2)*Q(3)*Q(3)*Q(3)+C(121)*Q(2)*Q(2)*Q(3)*Q(3)
     $   *Q(3)+C(122)*Q(0)*Q(3)*Q(3)*Q(3)*Q(3)+C(123)*Q(1)*Q(3)*Q(3)
     $   *Q(3)*Q(3)+C(124)*Q(2)*Q(3)*Q(3)*Q(3)*Q(3)+C(125)*Q(3)*Q(3)
     $   *Q(3)*Q(3)*Q(3)
      ENDIF
      END

      SUBROUTINE MP_ML5_0_0_1_EVAL_POLY(C,R,Q,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      COMPLEX*32 C(0:LOOPMAXCOEFS-1)
      INTEGER R
      COMPLEX*32 Q(0:3)
      COMPLEX*32 OUT

      OUT=C(0)
      IF (R.GE.1) THEN
        OUT=OUT+C(1)*Q(0)+C(2)*Q(1)+C(3)*Q(2)+C(4)*Q(3)
      ENDIF
      IF (R.GE.2) THEN
        OUT=OUT+C(5)*Q(0)*Q(0)+C(6)*Q(0)*Q(1)+C(7)*Q(1)*Q(1)+C(8)*Q(0)
     $   *Q(2)+C(9)*Q(1)*Q(2)+C(10)*Q(2)*Q(2)+C(11)*Q(0)*Q(3)+C(12)
     $   *Q(1)*Q(3)+C(13)*Q(2)*Q(3)+C(14)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.3) THEN
        OUT=OUT+C(15)*Q(0)*Q(0)*Q(0)+C(16)*Q(0)*Q(0)*Q(1)+C(17)*Q(0)
     $   *Q(1)*Q(1)+C(18)*Q(1)*Q(1)*Q(1)+C(19)*Q(0)*Q(0)*Q(2)+C(20)
     $   *Q(0)*Q(1)*Q(2)+C(21)*Q(1)*Q(1)*Q(2)+C(22)*Q(0)*Q(2)*Q(2)
     $   +C(23)*Q(1)*Q(2)*Q(2)+C(24)*Q(2)*Q(2)*Q(2)+C(25)*Q(0)*Q(0)
     $   *Q(3)+C(26)*Q(0)*Q(1)*Q(3)+C(27)*Q(1)*Q(1)*Q(3)+C(28)*Q(0)
     $   *Q(2)*Q(3)+C(29)*Q(1)*Q(2)*Q(3)+C(30)*Q(2)*Q(2)*Q(3)+C(31)
     $   *Q(0)*Q(3)*Q(3)+C(32)*Q(1)*Q(3)*Q(3)+C(33)*Q(2)*Q(3)*Q(3)
     $   +C(34)*Q(3)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.4) THEN
        OUT=OUT+C(35)*Q(0)*Q(0)*Q(0)*Q(0)+C(36)*Q(0)*Q(0)*Q(0)*Q(1)
     $   +C(37)*Q(0)*Q(0)*Q(1)*Q(1)+C(38)*Q(0)*Q(1)*Q(1)*Q(1)+C(39)
     $   *Q(1)*Q(1)*Q(1)*Q(1)+C(40)*Q(0)*Q(0)*Q(0)*Q(2)+C(41)*Q(0)*Q(0)
     $   *Q(1)*Q(2)+C(42)*Q(0)*Q(1)*Q(1)*Q(2)+C(43)*Q(1)*Q(1)*Q(1)*Q(2)
     $   +C(44)*Q(0)*Q(0)*Q(2)*Q(2)+C(45)*Q(0)*Q(1)*Q(2)*Q(2)+C(46)
     $   *Q(1)*Q(1)*Q(2)*Q(2)+C(47)*Q(0)*Q(2)*Q(2)*Q(2)+C(48)*Q(1)*Q(2)
     $   *Q(2)*Q(2)+C(49)*Q(2)*Q(2)*Q(2)*Q(2)+C(50)*Q(0)*Q(0)*Q(0)*Q(3)
     $   +C(51)*Q(0)*Q(0)*Q(1)*Q(3)+C(52)*Q(0)*Q(1)*Q(1)*Q(3)+C(53)
     $   *Q(1)*Q(1)*Q(1)*Q(3)+C(54)*Q(0)*Q(0)*Q(2)*Q(3)+C(55)*Q(0)*Q(1)
     $   *Q(2)*Q(3)+C(56)*Q(1)*Q(1)*Q(2)*Q(3)+C(57)*Q(0)*Q(2)*Q(2)*Q(3)
     $   +C(58)*Q(1)*Q(2)*Q(2)*Q(3)+C(59)*Q(2)*Q(2)*Q(2)*Q(3)+C(60)
     $   *Q(0)*Q(0)*Q(3)*Q(3)+C(61)*Q(0)*Q(1)*Q(3)*Q(3)+C(62)*Q(1)*Q(1)
     $   *Q(3)*Q(3)+C(63)*Q(0)*Q(2)*Q(3)*Q(3)+C(64)*Q(1)*Q(2)*Q(3)*Q(3)
        OUT=OUT+C(65)*Q(2)*Q(2)*Q(3)*Q(3)+C(66)*Q(0)*Q(3)*Q(3)*Q(3)
     $   +C(67)*Q(1)*Q(3)*Q(3)*Q(3)+C(68)*Q(2)*Q(3)*Q(3)*Q(3)+C(69)
     $   *Q(3)*Q(3)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.5) THEN
        OUT=OUT+C(70)*Q(0)*Q(0)*Q(0)*Q(0)*Q(0)+C(71)*Q(0)*Q(0)*Q(0)
     $   *Q(0)*Q(1)+C(72)*Q(0)*Q(0)*Q(0)*Q(1)*Q(1)+C(73)*Q(0)*Q(0)*Q(1)
     $   *Q(1)*Q(1)+C(74)*Q(0)*Q(1)*Q(1)*Q(1)*Q(1)+C(75)*Q(1)*Q(1)*Q(1)
     $   *Q(1)*Q(1)+C(76)*Q(0)*Q(0)*Q(0)*Q(0)*Q(2)+C(77)*Q(0)*Q(0)*Q(0)
     $   *Q(1)*Q(2)+C(78)*Q(0)*Q(0)*Q(1)*Q(1)*Q(2)+C(79)*Q(0)*Q(1)*Q(1)
     $   *Q(1)*Q(2)+C(80)*Q(1)*Q(1)*Q(1)*Q(1)*Q(2)+C(81)*Q(0)*Q(0)*Q(0)
     $   *Q(2)*Q(2)+C(82)*Q(0)*Q(0)*Q(1)*Q(2)*Q(2)+C(83)*Q(0)*Q(1)*Q(1)
     $   *Q(2)*Q(2)+C(84)*Q(1)*Q(1)*Q(1)*Q(2)*Q(2)+C(85)*Q(0)*Q(0)*Q(2)
     $   *Q(2)*Q(2)+C(86)*Q(0)*Q(1)*Q(2)*Q(2)*Q(2)+C(87)*Q(1)*Q(1)*Q(2)
     $   *Q(2)*Q(2)+C(88)*Q(0)*Q(2)*Q(2)*Q(2)*Q(2)+C(89)*Q(1)*Q(2)*Q(2)
     $   *Q(2)*Q(2)+C(90)*Q(2)*Q(2)*Q(2)*Q(2)*Q(2)+C(91)*Q(0)*Q(0)*Q(0)
     $   *Q(0)*Q(3)+C(92)*Q(0)*Q(0)*Q(0)*Q(1)*Q(3)+C(93)*Q(0)*Q(0)*Q(1)
     $   *Q(1)*Q(3)+C(94)*Q(0)*Q(1)*Q(1)*Q(1)*Q(3)+C(95)*Q(1)*Q(1)*Q(1)
     $   *Q(1)*Q(3)+C(96)*Q(0)*Q(0)*Q(0)*Q(2)*Q(3)+C(97)*Q(0)*Q(0)*Q(1)
     $   *Q(2)*Q(3)+C(98)*Q(0)*Q(1)*Q(1)*Q(2)*Q(3)+C(99)*Q(1)*Q(1)*Q(1)
     $   *Q(2)*Q(3)
        OUT=OUT+C(100)*Q(0)*Q(0)*Q(2)*Q(2)*Q(3)+C(101)*Q(0)*Q(1)*Q(2)
     $   *Q(2)*Q(3)+C(102)*Q(1)*Q(1)*Q(2)*Q(2)*Q(3)+C(103)*Q(0)*Q(2)
     $   *Q(2)*Q(2)*Q(3)+C(104)*Q(1)*Q(2)*Q(2)*Q(2)*Q(3)+C(105)*Q(2)
     $   *Q(2)*Q(2)*Q(2)*Q(3)+C(106)*Q(0)*Q(0)*Q(0)*Q(3)*Q(3)+C(107)
     $   *Q(0)*Q(0)*Q(1)*Q(3)*Q(3)+C(108)*Q(0)*Q(1)*Q(1)*Q(3)*Q(3)
     $   +C(109)*Q(1)*Q(1)*Q(1)*Q(3)*Q(3)+C(110)*Q(0)*Q(0)*Q(2)*Q(3)
     $   *Q(3)+C(111)*Q(0)*Q(1)*Q(2)*Q(3)*Q(3)+C(112)*Q(1)*Q(1)*Q(2)
     $   *Q(3)*Q(3)+C(113)*Q(0)*Q(2)*Q(2)*Q(3)*Q(3)+C(114)*Q(1)*Q(2)
     $   *Q(2)*Q(3)*Q(3)+C(115)*Q(2)*Q(2)*Q(2)*Q(3)*Q(3)+C(116)*Q(0)
     $   *Q(0)*Q(3)*Q(3)*Q(3)+C(117)*Q(0)*Q(1)*Q(3)*Q(3)*Q(3)+C(118)
     $   *Q(1)*Q(1)*Q(3)*Q(3)*Q(3)+C(119)*Q(0)*Q(2)*Q(3)*Q(3)*Q(3)
     $   +C(120)*Q(1)*Q(2)*Q(3)*Q(3)*Q(3)+C(121)*Q(2)*Q(2)*Q(3)*Q(3)
     $   *Q(3)+C(122)*Q(0)*Q(3)*Q(3)*Q(3)*Q(3)+C(123)*Q(1)*Q(3)*Q(3)
     $   *Q(3)*Q(3)+C(124)*Q(2)*Q(3)*Q(3)*Q(3)*Q(3)+C(125)*Q(3)*Q(3)
     $   *Q(3)*Q(3)*Q(3)
      ENDIF
      END

      SUBROUTINE ML5_0_0_1_ADD_COEFS(A,RA,B,RB)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      INTEGER I
      COMPLEX*16 A(0:LOOPMAXCOEFS-1),B(0:LOOPMAXCOEFS-1)
      INTEGER RA,RB

      DO I=0,NCOEF_R(RB)-1
        A(I)=A(I)+B(I)
      ENDDO
      END

      SUBROUTINE MP_ML5_0_0_1_ADD_COEFS(A,RA,B,RB)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      INTEGER I
      COMPLEX*32 A(0:LOOPMAXCOEFS-1),B(0:LOOPMAXCOEFS-1)
      INTEGER RA,RB

      DO I=0,NCOEF_R(RB)-1
        A(I)=A(I)+B(I)
      ENDDO
      END

      SUBROUTINE ML5_0_0_1_MERGE_WL(WL,R,LCUT_SIZE,CONST,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      INTEGER I,J
      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER R,LCUT_SIZE
      COMPLEX*16 CONST
      COMPLEX*16 OUT(0:LOOPMAXCOEFS-1)

      DO I=1,LCUT_SIZE
        DO J=0,NCOEF_R(R)-1
          OUT(J)=OUT(J)+WL(I,J,I)*CONST
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_0_0_1_MERGE_WL(WL,R,LCUT_SIZE,CONST,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      INTEGER I,J
      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER R,LCUT_SIZE
      COMPLEX*32 CONST
      COMPLEX*32 OUT(0:LOOPMAXCOEFS-1)

      DO I=1,LCUT_SIZE
        DO J=0,NCOEF_R(R)-1
          OUT(J)=OUT(J)+WL(I,J,I)*CONST
        ENDDO
      ENDDO
      END

      SUBROUTINE ML5_0_0_1_UPDATE_WL_0_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_0_0_1_UPDATE_WL_0_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE ML5_0_0_1_UPDATE_WL_3_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INTEGER I,J,K,L,M
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE
      INTEGER NEW_POSITION
      COMPLEX*16 UPDATER_COEF

C     Welcome to the computational heart of MadLoop...
      OUT(:,:,:)=(0.0D0,0.0D0)
      DO J=1,OUT_SIZE
        DO M=0,4
          DO K=1,IN_SIZE
            UPDATER_COEF = B(J,M,K)
            IF (UPDATER_COEF.EQ.(0.0D0,0.0D0)) CYCLE
            DO L=0,34
              NEW_POSITION = COMB_COEF_POS(L,M)
              DO I=1,LCUT_SIZE
                OUT(J,NEW_POSITION,I)=OUT(J,NEW_POSITION,I) + A(K,L,I)
     $           *UPDATER_COEF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END

      SUBROUTINE MP_ML5_0_0_1_UPDATE_WL_3_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INTEGER I,J,K,L,M
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE
      INTEGER NEW_POSITION
      COMPLEX*32 UPDATER_COEF

C     Welcome to the computational heart of MadLoop...
      OUT(:,:,:)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
      DO J=1,OUT_SIZE
        DO M=0,4
          DO K=1,IN_SIZE
            UPDATER_COEF = B(J,M,K)
            IF (UPDATER_COEF.EQ.CMPLX(0.0E0_16,0.0E0_16,KIND=16)) CYCLE
            DO L=0,34
              NEW_POSITION = COMB_COEF_POS(L,M)
              DO I=1,LCUT_SIZE
                OUT(J,NEW_POSITION,I)=OUT(J,NEW_POSITION,I) + A(K,L,I)
     $           *UPDATER_COEF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END

      SUBROUTINE ML5_0_0_1_UPDATE_WL_4_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INTEGER I,J,K,L,M
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE
      INTEGER NEW_POSITION
      COMPLEX*16 UPDATER_COEF

C     Welcome to the computational heart of MadLoop...
      OUT(:,:,:)=(0.0D0,0.0D0)
      DO J=1,OUT_SIZE
        DO M=0,4
          DO K=1,IN_SIZE
            UPDATER_COEF = B(J,M,K)
            IF (UPDATER_COEF.EQ.(0.0D0,0.0D0)) CYCLE
            DO L=0,69
              NEW_POSITION = COMB_COEF_POS(L,M)
              DO I=1,LCUT_SIZE
                OUT(J,NEW_POSITION,I)=OUT(J,NEW_POSITION,I) + A(K,L,I)
     $           *UPDATER_COEF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END

      SUBROUTINE MP_ML5_0_0_1_UPDATE_WL_4_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INTEGER I,J,K,L,M
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE
      INTEGER NEW_POSITION
      COMPLEX*32 UPDATER_COEF

C     Welcome to the computational heart of MadLoop...
      OUT(:,:,:)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
      DO J=1,OUT_SIZE
        DO M=0,4
          DO K=1,IN_SIZE
            UPDATER_COEF = B(J,M,K)
            IF (UPDATER_COEF.EQ.CMPLX(0.0E0_16,0.0E0_16,KIND=16)) CYCLE
            DO L=0,69
              NEW_POSITION = COMB_COEF_POS(L,M)
              DO I=1,LCUT_SIZE
                OUT(J,NEW_POSITION,I)=OUT(J,NEW_POSITION,I) + A(K,L,I)
     $           *UPDATER_COEF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END

      SUBROUTINE ML5_0_0_1_UPDATE_WL_1_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)*B(J,1,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)*B(J,2,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)*B(J,3,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_0_0_1_UPDATE_WL_1_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)*B(J,1,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)*B(J,2,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)*B(J,3,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE ML5_0_0_1_UPDATE_WL_2_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INTEGER I,J,K,L,M
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE
      INTEGER NEW_POSITION
      COMPLEX*16 UPDATER_COEF

C     Welcome to the computational heart of MadLoop...
      OUT(:,:,:)=(0.0D0,0.0D0)
      DO J=1,OUT_SIZE
        DO M=0,4
          DO K=1,IN_SIZE
            UPDATER_COEF = B(J,M,K)
            IF (UPDATER_COEF.EQ.(0.0D0,0.0D0)) CYCLE
            DO L=0,14
              NEW_POSITION = COMB_COEF_POS(L,M)
              DO I=1,LCUT_SIZE
                OUT(J,NEW_POSITION,I)=OUT(J,NEW_POSITION,I) + A(K,L,I)
     $           *UPDATER_COEF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END

      SUBROUTINE MP_ML5_0_0_1_UPDATE_WL_2_1(A,LCUT_SIZE,B,IN_SIZE
     $ ,OUT_SIZE,OUT)
      USE ML5_0_0_1_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INTEGER I,J,K,L,M
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE
      INTEGER NEW_POSITION
      COMPLEX*32 UPDATER_COEF

C     Welcome to the computational heart of MadLoop...
      OUT(:,:,:)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
      DO J=1,OUT_SIZE
        DO M=0,4
          DO K=1,IN_SIZE
            UPDATER_COEF = B(J,M,K)
            IF (UPDATER_COEF.EQ.CMPLX(0.0E0_16,0.0E0_16,KIND=16)) CYCLE
            DO L=0,14
              NEW_POSITION = COMB_COEF_POS(L,M)
              DO I=1,LCUT_SIZE
                OUT(J,NEW_POSITION,I)=OUT(J,NEW_POSITION,I) + A(K,L,I)
     $           *UPDATER_COEF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END
