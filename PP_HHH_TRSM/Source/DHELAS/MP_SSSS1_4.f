C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     1
C     
      SUBROUTINE MP_SSSS1_4(S1, S2, S3, COUP, M4, W4,S4)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 S3(*)
      REAL*16 W4
      COMPLEX*32 P4(0:3)
      COMPLEX*32 S2(*)
      COMPLEX*32 S1(*)
      COMPLEX*32 DENOM
      REAL*16 M4
      COMPLEX*32 COUP
      COMPLEX*32 S4(5)
      S4(1) = +S1(1)+S2(1)+S3(1)
      S4(2) = +S1(2)+S2(2)+S3(2)
      S4(3) = +S1(3)+S2(3)+S3(3)
      S4(4) = +S1(4)+S2(4)+S3(4)
      P4(0) = -S4(1)
      P4(1) = -S4(2)
      P4(2) = -S4(3)
      P4(3) = -S4(4)
      DENOM = COUP/(P4(0)**2-P4(1)**2-P4(2)**2-P4(3)**2 - M4 * (M4 -CI
     $ * W4))
      S4(5)= DENOM*CI * S3(5)*S2(5)*S1(5)
      END


