ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL*16 MP__MDL_SQRT__AS,MP__MDL_G__EXP__4,MP__MDL_G__EXP__2
     $ ,MP__MDL_GWCFT_UV_B_1EPS_,MP__MDL_GWCFT_UV_T_1EPS_
     $ ,MP__MDL_BWCFT_UV_1EPS_,MP__MDL_TWCFT_UV_1EPS_
     $ ,MP__MDL_G__EXP__3,MP__MDL_MU_R__EXP__2,MP__MDL_GWCFT_UV_B_FIN_
     $ ,MP__MDL_GWCFT_UV_T_FIN_,MP__MDL_BWCFT_UV_FIN_
     $ ,MP__MDL_TWCFT_UV_FIN_,MP__MDL_CKM33,MP__MDL_CONJG__CKM3X3
     $ ,MP__MDL_CKM22,MP__MDL_LHV,MP__MDL_CKM3X3,MP__MDL_CONJG__CKM22
     $ ,MP__MDL_NCOL,MP__MDL_CA,MP__MDL_TF,MP__MDL_CF
     $ ,MP__MDL_CTHETA1__EXP__2,MP__MDL_STHETA1
     $ ,MP__MDL_CTHETA2__EXP__2,MP__MDL_STHETA2
     $ ,MP__MDL_CTHETA3__EXP__2,MP__MDL_STHETA3,MP__MDL_R23
     $ ,MP__MDL_R32,MP__MDL_R22,MP__MDL_R21,MP__MDL_R31,MP__MDL_R11
     $ ,MP__MDL_R33,MP__MDL_R13,MP__MDL_R12,MP__MDL_RH1,MP__MDL_RH2
     $ ,MP__MDL_RH3,MP__MDL_MZ__EXP__2,MP__MDL_MZ__EXP__4
     $ ,MP__MDL_SQRT__2,MP__MDL_MH__EXP__2,MP__MDL_MH__EXP__4
     $ ,MP__MDL_MT__EXP__4,MP__MDL_MT__EXP__2,MP__MDL_MH__EXP__12
     $ ,MP__MDL_MH__EXP__10,MP__MDL_MH__EXP__8,MP__MDL_MH__EXP__6
     $ ,MP__MDL_NCOL__EXP__2,MP__MDL_MB__EXP__2,MP__MDL_RH1__EXP__2
     $ ,MP__MDL_RH2__EXP__2,MP__MDL_RH3__EXP__2,MP__MDL_AEW
     $ ,MP__MDL_SQRT__AEW,MP__MDL_EE,MP__MDL_MW__EXP__2,MP__MDL_SW2
     $ ,MP__MDL_CW,MP__MDL_SQRT__SW2,MP__MDL_SW,MP__MDL_G1,MP__MDL_GW
     $ ,MP__MDL_V,MP__MDL_V__EXP__2,MP__MDL_LAM,MP__MDL_YB,MP__MDL_YT
     $ ,MP__MDL_YTAU,MP__MDL_MUH,MP__MDL_EE__EXP__2
     $ ,MP__MDL_MW__EXP__12,MP__MDL_MW__EXP__10,MP__MDL_MW__EXP__8
     $ ,MP__MDL_MW__EXP__6,MP__MDL_MW__EXP__4,MP__MDL_AH
     $ ,MP__MDL_AXIALZUP,MP__MDL_AXIALZDOWN,MP__MDL_VECTORZUP
     $ ,MP__MDL_VECTORZDOWN,MP__MDL_VECTORAUP,MP__MDL_VECTORADOWN
     $ ,MP__MDL_VECTORWMDXU,MP__MDL_AXIALWMDXU,MP__MDL_VECTORWPUXD
     $ ,MP__MDL_AXIALWPUXD,MP__MDL_GW__EXP__2,MP__MDL_CW__EXP__2
     $ ,MP__MDL_SW__EXP__2,MP__MDL_YB__EXP__2,MP__MDL_YT__EXP__2
     $ ,MP__MDL_CTHETA1,MP__MDL_CTHETA2,MP__MDL_CTHETA3,MP__MDL_KAP111
     $ ,MP__MDL_KAP112,MP__MDL_KAP122,MP__MDL_KAP222,MP__MDL_KAP133
     $ ,MP__MDL_KAP113,MP__MDL_KAP123,MP__MDL_KAP333,MP__MDL_KAP233
     $ ,MP__MDL_KAP223,MP__MDL_KAP1111,MP__MDL_KAP1112,MP__MDL_KAP1122
     $ ,MP__MDL_KAP1222,MP__MDL_KAP2222,MP__MDL_KAP1113
     $ ,MP__MDL_KAP1133,MP__MDL_KAP1333,MP__MDL_KAP3333
     $ ,MP__MDL_KAP2223,MP__MDL_KAP2333,MP__MDL_KAP2233
     $ ,MP__MDL_KAP1233,MP__MDL_KAP1223,MP__MDL_KAP1123,MP__AEWM1
     $ ,MP__MDL_GF,MP__AS,MP__MDL_YMB,MP__MDL_YMT,MP__MDL_YMTAU

      COMMON/MP_PARAMS_R/ MP__MDL_SQRT__AS,MP__MDL_G__EXP__4
     $ ,MP__MDL_G__EXP__2,MP__MDL_GWCFT_UV_B_1EPS_
     $ ,MP__MDL_GWCFT_UV_T_1EPS_,MP__MDL_BWCFT_UV_1EPS_
     $ ,MP__MDL_TWCFT_UV_1EPS_,MP__MDL_G__EXP__3,MP__MDL_MU_R__EXP__2
     $ ,MP__MDL_GWCFT_UV_B_FIN_,MP__MDL_GWCFT_UV_T_FIN_
     $ ,MP__MDL_BWCFT_UV_FIN_,MP__MDL_TWCFT_UV_FIN_,MP__MDL_CKM33
     $ ,MP__MDL_CONJG__CKM3X3,MP__MDL_CKM22,MP__MDL_LHV,MP__MDL_CKM3X3
     $ ,MP__MDL_CONJG__CKM22,MP__MDL_NCOL,MP__MDL_CA,MP__MDL_TF
     $ ,MP__MDL_CF,MP__MDL_CTHETA1__EXP__2,MP__MDL_STHETA1
     $ ,MP__MDL_CTHETA2__EXP__2,MP__MDL_STHETA2
     $ ,MP__MDL_CTHETA3__EXP__2,MP__MDL_STHETA3,MP__MDL_R23
     $ ,MP__MDL_R32,MP__MDL_R22,MP__MDL_R21,MP__MDL_R31,MP__MDL_R11
     $ ,MP__MDL_R33,MP__MDL_R13,MP__MDL_R12,MP__MDL_RH1,MP__MDL_RH2
     $ ,MP__MDL_RH3,MP__MDL_MZ__EXP__2,MP__MDL_MZ__EXP__4
     $ ,MP__MDL_SQRT__2,MP__MDL_MH__EXP__2,MP__MDL_MH__EXP__4
     $ ,MP__MDL_MT__EXP__4,MP__MDL_MT__EXP__2,MP__MDL_MH__EXP__12
     $ ,MP__MDL_MH__EXP__10,MP__MDL_MH__EXP__8,MP__MDL_MH__EXP__6
     $ ,MP__MDL_NCOL__EXP__2,MP__MDL_MB__EXP__2,MP__MDL_RH1__EXP__2
     $ ,MP__MDL_RH2__EXP__2,MP__MDL_RH3__EXP__2,MP__MDL_AEW
     $ ,MP__MDL_SQRT__AEW,MP__MDL_EE,MP__MDL_MW__EXP__2,MP__MDL_SW2
     $ ,MP__MDL_CW,MP__MDL_SQRT__SW2,MP__MDL_SW,MP__MDL_G1,MP__MDL_GW
     $ ,MP__MDL_V,MP__MDL_V__EXP__2,MP__MDL_LAM,MP__MDL_YB,MP__MDL_YT
     $ ,MP__MDL_YTAU,MP__MDL_MUH,MP__MDL_EE__EXP__2
     $ ,MP__MDL_MW__EXP__12,MP__MDL_MW__EXP__10,MP__MDL_MW__EXP__8
     $ ,MP__MDL_MW__EXP__6,MP__MDL_MW__EXP__4,MP__MDL_AH
     $ ,MP__MDL_AXIALZUP,MP__MDL_AXIALZDOWN,MP__MDL_VECTORZUP
     $ ,MP__MDL_VECTORZDOWN,MP__MDL_VECTORAUP,MP__MDL_VECTORADOWN
     $ ,MP__MDL_VECTORWMDXU,MP__MDL_AXIALWMDXU,MP__MDL_VECTORWPUXD
     $ ,MP__MDL_AXIALWPUXD,MP__MDL_GW__EXP__2,MP__MDL_CW__EXP__2
     $ ,MP__MDL_SW__EXP__2,MP__MDL_YB__EXP__2,MP__MDL_YT__EXP__2
     $ ,MP__MDL_CTHETA1,MP__MDL_CTHETA2,MP__MDL_CTHETA3,MP__MDL_KAP111
     $ ,MP__MDL_KAP112,MP__MDL_KAP122,MP__MDL_KAP222,MP__MDL_KAP133
     $ ,MP__MDL_KAP113,MP__MDL_KAP123,MP__MDL_KAP333,MP__MDL_KAP233
     $ ,MP__MDL_KAP223,MP__MDL_KAP1111,MP__MDL_KAP1112,MP__MDL_KAP1122
     $ ,MP__MDL_KAP1222,MP__MDL_KAP2222,MP__MDL_KAP1113
     $ ,MP__MDL_KAP1133,MP__MDL_KAP1333,MP__MDL_KAP3333
     $ ,MP__MDL_KAP2223,MP__MDL_KAP2333,MP__MDL_KAP2233
     $ ,MP__MDL_KAP1233,MP__MDL_KAP1223,MP__MDL_KAP1123,MP__AEWM1
     $ ,MP__MDL_GF,MP__AS,MP__MDL_YMB,MP__MDL_YMT,MP__MDL_YMTAU


      COMPLEX*32 MP__MDL_COMPLEXI,MP__MDL_I1X33,MP__MDL_I2X33
     $ ,MP__MDL_I3X33,MP__MDL_I4X33,MP__MDL_VECTOR_TBGP
     $ ,MP__MDL_AXIAL_TBGP,MP__MDL_VECTOR_TBGM,MP__MDL_AXIAL_TBGM

      COMMON/MP_PARAMS_C/ MP__MDL_COMPLEXI,MP__MDL_I1X33,MP__MDL_I2X33
     $ ,MP__MDL_I3X33,MP__MDL_I4X33,MP__MDL_VECTOR_TBGP
     $ ,MP__MDL_AXIAL_TBGP,MP__MDL_VECTOR_TBGM,MP__MDL_AXIAL_TBGM

