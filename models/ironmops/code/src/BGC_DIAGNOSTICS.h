C$Header: /Users/ikriest/CVS/mops/BGC_DIAGNOSTICS.h,v 1.3 2018/03/12 06:44:38 ikriest Exp $
C$Name: mops-2_3 $

! arrays for diagnostics
      real*8 f1_out(bgc_ktotal),
     &       f2_out(bgc_ktotal),
     &       f3_out(bgc_ktotal),
     &       f4_out(bgc_ktotal),
     &       f5_out(bgc_ktotal),
     &       f6_out(bgc_ktotal),
     &       f7_out(bgc_ktotal),
     &       f8_out(bgc_ktotal),
     &       f9_out(bgc_ktotal),
     &       f10_out(bgc_ktotal),
     &       f11_out(bgc_ktotal),
     &       f12_out(bgc_ktotal)

      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
     &                f8_out,f9_out,f10_out,f11_out,f12_out 

