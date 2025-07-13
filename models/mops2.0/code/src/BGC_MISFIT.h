C$Header: /Users/ikriest/CVS/mops/BGC_MISFIT.h,v 1.1 2015/11/17 14:18:51 ikriest Exp $
C$Name: mops-2_0 $

! arrays for diagnostics
      real*8 m1_out(bgc_ktotal),
     &       m2_out(bgc_ktotal),
     &       m3_out(bgc_ktotal)

      COMMON/COSTVARS/m1_out,m2_out,m3_out

