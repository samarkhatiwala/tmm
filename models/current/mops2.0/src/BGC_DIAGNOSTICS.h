C$Header: /Users/ikriest/CVS/mops/BGC_DIAGNOSTICS.h,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

! arrays for diagnostics
      real*8 f1_out(bgc_ktotal),
     &       f2_out(bgc_ktotal),
     &       f3_out(bgc_ktotal),
     &       f4_out(bgc_ktotal),
     &       f5_out(bgc_ktotal),
     &       f6_out(bgc_ktotal),
     &       f7_out(bgc_ktotal)

      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out

