C$Header: /Users/ikriest/CVS/mops/BGC_CONTROL.h,v 1.1.1.1 2015/06/03 17:02:09 ikriest Exp $
C$Name: mops-1_2 $

C EVERYTHING RELATED TO TIMESTEPPING AND GEOMETRY

c the number of biogeochemical time steps per ocean time step 
      integer bgc_timesteps,bgc_kmax,bgc_keuph

c the time step length of biogeochemistry in days. 
      real*8 bgc_dt

      COMMON/BGCICONTROL/bgc_timesteps,bgc_kmax,bgc_keuph
      COMMON/BGCCONTROL/bgc_dt

      


            
