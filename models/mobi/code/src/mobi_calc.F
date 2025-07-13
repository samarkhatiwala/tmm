#include "CPP_OPTIONS_CHECK.h"

#if defined O_carbon || defined O_mobi_alk || defined O_mobi_o2 || defined O_mobi
#define BGC
#endif

#if defined BGC || (defined O_PaTh && defined O_PaTh_vflux)
#define NEEDEMP
#endif 

      SUBROUTINE MOBI_CALC(lSize, nx_loc, 
     &     day_loc, relyr_loc, 
     &     TEMP, SALT, TR_surf_glob,
     &     emp_loc, emp_glob, 
#if defined BGC     
     &     winds_loc, aice_loc, hice_loc, hsno_loc,
#endif     
#if defined O_carbon
# if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#  if defined O_carbon_co2_2d
     &     co2atm, 
#  else
     &     co2atm, 
#  endif      
# endif     
# if defined O_c14ccn_data
     &     dc14ccnnhatm, dc14ccneqatm, dc14ccnshatm,
# endif      
# if defined O_c13ccn_data || defined O_carbon_13_coupled
     &     c13o2atm,
# endif
#endif ! O_carbon
#if defined O_mobi
     &     swr_loc,
     &     globdisch, localdisch,
# if defined O_mobi_iron
     &     fe_adep_loc, 
# endif
# if defined O_mobi_silicon
     &     si_dep_loc, si_dep_glob, 
# endif
#endif ! O_mobi
#if defined O_PaTh
     &     dust_adep_loc, 
#endif
#if defined O_carbon
     &     gasexfluxloc, totfluxloc, 
# if defined O_carbon_13_coupled
     &     c13gasexfluxloc,
# endif
#endif
#if defined O_sed
     &     timeToRunSedModel,
     &     globwflx, totlocalwflx,
#endif
     &     debugFlag)

      implicit none

#include "size.h"

      integer lSize, nx_loc
      real day_loc, relyr_loc
      real TEMP(lSize), SALT(lSize), TR_surf_glob(nsrc)
      real emp_loc(nx_loc), emp_glob
#if defined BGC
      real winds_loc(nx_loc), aice_loc(nx_loc), hice_loc(nx_loc), 
     &     hsno_loc(nx_loc)
#endif
#if defined O_carbon
# if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#  if defined O_carbon_co2_2d
      real co2atm(nx_loc)
#  else
      real co2atm
#  endif      
# endif      
# if defined O_c14ccn_data
      real dc14ccnnhatm, dc14ccneqatm, dc14ccnshatm
# endif      
# if defined O_c13ccn_data || defined O_carbon_13_coupled
      real c13o2atm
# endif
#endif ! O_carbon

#if defined O_mobi
	  real swr_loc(nx_loc)
      real globdisch, localdisch(nx_loc)
# if defined O_mobi_iron
      real fe_adep_loc(nx_loc)
# endif
# if defined O_mobi_silicon
      real si_dep_loc(nx_loc), si_dep_glob
# endif
#endif ! O_mobi

#if defined O_PaTh
      real dust_adep_loc(nx_loc)
#endif

#if defined O_carbon
      real gasexfluxloc(nx_loc), totfluxloc(nx_loc)
# if defined O_carbon_13_coupled
      real c13gasexfluxloc(nx_loc)
# endif
#endif
#if defined O_sed
      integer timeToRunSedModel
      real globwflx, totlocalwflx
#endif
      integer debugFlag

      integer i, k, j, ii, ip, n
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW

      integer joff, is, ie, js, je, itr

#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "csbc.h"
#include "grdvar.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "tmngr.h"
#if defined O_mobi || defined O_carbon_14
#include "mobi.h"
#include "calendar.h"
#endif
#if defined O_carbon_14
#include "diaga.h"
      real rrc14std
#endif
#if defined O_ice
#  if defined O_ice_cpts
#include "cpts.h"
#  endif
#include "ice.h"
#endif
#if defined O_mobi
# if defined O_embm
#include "atm.h"
#  if defined O_carbon_13
#include "cembm.h"
#  endif
# endif
#endif
#if defined O_PaTh
#include "protac_thor.h"
#endif
#if defined O_sed
#include "sed.h"
#endif

      real src(imt,km,jsmw:jemw,nsrc)
      common/mobicomm/src

CSPK  Set some parameters
	  dayoyr = day_loc ! relative day number referenced to the beginning of the current year
	  relyr = relyr_loc ! number of years (and fractional years) of model

	  j = 1

#if defined BGC
CSPK Copy T/S from TMM to UVic
CSPK Other tracers should already have been copied by a previous call to S/R mobi_copy_data
      ii=1
      do i=1,nx_loc
        do k=1,kmt(i,j)
	      t(i,k,j,itemp,taum1) = TEMP(ii) ! !degree Celsius
	      t(i,k,j,isalt,taum1) = SALT(ii)       
          ii=ii+1
        enddo
      enddo
#endif

CSPK set SBC arrays
CSPK note: we set tau, taum1 and taup1 to all point to the same location in S/R mobi_ini
	   do n=1,nt
		 if (trsbcindex(n) .ne. 0) then
		   call set_sbc(0, j, j, 1, nx_loc, trsbcindex(n), n, 
     &                  taup1, .false.)
     
		 endif
	   enddo                      

#if defined O_sed
!     S/R gosbc needs the global weathering flux (stored in variable weathflx and 
!     calculated in S/R sed in the online model (but in the offline model that is 
!     just the local flux). Variable weathflx is subsequently overwritten by S/R sed.
      weathflx = globwflx
#endif      
       
CSPK The sequence here is based on the time stepping in UVic_ESCM.F

CSPK copy over winds and other data needed subsequently by S/R gasbc and S/R tracer
#if defined BGC
	  sbc(1:nx_loc,j,iws) = winds_loc(1:nx_loc)
	  aice(1:nx_loc,j,2) = aice_loc(1:nx_loc)
	  hice(1:nx_loc,j,2) = hice_loc(1:nx_loc)
	  hsno(1:nx_loc,j,2) = hsno_loc(1:nx_loc)
#endif

#if defined O_mobi
	  dnswr(1:nx_loc,j) = swr_loc(1:nx_loc)
#endif

#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
CSPK set atmospheric CO2 here as the call to co2ccndata in gasbc is to a stub routine where 
C    it would otherwise be read from file
# if defined O_carbon_co2_2d
        at(1:nx_loc,j,2,ico2) = co2atm(1:nx_loc)
# endif
      co2ccn = co2atm ! we always set this here because gasbc seems to use co2ccn even when O_carbon_co2_2d is defined
#endif

#if defined O_c13ccn_data || defined O_carbon_13_coupled
      c13ccn = c13o2atm
      dc13ccn = 1000.*(c13ccn/(co2ccn-c13ccn)/rc13std - 1.)        
#endif

#if defined O_c14ccn_data
CSPK set atmospheric dc14 here as the call to c14data in gasbc is to a stub routine where 
C    it would otherwise be read from file
      dc14ccnn = dc14ccnnhatm
      dc14ccne = dc14ccneqatm
      dc14ccns = dc14ccnshatm
#endif      

#if defined NEEDEMP
# if !defined O_constant_flux_reference
CSPK set global averages here if necessary as the calls in gasbc are to stub routines
      itr = 0
      do n=3,nt
C        if (itrc(n) .ne. 0) then
         itr = itr + 1
		 gaost(n) = TR_surf_glob(itr)
C        endif
      enddo
# endif
#endif

#if defined O_mtlm && defined O_carbon
CSPK These arrays are first set to zero in S/R glsbc and then updated in S/R mtlm. They 
C    are used in S/R gasbc but as we're not using the land model they should not be 
C    referenced anywhere. Still, we set them to zero anyway to be safe.
	  sbc(1:nx_loc,j,inpp) = 0.0
	  sbc(1:nx_loc,j,isr) = 0.0
	  sbc(1:nx_loc,j,iburn) = 0.0
#endif

!-----------------------------------------------------------------------
!       get the atmospheric S.B.C.
!-----------------------------------------------------------------------

        call gasbc (1, nx_loc, 1, jmt)

#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
CSPK S/R gasbc modifies co2ccn. This should not have any effect if co2emit and global mean 
C    dic flux calculated by a call to stub routine areaavg are exactly zero as we want them 
C    to be. But to be safe, we reset its value again here.
# if defined O_carbon_co2_2d
        at(1:nx_loc,j,2,ico2) = co2atm(1:nx_loc)
# endif
      co2ccn = co2atm ! we always set this here because gasbc seems to use co2ccn even when O_carbon_co2_2d is defined
#endif

CSPK S/R gasbc modifies c13ccn. We reset it here. 
#if defined O_c13ccn_data || defined O_carbon_13_coupled
	  c13ccn = c13o2atm
	  dc13ccn = 1000.*(c13ccn/(co2ccn-c13ccn)/rc13std - 1.)        
#endif

#if defined O_carbon
CSPK S/R gasbc computes the air-sea gas exchange flux in umol cm-2 s-1. We convert it here 
C    to mol m-2 s-1.
      gasexfluxloc(1:nx_loc) = sbc(1:nx_loc,j,idicflx)*1.e-2
#endif

#if defined O_carbon_13_coupled
CSPK S/R gasbc computes the air-sea gas exchange flux in umol cm-2 s-1. We convert it here 
C    to mol m-2 s-1.
      c13gasexfluxloc(1:nx_loc) = sbc(1:nx_loc,j,idic13flx)*1.e-2
#endif      

#if defined O_mobi
CSPK  discharge is used in S/R gosbc and tracer. Units are [g/cm^2/s]
      globaldisch = globdisch
      disch(1:nx_loc,1) = localdisch(1:nx_loc)
# if defined O_mobi_iron
C Note: we do this after call to S/R gasbc because sbc(1,1,idfeadep) is set to zero there (before being read from file)
      sbc(1:nx_loc,j,idfeadep)=fe_adep_loc(1:nx_loc) ! note input is in molFe/m^2/s; it is converted to mmolFe/m^2/s in tracer.F so we don't do it here
# endif

# if defined O_mobi_silicon
C Note: we do this after call to S/R gasbc because sbc(1,1,isilflx) is set to zero there (before being read from file)
      sbc(1:nx_loc,j,isilflx)=si_dep_loc(1:nx_loc)*100.0*sildustfluxfac ! input si flux is otherwise scaled in data.F
      sildustflux=si_dep_glob*100.0*sildustfluxfac
# endif
#endif ! O_mobi

#if defined O_PaTh
C Note: we do this after call to S/R gasbc because sbc(1,1,idustadep) is set to zero there (before being read from file)
      sbc(1:nx_loc,j,idustadep)=dust_adep_loc(1:nx_loc) ! kg/m^2/s
#endif

#if defined NEEDEMP
C To be safe this field is initialized to 0 in mobi_ini
C Salt flux for virtual flux calculation in gosbc; emp_loc and emp_glob should be in cm/s. 
C In UVic, the global freshwater flux at every time step is subtracted from the local 
C freshwater flux at that time step before being applied for the virtual flux calculation. 
C UVic online does that in gosbc; we do that here instead. 
      sbc(1:nx_loc,j,isflx) = (emp_loc(1:nx_loc)-emp_glob)*socn
#endif
      
C#if defined O_mtlm
C!-----------------------------------------------------------------------
C!       get land S.B.C.s
C!-----------------------------------------------------------------------
C
C       call glsbc (1, imt, 1, jmt)
C
C!----------------------------------------------------------------------
C!       call the land-surface and vegetation  model once for each time
C!       step until one segment of "segtim" days is complete.
C!-----------------------------------------------------------------------
C
C       do loop=1,ntspls
C         call mtlm (1, imt, 1, jmt)
C       enddo
C
C#endif

C
C!-----------------------------------------------------------------------
C!       get ocean S.B.C.s
C!-----------------------------------------------------------------------

      call gosbc (1, nx_loc, 1, jmt)

C
C!-----------------------------------------------------------------------
C!       call the ocean model once for each time step until one
C!       segment of "segtim" days is complete. hold ocean S.B.C. fixed
C!       during each segment and predict average S.B.C. for atmos
C!-----------------------------------------------------------------------
C
C       do loop=1,ntspos
C         call mom

CSPK S/R mom calls S/R setvbc, S/R tracer and S/R sed

       js = 1
       je = 1
       is = 1
       ie = nx_loc
       joff = 0 ! just to be safe; for the standard UVic options joff=0

CSPK S/R setvbc copies fluxes calculated and stored in the sbc arrays to the corresponding 
C      stf arrays. These fluxes are converted to and applied as a body force term below.
       call setvbc(joff, js, je, is, ie)
       call tracer(joff, js, je, is, ie)

#if defined O_sed
       if (timeToRunSedModel .eq. 1) then
C		write(*,*)'SPK: calling sed ...'
		call sed (1, nx_loc, 1, jmt)
C		write(*,*)'SPK: done sed'		
# if defined O_sed_weath_diag
        totlocalwflx = weathflx ! this is the local total weathering flux; we compute a global flux in external_forcing
# endif		
       endif
#endif

CSPK S/R sed computes DIC and alkalinity fluxes from the sediments to the ocean, converts 
C    them to body force terms, and stores them in sbc(...,ibdicfx) and sbc(...,ialkflx). 
C    These are applied in a subsequent time step in S/R tracer.
       
CSPK Add surface fluxes (calculated in setvbc.F) as a body force term in surface layer. Be very careful of tracer indexing!
      do j=js,je
        do i=is,ie
		  do n=3,nt
C	  if (itrc(n) .ne. 0) then
		   src(i,1,j,itrc(n)) = src(i,1,j,itrc(n))+stf(i,j,n)*dztr(1)
C	  endif
		  enddo
		enddo
	  enddo	  

#if defined O_carbon
CSPK Virtual flux and weathering flux (if using sediments) is added in S/R gosbc above to the 
C    air-sea gas exchange flux computed in S/R gasbc. Units are umol cm-2 s-1. We convert here 
C    to mol m-2 s-1. Note that 1 umol/cm^3 is numerically equal to 1 mol/m^3.
      totfluxloc(1:nx_loc) = sbc(1:nx_loc,j,idicflx)*1.e-2
#endif
		
      END
      
