C$Header: /Users/ikriest/CVS/mops/CAR_PARAMS.h,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

! RELATED ONLY TO CARBON CHEMISTRY
! This common block passes 
! parameters for air-sea gas exchange from bgc_ini to car_surface_pco2 (CO2PARAMS)
! T and S dependent constants from car_coeffs to car_surface_pco2 (CO2COEFFS)
! Later: 
! T and S dependent constants from car_coeffs to car_carbonate (CO2COEFFS)

      INTEGER car_ktotal
      
      PARAMETER(car_ktotal=100)

      real*8 scar1,scar2,scar3,scar4,phlo,phhi,
     &       convert_mol_to_mmol,rho0,permil,permeg

      COMMON/CO2PARAMS/scar1,scar2,scar3,scar4,phlo,phhi,
     &       convert_mol_to_mmol,rho0,permil,permeg
      
      real*8 ak0,ak1,ak2,ak1p,ak2p,ak3p,
     &       aks,akb,akw,aksi,akf,ff,bt,st,ft,sph,emp,dicgave,alkgave,
     &       pco2atm,co2airseaflux

      COMMON/CO2COEFFS/ak0(car_ktotal),ak1(car_ktotal),ak2(car_ktotal),
     &       ak1p(car_ktotal),ak2p(car_ktotal),ak3p(car_ktotal),
     &       aks(car_ktotal),akb(car_ktotal),akw(car_ktotal),
     &       aksi(car_ktotal),akf(car_ktotal),
     &       ff(car_ktotal),bt(car_ktotal),st(car_ktotal),
     &       ft(car_ktotal),sph,emp,dicgave,alkgave,
     &       pco2atm,co2airseaflux

