// $Header: /u/gcmpack/MITgcm/pkg/dic/DIC_OPTIONS.h,v 1.12 2014/12/05 01:43:40 jmc Exp $
// $Name:  $

#ifndef DIC_OPTIONS_H
#define DIC_OPTIONS_H
#include "PACKAGES_CONFIG.h"
//#include "CPP_OPTIONS.h"

#ifdef ALLOW_DIC
//     Package-specific Options & Macros go here

#define DIC_BIOTIC
#define ALLOW_O2
#define ALLOW_FE
#undef READ_PAR
#define MINFE
#define DIC_NO_NEG
// these all need to be defined for coupling to atmospheric model:
#undef USE_QSW
#undef USE_QSW_UNDERICE
#undef USE_ATMOSCO2
#undef USE_PLOAD

// use surface salinity forcing (scaled by mean surf value) for DIC & ALK forcing
#undef ALLOW_OLD_VIRTUALFLUX

// put back bugs related to Water-Vapour in carbonate chemistry & air-sea fluxes
#undef WATERVAP_BUG

// dissolution only below saturation horizon following method by Karsten Friis
#undef CAR_DISS

// Include self-shading effect by phytoplankton
#undef LIGHT_CHL
// Include iron sediment source using DOP flux
#define SEDFE

#endif /* ALLOW_DIC */
#endif /* DIC_OPTIONS_H */

//EH3 ;;; Local Variables: ***
//EH3 ;;; mode:fortran ***
//EH3 ;;; End: ***
