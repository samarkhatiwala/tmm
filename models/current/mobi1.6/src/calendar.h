!====================== include file "calendar.h" ======================

!                       calendar specification arrays

!-----------------------------------------------------------------------
!     eqyear  = true to select a calendar in which each year
!               has the same number of days (i.e., no leap years)
!               false selects a Julian calendar

!     eqmon   = true to force all months to have the same number of days
!               false => the usual 31, 28, 31, 30, ..., days per month.
!               only used when eqyear = true

!     dayname = character names of days

!     monname = character names of months

!     monlen  = the length of each month (in days) when eqmon is true

!     yrlen   = the length of a typical (non-leap) year in days

!     daypm   = array of month lengths in days   (non-leap)

!     msum    = array of cumulative days preceding each month
!               (again, non-leap)

!     daylen  = the length of a day in seconds
!-----------------------------------------------------------------------

      character(10) :: dayname
      character(12) :: monname
      common /calen_c/ dayname(7), monname(12)

      logical eqyear, eqmon
      common /calen_l/ eqyear, eqmon

      integer daypm, msum, yrlen, monlen
      common /calen_i/ daypm(12), msum(12), yrlen, monlen

      real daylen
      common /calen_r/ daylen
