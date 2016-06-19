!======================= include file "tmngr.h" ========================

!                       time manager variables

!-----------------------------------------------------------------------
!     time manager inputs:
!-----------------------------------------------------------------------

!     how to choose a reference time:

!     refrun  = (true,false) to base calculation for diagnostic switches
!              on (the start of each job, other reference time)
!              example:
!              suppose each job submission integrates
!              for one month but the number of days per month changes.
!              setting "refrun" = true and setting
!              "timavgint" = (days in month)/3 will give 3 averaging
!              periods per month of approximately 10 days each. the
!              only restriction is that "timavgint"is an integral number
!              of time steps (if not then "timavgint" is reset to insure
!              this condition. other diagnostic switches do not have
!              this restriction).

!     refinit = (true, false) for basing calculation of logical switches
!              on (initial conditions, other reference time)
!              example: if term balances are desired every 20 days
!              (trmbint=20.0) and refinit = true, then they
!              will be done every 20 days starting from initial
!              condition time.

!     refuser = (true, false) to base calculations of logical switches
!              on (user-chosen reference time, other reference time)
!              if refuser = true, the user must also supply values for
!              ryear, rmonth, rday, rhour, rmin, rsec (integer)
!              example: if term balances are desired every 20 days
!              (trmbint=20.0) and refuser = true, then they will be done
!              every 20 days counting from reference time, ignoring the
!              initial condition time. for comparing diagnostics from
!              various experiments with different initial condition
!              times, refuser = true will be more appropriate. setting
!              refuser = true and choosing the reference time to be
!              the initial condition time is the same as refinit = true.

!     summary of how to choose the time for referencing calculations
!     of logical switches

!     refrun  = T ==>  referenced to the start of each run
!     refinit = T ==>  referenced to initial condition time given by:
!                     year0, month0, day0, hour0, min0, sec0
!     refuser = T ==>  referenced to user specified reference time so
!                     must set: ryear, rmonth, rday, rhour, rmin, rsec

!-----------------------------------------------------------------------

!     time variable arrays

!     arrays "iday" and "msday" contain the primary internal
!     representation of all times within the time manager. they are
!     referenced by using a subscript to indicate which time.

!     iday    = integer days (since Dec 31, 1899 when specifying a date)
!     msday   = non-negative integer milliseconds after midnight

!     it is desirable to have time information expanded to include the
!     following secondary time fields:

!     year       =
!     month      =
!     day        =
!     hour       =
!     minute     =
!     second     =
!     tstamp     = 32 character date and time stamp m/d/y h:m:s
!     dayofyear  = integer day of the year (1..yrlen)
!     dayofweek  = 1=sun - 7=sat
!     daysinmon  = days in the month
!     daysinyear = days in the year

!     those times for which primary and secondary information is
!     maintained by the time manager are called "full times". those for
!     which only primary information is kept are called "short times"

!     indices to  "full times" (including year, month ,day, etc).

!     itime     = simulation time corresponding to "itt"
!     initial   = time of the initial conditions
!     irunstart = time of the start of the run
!     iuser     = user defined reference time
!     iref      = one of the three above selected by logicals
!                 (refinit, refrun, refuser)

!     indices to  "short times". ("iday", "msday" only)

!     isunday    = time of a sunday for week and two week switches
!     ihalfstep  = dt/2 beyond itime
!     imodeltime = time since initial conditions
!     iruntime   = time since run start
!     iusertime  = time since user specified reference time
!     idt        = integer days and milliseconds of dt
!     idtd2      = integer days and milliseconds of dt/2

!     ireftime   = time used locally in alarm function

!     for any time index (short or full) the internal representation
!     may be converted to either real days or real seconds using
!     the functions:
!                  realdays(index)
!                  realsecs(index)

!     dayoyr   = relative day number referenced to the beginning
!               of the current year.  (real)
!     relyr    = number of years (and fractional years) of model
!               integration (for time tau+1 {itt}) relative to
!               initial condition
!     prelyr   = relyr for previous time step

!     stamp    = 32 character date and time for current model timestep
!     pstamp   = 32 character date and time for previous model timestep
!     calendar = calendar name

!     itt      = current time step counter (from initial cond.)
!     itt0     = time step at start of current run

!               variables used for initialization

!     irstdy   = integer number of days at start of run
!     msrsdy   = fractional day in millisec at start of run

!     year0    = year of initial conditions
!     month0   = month of initial conditions
!     day0     = day of initial conditions
!     hour0    = hour of initial conditions
!     min0     = minute of initial conditions
!     sec0     = second of initial conditions

!     ryear    = year of user specified reference time
!     rmonth   = month of user specified reference time
!     rday     = day of user specified reference time
!     rhour    = hour of user specified reference time
!     rmin     = minute of user specified reference time
!     rsec     = second of user specified reference time

!-----------------------------------------------------------------------

      integer ntimes, nfulltimes
      parameter (ntimes = 100, nfulltimes = 30)

      character(32) :: tstamp, stamp, pstamp
      character(120) :: calendar
      common /tmngrc/ tstamp(nfulltimes), stamp, pstamp, calendar

      integer nextfulltime, nexttime
      integer initial, iref, irunstart, itime, iuser
      integer iruntime, imodeltime, ireftime, iusertime
      integer ihalfstep, isunday
      integer itemptime,itemptime2,itmptime,itmptime2,itmptime3
      integer idt, idtd2
      integer iday, msday, year, month, day, hour, minute, second
      integer dayofyear, dayofweek, daysinmon, daysinyear
      integer itt0, itt, irstdy, msrsdy
      integer year0, month0, day0, hour0, min0, sec0
      integer ryear, rmonth, rday, rhour, rmin, rsec

      common /tmngri_i/ nextfulltime, nexttime
      common /tmngri_i/ initial, iref, irunstart, itime, iuser
      common /tmngri_i/ iruntime, imodeltime, ireftime, iusertime
      common /tmngri_i/ ihalfstep, isunday
      common /tmngri_i/ itemptime, itemptime2, itmptime, itmptime2
      common /tmngri_i/ itmptime3
      common /tmngri_i/ idt, idtd2
      common /tmngri_i/ iday(ntimes), msday(ntimes)
      common /tmngri_i/ year(nfulltimes), month(nfulltimes)
      common /tmngri_i/ day(nfulltimes), hour(nfulltimes)
      common /tmngri_i/ minute(nfulltimes), second(nfulltimes)
      common /tmngri_i/ dayofyear(nfulltimes), dayofweek(nfulltimes)
      common /tmngri_i/ daysinmon(nfulltimes), daysinyear(nfulltimes)
      common /tmngri_i/ itt0, itt, irstdy, msrsdy
      common /tmngri_i/ year0, month0, day0, hour0, min0, sec0
      common /tmngri_i/ ryear, rmonth, rday, rhour, rmin, rsec

      logical refrun, refinit, refuser
      common /tmngr_l/ refrun, refinit, refuser

      real dayoyr, relyr, prelyr
      common /tmngr_r/ dayoyr, relyr, prelyr

!-----------------------------------------------------------------------
!     acceleration parameters for long term forcing
!-----------------------------------------------------------------------
!     an acceleration factor can be used to accelerate long term
!     forcing. forcing and output at time intervals of less than one
!     year are not affected (eg. seasonal cycle of insolation).
!     averaging intervals are also not changed. only longer term
!     forcings are accelerated (eg. CO2 or orbital changes).

!     runlen, runstep and most output intervals will will be divided by
!     this factor (accel). output intervals are divided by accel if the
!     interval divided by accel is larger or equal to a year. both
!     runlen and runstep divided accel must still be evenly divisible
!     by the coupling time (segtim). output intervals divided by accel
!     should still be evenly divisible by the coupling time and the
!     averaging period.

!     times in various output files may appear to be inconsistent.
!     for example the year may switch part way through a seasonal cycle.
!     be cautious when looking at output averaged over less than a year.
!     it is probably best to use an integer value for accel that will
!     leave the runlen and runstep as integer multiples of a year. the
!     initial year can be set to offset any existing time in a restart.
!     choosing an inappropriate acceleration factor may mess up output
!     intervals and averaging periods. think before you use.

!     accel     = factor for accelerating forcing transient data
!     accel_yr0 = relyr to start acceleration

      real accel, accel_yr0
      common /tmngr_r/ accel, accel_yr0
