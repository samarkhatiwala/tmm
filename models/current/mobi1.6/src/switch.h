!====================== include file "switch.h" ========================

!     all time dependent decisions are made by time manager "tmngr.F"
!     and communicated elsewhere to the model via logical switches.

!     inputs: (defaulted in "blkdta.F", optionally reset via namelist)

!     runlen  = integration period (see rununits). note "runlen" should
!               be an integral number of density time steps. if not,
!               then "runlen" is automatically adjusted to insure this.
!               fractional days are supported but not fractional months
!               or years.
!     rununits= units of "runlen". may be "days", "months", or "years".
!               tmngr will convert "runlen" which is in "rununits"
!               to "rundays" in units of days.

!     segtim  = the integration time "runlen" is broken into a number of
!               segments each of length "segtim" days. updated surface
!               boundary conditions are applied to MOM every "segtim"
!               days. this is useful when coupling to atmospheric models
!               in which case both models exchange surface boundary
!               conditions every "segtim" days where "segtim"
!               is 1/(coupling frequency). without an atmospheric model,
!               when getting surface boundary conditions from data,
!               "segtim" is set to the time step (in days) by mom.F. in
!               either case, "runlen" (in days) should be an integral
!               number of "segtim".

!     nmix    = number of time steps between mixing timesteps. used
!               to damp timestep splitting due to centred leapfrog.

!     init    = (true,false)  indicates that this run is a
!               (start from initial conditions, restart)

!     restrt  = (true,false) = (do,don`t) write a restart at the end
!               of the run

!     eb      = (true,false) configures for the use of a
!               (euler backward,forward) type mixing timestep

!     init_time     = (true,false) sets restarts to initial time

!     init_time_in  = (true,false) sets input restart to initial time

!     init_time_out = (true,false) sets output restart to initial time

!-----------------------------------------------------------------------
!     inputs to tmngr.F: diagnostic intervals
!-----------------------------------------------------------------------

!     note: switches are used to control the interval between doing
!           diagnostics. units for all switches are in days.
!           setting a switch < 0.0 disables whatever the switch is
!           controlling. setting it = 0.0 causes the diagnostic to be
!           done every time step, and setting it > 0.0 causes the
!           diagnostic to be done repeatedly on the specified interval.

!     cmixint = number of days between writing estimated mixing coeffs
!               on faces of T cells and U cells

!     crossint = number of days between writing diapycnal and isopycnal
!               components of flow

!     fctint = number of days between writing difference between
!              FCT and leapfrog advection

!     densityint = number of days between writing density

!     exconvint = number of days between writing temperature rate of
!                 change due to explicit convection

!     glenint =  number of days between global energetics integrals.

!     trmbint =  number of days between momentum and tracer term
!                balances (global and regional).

!     itrmb   = (true,false) = (do,don`t) write regional mask info for
!               the term balance diagnostic. Typically set true
!               at the beginning of a run; otherwise false since it is
!               not necessary to keep writing a time independent field
!               particularly when it may be a significant part of the
!               time dependent part of the diagnostic.

!     gyreint =  number of days between calculation of tracer northward
!                transport.
!     igyre   = (true,false) = (do,don`t) write regional mask info for
!               the gyre diagnostic. Typically set true
!               at the beginning of a run; otherwise false since it is
!               not necessary to keep writing a time independent field
!               particularly when it may be a significant part of the
!               time dependent part of the diagnostic.

!     vmsfint =  number of days between calculation of vertical and
!                meridional stream function.

!     tyzint  =  number of days between calculation of zonally averaged
!                tracer components.

!     prxzint =  number of days between printouts of x-z data.

!     extint  =  number of days between printouts of external mode.

!     dspint  =  number of days between surface pressure calculation.
!                Note: only when "diagnostic_surface_height" is enabled.
!     dspper  = averaging period for "diagnostic_surface_height"

!     tavgint = number of days between regional tracer averages (under
!               horizontal regions).

!     itavg   = (true,false) = (do,don`t) write regional mask info for
!               the tracer average diagnostic. Typically set true
!               at the beginning of a run; otherwise false since it is
!               not necessary to keep writing a time independent field
!               particularly when it may be a significant part of the
!               time dependent part of the diagnostic.

!     tmbint  = number of days over which tracer equation in averaged
!               in depth and longitude to determine the meridional
!               balance among storage, divergence, dissipation and
!               forcing.
!     tmbper  = averaging period for "meridional_tracer_balance"

!     itmb    = (true,false) = (do,don`t) write "msktmb" for tracer
!               the meridional balance diagnostic. Typically set true
!               at the beginning of a run; otherwise false since it is
!               not necessary to keep writing a time independent field
!               particularly when it may be a significant part of the
!               time dependent part of the diagnostic.

!     tsiint  = number of days between printing of time step integrals.
!     tsiper  = averaging period for "time_step_monitor"

!     stabint = number of days between sampling for various stability
!               criteria.

!     timavgint= interval (days) for writing time mean data from
!               the "averaging" grid (only when "time_averages" is
!               enabled). if "timavgint" is not an integral number of
!               density time steps,"timavgint" is automatically adjusted
!               to insure this. if the number of days to integrate is
!               not an integral number of "timavgint" then the last
!               averaging period will be less than "timavgint" days.this
!               may lead to one more averaging period than expected.
!               see "iounit.h" for more details.
!     timavgper= averaging period for "time_averages"

!     xbtint  = averaging period (days) for writing XBT data (only when
!               "xbts" is enabled). if "xbtint" is not an integral
!               number of density time steps, "xbtint" is automatically
!               adjusted to insure this. if the number of days to
!               integrate is not an integral number of "xbtint" then the
!               last averaging period will be less than "xbtint" days.
!               this may lead to one more averaging period than
!               expected. see "iounit.h" for more details.
!     xbtper  = averaging period for "xbts"

!     zmbcint = number of days between calculation of zonal mean
!               surface boundary conditions (and related  quantities)

!     restint = number of days between saving restarts

!     tbtint  = averaging period (days) for writing term balances
!     tbtper  = averaging period for tracer term balances

!-----------------------------------------------------------------------
!     outputs from tmngr.F: logical switches
!-----------------------------------------------------------------------

!     rundays = integration time in days (from "runlen")

!     the following are logical counterparts to the above switches are
!     set within "tmngr" every time step. logical switches control all
!     decisions about when to do things in MOM.

!     cmixts  = (false,true) = (don`t, do) do write estimated mixing
!               coefficients on this time step.
!               based on "cmixint".

!     crossts  = (false,true) = (don`t, do) write diapycnal and
!               isopycnal components of flow on this time step.
!               based on "crossint".

!     fctts    = (false,true) = (don`t, do) write difference between
!               FCT and leapfrog advection on this time step.
!               based on "fctint".

!     densityts  = (false,true) = (don`t, do) write density on this time
!               step. based on "densityint".

!     exconvts  = (false,true) = (don`t, do) do write temperature change
!               due to explicit convection on this time step.
!               based on "exconvint".

!     glents  = (false,true) = (don`t, do) do calculation of global
!               energy integrals on this time step. based on "glenint".

!     trmbts  = (false,true) = (don`t, do) do calculation of momentum &
!               tracer term balance on this timestep. based on "trmbint"

!     gyrets  = (false,true) = (don`t, do) do calculation of tracer
!               northward transport on this timestep. based on "gyreint"

!     vmsfts  = (false,true) = (don`t, do) do calculation of vertical
!               and meridional stream function on this time step.
!               based on "vmsfint"

!     tyzts   = (false,true) = (don`t, do) do calculation of zonally
!               averaged tracer components on this time step.
!               based on "tyzint"

!     prxzts  = (false,true) = (don`t, do) do printouts of x-z data
!               on this time step. based on "prxzint"

!     extts  = (false,true) = (don`t, do) do printout of external mode
!               on this time step. based on "extint"

!     dspts  = (false,true) = (don`t, do) do calculation of diagnostic
!              surface pressure on this time step. based on "dspint"

!     stabts  = (false,true) = (don`t, do) test for stability on this
!               time step. based on "stabint"

!     tavgts  = (false,true) = (don`t do) do tracer averages on this
!               time step. based on "tavgint"

!     tmbts   = (false,true) = (don`t, do) write out tracer meridional .
!               balance on this time step. based on "tmbint"

!     tsits   = (false,true) = (don`t, do) print time step integrals
!               on this time step. based on "tsiint"

!     zmbcts  = (false,true) = (don`t, do) print zonal mean boundary
!               conditions on this time step.  based on "zmbcint"

!     timavgts = (false,true) = (don`t, do) write time mean data
!               on this time step. based on "timavgint"

!     xbtts   = (false,true) = (don`t, do) write averaged XBT data on
!               this time step based on "xbtint"

!     restts  = (false,true) = (don`t, do) save a restart on this time
!               step based on "restint"

!     tbtts   = (false,true) = (don`t, do) write averaged tracer term
!               balance data on this time step based on "tbint"

!     leapfrog= (false,true) on a (mixing, normal leapfrog) time step
!                based on "nmix"

!     euler1  = true on the 1st pass of an euler backward time step
!               otherwise false. (applies when "eb" = true)
!     euler2  = true on the 2nd pass of an euler backward time step
!               otherwise false. (applies when "eb" = true)
!     forward = true on a forward time step. otherwise false
!                (applies when "eb" = false)

!     the following logical switches are based on the model time step.

!     first   = (true,false) =  when it`s (the first, not the first)
!                               time step of a run
!     eots    = end of a time step. always true except for first
!               pass of an euler backward time step
!     eorun   = last time step of a run. always false except during the
!               last time step of the run.

!     eoday   = true when within 1/2 time step of the end of a day
!               else ... false
!     eoweek  = true when within 1/2 time step of the end of a 7 day
!               week (referenced to the start of a year) else ...false
!     eo2wks  = true when within 1/2 time step of the end of two weeks
!               (referenced to the start of a year) else ... false
!     midmon  = true when within 1/2 time step of the middle of a month
!               else ... false
!     eomon   = true when within 1/2 time step of the end of a month
!               else ... false
!     eoyear  = true when within 1/2 time step of the end of a year
!               else ... false
!     osegs   = true on the 1st time step of an ocean segment in mom.F
!               otherwise false.
!     osege  =  true on the last time step of an ocean segment in mom.F
!               otherwise false.

      character(len=8) :: rununits
      common /switc_c/ rununits

      integer nmix, ieoday,ieoweek,ieo2wks
      integer ieomon,imidmon,ieoyear,ieorun
      common /switc_i/ nmix, ieoday,ieoweek,ieo2wks
      common /switc_i/ ieomon,imidmon,ieoyear,ieorun

      logical eb, leapfrog, euler1, euler2, forward, eots
      logical init, first, restrt, itavg, itmb, itrmb, igyre
      logical eoday, eoweek, eo2wks, eomon, midmon, eoyear, eorun
      common /switc_l/ eb, leapfrog, euler1, euler2, forward, eots
      common /switc_l/ init, first, restrt
      common /switc_l/ itavg, itmb, itrmb, igyre
      common /switc_l/  eoday, eoweek, eo2wks
      common /switc_l/  eomon, midmon, eoyear, eorun
      logical init_time, init_time_in, init_time_out
      common /switc_l/  init_time, init_time_in, init_time_out

      real runlen, rundays
      common /switc_r/ runlen, rundays

!-----------------------------------------------------------------------

!     S W I T C H E S    B A S E D    O N    A N    I N T E R V A L

!     each interval switch needs three variables in common. The
!     following naming convention is used.

!         1) an interval (real) for diagnostic output (e.g,.  glenint)
!         2) a switch (logical) for the interval (e.g.,  glents )

!     the third is an internal variable needed by the time manager
!     to support calculation of the logical switch

!         3) an index (integer)                       (e.g., iglenint)

!     the user must specify the interval [e.g., glenint] for diagnostic
!     output in units of days. tmngr sets the corresponding logical
!     switch [e.g., glents] every time step. It is set to true when
!     within half a time step of the requested interval, otherwise it is
!     false. All decisions relating to the interval [e.g., glenint]
!     are based on the logical switch [e.g., glents].

!     internal time structures

!     The switch index [e.g., iglenint] is used to subsrcipt into
!     internal arrays maintained by tmngr.F. The switch index is
!     allocated on the first call to function "alarm".
!     The array entry [e.g., iinterval(iglenint)] is a time index to the
!     internal representation of the interval [e.g., glenint].
!     The array entry [e.g., ialarm(iglenint)] is a time index to the
!     next time the alarm will be true.
!-----------------------------------------------------------------------

      integer          itavgint,  iglenint,  itrmbint, iprxzint
      common /switc_i/ itavgint,  iglenint,  itrmbint, iprxzint
      logical          tavgts,    glents,    trmbts,   prxzts
      common /switc_l/ tavgts,    glents,    trmbts,   prxzts
      real             tavgint,   glenint,   trmbint,  prxzint
      common /switc_r/ tavgint,   glenint,   trmbint,  prxzint

      integer          iextint,  iexconvint, icmixint
      common /switc_i/ iextint,  iexconvint, icmixint
      logical          extts,    exconvts,   cmixts
      common /switc_l/ extts,    exconvts,   cmixts
      real             extint,   exconvint,  cmixint
      common /switc_r/ extint,   exconvint,  cmixint

      integer          ivmsfint, igyreint, ityzint, ifctint
      common /switc_i/ ivmsfint, igyreint, ityzint, ifctint
      logical          vmsfts,   gyrets,   tyzts,   fctts
      common /switc_l/ vmsfts,   gyrets,   tyzts,   fctts
      real             vmsfint,  gyreint,  tyzint,  fctint
      common /switc_r/ vmsfint,  gyreint,  tyzint,  fctint

      integer          istabint, izmbcint, icrossint, idensityint
      common /switc_i/ istabint, izmbcint, icrossint, idensityint
      logical          stabts,   zmbcts,   crossts,   densityts
      common /switc_l/ stabts,   zmbcts,   crossts,   densityts
      real             stabint,  zmbcint,  crossint,  densityint
      common /switc_r/ stabint,  zmbcint,  crossint,  densityint

      integer          iosegs, iosege
      common /switc_i/ iosegs, iosege
      logical          osegs,  osege
      common /switc_l/ osegs,  osege
      real             segtim
      common /switc_r/ segtim

      integer          irestint
      common /switc_i/ irestint
      logical          restts
      common /switcl/  restts
      real             restint
      common /switcr/  restint

!-----------------------------------------------------------------------

!     S W I T C H E S    B A S E D    O N    A N    I N T E R V A L

!              A N D   A V E R A G I N G   P E R I O D

!     each averaging period switch needs five variables in common. The
!     following naming convention is used.

!         1) an interval (real) for diagnostic output    (e.g. xbtint  )
!         2) a switch (logical) for the interval         (e.g. xbtts   )
!         3) an averaging period (real)                  (e.g. xbtper  )
!         4) a switch (logical) for accumulating         (e.g. xbtperts)

!     the third is an internal variable needed by the time manager
!     to support calculation of the logical switches

!         5) an index (integer)                         (e.g. ixbtint  )

!     The user must specify the interval [e.g., xbtint] for diagnostic
!     output in units of days and the averaging period [e.g., xbtper]
!     in units of days. The averaging period may be less than or equal
!     to the interval. For example, if the interval is 30.0 days and the
!     averaging period is 5.0 days, results will be averaged over all
!     time steps within days 26, 27, 28, 29, and 30.  An averaging period
!     of 0.0 days averages over the last time step of the interval (as
!     does xbtper = dt), and an averaging period less than zero turns
!     the switches off for all time steps.

!     The logical switch for writing output at the specified interval
!     [e.g., xbtts] is set to true on the last time step of the
!     averaging period. The logical switch for accumulating results
!     [e.g., xbtperts] is true for all time steps within the averaging
!     period, otherwise it is false.

!     internal time structures

!     The index [e.g., ixbtint] is allocated on the first call to
!     function "avg_alarm". The array element iperiod(ixbtint) is an
!     index to the time structure for the internal representation of
!     "xbtper", and ilastsw(ixbtint) is the index of the switch that
!     flags the last time step of the accumulation period.
!     Depending on use,  ilastsw(ixbtint) may either be the index
!     of another "named" switch or the index of a new switch
!     allocated on the first time step.
!     In the latter case, iinterval(ilastsw(ixbtint)) is the index of
!     the time structure where "xbtint" is stored in internal form,
!     and ialarm(ilastsw(ixbtint)) is the index of the time when an
!     accumulation period will next end.
!     The variable nextts(ixbtint) is true whenever the next
!     time step will begin the accumulation period.

!-----------------------------------------------------------------------

      integer          ixbtint,   idspint,  itmbint,  itimavgint
      common /switc_i/ ixbtint,   idspint,  itmbint,  itimavgint

      logical          xbtts,     dspts,    tmbts,    timavgts
      logical          xbtperts,  dspperts, tmbperts, timavgperts
      common /switc_l/ xbtts,     dspts,    tmbts,    timavgts
      common /switc_l/ xbtperts,  dspperts, tmbperts, timavgperts

      real             xbtint,    dspint,   tmbint,   timavgint
      real             xbtper,    dspper,   tmbper,   timavgper
      common /switc_r/ xbtint,    dspint,   tmbint,   timavgint
      common /switc_r/ xbtper,    dspper,   tmbper,   timavgper

      integer          itbtint,   itsiint
      common /switc_i/ itbtint,   itsiint

      logical          tbtts,     tsits
      logical          tbtperts,  tsiperts
      common /switc_l/ tbtts,     tsits
      common /switc_l/ tbtperts,  tsiperts

      real             tbtint,    tsiint
      real             tbtper,    tsiper
      common /switc_r/ tbtint,    tsiint
      common /switc_r/ tbtper,    tsiper

!-----------------------------------------------------------------------

!                 S W I T C H E S    B A S E D    O N

!         C A L E N D A R   O R    P R E V I O U S    S W I T C H

!               A N D   A V E R A G I N G    P E R I O D

!     the following logical switches are based on any calendar or
!     interval switch and an averaging period (in days). The  averaging
!     period must be less than or equal to the interval. The last
!     time step of the averaging period is at the end of the interval.
!     If the averaging period is set to zero, the averaging period
!     consists only of the last time period of the interval.  If
!     the averaging period is less than zero, these switches are always
!     false.

!     each averaging period switch needs four variables in common. For
!     example, if the averaging period is before the end of each month
!     then the calendar switch (eomon), and index (ieomon) are presumed
!     to exist in common and need not be added.

!     Additionally, four items are needed.

!       1) an averaging period (real)                  (e.g.  testper  )
!       2) a switch (logical) for accumulating results (e.g.  testperts)
!       3) a switch (logical) for the end of interval  (e.g.  testts   )

!     the fourth is an internal variable needed by the time manager
!     to support calculation of the logical switch

!       4) an index (integer)                          (e.g.  itestper )

!     Suppose it is required to produce averages over all time steps
!     during the last 5 days of each month. Then "testper" = 5.0 and
!     the following will calculate the accumulating switch.

!      testts = avg_alarm(itestper, ihalfstep, 0, testper, iref, ieomon)
!      testperts = on(itestper)

!     Note the use of "ieomon" to key off the months.  The switch
!     "testts" will be true whenever "eomon" is true.
!     Also note that when an averaging switch is keyed off another
!     switch, the switch interval argument is not used, but is
!     retained for consistency with the form of other averaging
!     switches.
!-----------------------------------------------------------------------

      integer maxsw
      parameter (maxsw=100)

      integer itestper, nsw, ialarm, iinterval, iperiod, ilastsw
      common /switc_i/ itestper, nsw, ialarm(maxsw), iinterval(maxsw)
      common /switc_i/ iperiod(maxsw), ilastsw(maxsw)

      logical testts, testperts, on, lastts, nextts
      common /switc_l/ testts, testperts, on(maxsw), lastts(maxsw)
      common /switc_l/ nextts(maxsw)

      real testint, testper
      common /switc_r/ testint, testper
