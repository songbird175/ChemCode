      program smog

c  program smog
c
c  simulation of photochemical reactions of vocs, nox and hcho
c  and the formation of ozone
c
c  matt fraser, scott hersey
c  april 2002, feb. 2016
c
c  coded in fortran (f77) 
c
c
      implicit none
      real p_no2(24), p_rcho(24)  ! photolysis rates
      real dep(10)                ! deposition rates
      real E(10)                  ! emission rates
      real vent                   ! ventilation rate
      real rk(12)                 ! chemical reaction rate constants
      real o2,m                   ! contant value for o2,m
      real c(10)                  ! concentration of accumulation species
      integer iday                ! day of simulation
      real etime, tinc            ! time markers 
      real treport, trinc         ! time markers
      integer j1, j2              ! interpolation indexes
      real o, oh                  ! pssa values for o and oh radicals
      integer n                   ! number of accumulation species
      integer i                   ! index
      real tin, tout              ! times for integration




c  ARRAY ASSIGNMENTS FOR SPECIES
c  1 = no2
c  2 = no
c  3 = o3
c  4 = rh
c  5 = rcho
c  6 = hno3
c  7 = pan
c  8 = rcoo2
c  9 = ho2
c 10 = ro2
c  o and oh determined by pssa

c  photolysis rates per minute

      data p_no2 / 0.0,0.0,0.0,0.0,0.010,0.136,
     *             0.292,0.392,0.464,0.503,0.523,0.529,
     *             0.519,0.493,0.450,0.379,0.260,0.104,
     *             0.0,0.0,0.0,0.0,0.0,0.0 /

      data p_rcho / 0.0,0.0,0.0,0.0,0.0,0.21e-3,
     *              0.67e-3,1.20e-3,1.64e-3,1.96e-3,
     *              2.14e-3,2.20e-3,2.13e-3,1.92e-3,
     *              1.59e-3,1.12e-3,0.60e-3,0.16e-3,
     *              0.0,0.0,0.0,0.0,0.0,0.0 /
 
c reaction rate constants (units of ppm min)
      data rk / 0.0,2.183e-5,26.59,3.775e3,2.341e4,0.0,
     *          1.214e4,1.127e4,3.8e3,1.613e4,2.07e3,2.143e-2 /

c concentrations of species that are constant in ppm
      o2 = 2.1e5
      m  = 1.0e6

c initial concentrations of species in ppm (arbitrary)
      data c / 0.010,0.010,0.010,0.050,0.010,
     *         0.001,0.001,0.0,0.0,0.0 /

c deposition rates dependant on deposition velocity
c and thus mixing height and varies by how sticky compound is
c
c no2 = 0.42 cm/s = 1.1% per hour
c o3 = 2.5 cm/s = 6.4% per hour 
c rcho = 0.42 cm/s = 1.1% per hour
c hno3 = 2.5 cm/s = 6.4% per hour
c pan = 1.7 cm/s = 4.4% per hour

c data for deposition given in fraction of concentration per minute

      data dep / 0.18e-3,0.0,1.1e-3,0.0,0.18e-3,1.1e-3,0.73e-3,0.0,
     *           0.0,0.0 / 

c emission rates given in ppm/min
c base case emissions

      data E / 6.2e-6,55.8e-6,0.0,125.0e-6,3.5e-6,0.0,
     *        0.0,0.0,0.0,0.0 /

c loss rate through ventilation in fraction per minute

      data vent / 0.0007 /

c number of accumulation species 

      data n / 10 /

c time factors in hours
      iday = 1
      etime = 0.0
      tinc = 0.1
      treport = 0.3
      trinc = 0.3

c
c output file = smog.out
c
      open (unit=1,file='smog.out',status='new')
      write (1,99)
 99   format('time(h)    NO2       NO        O3        RH        '
     *       'RCHO      HNO3      PAN       RCOO2     HO2       '               
     *       'RO2       O         OH        j(NO2)    j(RCHO)')

 10   continue

c set photolysis rates by interpolation

      j1 = aint(etime)
      j2 = j1+1
      if (j1 .eq. 0) then
         rk(1) = 0.0
         rk(6) = 0.0
      else
         rk(1) = (etime-1.0*j1)*p_no2(j2) 
     *           + (1.0*j2-etime)*p_no2(j1)
         rk(6) = (etime-1.0*j1)*p_rcho(j2) 
     *           + (1.0*j2-etime)*p_rcho(j1)
      endif

c calculate a step for solution
      tin = 60.0*etime     ! convert from hours (etime) to min (tin)
      tout = 60.0*(etime + tinc)
      call hybrid(n,c,rk,tin,tout,E,dep,vent,etime,o,oh,m,o2)
      etime = tout/60.0

c update time and write to output if needed
      if (etime .ge. treport) then
         treport = treport + trinc
         write (1,98) etime,(c(i),i=1,10),o,oh,rk(1),rk(6)
      endif

      if (etime .ge. 24.0) then
         write (*,97) iday
         write (1,97) iday
         iday = iday + 1
         etime = etime - 24.0
         treport = treport - 24.0
      endif

      if (iday .ge. 8) then
          write (1,96) (E(i),i=1,10)
          stop
      endif

      goto 10  
     
 98   format (1x,f4.1,5x,14(1x,0pe9.3))
 97   format (/'end day',i2/)
 96   format (/'Emissions:',10(1x,0pe9.3))
  
      end


      subroutine hybrid(n,c,rk,tin,tout,E,dep,vent,etime,o,oh,m,o2)

c subroutine for integration of stiff set of coupled differential equations
c using a hybrid integration technique
c
c the general for of the equations is assumed to be:
c 
c   dc/dt = g(c) = q(c) - l(c)*c
c
c where q(c) is formation and emission and l(c) is first order loss
c
c must call a subroutine diffr(c,rk,q,l,g) to provide the g, q, and l 
c functions from the concentrations (c) and the rate constants (rk)
c time step over that interval (tin to tout) should be small so that the
c rate constants are constant (i.e. photolysis does not change).
c
c in calling hybrid, n = number of equations
c                    c = vector of concentrations
c                        at time tin on input
c                        at time tout on output
c                    rk = rate constants
c                    tin = time at beginning of integration
c                    tout = time at end of integration

      implicit none

c declare time values
      real tin, tout, tnow, dtovr2
      real dtmin, xxx, dt
  
c declare work variables
      real y2, y3, tau2, sumtau
   
c declare integration control parameters
      real e1,e2,e3,e4,e5,e6

c declare concentration variables
      real c(10),c2(10),c3(10)

c rate constants
      real rk(12)

c formation rate (rxn + emission) for each species
      real q(10), q2(10)

c sink (first order rate constants) for each species
      real l(10), l2(10)

c net formation (g = q - l*c)
      real g(10), g2(10)

c common blocks
      real o,oh,m,o2
      real E(10),dep(10),vent,etime

c flags for stiff species, fast convergence and convergence
      logical stiff(10), fascon, convrg
      integer ifast, itmax, j, k, n

c *time constants to set step size*  DO NOT CHANGE
c e1 = scaling factor to determine initial step size
c e2 = determines if individual ode is stiff or normal
c      (if normal -> EULER'S METHOD, if stiff -> ASYMPTOTIC)
c e3 = convergence criteria
c e4 = if convergence is not achieved after a few steps, reduce
c      the step size by this factor
c e5 = if convergece is reached rapidly, increase the step size by
c      10% using this factor
c e6 = for species with low concentrations, automatically assume 
c      convergence.  If concentration below e6 -> assume convergence
c
c dtmin = minimum time step allowed before the integration halts with error
c ifast = number of corrections for fast convergence (if .le. ifast then
c         consider fast convergence
c itmax = if integration not completed by itmax, then reduce time step and
c         try again
 
      data  e1,  e2,    e3,  e4,  e5,     e6,  dtmin,ifast,itmax 
     * / 0.001, 1.0, 0.001, 0.5, 1.1, 1.0e-9, 1.0e-5,    2, 4 / 

c determine q, l, and g by subroutine
      call diffr(c,rk,q,l,g,E,dep,vent,etime,o,oh,m,o2)

c set the time step with a first estimate
      dt = 1.0
      xxx = 1.0
      do 10 j= 1, n
        if(g(j) .ne. 0.0) xxx=abs(c(j)/g(j))
        if(xxx .gt. 0.0) dt=min(dt,xxx)
 10   continue

      dt = e1*dt
      if (dt .lt. dtmin) dt=dtmin
      tnow = tin

c check step size will not put integration past tout
 20   if (tnow+dt .gt. tout) dt = tout - tnow    
       dtovr2 = dt*0.5

c check for stiffness
      do 30 j = 1, n
         stiff(j) = (dt*l(j)) .gt. e2
 30   continue

c predict concentration at t+dt
      do 40 j = 1,n
         if(stiff(j)) then
          tau2 = 2.0/l(j)
           c2(j) = (c(j)*(tau2-dt) + tau2*dt*q(j))/(tau2+dt)
         else
           c2(j) = c(j) + g(j)*dt
         endif

c reset any negative concentrations
      if (c2(j) .lt. 0.0) c2(j) = 0.0
 40   continue

c iteratively correct the prediction for t+dt
      k = 1
      fascon = .true.

c compute q, l and g using lastest concentration at t+dt
 50   call diffr(c2,rk,q2,l2,g2,E,dep,vent,etime,o,oh,m,o2)

c corrected predictions for c(t+dt)
      do 60 j =1,n
         if (stiff(j)) then
            sumtau = 1.0/l(j) + 1.0/l2(j)
            c3(j) = (dtovr2*sumtau*(q(j)+q2(j)) + c(j)*(sumtau-dt))
     *              / (sumtau+dt)
          else
            c3(j)=c(j) + (g(j)+g2(j))*dtovr2
          endif

c reset any negative concentrations
      if (c3(j) .lt. 0.0) c3(j) = 0.0
 60   continue

c test for convergence
      convrg = .true.
      do 70 j=1,n
         y2=c2(j)
         y3=c3(j)
         if (y3 .gt. e6 .and. abs(y3-y2) .gt. e3*min(y2,y3)) then
            convrg = .false.
            goto 80
        endif
 70   continue

 80   if (.not. convrg) then

c if it has converged, assign concentration values for next step
c if too many iterations, decrease step size
        
        if (k .lt. itmax) then
            do 90 j = 1, n
            c2(j) = c3(j)
 90   continue        
            fascon = k .lt. ifast
            k = k + 1
            goto 50
         else
c if iteration did not converge in itmax steps, reduce step size

        dt = e4*dt
        if (dt .ge. dtmin) go to 20
        write (*,92)
        stop
 92     format(/' stepsize below minimum allowed in hybrid solution')
        endif
       
      else
       
c load c values at time d+dt into c array and assign q, l, g values
c for next step
       
      do 94 j= 1, n
       c(j) = c3(j)
 94   continue

      tnow = tnow + dt

c if iteration converged in less than ifast tries, increase step size

      if (fascon) dt = e5*dt

      if (tnow .ge. tout) then
         return
      else
         do 98 j = 1, n
          q(j) = q2(j)
          l(j) = l2(j)
          g(j) = g2(j)
 98   continue
      go to 20
         endif
      endif
        
      end



      subroutine diffr(c,rk,q,l,g,E,dep,vent,etime,o,oh,m,o2)

c evaluation of the generation (q), loss (l) and net growth (g)
c where g = dc/dt for each of the n species given in concentration
c array c.  rk is the array of reaction rates, and emissions, ventilation, 
c and deposition given in this subroutine

      implicit none
      real dadj
      integer i
      real c(10)
      real rk(12)
      real g(10)
      real l(10)
      real q(10)
      real o,oh,o2,m
      real E(10),dep(10),vent,etime

c use pssa to find radical concentrations
       
      o = rk(1)*c(1)/(rk(2)*m*o2)
      oh = rk(7)*c(9)*c(2)/(rk(4)*c(4)+rk(5)*c(5)+rk(10)*c(1))

c determine source rates for each species

      q(1) = E(1)+rk(3)*c(2)*c(3) + rk(7)*c(9)*C(2) +
     *       rk(8)*c(10)*c(2) + rk(9)*c(8)*c(2) + 
     *       rk(12)*c(7)
      q(2) = E(2) + rk(1)*c(1)
      q(3) = rk(2)*o*o2*m
      q(4) = E(4)
      q(5) = E(5)+rk(8)*c(10)*c(2)
      q(6) = rk(10)*oh*c(1)
      q(7) = rk(11)*c(8)*c(1)
      q(8) = rk(5)*c(5)*oh + rk(12)*c(7)
      q(9) = rk(6)*c(5) + rk(8)*c(10)*c(2)
      q(10)= rk(4)*c(4)*oh + rk(6)*c(5) + rk(9)*c(8)*C(2)

c at nighttime decrease deposition under stagnant conditions
c  2000 hours to 0800 hours
        
      if (etime .lt. 8.0 .or. etime .gt. 20.0) then
         dadj = 0.1
      else
         dadj = 1.0
      endif

c determine loss rate for each species

      l(1) = rk(1) + rk(10)*oh + rk(11)*c(8) + 
     *       vent + dadj*dep(1)
      l(2) = rk(3)*c(3) + rk(7)*c(9) + rk(8)*c(10) + 
     *       rk(9)*c(8) + vent + dadj*dep(2)
      l(3) = rk(3)*c(2) + vent + dadj*dep(3)
      l(4) = rk(4)*oh + vent + dadj*dep(4)
      l(5) = rk(5)*oh + rk(6) + vent + dadj*dep(5)
      l(6) = vent + dadj*dep(6)
      l(7) = rk(12) + vent + dadj*dep(7)
      l(8) = rk(9)*c(2) + rk(11)*c(1) + vent + dadj*dep(8)
      l(9) = rk(7)*c(2) + vent + dadj*dep(9)
      l(10) = rk(8)*c(2) + vent + dadj*dep(10)

c determine net loss rates

      do 200 i = 1,10

          g(i) = q(i) - l(i)*c(i)

 200  continue

      return
      end



