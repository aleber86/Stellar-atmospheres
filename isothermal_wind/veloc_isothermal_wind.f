      program ecvel
c     ****************************************************************
c     Solutions of the velocity equation of an isothermal wind driven
c     only by gradient of pressure and gradient of gravitational potential.
c     This solutions are known as "Parker Solutions" (Parker E.N 1958).
      implicit none
      real*8 v(10000), rad(10000), C, vel, vel2
      real*8 consra, consve, funra1, funra2
      integer i, j, n, code
      logical debug, verb, NOTNAN
      parameter (n = 10000)
      external vel, vel2, consra, consve, funra1, funra2, NOTNAN
      
c     -------------------------------------------------------------
c     Initial conditions:
      data debug/.false./,code/0/,v(2000)/1.d0/,rad(2000)/1.d0/,
     &     verb/.false./
c     verb, debug = .true. Debug only 
c     -------------------------------------------------------------
      
      open(50, file='v-vs-r.dat')

      do j=0,100
         C = -3.0d0+(-1)**j*(8.d-2)*j

         call conver(vel, consra, v, rad, C, n, code, debug, verb)
         if(code.eq.1) goto 100
         call escrit(rad, v, n, 50, j, NOTNAN)
c     ------------------------------------------------------------
c     Padding  4*k
         write(50,*)
         write(50,*)
c     ------------------------------------------------------------

         call conver(vel2, consra, v, rad, C, n, code, debug, verb)
         if(code.eq.1) goto 100
         call escrit(rad, v, n, 50, j, NOTNAN )

c     ------------------------------------------------------------
c     Padding  4*k+1
         write(50,*)
         write(50,*)
c     ------------------------------------------------------------
         
         call conver(funra1, consve, rad, v, C, n, code, debug, verb)
         if(code.eq.1) goto 100
         call escrit(rad, v, n, 50, j, NOTNAN)
c     ------------------------------------------------------------
c     Padding 4*k+2
         write(50,*)
         write(50,*)
c     ------------------------------------------------------------
         call conver(funra2, consve, rad, v, C, n, code, debug, verb)
         if(code.eq.1) goto 100
         call escrit(rad, v, n, 50, j, NOTNAN)

c     ------------------------------------------------------------
c     Padding 4*k+3
         write(50,*)
         write(50,*)
      end do

c     ------------------------------------------------------------
c     ************************************************************
c     ------------------------------------------------------------
      

 100  continue
      if(code .eq. 1) write(*,*) 'ITERATION ERROR'
      
      close(50)
      end

      real*8 function vel(varg,cons)
c     *************************************************************
c     varg: square velocity
c     cons: R constant
c     *************************************************************
      implicit none
      real*8 varg, cons

      vel = dlog(varg) + cons
      return
      end

      real*8 function vel2(varg,cons)
c     *************************************************************
c     varg: square velocity
c     cons: R constant
c     *************************************************************
      implicit none
      real*8 varg, cons

      vel2 = dexp(varg - cons)
      return
      end
      
      real*8 function funra1(radio, cons)
      implicit none
      real*8 radio, cons
      
      funra1 = dexp((cons/4.d0)-(1.d0/radio))
      return
      end

      real*8 function funra2(radio, cons)
      implicit none
      real*8 radio, cons
      funra2 = 4.d0/(cons-4.d0*dlog(radio))
      return
      end
      
      real*8 function consra(rad, cons)
c     This function calculates the R constant per
c     r(i). Value needed for every v(i)^2

      real*8 rad, cons
      consra = 4.d0*((1.d0/rad)+dlog(rad))+cons
      
      return
      end

      real*8 function consve(velo, cons)
c     This function calculates the R constant per
c     v(i)^2. Value needed for every r(i).
      
      real*8 velo, cons
      consve = velo - dlog(velo) - cons
      
      return
      end


      
      subroutine conver(func, funcon, arg1, arg2, const,
     & nmax, code, debug, verb)
c     *************************************************************
c     func: tracendental function of velocity or radius
c     funcon: function that determines R constant
c     arg1: vector vel or radius (dependent of arg2)
c     arg2: vector vel or radius (indep.)
c     pasin: iteration step inside critical point
c     pasou: iteration step outside critical point
c     nmax: dim vector vel and rad
c     error: toleration of difference
c     
c    
c     result: calculated value in present step
c     diff: diference in present step
c     maxit: maximum iteration steps
c     debug: .true. print debug values
c     verb: .true. print verbose (debug too)
c     *************************************************************
      implicit none
      integer nmax, i, j, k, l, maxit, paso, code, ultimo
      integer kk
      real*8 func, arg1(nmax), const, diff, result, R
      real*8 pasin, pasou, error, arg2(nmax), resdeb
      real*8 step, funcon
      logical debug, verb
      parameter(pasin = 5.d-4, pasou = 5.d-3, error = 1.d-15)
      parameter(maxit = 1e6)
      
      do j=1,2
         if(j.eq.1) then
c     Condition inside critical point
            paso = 1
            step = pasou
            l= 2001
            ultimo = nmax
            
         else
c     Condition outside critical point 
            paso = -1
            step = -pasin
            l = 1999
            ultimo = 1
         end if
         
      do i = l, ultimo, paso
         arg1(i) = arg1(i-paso)
         arg2(i) = arg2(i-paso) + step

         R = funcon(arg2(i), const)
         diff = 100.d0
         k = 0
         
         do while(k.le.maxit .and. diff.gt.error)
            result = func(arg1(i), R)
            diff = dabs(result-arg1(i))
            arg1(i) = result
c     -------------------------------------------------------------
c     DEBUG: 
            if (debug) then
               if(j==1)then
                  write(*,*)'Inside critical point'
               else
                  write(*,*)'Outside critial point'
               end if
               
                  
               write(*,*)'Iteration step:',k, 'R:',R,
     &              '|arg-iterado|:',diff, 'arg:',arg1(i),
     &              'res. step:', result
            end if
c     -------------------------------------------------------------
c     BREAK CONDITION
            k = k + 1
            
         end do
c     -------------------------------------------------------------
c     VERBOSE:
         if(verb .and. k .le. maxit) then
            call system('echo "ITERATION \033[5;32;7m OK \33[0m" ')
            write(*,*)'Present step: ', k
         else if(verb .and. k .gt. maxit) then
            call system('echo "ITERATION \033[5;37;7m ERROR \33[0m" ')
         end if
         
c     -------------------------------------------------------------
c     DEBUG:
         if(debug .and. k .le. maxit) write(*,*)'Step OK'
c     -------------------------------------------------------------
c     FLUX CONTROL
c     -------------------------------------------------------------
         if (k.gt.maxit) then
            code = 1
            if(debug) write(*,*)'Iter. step:',k,'|arg-iter|:',
     & diff, 'Tolerance:',error, 'C:', const
            return
         end if
c     -------------------------------------------------------------
c     *************************************************************
c     -------------------------------------------------------------
      end do
      end do      
      return
      end

      logical function NOTNAN(value)
      real*8 value
c     Function determines if a value is NOT NaN
      NOTNAN = .false.
      if(value .gt. 0.d0) then
         NOTNAN = .true.
      else if(value .le.0.d0) then
         NOTNAN = .true.
      else
         NOTNAN = .false.
      end if
      end

      subroutine escrit(radio, velo, dim, unit, index, NTNAN)
c     *************************************************************
c     unit: number of file unit
c     index: index for C value. 0 for C = -3
c     NTNAN: Not NaN function
c     *************************************************************
      implicit none
      integer dim, i, unit, index
      real*8 radio(dim), velo(dim)
      logical NTNAN

      do i=1,dim
c     Critical point is contained in transonic and supersoinc-subsonic
c     solution.
         if(index.ne.0 .and. i.eq.2000) then
            continue
         else
            if((radio(i).le.40).and.
     &       (NTNAN(radio(i)).and.NTNAN(velo(i)))) then
               write(unit,*)radio(i), dsqrt(velo(i))
            end if
            
         end if
         
      end do
      return
      end
      
