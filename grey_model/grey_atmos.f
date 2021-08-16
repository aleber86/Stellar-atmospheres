      program modatm
c     *******************************************************
c     Numerical model for a simplified grey atmosphere. 
c     #Parameters:
c     -Surface gravity.
c     -Effective temeprature.
c     -Element abundance (full hydorgen composition).
c     (Tested with data from Kurucz(1979) model.
c      Tau(lambda) and Kappa(lambda); lambda = 500 )
c     *******************************************************
c     Variable names and references
c     *******************************************************
c     Input data:
c     Tau:   Optical depth
c     Kappa: Opacity
c     Teff:  Effective temperature
c     G:     Surface gravity
c     *******************************************************
c     Output data:
c     P:     Pressure per shell
c     T:     Temperature per shell
c     N:     Number of atoms per volume
c     Ne:    Number of electrons per volume
c     Flujo: Energy flux
c     Lumi:  Surface luminosity
c     J:     Mean intensity
c     S:     Font function
c     ********************************************************
c     (International System of Units)
c     ********************************************************
      implicit none
      character*14 archiv
      integer i, n, code
      parameter(n = 30, archiv = 'grey_atmos.dat')
      real*8 temp(n), tau(n), kappa(n), P(n)
      real*8 Nc(n), Nec(n) ,Pe(n), J(n), S(n)
      real*8 tmpfun, press, nfun, presse, intmed, funfue
      real*8 flujo, lumi
      real*8 pi, sigma, kb, teff, g
    
      external tmpfun, press, nfun, presse

c     ******************************************************
c     Constants definition block
c     ******************************************************      
      COMMON /C1/ pi, kb, sigma
      COMMON /C2/ teff
      COMMON /C3/ g
      data g/27542.28d0/, teff/5770.0d0/, kb/1.3807d-16/,
     & sigma/5.6704d-5/
      pi = 4.0d0*datan(1.0d0)
c     ******************************************************
      code = 0
c     ******************************************************
c     Flux and luminosity
c     sigma = 5.6704*10^-5 erg/(cm^2*sec*K^4)
c     ******************************************************
      flujo = sigma * teff**4
      lumi = 4.0d0*pi*(6.96d10)**2 * flujo

      
      open(50, file='temp_tau.dat', status='old')
      read(50, *)
      read(50, *)
      read(50, *)
      
c     Read data from file. tau and kappa
c     Temp(tau) calculation
      do i = 1,n
         read(50,*)tau(i), kappa(i)
c     Entry values (tau and kappa) I: log10(I)
         kappa(i) = 10.d0**kappa(i)
         if (i .eq. 1) then
            tau(i)=0.d0
         else
            tau(i) = 10.d0**tau(i)
         end if
         
         temp(i) = tmpfun(tau(i))
      end do
      close(50)
      

c     ##############################
c     DEBUG function temp(tau)
c      do i = 1,n
c         write(*,*)temp(i)
c      end do
c     #############################

c     #############################
c     Pressure function
c     #############################
      call iter(press, tau, kappa, P, n, code)

c***************************************************************
c     *ERROR OUTPUT
c***************************************************************      
      if (code .eq. 1) then
         call system('echo "\033[5;31;7m  Iteration ERROR. 
     &\033[0m"')
         
         write(*,*)"Verify pressure convergence"
         
         GOTO 100
      end if
c***************************************************************
c###############################################################
c***************************************************************      
      write(*,*)
      write(*,62)
      write(*,*) "Flux [erg/(cm^2*s)]: ", flujo, "Lum [erg/s]: ", lumi
      write(*,62)
      write(*,*)
      write(*,*)'N:  Numerical density H'
      write(*,*)'Ne: Electron density'
      write(*,*)'Pe: Electron pressure'
      write(*,*)'J : Mean intensity'
      write(*,*)'S : Font function'
      write(*,*)
      
      write(*,61)"N [cm^-3]   ", "Ne [cm^-3]   ", "Pe [dyn/cm^2]",
     & "J [erg/(cm^2*s)]" ,"S[erg/(cm^2*s)]"
      do i=1,n
         Nc(i) = nfun(P(i),Temp(i))     
         call saha(Nc(i),Nec(i),Temp(i))
         Pe(i) = presse(Nec(i),Temp(i))
         J(i) = intmed(Temp(i))
         S(i) = funfue(flujo, tau(i))
       
         write(*,*)i, Nc(i), Nec(i), Pe(i) ,J(i),
     &    S(i)
      end do

     
      
      open(51, file = archiv)
c*****************************************************************
c     Header of file
      write(51,*)"#Temperatura y presion del modelo Teff = 5770 K"
      write(51,60)"#Tau", "Kappa","Temp", "Presion", "N", "Ne"
     &     ,"Pe", "J", "S", "|J-S|"
c******************************************************************

      write(*,64)'Pressure value per shell [dyn/cm^2]'
      do i=1,n

c     Imprime en la pantalla los datos de presion por capa
         write(*,63)i, P(i)
         
c     ESCRITURA EN EL ARCHIVO DE SALIDA U = 51
         
         write(51,*)tau(i),kappa(i) ,temp(i), P(i), Nc(i),
     &           Nec(i), Pe(i), J(i), S(i), abs(J(i)-S(i))
                  
      end do
      close(51)
      
 60   format(8(A,10x))
 61   format(9x,5(9x,A,3x))
 62   format(130('*'))
 63   format("P(",I2")= ", g26.20)
 64   format(2/,44('*'),A,44('*'),/)
 

      write(*,*)
      write(*,*)'**** Results written in : ', archiv,' ****'

c------------------------------------------------------------
c     Error output 
 100  continue
c------------------------------------------------------------      
      end

      real*8 function tmpfun(tau)
      implicit none
c     Temperature law (tau dependent)
      real*8 tau, teff
      COMMON /C2/ teff
      tmpfun = teff*(3.d0/4.d0*tau + 1.d0/2.d0)**(1.d0/4.d0)
      return
      end

      

      
      real*8 function press(Pin, Pf, taui, tauf ,kappai, kappaf)
      implicit none
c     Pressure function (trapezoidal rule)
c     P_1 = F(P1,P2,tau_i,tau_f,kappa_i, kappaf)
      real*8 Pin, Pf, taui, tauf, kappai, kappaf, g

      COMMON /C3/ g
     
     
      press = ( dsqrt(Pin)**3.d0 + 3.d0/2.d0*g*(dsqrt(Pf)/kappaf
     &     + dsqrt(Pin)/kappai)*(tauf - taui)/2.d0)**(2.d0/3.d0)
      return
      end

      real*8 function nfun(P,T)
      implicit none
c     Number of H atoms
c     (State equation - ideal gas)
      real*8 P, T, pi, kb, sigma
      COMMON /C1/ pi, kb, sigma
  
      nfun = P /( T * kb)
      return
      end

      real*8 function presse(Nel, T)
      implicit none
c     Elec. Press.
c     (Sate equation - ideal gas)
      real*8 Nel, T, kb, pi, sigma
      COMMON /C1/ pi, kb, sigma
      
      presse = Nel*kb*T
      return
      end

      real*8 function intmed(temp)
      real*8 temp, pi, kb, sigma
      COMMON /C1/ pi, kb, sigma
      
      intmed = sigma*(temp**4)/pi
      return
      end

      real*8 function funfue(flux, tau)
      real*8 tau, flux, pi, kb, sigma
      COMMON /C1/ pi, kb, sigma

      funfue = 3.d0/4.d0*(tau + 2.d0/3.d0)*flux/pi
      return
      end
      
      
      subroutine saha(N, Ne, T)
c     Electron density 
      implicit none
      real*8 Ne, T, N, pi, kb, sigma
      real*8 phi, disc, respos, resneg
      COMMON /C1/ pi, kb, sigma
      
      phi(T) = 0.6665d0*2.0d0*T**(5.d0/2.d0)
     & *10.0d0**(-5040.d0 * 13.6d0/T)
    
c     ------------------------------------------------
c     Quadratic equation
      disc = dsqrt(4.d0*phi(T)**2.d0 + 4.d0*T*kb*N*phi(T))

      respos = (-2.d0 * phi(T) + disc)/(2.d0*kb*T)
      resneg = (-2.d0 * phi(T) - disc)/(2.d0*kb*T)
c     -------------------------------------------------

c     Condition for phisical solution.
c     
      if (respos .gt. 0.0d0 .and. respos .le. N) then
         Ne = respos
c     --------------------------------------------------
c     ****DEBUG         
c     write(*,*)'POS'
c     ---------------------------------------------------         
      else if (resneg .gt. 0.0d0 .and. resneg .le. N) then
c     ---------------------------------------------------
c     ****DEBUG         
c     write(*,*)'NEG'
c     ---------------------------------------------------         
         Ne = resneg
      end if
      return
      end

      
      subroutine iter(func, tau, kappa, P, n, code)
      implicit none
c     Iterative method.
c     *********************************************
c     User defined values:
c     nmax:  max number of iteration.
c     diff: |P_f - P_i|<error (min. cut condition)
c     code: 0 OK, 1 ERROR ITER.
c     P(1): first value (known condition)
c     *********************************************
      integer i, n, nmax, code, k
      real*8 func, tau(n), kappa(n), P(n)
      real*8 error, diff, result

      error = 1.0d-5
      nmax = 10**4
      P(1) = 1.0d3
     
      do k = 2,n
         P(k) = P(k-1)
         diff = 100.d0
         i= 0
         do while (i .le. nmax .and. diff .gt. error)

          result = func(P(k-1), P(k), tau(k-1),tau(k),
     &        kappa(k-1),kappa(k))
c    ----DEBUG--------------------------------------------         
         write(*,*)"i: ", i, abs(result-P(k)), result,
     &        P(k), "k:", k
c    -------------------------------------------------------
         diff = abs(result - P(k))
         P(k) = result
         i=i+1
      end do
      
      write(*,*)"Index number of iteration: ", k
      
c     ***********************************************************      
c     ERROR CONDITION. BREAK.
c     **********************************************************      
      if (i>nmax) then
         write(*,*)'Number of max. it:',nmax,'(namx)'
         write(*,*)'Iteration number:',i, '(i)'
         write(*,*)'Max difference:',error, '(error)'
         code = 1
         return
      else
c     **********************************************************
c     *OK*
c     **********************************************************
c     ****DEBUG         
c         write(*,65)diff
 65   format(g50.30e3)

c     **********************************************************
         code = 0
         call system('echo "ITERATION \033[5;32;7m OK \33[0m"')
      end if
      end do
      
      return
     
      end
