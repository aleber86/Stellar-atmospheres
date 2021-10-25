c     Program calculates bound-free and free-free
c     opacity coefficient of HI.
c     Values set for spectral types B5V and G2V
c     (15200 K and 5800 K).
c     B5 type stars are expected to be the upper limit for
c     Balmer Jump. For G2 types there is no Balmer Jump and
c     the spectrum takes the "Planckian" form.
      
c     ***************************************************
c     lambda: wavelength
c     n(1-10):  lambda wavelength for (Lyman, Balmer, Paschen...)
c     R: Rydberg's constant [cm^-1]
c     alfa0: constant for the opacity coefficient
c     kb: Boltzmann's constant [eV K^-1]
c     h: Planck's constant [eV s]
c     c: light speed [Ang s^-1]
c     XI: Ionization energy for HI 
c     kbf: Opacity coefficient bound-free
c     kff: Opacity coefficient free-free
c     ***************************************************

      program coefop
      implicit none
      real*8 lambda, n(10), temp(2), R, alfa0
      real*8 XI, kb, kbf, lamini, step, Z, kff
      real*8 h, c
      integer u, i, k, level
      character*3 file(2)
      common /const1/ R, alfa0, XI, kb
      common /const2/ h, c
      
      data file(1)/'B5V'/, file(2)/'G2V'/, temp(1)/15200.d0/,
     &     temp(2)/5800.d0/, n(1)/912.d0/, n(2)/3647.d0/, n(3)/8206.d0/,
     &     n(4)/14588.d0/, n(5)/22790.d0/, n(6)/32820.d0/,
     &     n(7)/44676.d0/, n(8)/58353.d0/, n(9)/73853.d0/,
     &     n(10)/91176.d0/, R/1.09677581d5/, alfa0 /1.0406d-26/,
     &     XI/13.6d0/, kb/8.61733d-5/, lamini/100.d0/,
     &     step/10.d0/, Z/1.d0/, h/4.135667d-15/, c/2.99792458d18/
     
      
      do i = 1,2
         u = 50+i
         open(u, file='data'//file(i)//'.dat')
         
         write(u,*)'#Opacity coefficient HI ('
     &        //file(i)//')'
         write(u,60)'#','1/Lambda [mu m^-1]','k_bf [cm^2]',
     &    'k_ff [cm^2]'
         
         lambda = lamini
         do level = 1,10
            do while(lambda.le.n(level))
            
               write(u,*)1.0d4/(lambda), kbf(lambda, level,
     &              temp(i), Z), kff(lambda, temp(i))
            
               lambda = lambda + step
            end do
         end do

         close(u)

      end do
      
 60   format(x,A,x,A,14x,2(A,15x))
      end

      real*8 function kbf(llam ,ni, teff, Zatom)
c     Bound-free opacity coefficient .
c     *************************************************************
c     gbf: Gaunt bound-free function (quantum mechanics correction)
c     Xn: Energy level function
c     llam: wavelenght 
c     
c     *************************************************************
      implicit none
      integer ni, index, arg, ar
      real*8 teff, sum, gbf, Zatom, R, alfa0, kb
      real*8 Xn, XI, llam, ll
      common /const1/ R, alfa0, XI, kb
      
      gbf(ll, ar)= 1.d0 - 0.3456d0*(Zatom**2/(R*1.d-8*ll))**(1.d0/3.d0)*
     &     (Zatom**2*ll*1.d-8*R/(1.d0*ar)**2 - .5d0)

      Xn(arg) = XI * (1.d0 - 1.d0/(1.d0*arg)**2) 

     
         sum = 0.d0
         do index = ni, ni + 2
            sum = sum + (1.d0/(index*1.d0)**3) * gbf(llam, index) *
     &           exp(-Xn(index)/(kb * Teff))
            
         end do
      
         kbf = alfa0 * llam**3 * (sum + kb*Teff/(2.d0*XI)*   
     &        (dexp(-Xn(ni+3)/(kb * Teff)) - dexp(-XI/(kb * Teff))))
        

      
      return
      end

      real*8 function kff(llam, teff)
c     Free-free opacity coefficient .
c     *************************************************************
c     gff: Gaunt free-free function (quantum mechanics correction)
c     llam: Wavelenght 
c     teff: Surface temperature
c     *************************************************************
      implicit none
      real*8 teff, sum, gff, R, alfa0, kb
      real*8 XI, llam, ll, h, c, T
      common /const1/ R, alfa0, XI, kb
      common /const2/ h, c
      
      gff(ll, T)= 1.d0 + 0.3456d0*(1.d0/(R*1.d-8*ll))**(1.d0/3.d0)*
     &     (kb*T/(h*c/(ll*1.d-8)) + .5d0)

      kff = 1.0406d-40 * gff(llam, teff)*llam**3*
     &     teff*(10.d0)**(1.098d-7/teff)
      
      
      return
      end
