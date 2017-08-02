subroutine nbicur()
! Calculates the neutral beam current contribution for each flux surface
! Follows method in nbeams program - diffuse beam model
  
  use param

  implicit none

  double precision :: rin, zin
  double precision :: tau_s, v_beam, xi_b, coolog
  double precision :: fid, fidf, nfd
  double precision ::  sig_eff
  !double precision :: J_f, J_nb
  double precision :: rho
  double precision :: zeff, zni
  integer :: l, i, j
  double precision :: psi,te, ne
  double precision :: tempe, densi, dense
  double precision :: v_c, y_c, Z_hat, I_func, E_c

  integer :: con, rb, zb, Z_b

  double precision, dimension(ncon) :: kaps, kapps, dels, delps, volps, lambdas
  double precision :: shift, elong
  double precision, dimension(ncon) :: dpdrs 

  double precision :: rr, zz, Ru
  integer :: k, ik,err
  double precision :: rat, del, delp, kap, kapp, dpdr, volp, lambda, D0,D1, beam_atten

  double precision :: rrange, zrange, beam_int
  double precision, dimension(nr,nz) :: h, n_f, h_tot, j_f, J_nb

  !See no. of cells in beam path
  integer :: count

  sig_z = sig_r
  rrange = 4.*sig_r
  zrange = 4.*sig_z
  write(nw,*) 'Calculating Neutral beam current'
  kaps = 0.
  kapps = 0.
  dels = 0.
  delps = 0.
  h=0.
  h_tot=0.
  n_f=0.
  count=0

  Z_b = 1
 
  !V beam (E_b in keV so change to J)
  v_beam = sqrt(2.*E_b*1000*eq/(2.*mp))



  
  ! Calculate shafranov shift and elongation for each flux surface 
 do con = 1, ncon

     dels(con) = shift(con,0)
     delps(con) = shift(con,1)

     kaps(con) = elong(con, 0)
     kapps(con) = elong(con,1)

  end do

  !Calculate dpsi/drho and dvol/drho for each flux surface
  call dpsidrho(dpdrs)
  
  call dVdrho(volps)

  write(nw,*) 'calculated dVdrho'
  call lam(lambdas)
 
  write(nw,*) 'calculated lambdas'
    !write(nw,*) lambdas
!  if (icont .gt. -3) then
     allocate ( nbph(nr,nz))

     !
     do i = 1,nr
        do j=1,nz

        !outside plasma boundary
           if (ixout(i,j) .le. 0) then

              nbph(i,j) = 0.
           else

              !r,z, and psi
              rr = r(i)
              zz = z(j)
              psi = umax - u(i,j)


              !Values for current contribution
              !taken from NEUTRAL-BEAM-HEATING APPLICATIONS AND DEVELOPMENT
              !Menon M., Oak ridge national lab
              ne = dense(psi,0)
              te = tempe(psi,0)
              coolog=log(sqrt(ne*1.0d-6)/te)
              coolog=24.-coolog
              tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)

              E_c = 1.2 * Z_b**(4./3.)  * 2*mp * te /((2.5*mp)**(2./3.) * me**(1./3.))
              v_c = sqrt(2.*E_c*eq/(2.*mp))
              y_c = v_c/v_beam

              zeff=zm
              if (imp.eq.1) then
                 if (ne.gt.0.) then
                    zeff=0.
                    do l=1,nimp+1
                       zni=densi(psi,l,0)
                       zeff=zeff+(zni*iz(l)**2)/ne
                    end do
                 end if
              end if
              Z_hat = 4.*zeff/(5.*A_beam)
              
              !Ru = 
              ! if its within beam height and greater than min R reached by nb
              if (abs(zz - Z_beam) .le. zrange .and. rr .ge. (R_t-rrange)) then
                 !if (
                 ik = 0
                 count = count+1
                 !write(nw,*) '2'
                 do k=1,ncon
                    if (psi .ge. psiv(k)) exit
                    ik = ik+1
                 end do
                 if (ik .eq. 0) then
                    ik=1
                    rat = 0.
                 else if (ik .lt. ncon) then
                    rat = (psi - psiv(ik)) / (psiv(ik+1) - psiv(ik))
                 else
                    write(nw,*) 'cannot interp psi in nbicur'
                    stop
                 end if
               !  write(nw,*) 'Starting lin interp'
                 !Interpolate flux grid values of del/kap to equilibrium grid values
                 del = dels(ik) + rat * (dels(ik+1) - dels(ik))
                 delp = delps(ik) + rat * (delps(ik+1) - delps(ik))
                 kap = kaps(ik) + rat * (kaps(ik+1) - kaps(ik))
                 kapp = kapps(ik) + rat * (kapps(ik+1) - kapps(ik))
                 volp = volps(ik) + rat * (volps(ik+1) - volps(ik))
                 dpdr = dpdrs(ik) + rat * (dpdrs(ik+1) - dpdrs(ik))
                 lambda = lambdas(ik) + rat * (lambdas(ik+1) - lambdas(ik))

                 kapp = kapp * dpdr
                 delp = delp * dpdr
                 volp = volp * dpdr
                 
                 rho = sqrt( (rr - (rcen + del))**2 + (zz/kap)**2)

                 
             !    write(nw,*) rho, volp, kap, kapp, del, delp

                 call dep(rho, volp, kap, kapp, del, delp, rr, zz, lambda, h(i,j))


                 h_tot(i,j) = beam_int(i,j) * h(i,j) * exp ( - ((zz - Z_beam)/sig_z)**2) &
                      / (sqrt(pi)*sig_z * (1-exp(-(rrange/sig_r)**2)) )  

                 !!!!!Need to look at negative deposition values
                 !if(h_tot(i,j).lt.0) write(nw,*) h_tot(i,j), rr,zz
                 
                 
                 !!!! NEED This later for n_f
                 !n_f(i,j) = beam_int(i,j) * h(i,j) * exp ( - ((zz - Z_beam)/sig_z)**2) * (I_0/ (eq* vol)) / (sqrt(pi)*sig_z)

                 n_f(i,j) = h_tot(i,j)*I_0/(eq*vol)

                 xi_b = R_t/(rho+rcen)

                 J_f(i,j) = eq * Z_b * n_f(i,j) * tau_s * xi_b * v_beam * I_func(y_c, Z_hat)

                 if (j.eq.nsym) then
                    print*, 'J_F', J_f(i,j)
                    print*, 'tau_s ', tau_s
                    print*, 'xi_b ', xi_b
                    print*, 'v_c ', v_c
                    print*, 'y_c ', y_c
                    print*, 'E_c ', E_c
                    print*, 'I_func ', I_func(y_c, Z_hat)
                    print*, 'Z_hat ', Z_hat
                    print*, ' '
                 end if
                 
                 !if (n_f(i,j) .ne.0) write(nw,*) n_f(i,j), rr, zz
                 

              end if

           end if
        
              

        end do
     end do
     
!  end if
  write(nw,*) 'No. of cells in beam = ', count
  write(nw,*) 'exiting loop'

  !te = tempe(psi,0)
!!! Look at E_b units

 ! write(nw,*) J_f(:,nsym)

  !zeff...
!  zeff=zm
!  if (imp.eq.1) then
!     if (ne.gt.0.) then
!        zeff=0.
!        do l=1,nimp+1
!           zni=densi(psi,l,0)
!           zeff=zeff+(zni*iz(l)**2)/ne
!        end do
!     else
!        zni=densi(psi,1,0)
!        if (zni.gt.0.) then
!           write(nw,*)'error***problem in nbicur, ne=0'
!           write(nw,*)'cannot evaluate zeff'
!           write(nw,*)'psi=',psi,' ne=',ne,' ni=',zni,' umax=',umax
!           stop
!        else
!           !  zero density so no current-> can arbitrarily set zeff=zm
!           zeff=zm
!        end if
!     end if
!  end if

!  write(nw,*) 'calculated zeff'
 
 ! Z_hat = 4.*zeff/(5.*A_beam)

  ! critcal velocity
!  v_c = 4.5d6 * ( (te/1.d3) / 10)**0.5

!  y_c = v_c/v_beam


  !fid = FIDF()
!  fid  = 0
!  nfd = 0
  !
!  rho = ( (rin - r0) **2 + (zin - Z_beam)**2) ** 0.5

!  xi_b = R_t / (r0 + rho)

!  tau_s = 0.


!  J_f = eq * Z_beam * tau_s *xi_b * v_beam * fid * nfd 

  write(nw,*) 'Completed nb contribution'
  
end subroutine nbicur







! Beam profile for r < r_b inside the beam section
subroutine beam_prof(R_b)

  use param
  implicit none

  double precision :: J_b, C,  R_b, r_beam
  ! Position in beam
  r_beam = (Z_beam**2 + (R_t - R_b)**2) **0.5


  C = 1./(pi* sig_r * sig_r * (1 - exp(-r_b*r_b/(sig_r*sig_r))))

  ! Gaussian Profile
  J_b = C * exp(- r_beam*r_beam /( sig_r * sig_r))

end subroutine beam_prof



!Ion distribution function
function I_func(y_c, Z_hat) 

  use param
  implicit none
  
  double precision, intent(in) :: y_c, Z_hat
  double precision :: I_func, y
  integer :: i,simfac

  I_func = 0.

  

  do i= 0,1000

     y = 1.*i/1000

     if (i.eq.0 .or. i.eq.1000) then
        simfac=1
     else if (mod(i,2) .eq. 0) then
        simfac=2
     else
        simfac=4
     end if
     
     
     I_func =  I_func + simfac* (y**3 / (y**3 + y_c**3)) ** (Z_hat/3. + 1)

  end do

  I_func = I_func * 1./(1000 * 3) * (1+y_c)**(Z_hat/3.) 
end function I_func

     
subroutine dep(rho, volp, kap, kapp, del, delp, rr, zz, lam, h)
!Deposition profile due to beamlet section independant of tangency radius
  use param
  implicit none

  double precision :: rho, volp, kap, kapp, del, delp, rr, zz, lam
  double precision :: h
  double precision :: D0, D1

  
  h  = 2*rho*vol/(volp*lam) * rr

  h = h * ( (1+(kapp* Z_beam**2/(rho*kap**3)))/sqrt(rho**2 - (zz**2/kap**2)) + delp/rho)





end subroutine dep




  

function beam_atten(R_index,Z_index, R_tan, id)

    ! calculates beam attentuation for a given radius and tangency radius

    
    use param
    implicit none

    double precision :: D,rr, beam_atten

    double precision :: psi, dense, inv_lam, Ru, R_tan
    integer :: i,id, R_index, Z_index, sim_fac

    
      

    D = 0.

    !!! NEED to fix simpsons 1/3 for even number of terms
    ! Beam line on the first half of the injection
    if (id .eq. 0) then 
       do i = R_index, nr

          Ru = rcen + (amin**2-z(Z_index)**2/elon**2)**0.5

          if (r(i) .gt. Ru) exit

          rr = r(i)
          psi = umax - u(i,Z_index)

          inv_lam = dense(psi,0)/(E_b*2.8e17)

          ! Simpson's 1/3 rule 
          if ( i .eq. R_index .or. i .eq. nr) then
             sim_fac = 1
          else if ( mod(i-R_index,2) .eq. 0) then
             sim_fac = 2
          else
             sim_fac = 4
          end if
          
          

          D = D + sim_fac *(rr * inv_lam/sqrt(rr**2 - R_tan**2))
       end do


    ! Beam line on the second half of the injection   
    else if (id .eq. 1) then

       do i = R_index,1, -1

          rr = r(i)

          if (rr .lt. R_tan) exit
          
          psi = umax - u(i,Z_index)

          inv_lam = dense(psi,0)/(E_b*2.8e17)

          if ( i .eq. R_index .or. i .eq. 1) then
             sim_fac = 1
          else if (mod(i-R_index,2) .eq. 0) then
             sim_fac = 2
          else
             sim_fac = 4
          end if
          
             
          D = D + sim_fac*(rr * inv_lam/sqrt(rr**2 - R_tan**2))
      
       end do

    end if
    

    beam_atten = D*dr/3.
    !write(nw,*) beam_atten,  rr, z(Z_index)
  end function beam_atten


  

  function beam_int (R_index, Z_index)

    ! Calculates the integral of the beam attenuation and power distribution
    ! using the trapezoid method. Assumes gaussian profile and accounts for
    ! varying tangency radii for different parts of the gaussian

    
    use param
    implicit none

    integer :: R_index, Z_index, r_in
    double precision :: C, beam_int, range, x, beam_atten, R_tan
    double precision :: D0, D1

         
    integer :: i, nit, gam_d, sim_fac

    
!!! NORMALISATION ???!!!

    C = 1./ (sqrt(pi) * sig_r)
    beam_int = 0.
    range = 4.* sig_r
    nit = int(range / dr)


        
    do i = -nit, nit

       r_in = R_index + i
       R_tan = R_t + i*dr
       ! Beam can't hit below R_tan
       if (r(r_in) .lt. (R_tan) ) cycle

       ! How far into the gaussian the beam goes
       x = r(R_index) - R_t + range

       ! Exit if the NB doesn't go all the way through the gaussian
       if ( ( (i+nit) * dr) .gt. x) exit


       ! If the beam goes through the flux surface twice
       gam_d = 1
!       if (r(i) .lt. R_t) then
!          gam_d = 0
!       end if


       ! Simpsons 1/3 factor for integration
       if (abs(i)-nit .eq. 0) then
          sim_fac = 1
       else if (mod(i+nit,2) .eq. 0) then
          sim_fac = 2
       else
          sim_fac = 4
       end if
       
       
       D0 = beam_atten(r_in, Z_index, R_tan,0)
       D1 = beam_atten(r_in,Z_index,R_tan,1)
       
       !Integration for the beam shape
       beam_int = beam_int + sim_fac/sqrt(r(r_in)**2 - R_tan**2)*(exp(- ((i*dr/sig_r)**2)) * &
            (exp(-D0) + gam_d*(exp(-D0 - 2*D1)) ) )   

    end do

    beam_int = beam_int * dr * C /3.

   
  end function beam_int
  
          
    

    
