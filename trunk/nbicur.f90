subroutine nbicur()
! Calculates the neutral beam current contribution for each flux surface
! Follows method in nbeams program - diffuse beam model
  
  use param

  implicit none

  double precision :: rin, zin
  double precision :: tau_s, v_beam, xi_b
  double precision :: fid, fidf, nfd
  double precision ::  sig_eff
  double precision :: J_f, J_nb
  double precision :: rho
  double precision :: zeff, zni
  integer :: l, i, j
  double precision :: psi,te, ne
  double precision :: tempe, densi, dense
  double precision :: v_c

  integer :: con, rb, zb

  double precision, dimension(ncon) :: kaps, kapps, dels, delps, volps, lambdas
  double precision :: shift, elong
  double precision, dimension(ncon) :: dpdrs 

  double precision :: rr, zz
  integer :: k, ik
  double precision :: rat, del, delp, kap, kapp, dpdr, volp, lambda, D0,D1, beam_atten

  double precision :: rrange, zrange, beam_int
  double precision, dimension(nr,nz) :: h, n_f

  !See no. of cells in beam path
  integer :: count


  rrange = 4.*sig_r
  zrange = 4.*sig_z
  write(nw,*) 'Calculating Neutral beam current'
  kaps = 0.
  kapps = 0.
  dels = 0.
  delps = 0.
  h=0.
  n_f=0.
  count=0
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
  if (icont .gt. -3) then
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
              ne = dense(psi,0)
              !write(nw,*) '1'
              !write(nw,*) zz, Z_beam, sig_z, rr, R_t, sig_r

              !!!!!! CHECK LIMITS WHEN DONE !!!!!
              ! if its within beam height and greater than min R reached by nb
              if (abs(zz - Z_beam) .le. zrange .and. rr .ge. (R_t-rrange) ) then
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
                 
                 rho = sqrt( (rr - (r0 + del))**2 + (zz/kap)**2)

                 !call R_pm(rho, del, kap, zz, R_p, R_m)
              !   write(nw,*) 'Finshed Lin interp'
             !    write(nw,*) rho, volp, kap, kapp, del, delp


                 call dep(rho, volp, kap, kapp, del, delp, rr, zz, lambda, h(i,j))

            !     write(nw,*) 'Dep is ', h(i,j)
                 
                 n_f(i,j) = beam_int(i,j) * h(i,j) * exp ( - ((zz - Z_beam)/sig_z)**2) * (I_0/ (eq* vol)) / (sqrt(pi)*sig_z)  

                 if (n_f(i,j) .ne.0) write(nw,*) n_f(i,j), rr, zz
                 

              end if

           end if
        
              

        end do
     end do
     
  end if
           
  write(nw,*) 'No. of cells in beam = ', count
  write(nw,*) 'exiting loop'

  !te = tempe(psi,0)
!!! Look at E_b units
  v_beam = (2.*E_b/(2.*mp))**0.5


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





! fast ion distribution function
function FIDF()

  use param
  implicit none

  double precision, dimension(1000) :: y, func, d
  
  double precision :: abserr
  double precision   :: int, fidf, I_func
  integer            :: err, size, i, n
  integer :: ierr, neval

!  do i=1,1000
!     y(i) = 1.*i/1000
!     func(i) = I_func(y(i), Z_hat, y_c)
!  end do

  !call qag (I_func, 0, 1, 0.0, 0.001, 6, int, abserr, neval, ierr)


  fidf = (1+y_c)**Z_hat
end function FIDF
  


!Ion distribution function
function I_func(y) 

  use param
  implicit none
  
  double precision, intent(in) :: y
  double precision :: I_func

  I_func = (y**3 / (y**3 + y_c**3)) ** (Z_hat + 1)

end function I_func


subroutine R_pm(rho, del, kap,zz, R_p, R_m)

  use param
  implicit none

  double precision :: rho, del, kap, R_p, R_m, zz
  double precision :: root

  root = sqrt(rho*rho - (zz*zz)/(kap*kap))
  
  R_p = r0 + del + root


  R_m = r0 + del - root

end subroutine R_pm


     
subroutine dep(rho, volp, kap, kapp, del, delp, rr, zz, lam, h)

  use param
  implicit none

  double precision :: rho, volp, kap, kapp, del, delp, rr, zz, lam
  double precision :: h
  double precision :: D0, D1

  
  h  = 2*rho*vol/(volp*lam) * rr/sqrt(rr**2 - R_t**2)

  h = h * ( (1+(kapp* Z_beam**2/(rho*kap**3)))/sqrt(rho**2 - (Z_beam**2/kap**2)) + delp/rho)


     
  !h = h * (exp(-D0) +  gam_d * exp(-D0+ 2*D1) )


end subroutine dep




  

  function beam_atten(R_index,Z_index, id)

    ! calculates beam attentuation for a given radius
    use param
    implicit none

    double precision :: D,rr, beam_atten

    double precision :: psi, dense, inv_lam, Ru
    integer :: i,id, R_index, Z_index, sim_fac

    !!!!! NOT CORRECT Ru !!!!!
    Ru = 2.
    D = 0.

    ! Beam line on the first half of the injection
    if (id .eq. 0) then 
       do i = R_index, nr

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
          
          

          D = D + sim_fac *(rr * inv_lam/sqrt(rr**2 - R_t**2))
       end do


    ! Beam line on the second half of the injection   
    else if (id .eq. 1) then

       do i = R_index,1, -1

          rr = r(i)

          if (rr .lt. R_t) exit
          
          psi = umax - u(i,Z_index)

          inv_lam = dense(psi,0)/(E_b*2.8e17)

          if ( i .eq. R_index .or. i .eq. 1) then
             sim_fac = 1
          else if (mod(i-R_index,2) .eq. 0) then
             sim_fac = 2
          else
             sim_fac = 4
          end if
          
             
          D = D + sim_fac*(rr * inv_lam/sqrt(rr**2 - R_t**2))
      
       end do

    end if
    
       

    beam_atten = D*dr/3.

  end function beam_atten


  

  function beam_int (R_index, Z_index)

    ! Calculates the integral of the beam attenuation and power distribution
    ! using the trapezoid method. Assumes gaussian profile

    
    use param
    implicit none

    integer :: R_index, Z_index, r_in
    double precision :: C, beam_int, range, x, beam_atten

         
    integer :: i, nit, gam_d, sim_fac

    
!!! NORMALISATION ???!!!

    C = 1./ (sqrt(pi) * sig_r)
    beam_int = 0.
    range = 4.* sig_r
    nit = int(range / dr)


        
    do i = -nit, nit

       r_in = R_index + i

       ! Beam can't hit below R_t-range
       if (r(r_in) .lt. (R_t - range) ) cycle

       ! How far into the gaussian the beam goes
       x = r(R_index) - R_t + range

       ! Exit if the NB doesn't go all the way through the gaussian
       if ( ( (i+nit) * dr) .gt. x) exit


       ! If the beam goes through the flux surface twice
       gam_d = 1
       if (r(i) .lt. R_t) then
          gam_d = 0
       end if


       ! Simpsons 1/3 factor for integration
       if (abs(i)-nit .eq. 0) then
          sim_fac = 1
       else if (mod(i+nit,2) .eq. 0) then
          sim_fac = 2
       else
          sim_fac = 4
       end if
       
       

       
       !Integration 
       beam_int = beam_int + sim_fac*(exp(- ((i*dr/sig_r)**2)) * (beam_atten(r_in, Z_index,0) &
            +gam_d*(beam_atten(i, Z_index, 1))) ) 



    end do

    beam_int = beam_int * dr * C /3.

  end function beam_int
  
          
    

    
