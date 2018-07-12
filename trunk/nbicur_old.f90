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

  double precision :: rrange, zrange, beam_int, elec_return
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

                 ! Deposition for pencil beam
                 call dep(rho, volp, kap, kapp, del, delp, rr, zz, lambda, h(i,j))


                 !Account for beam shape
                 h_tot(i,j) = beam_int(i,j) * h(i,j) * exp ( - ((zz - Z_beam)/sig_z)**2) &
                      / (sqrt(pi)*sig_z * (1-exp(-(rrange/sig_r)**2)) )

                 !!!!!Need to look at negative deposition values
                 !if(h_tot(i,j).lt.0) write(nw,*) h_tot(i,j), rr,zz


                 !!!! NEED This later for n_f
                 !n_f(i,j) = beam_int(i,j) * h(i,j) * exp ( - ((zz - Z_beam)/sig_z)**2) * (I_0/ (eq* vol)) / (sqrt(pi)*sig_z)

                 n_f(i,j) = h_tot(i,j)*I_0/(eq*vol)

                 xi_b = R_t/(rho+rcen)

                 J_f(i,j) = eq * Z_b * n_f(i,j) * tau_s * xi_b * v_beam * I_func(y_c, Z_hat)

                 nbph(i,j) = elec_return(J_f(i,j), zeff, Z_b,u(i,j))

                 if (j.eq.nsym .and. rr.gt. 1.1) then

                    print*, 'NB', nbph(i,j)
                    print*, 'Beam atten', beam_int(i,j)
                    print*, 'tau_s ', tau_s
                    print*, 'xi_b ', xi_b
                    print*, 'n_f', n_f(i,j)
                    print*, 'I_func ', I_func(y_c, Z_hat)

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



!  write(nw,*) 'calculated zeff'

 ! Z_hat = 4.*zeff/(5.*A_beam)

  ! critcal velocity
!  v_c = 4.5d6 * ( (te/1.d3) / 10)**0.5

!  y_c = v_c/v_beam


  write(nw,*) 'Completed nb contribution'

end subroutine nbicur



!Do it via flux surfaces

subroutine nbicur2()
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
  double precision :: psi,te, ne, dpsi, shine, depsum
  double precision :: tempe, densi, dense
  double precision :: v_c, y_c, Z_hat, I_func, E_c

  integer :: con, rb, zb, Z_b

  double precision, dimension(ncon) ::  volps, lambdas
  double precision :: shift, elong, kaps, kapps, dels, delps
  double precision, dimension(ncon) :: dpdrs

  double precision :: rr, zz, Ru

  double precision :: depp, depm, dep2

  integer :: k, ik,err
  double precision :: rat, del, delp, kap, kapp, dpdr, volp, lambda

  double precision :: rrange, zrange, beam_int, elec_return, uval
  double precision, dimension(ncon) :: n_f, h_tot, j_f, dep
  double precision, dimension(nr) :: J_tot


  double precision :: fprof, bphi, bth, fsi, B_tot, bp
  !See no. of cells in beam path
  integer :: count


  write(nw,*) 'Calculating Neutral beam current'



  n_f=0.
  count=0
  dep = 0.
  Z_b = 1

  !V beam (E_b in keV so change to J)
  v_beam = sqrt(2.*E_b*1000*eq/(2.*mp))

 !Calculate dpsi/drho and dvol/drho for each flux surface
  call dpsidrho(dpdrs)

  call dVdrho(volps)

  write(nw,*) 'calculated dVdrho'
  call lam(lambdas)


  write(nw,*) 'calculated lambdas'
    !write(nw,*) lambdas
!  if (icont .gt. -3) then


  J_nb(1) = 0
  ! Calculate shafranov shift and elongation for each flux surface
  do con = 2, ncon

     psi = psiv(con)
     rr = maxval(rpts(con,:))

     if (rr .lt. R_t) cycle
     lambda = lambdas(con)
     volp = volps(con)
     dpdr = dpdrs(con)

     kap = elong(con,0)
     kapp = elong(con,1) * dpdr
     del = shift(con,0)
     delp = shift(con,1) * dpdr

     rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)

     depp =  dep2(con, volp, lambda, kap, kapp, del, delp, 1)

     if ( rho**2 .ge. (Z_beam/kap)**2) then
        depm =  dep2(con, volp, lambda, kap, kapp, del, delp, -1)
     else
        depm = 0.
     end if

     dep(con) = depp + depm


     !write(nw,*) 'dep found', con, depp, depm


     ne = dense(psi,0)
     te = tempe(psi,0)

     coolog=log(sqrt(ne*1.0d-6)/te)
     coolog=24.-coolog
     tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)



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

     E_c = 1.2 * Z_b**(4./3.)  * 2*mp * te /((2.5*mp)**(2./3.) * me**(1./3.))
     v_c = sqrt(2.*E_c*eq/(2.*mp))
     y_c = v_c/v_beam


     n_f(con) = I_0/eq * dep(con) / vol

     xi_b = R_t/(rho+rcen)


     J_f(con) = eq * Z_b * n_f(con) * tau_s * xi_b * v_beam * I_func(y_c, Z_hat)

     uval = umax - psi

     J_nb(con) = elec_return(J_f(con), zeff, Z_b, uval) * bdl(con) / bsqav(con)

  end do

  J_tot = 0.



  depsum = 0.
!!! Need to see what happens at i=ncon (core)
  do i=2,ncon-1
     dpsi = (psiv(i-1) + psiv(i+1))/2.
     print*, i, dep(i), dpsi, volps(i)
     depsum = depsum + dep(i)*volps(i)*dpsi

  end do

  print*, 'depsum is ', depsum
  shine = 1 - depsum/vol

  dep = dep/(1.-shine)
  J_nb = J_nb/(1.-shine)



  ! do i=1,nr
 !    if (ixout(i,nsym) .ne. 1) cycle
 !    psi = umax - u(i,nsym)
 !    do j=1,ncon
 !       if (psi .ge. psiv(j)) then

 !          rat = (psi - psiv(j))/(psiv(j-1) - psiv(j))
 !          exit
 !       end if
 !    end do


 !       rr = r(i)
 !       zz = z(nsym)
 !       bth = bp(rr,zz)

 !       fsi = fprof(psi,2)
 !       bphi = fsi/rr
 !       B_tot = sqrt( bth*bth+bphi*bphi )

 !       J_tot(i) =  (J_nb(j) + rat*(J_nb(j-1) - J_nb(j)) ) * B_tot


  !   end do





!  end if
  write(nw,*) 'No. of cells in beam2 = ', count
  write(nw,*) 'exiting loop2'
  write(nw,*) 'Completed nb2 contribution'

end subroutine nbicur2






subroutine deposition(rho, con, rr, zz, dpdrs, volp,lambda, h)

  use param
  implicit none

  double precision :: rr, zz, dpdrs, volp, lambda, rho
  integer :: con

  double precision :: del, delp, kap, kapp, shift, elong

  double precision :: psi, ne, te, dense, tempe

  double precision :: coolog, tau_s

  double precision :: E_c, v_c, y_c, Z_b, v_beam, Z_hat

  double precision :: zeff, zni, densi
  integer :: l, gam_d

  double precision :: h, attenuation, D0, D1



  del = shift(con,0)
  delp = shift(con,1)

  kap = elong(con, 0)
  kapp = elong(con,1)

  psi = psiv(con)


  ne = dense(psi,0)
  te = tempe(psi,0)
  coolog=log(sqrt(ne*1.0d-6)/te)
  coolog=24.-coolog
  tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)

  v_beam = sqrt(2.*E_b*1000*eq/(2.*mp))


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

  kapp = kapp * dpdrs
  delp = delp * dpdrs
  volp = volp * dpdrs




  !    write(nw,*) rho, volp, kap, kapp, del, delp

  ! Deposition for pencil beam
  call dep(rho, volp, kap, kapp, del, delp, rr, zz, lambda, h)

  D0 = attenuation(rr,0)
  if (rpts(con,1) .ge. R_t) then
     gam_d = 1

     D1 = attenuation(rr,1)

  else
     gam_d = 0
     D1 = 0.
  end if



  h = h *( exp(-D0) + gam_d*(exp(-(D0+2*D1)))  )


end subroutine deposition







function dep2(con, volp, lambda, kap, kapp, del, delp, id)
  !calculates deposition profile from beams module
  !performs two integrals, over R for a fixed Z, and then integrates that
  !for all Z in the beam path

  use param
  implicit none

  double precision :: volp, lambda, kap, kapp, del, delp, dep2
  integer :: con, id, sgn

  double precision :: rr, rho, Zu, Rl, Ru, r_beam, rpm

  double precision :: wi, zti, zxi, rti, rxi

  double precision :: zterm, zsum, rterm, rsum, expterm, gauss

  double precision :: D0, D1, attenuation

  integer :: i, j, dpass


  !rho for each flux surface (outboard midplane for each flux sur
  rr = maxval(rpts(con,:))
  rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)


  r_beam = sig_r * 2

  Zu = min(rho*kap, r_beam)

  wi = pi/6
  zterm = 0.
  zsum = 0.


  !Z integration
  do i=1,6

     !Gaussian terms
     zti = cos((2*i -1)*pi/(2*6))
     zxi = (Zu/2) * (1+zti)

     !Calc R+ or R-
     if (id .eq. 1) then
        rpm = rcen + del + sqrt(rho**2 - (zxi/kap)**2)
        sgn = 1
     else if (id .eq. -1) then

        rpm = rcen + del - sqrt(rho**2 - (zxi/kap)**2)
        sgn = -1
     end if

     !Checks to see if R+- is below the tangency radius and if so cycles
     if (rpm .lt. (R_t - sqrt(r_beam**2-zxi**2) ) ) cycle


     !Terms in h(p) that involve Z
     zterm = (1 + kapp*zxi**2/(rho*kap**3))/sqrt(rho**2 - zxi**2/kap**2) + sgn*delp/rho

     !Multiplied by sqrt(1-ti^2)
     zterm = zterm * sqrt(1 - zti**2)


     !Limits for R integration
     Rl = R_t - sqrt(r_beam**2 - zxi**2)
     Ru = min(rpm, R_t + sqrt(r_beam**2 - zxi**2))


     !R integration
     rsum = 0.
     rterm = 0.
     do j=1,6

        !Gaussian terms for R
        rti = cos( (2*j - 1) * pi / (2*6) )
        rxi = (Ru-Rl)/(2.*rti)  + (Ru+Rl)/2

        !Checks to see if R+- is within integral
        if (rpm .lt. rxi ) cycle

        !Calc attenuation and checks for double pass
        D0 = attenuation(rpm,rxi,0)
        dpass = 0
        if (rxi .ge. R_t) then
           dpass = 1

           D1 = attenuation(rpm,rxi,1)
        end if

        !Total atten
        expterm = exp(-D0) + dpass*exp(-(D0+2*D1))

        !R terms in h(p)
        rterm = rpm/sqrt(rpm*2 - rxi**2) * expterm

!!! Check this again if numbers are weird
        !Gaussian for beam shape
        gauss = exp( - (zxi**2 + (R_t - rxi)**2)/sig_r**2) / (pi*sig_r**2 * (1 - exp(-r_beam**2/sig_r**2)) )




        !R term
        rterm = rterm * wi * gauss * sqrt(1 - rti**2)

        rsum = rsum + rterm
     end do

     !multiply by limits (b-a/2)
     rsum = rsum * (Ru-Rl)/2.

     !Multiplies R and Z contribution of integral and adds to Z total
     zsum = zsum + zterm*rsum

  end do

  zsum = zsum * wi * 2 *Zu/2.


  dep2 = 2 * rho * vol * zsum / (volp*lambda)

end function dep2




function attenuation(rpm, rb, id)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, attenuation
  integer :: id, i, con, simfac

  double precision :: Ru, D, dense, invlam

  double precision :: rmin, rmax, r_int, rat, dr_nb

  D = 0.



  if (id .eq. 0) then

     Ru = rcen + (amin**2-Z_beam**2/elon**2)**0.5


     dr_nb = (Ru-rpm)/100

     do i= 0,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rpm + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b)




        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))


     end do

  else if (id .eq. 1) then


     dr_nb = (rpm-rb)/100

     !stop atten go to inf
     do i= 1,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rb + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b)

        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))

     end do


     end if





     attenuation = D*dr_nb/3.


end function attenuation





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

  !Account for inner or outer part of flux surface
  if (rr .lt. r0) then
     h = h * ( (1+(kapp* Z_beam**2/(rho*kap**3)))/sqrt(rho**2 - (zz**2/kap**2)) + delp/rho)

  else

     h = h * ( (1+(kapp* Z_beam**2/(rho*kap**3)))/sqrt(rho**2 + (zz**2/kap**2)) + delp/rho)
  end if



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

          Ru = rcen + (amin**2-Z_beam**2/elon**2)**0.5

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


  function elec_return(J_f, zeff, Z_b,uval)
    !Calculates the current after the electron return
    ! current is accounted for

    use param
    implicit none

    double precision :: J_f, zeff, epsi, elec_return, uval
    double precision :: G, sqrt, J_nb
    integer :: Z_b

    call argen(uval, epsi,0)

    G = ( 1.55+(0.85/zeff) ) * sqrt(epsi)  - ( 0.2 + (1.55/zeff) )*epsi



    J_nb = (1 - ( Z_b/zeff * (1 - G) ) )


    J_nb = J_nb *J_f


    elec_return = J_nb

  end function elec_return


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

                 ! Deposition for pencil beam
                 call dep(rho, volp, kap, kapp, del, delp, rr, zz, lambda, h(i,j))


                 !Account for beam shape
                 h_tot(i,j) = beam_int(i,j) * h(i,j) * exp ( - ((zz - Z_beam)/sig_z)**2) &
                      / (sqrt(pi)*sig_z * (1-exp(-(rrange/sig_r)**2)) )

                 !!!!!Need to look at negative deposition values
                 !if(h_tot(i,j).lt.0) write(nw,*) h_tot(i,j), rr,zz


                 !!!! NEED This later for n_f
                 !n_f(i,j) = beam_int(i,j) * h(i,j) * exp ( - ((zz - Z_beam)/sig_z)**2) * (I_0/ (eq* vol)) / (sqrt(pi)*sig_z)

                 n_f(i,j) = h_tot(i,j)*I_0/(eq*vol)

                 xi_b = R_t/(rho+rcen)

                 J_f(i,j) = eq * Z_b * n_f(i,j) * tau_s * xi_b * v_beam * I_func(y_c, Z_hat)

                 nbph(i,j) = elec_return(J_f(i,j), zeff, Z_b,u(i,j))

                 if (j.eq.nsym .and. rr.gt. 1.1) then

                    print*, 'NB', nbph(i,j)
                    print*, 'Beam atten', beam_int(i,j)
                    print*, 'tau_s ', tau_s
                    print*, 'xi_b ', xi_b
                    print*, 'n_f', n_f(i,j)
                    print*, 'I_func ', I_func(y_c, Z_hat)

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



!  write(nw,*) 'calculated zeff'

 ! Z_hat = 4.*zeff/(5.*A_beam)

  ! critcal velocity
!  v_c = 4.5d6 * ( (te/1.d3) / 10)**0.5

!  y_c = v_c/v_beam


  write(nw,*) 'Completed nb contribution'

end subroutine nbicur



!Do it via flux surfaces

subroutine nbicur2()
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
  double precision :: psi,te, ne, dpsi, shine, depsum
  double precision :: tempe, densi, dense
  double precision :: v_c, y_c, Z_hat, I_func, E_c

  integer :: con, rb, zb, Z_b

  double precision, dimension(ncon) ::  volps, lambdas
  double precision :: shift, elong, kaps, kapps, dels, delps
  double precision, dimension(ncon) :: dpdrs

  double precision :: rr, zz, Ru

  double precision :: depp, depm, dep2

  integer :: k, ik,err
  double precision :: rat, del, delp, kap, kapp, dpdr, volp, lambda

  double precision :: rrange, zrange, beam_int, elec_return, uval
  double precision, dimension(ncon) :: n_f, h_tot, j_f, dep
  double precision, dimension(nr) :: J_tot


  double precision :: fprof, bphi, bth, fsi, B_tot, bp
  !See no. of cells in beam path
  integer :: count


  write(nw,*) 'Calculating Neutral beam current'



  n_f=0.
  count=0
  dep = 0.
  Z_b = 1

  !V beam (E_b in keV so change to J)
  v_beam = sqrt(2.*E_b*1000*eq/(2.*mp))

 !Calculate dpsi/drho and dvol/drho for each flux surface
  call dpsidrho(dpdrs)

  call dVdrho(volps)

  write(nw,*) 'calculated dVdrho'
  call lam(lambdas)


  write(nw,*) 'calculated lambdas'
    !write(nw,*) lambdas
!  if (icont .gt. -3) then


  J_nb(1) = 0
  ! Calculate shafranov shift and elongation for each flux surface
  do con = 2, ncon

     psi = psiv(con)
     rr = maxval(rpts(con,:))

     if (rr .lt. R_t) cycle
     lambda = lambdas(con)
     volp = volps(con)
     dpdr = dpdrs(con)

     kap = elong(con,0)
     kapp = elong(con,1) * dpdr
     del = shift(con,0)
     delp = shift(con,1) * dpdr

     rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)

     depp =  dep2(con, volp, lambda, kap, kapp, del, delp, 1)

     if ( rho**2 .ge. (Z_beam/kap)**2) then
        depm =  dep2(con, volp, lambda, kap, kapp, del, delp, -1)
     else
        depm = 0.
     end if

     dep(con) = depp + depm


     !write(nw,*) 'dep found', con, depp, depm


     ne = dense(psi,0)
     te = tempe(psi,0)

     coolog=log(sqrt(ne*1.0d-6)/te)
     coolog=24.-coolog
     tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)



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

     E_c = 1.2 * Z_b**(4./3.)  * 2*mp * te /((2.5*mp)**(2./3.) * me**(1./3.))
     v_c = sqrt(2.*E_c*eq/(2.*mp))
     y_c = v_c/v_beam


     n_f(con) = I_0/eq * dep(con) / vol

     xi_b = R_t/(rho+rcen)


     J_f(con) = eq * Z_b * n_f(con) * tau_s * xi_b * v_beam * I_func(y_c, Z_hat)

     uval = umax - psi

     J_nb(con) = elec_return(J_f(con), zeff, Z_b, uval) * bdl(con) / bsqav(con)

  end do

  J_tot = 0.



  depsum = 0.
!!! Need to see what happens at i=ncon (core)
  do i=2,ncon-1
     dpsi = (psiv(i-1) + psiv(i+1))/2.
     print*, i, dep(i), dpsi, volps(i)
     depsum = depsum + dep(i)*volps(i)*dpsi

  end do

  print*, 'depsum is ', depsum
  shine = 1 - depsum/vol

  dep = dep/(1.-shine)
  J_nb = J_nb/(1.-shine)



  ! do i=1,nr
 !    if (ixout(i,nsym) .ne. 1) cycle
 !    psi = umax - u(i,nsym)
 !    do j=1,ncon
 !       if (psi .ge. psiv(j)) then

 !          rat = (psi - psiv(j))/(psiv(j-1) - psiv(j))
 !          exit
 !       end if
 !    end do


 !       rr = r(i)
 !       zz = z(nsym)
 !       bth = bp(rr,zz)

 !       fsi = fprof(psi,2)
 !       bphi = fsi/rr
 !       B_tot = sqrt( bth*bth+bphi*bphi )

 !       J_tot(i) =  (J_nb(j) + rat*(J_nb(j-1) - J_nb(j)) ) * B_tot


  !   end do





!  end if
  write(nw,*) 'No. of cells in beam2 = ', count
  write(nw,*) 'exiting loop2'
  write(nw,*) 'Completed nb2 contribution'

end subroutine nbicur2






subroutine deposition(rho, con, rr, zz, dpdrs, volp,lambda, h)

  use param
  implicit none

  double precision :: rr, zz, dpdrs, volp, lambda, rho
  integer :: con

  double precision :: del, delp, kap, kapp, shift, elong

  double precision :: psi, ne, te, dense, tempe

  double precision :: coolog, tau_s

  double precision :: E_c, v_c, y_c, Z_b, v_beam, Z_hat

  double precision :: zeff, zni, densi
  integer :: l, gam_d

  double precision :: h, attenuation, D0, D1



  del = shift(con,0)
  delp = shift(con,1)

  kap = elong(con, 0)
  kapp = elong(con,1)

  psi = psiv(con)


  ne = dense(psi,0)
  te = tempe(psi,0)
  coolog=log(sqrt(ne*1.0d-6)/te)
  coolog=24.-coolog
  tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)

  v_beam = sqrt(2.*E_b*1000*eq/(2.*mp))


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

  kapp = kapp * dpdrs
  delp = delp * dpdrs
  volp = volp * dpdrs




  !    write(nw,*) rho, volp, kap, kapp, del, delp

  ! Deposition for pencil beam
  call dep(rho, volp, kap, kapp, del, delp, rr, zz, lambda, h)

  D0 = attenuation(rr,0)
  if (rpts(con,1) .ge. R_t) then
     gam_d = 1

     D1 = attenuation(rr,1)

  else
     gam_d = 0
     D1 = 0.
  end if



  h = h *( exp(-D0) + gam_d*(exp(-(D0+2*D1)))  )


end subroutine deposition







function dep2(con, volp, lambda, kap, kapp, del, delp, id)
  !calculates deposition profile from beams module
  !performs two integrals, over R for a fixed Z, and then integrates that
  !for all Z in the beam path

  use param
  implicit none

  double precision :: volp, lambda, kap, kapp, del, delp, dep2
  integer :: con, id, sgn

  double precision :: rr, rho, Zu, Rl, Ru, r_beam, rpm

  double precision :: wi, zti, zxi, rti, rxi

  double precision :: zterm, zsum, rterm, rsum, expterm, gauss

  double precision :: D0, D1, attenuation

  integer :: i, j, dpass


  !rho for each flux surface (outboard midplane for each flux sur
  rr = maxval(rpts(con,:))
  rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)


  r_beam = sig_r * 2

  Zu = min(rho*kap, r_beam)

  wi = pi/6
  zterm = 0.
  zsum = 0.


  !Z integration
  do i=1,6

     !Gaussian terms
     zti = cos((2*i -1)*pi/(2*6))
     zxi = (Zu/2) * (1+zti)

     !Calc R+ or R-
     if (id .eq. 1) then
        rpm = rcen + del + sqrt(rho**2 - (zxi/kap)**2)
        sgn = 1
     else if (id .eq. -1) then

        rpm = rcen + del - sqrt(rho**2 - (zxi/kap)**2)
        sgn = -1
     end if

     !Checks to see if R+- is below the tangency radius and if so cycles
     if (rpm .lt. (R_t - sqrt(r_beam**2-zxi**2) ) ) cycle


     !Terms in h(p) that involve Z
     zterm = (1 + kapp*zxi**2/(rho*kap**3))/sqrt(rho**2 - zxi**2/kap**2) + sgn*delp/rho

     !Multiplied by sqrt(1-ti^2)
     zterm = zterm * sqrt(1 - zti**2)


     !Limits for R integration
     Rl = R_t - sqrt(r_beam**2 - zxi**2)
     Ru = min(rpm, R_t + sqrt(r_beam**2 - zxi**2))


     !R integration
     rsum = 0.
     rterm = 0.
     do j=1,6

        !Gaussian terms for R
        rti = cos( (2*j - 1) * pi / (2*6) )
        rxi = (Ru-Rl)/(2.*rti)  + (Ru+Rl)/2

        !Checks to see if R+- is within integral
        if (rpm .lt. rxi ) cycle

        !Calc attenuation and checks for double pass
        D0 = attenuation(rpm,rxi,0)
        dpass = 0
        if (rxi .ge. R_t) then
           dpass = 1

           D1 = attenuation(rpm,rxi,1)
        end if

        !Total atten
        expterm = exp(-D0) + dpass*exp(-(D0+2*D1))

        !R terms in h(p)
        rterm = rpm/sqrt(rpm*2 - rxi**2) * expterm

!!! Check this again if numbers are weird
        !Gaussian for beam shape
        gauss = exp( - (zxi**2 + (R_t - rxi)**2)/sig_r**2) / (pi*sig_r**2 * (1 - exp(-r_beam**2/sig_r**2)) )




        !R term
        rterm = rterm * wi * gauss * sqrt(1 - rti**2)

        rsum = rsum + rterm
     end do

     !multiply by limits (b-a/2)
     rsum = rsum * (Ru-Rl)/2.

     !Multiplies R and Z contribution of integral and adds to Z total
     zsum = zsum + zterm*rsum

  end do

  zsum = zsum * wi * 2 *Zu/2.


  dep2 = 2 * rho * vol * zsum / (volp*lambda)

end function dep2




function attenuation(rpm, rb, id)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, attenuation
  integer :: id, i, con, simfac

  double precision :: Ru, D, dense, invlam

  double precision :: rmin, rmax, r_int, rat, dr_nb

  D = 0.



  if (id .eq. 0) then

     Ru = rcen + (amin**2-Z_beam**2/elon**2)**0.5


     dr_nb = (Ru-rpm)/100

     do i= 0,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rpm + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b)




        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))


     end do

  else if (id .eq. 1) then


     dr_nb = (rpm-rb)/100

     !stop atten go to inf
     do i= 1,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rb + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b)

        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))

     end do


     end if





     attenuation = D*dr_nb/3.


end function attenuation





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

  !Account for inner or outer part of flux surface
  if (rr .lt. r0) then
     h = h * ( (1+(kapp* Z_beam**2/(rho*kap**3)))/sqrt(rho**2 - (zz**2/kap**2)) + delp/rho)

  else

     h = h * ( (1+(kapp* Z_beam**2/(rho*kap**3)))/sqrt(rho**2 + (zz**2/kap**2)) + delp/rho)
  end if



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

          Ru = rcen + (amin**2-Z_beam**2/elon**2)**0.5

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


  function elec_return(J_f, zeff, Z_b,uval)
    !Calculates the current after the electron return
    ! current is accounted for

    use param
    implicit none

    double precision :: J_f, zeff, epsi, elec_return, uval
    double precision :: G, sqrt, J_nb
    integer :: Z_b

    call argen(uval, epsi,0)

    G = ( 1.55+(0.85/zeff) ) * sqrt(epsi)  - ( 0.2 + (1.55/zeff) )*epsi



    J_nb = (1 - ( Z_b/zeff * (1 - G) ) )


    J_nb = J_nb *J_f


    elec_return = J_nb

  end function elec_return




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

                 ! Deposition for pencil beam
                 call dep(rho, volp, kap, kapp, del, delp, rr, zz, lambda, h(i,j))


                 !Account for beam shape
                 h_tot(i,j) = beam_int(i,j) * h(i,j) * exp ( - ((zz - Z_beam)/sig_z)**2) &
                      / (sqrt(pi)*sig_z * (1-exp(-(rrange/sig_r)**2)) )

                 !!!!!Need to look at negative deposition values
                 !if(h_tot(i,j).lt.0) write(nw,*) h_tot(i,j), rr,zz


                 !!!! NEED This later for n_f
                 !n_f(i,j) = beam_int(i,j) * h(i,j) * exp ( - ((zz - Z_beam)/sig_z)**2) * (I_0/ (eq* vol)) / (sqrt(pi)*sig_z)

                 n_f(i,j) = h_tot(i,j)*I_0/(eq*vol)

                 xi_b = R_t/(rho+rcen)

                 J_f(i,j) = eq * Z_b * n_f(i,j) * tau_s * xi_b * v_beam * I_func(y_c, Z_hat)

                 nbph(i,j) = elec_return(J_f(i,j), zeff, Z_b,u(i,j))

                 if (j.eq.nsym .and. rr.gt. 1.1) then

                    print*, 'NB', nbph(i,j)
                    print*, 'Beam atten', beam_int(i,j)
                    print*, 'tau_s ', tau_s
                    print*, 'xi_b ', xi_b
                    print*, 'n_f', n_f(i,j)
                    print*, 'I_func ', I_func(y_c, Z_hat)

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



!  write(nw,*) 'calculated zeff'

 ! Z_hat = 4.*zeff/(5.*A_beam)

  ! critcal velocity
!  v_c = 4.5d6 * ( (te/1.d3) / 10)**0.5

!  y_c = v_c/v_beam


  write(nw,*) 'Completed nb contribution'

end subroutine nbicur



!Do it via flux surfaces

subroutine nbicur2()
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
  double precision :: psi,te, ne, dpsi, shine, depsum
  double precision :: tempe, densi, dense
  double precision :: v_c, y_c, Z_hat, I_func, E_c

  integer :: con, rb, zb, Z_b

  double precision, dimension(ncon) ::  volps, lambdas
  double precision :: shift, elong, kaps, kapps, dels, delps
  double precision, dimension(ncon) :: dpdrs

  double precision :: rr, zz, Ru

  double precision :: depp, depm, dep2

  integer :: k, ik,err
  double precision :: rat, del, delp, kap, kapp, dpdr, volp, lambda

  double precision :: rrange, zrange, beam_int, elec_return, uval
  double precision, dimension(ncon) :: n_f, h_tot, j_f, dep
  double precision, dimension(nr) :: J_tot


  double precision :: fprof, bphi, bth, fsi, B_tot, bp
  !See no. of cells in beam path
  integer :: count


  write(nw,*) 'Calculating Neutral beam current'



  n_f=0.
  count=0
  dep = 0.
  Z_b = 1

  !V beam (E_b in keV so change to J)
  v_beam = sqrt(2.*E_b*1000*eq/(2.*mp))

 !Calculate dpsi/drho and dvol/drho for each flux surface
  call dpsidrho(dpdrs)

  call dVdrho(volps)

  write(nw,*) 'calculated dVdrho'
  call lam(lambdas)


  write(nw,*) 'calculated lambdas'
    !write(nw,*) lambdas
!  if (icont .gt. -3) then


  J_nb(1) = 0
  ! Calculate shafranov shift and elongation for each flux surface
  do con = 2, ncon

     psi = psiv(con)
     rr = maxval(rpts(con,:))

     if (rr .lt. R_t) cycle
     lambda = lambdas(con)
     volp = volps(con)
     dpdr = dpdrs(con)

     kap = elong(con,0)
     kapp = elong(con,1) * dpdr
     del = shift(con,0)
     delp = shift(con,1) * dpdr

     rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)

     depp =  dep2(con, volp, lambda, kap, kapp, del, delp, 1)

     if ( rho**2 .ge. (Z_beam/kap)**2) then
        depm =  dep2(con, volp, lambda, kap, kapp, del, delp, -1)
     else
        depm = 0.
     end if

     dep(con) = depp + depm


     !write(nw,*) 'dep found', con, depp, depm


     ne = dense(psi,0)
     te = tempe(psi,0)

     coolog=log(sqrt(ne*1.0d-6)/te)
     coolog=24.-coolog
     tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)



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

     E_c = 1.2 * Z_b**(4./3.)  * 2*mp * te /((2.5*mp)**(2./3.) * me**(1./3.))
     v_c = sqrt(2.*E_c*eq/(2.*mp))
     y_c = v_c/v_beam


     n_f(con) = I_0/eq * dep(con) / vol

     xi_b = R_t/(rho+rcen)


     J_f(con) = eq * Z_b * n_f(con) * tau_s * xi_b * v_beam * I_func(y_c, Z_hat)

     uval = umax - psi

     J_nb(con) = elec_return(J_f(con), zeff, Z_b, uval) * bdl(con) / bsqav(con)

  end do

  J_tot = 0.



  depsum = 0.
!!! Need to see what happens at i=ncon (core)
  do i=2,ncon-1
     dpsi = (psiv(i-1) + psiv(i+1))/2.
     print*, i, dep(i), dpsi, volps(i)
     depsum = depsum + dep(i)*volps(i)*dpsi

  end do

  print*, 'depsum is ', depsum
  shine = 1 - depsum/vol

  dep = dep/(1.-shine)
  J_nb = J_nb/(1.-shine)



  ! do i=1,nr
 !    if (ixout(i,nsym) .ne. 1) cycle
 !    psi = umax - u(i,nsym)
 !    do j=1,ncon
 !       if (psi .ge. psiv(j)) then

 !          rat = (psi - psiv(j))/(psiv(j-1) - psiv(j))
 !          exit
 !       end if
 !    end do


 !       rr = r(i)
 !       zz = z(nsym)
 !       bth = bp(rr,zz)

 !       fsi = fprof(psi,2)
 !       bphi = fsi/rr
 !       B_tot = sqrt( bth*bth+bphi*bphi )

 !       J_tot(i) =  (J_nb(j) + rat*(J_nb(j-1) - J_nb(j)) ) * B_tot


  !   end do





!  end if
  write(nw,*) 'No. of cells in beam2 = ', count
  write(nw,*) 'exiting loop2'
  write(nw,*) 'Completed nb2 contribution'

end subroutine nbicur2






subroutine deposition(rho, con, rr, zz, dpdrs, volp,lambda, h)

  use param
  implicit none

  double precision :: rr, zz, dpdrs, volp, lambda, rho
  integer :: con

  double precision :: del, delp, kap, kapp, shift, elong

  double precision :: psi, ne, te, dense, tempe

  double precision :: coolog, tau_s

  double precision :: E_c, v_c, y_c, Z_b, v_beam, Z_hat

  double precision :: zeff, zni, densi
  integer :: l, gam_d

  double precision :: h, attenuation, D0, D1



  del = shift(con,0)
  delp = shift(con,1)

  kap = elong(con, 0)
  kapp = elong(con,1)

  psi = psiv(con)


  ne = dense(psi,0)
  te = tempe(psi,0)
  coolog=log(sqrt(ne*1.0d-6)/te)
  coolog=24.-coolog
  tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)

  v_beam = sqrt(2.*E_b*1000*eq/(2.*mp))


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

  kapp = kapp * dpdrs
  delp = delp * dpdrs
  volp = volp * dpdrs




  !    write(nw,*) rho, volp, kap, kapp, del, delp

  ! Deposition for pencil beam
  call dep(rho, volp, kap, kapp, del, delp, rr, zz, lambda, h)

  D0 = attenuation(rr,0)
  if (rpts(con,1) .ge. R_t) then
     gam_d = 1

     D1 = attenuation(rr,1)

  else
     gam_d = 0
     D1 = 0.
  end if



  h = h *( exp(-D0) + gam_d*(exp(-(D0+2*D1)))  )


end subroutine deposition







function dep2(con, volp, lambda, kap, kapp, del, delp, id)
  !calculates deposition profile from beams module
  !performs two integrals, over R for a fixed Z, and then integrates that
  !for all Z in the beam path

  use param
  implicit none

  double precision :: volp, lambda, kap, kapp, del, delp, dep2
  integer :: con, id, sgn

  double precision :: rr, rho, Zu, Rl, Ru, r_beam, rpm

  double precision :: wi, zti, zxi, rti, rxi

  double precision :: zterm, zsum, rterm, rsum, expterm, gauss

  double precision :: D0, D1, attenuation

  integer :: i, j, dpass


  !rho for each flux surface (outboard midplane for each flux sur
  rr = maxval(rpts(con,:))
  rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)


  r_beam = sig_r * 2

  Zu = min(rho*kap, r_beam)

  wi = pi/6
  zterm = 0.
  zsum = 0.


  !Z integration
  do i=1,6

     !Gaussian terms
     zti = cos((2*i -1)*pi/(2*6))
     zxi = (Zu/2) * (1+zti)

     !Calc R+ or R-
     if (id .eq. 1) then
        rpm = rcen + del + sqrt(rho**2 - (zxi/kap)**2)
        sgn = 1
     else if (id .eq. -1) then

        rpm = rcen + del - sqrt(rho**2 - (zxi/kap)**2)
        sgn = -1
     end if

     !Checks to see if R+- is below the tangency radius and if so cycles
     if (rpm .lt. (R_t - sqrt(r_beam**2-zxi**2) ) ) cycle


     !Terms in h(p) that involve Z
     zterm = (1 + kapp*zxi**2/(rho*kap**3))/sqrt(rho**2 - zxi**2/kap**2) + sgn*delp/rho

     !Multiplied by sqrt(1-ti^2)
     zterm = zterm * sqrt(1 - zti**2)


     !Limits for R integration
     Rl = R_t - sqrt(r_beam**2 - zxi**2)
     Ru = min(rpm, R_t + sqrt(r_beam**2 - zxi**2))


     !R integration
     rsum = 0.
     rterm = 0.
     do j=1,6

        !Gaussian terms for R
        rti = cos( (2*j - 1) * pi / (2*6) )
        rxi = (Ru-Rl)/(2.*rti)  + (Ru+Rl)/2

        !Checks to see if R+- is within integral
        if (rpm .lt. rxi ) cycle

        !Calc attenuation and checks for double pass
        D0 = attenuation(rpm,rxi,0)
        dpass = 0
        if (rxi .ge. R_t) then
           dpass = 1

           D1 = attenuation(rpm,rxi,1)
        end if

        !Total atten
        expterm = exp(-D0) + dpass*exp(-(D0+2*D1))

        !R terms in h(p)
        rterm = rpm/sqrt(rpm*2 - rxi**2) * expterm

!!! Check this again if numbers are weird
        !Gaussian for beam shape
        gauss = exp( - (zxi**2 + (R_t - rxi)**2)/sig_r**2) / (pi*sig_r**2 * (1 - exp(-r_beam**2/sig_r**2)) )




        !R term
        rterm = rterm * wi * gauss * sqrt(1 - rti**2)

        rsum = rsum + rterm
     end do

     !multiply by limits (b-a/2)
     rsum = rsum * (Ru-Rl)/2.

     !Multiplies R and Z contribution of integral and adds to Z total
     zsum = zsum + zterm*rsum

  end do

  zsum = zsum * wi * 2 *Zu/2.


  dep2 = 2 * rho * vol * zsum / (volp*lambda)

end function dep2




function attenuation(rpm, rb, id)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, attenuation
  integer :: id, i, con, simfac

  double precision :: Ru, D, dense, invlam

  double precision :: rmin, rmax, r_int, rat, dr_nb

  D = 0.



  if (id .eq. 0) then

     Ru = rcen + (amin**2-Z_beam**2/elon**2)**0.5


     dr_nb = (Ru-rpm)/100

     do i= 0,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rpm + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b)




        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))


     end do

  else if (id .eq. 1) then


     dr_nb = (rpm-rb)/100

     !stop atten go to inf
     do i= 1,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rb + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b)

        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))

     end do


     end if





     attenuation = D*dr_nb/3.


end function attenuation





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

  !Account for inner or outer part of flux surface
  if (rr .lt. r0) then
     h = h * ( (1+(kapp* Z_beam**2/(rho*kap**3)))/sqrt(rho**2 - (zz**2/kap**2)) + delp/rho)

  else

     h = h * ( (1+(kapp* Z_beam**2/(rho*kap**3)))/sqrt(rho**2 + (zz**2/kap**2)) + delp/rho)
  end if



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

          Ru = rcen + (amin**2-Z_beam**2/elon**2)**0.5

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


  function elec_return(J_f, zeff, Z_b,uval)
    !Calculates the current after the electron return
    ! current is accounted for

    use param
    implicit none

    double precision :: J_f, zeff, epsi, elec_return, uval
    double precision :: G, sqrt, J_nb
    integer :: Z_b

    call argen(uval, epsi,0)

    G = ( 1.55+(0.85/zeff) ) * sqrt(epsi)  - ( 0.2 + (1.55/zeff) )*epsi



    J_nb = (1 - ( Z_b/zeff * (1 - G) ) )


    J_nb = J_nb *J_f


    elec_return = J_nb

  end function elec_return


subroutine nbicur2()
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
  double precision :: psi,te, ne, dpsi, shine, depsum
  double precision :: tempe, densi, dense, avte, avne
  double precision :: v_c, y_c, Z_hat, I_func, E_c, I_0

  integer :: con, rb, zb, Z_b

  double precision, dimension(ncon) ::  volps, lambdas
  double precision :: shift, elong, kaps, kapps, dels, delps
  double precision, dimension(ncon) :: dpdrs

  double precision :: rr, zz, Ru, rhomax

  double precision :: depp, depm, dep2

  integer :: k, ik,err, beam, ecomp
  double precision :: rat, del, delp, kap, kapp, dpdr, volp, lambda

  double precision :: rrange, zrange, beam_int, elec_return, uval
  double precision, dimension(ncon) :: n_f, h_tot, j_f, dep, J_nbeam, rhos, rnormnb
  double precision, dimension(nr) :: J_tot


  double precision :: fprof, bphi, bth, fsi, B_tot, bp
  !See no. of cells in beam path
  integer :: count


  write(nw,*) 'Calculating Neutral beam current'

  vol = 0.
  do i=1,nr

     rr = r(i)
     do j=1,nz
        if (ixout(i,j) .le. 0) cycle

        vol = vol+rr*dr*dz

     end do
  end do

  vol = vol*2*pi

  !print*, vol
  n_f=0.
  count=0
  dep = 0.
  Z_b = 1


  !V beam (E_b in keV so change to J)

  !Calculate dpsi/drho and dvol/drho for each flux surface
  call dpsidrho(dpdrs, rnormnb)

  call dVdrho(volps)

  write(nw,*) 'calculated dVdrho'



  write(nw,*) 'calculated lambdas'
  !write(nw,*) lambdas
  !  if (icont .gt. -3) then

  J_nb= 0.
  J_nbeam = 0.


  do beam=1,1

     do ecomp = 1,3

        v_beam = sqrt(2.*E_b(beam)*1000*eq/(2.*mp * ecomp))
        call lam(lambdas, beam, ecomp)
        ! Calculate shafranov shift and elongation for each flux surface
        do con = 2, ncon

           psi = psiv(con)
           rr = maxval(rpts(con,:))

           if (rr .lt. R_t(beam)) cycle
           lambda = lambdas(con)
           volp = volps(con)
           dpdr = dpdrs(con)

           kap = elong(con,0)
           kapp = elong(con,1) * dpdr
           del = shift(con,0)
           delp = shift(con,1) * dpdr

           rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)
           rhos(con) = rho

           depp =  dep2(con, volp, lambda, kap, kapp, del, delp, 1, beam, ecomp)


           xi_b = R_t(beam)/(rho+rcen)

           !Pitch angle is different for inner and outer sections
           depp = depp * xi_b

           if ( rho**2 .ge. (Z_beam/kap)**2) then
              depm =  dep2(con, volp, lambda, kap, kapp, del, delp, -1, beam, ecomp)


              !Inner pitch angle
              xi_b = R_t(beam)/(rcen- rho)

              depm = depm * xi_b
           else
              depm = 0.
           end if

           dep(con) = depp + depm


           !write(nw,*) 'dep found', con, depp, depm


           ne = dense(psi,0)
           te = tempe(psi,0)

           coolog=log(sqrt(ne*1.0d-6)/te)
           coolog=24.-coolog
           tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)



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


           E_c = 1.2 * Z_b**(4./3.)  * 2*mp * te /((2.5*mp)**(2./3.) * me**(1./3.))
           v_c = sqrt(2.*E_c*eq/(2.*mp))
           y_c = v_c/v_beam

           I_0 = P_beam(beam)*P_frac(ecomp) * 1000 /(E_b(beam)/ecomp)

           n_f(con) = I_0/eq * dep(con) / vol



           J_f(con) = eq * Z_b * n_f(con) * tau_s  * v_beam * I_func(y_c, Z_hat)

           uval = umax - psi

           J_nbeam(con) = elec_return(J_f(con), zeff, Z_b, uval) * bdl(con) / bsqav(con)

           write(nw,*) con, beam, ecomp
           print*, dep(con), n_f(con)

           !dense(psi,0)*1e-14*pi*amin/cur



           !if (con .ge. 10 .and. beam .eq.1) then
           !   print*, 'Flux surface ', con
           !   print*, J_nbeam(con), dep(con)
           !end if
!!!!! Not sure if this is correct !!!!!!!!
           !J_nbeam(con) = J_nbeam(con)/rho


           !if(con.eq.49 .or. con .eq. 12) then
           !   print*, 'Dep and N_f ',dep(con), n_f(con)
           !end if

           !print*, con , J_nb(con)/rho, J_nb(con)



        end do




        J_tot = 0.



        depsum = 0.
!!! Need to see what happens at i=ncon (core)
        do i=2,ncon-1
           dpsi = (psiv(i-1) - psiv(i+1))/2.

           depsum = depsum + dep(i)*volps(i)*dpsi
           !print*, 'Dep sum', dep(i), volps(i), dpsi, depsum
           !print*, 'Dep is ', i, dep(i), J_nbeam(i)
        end do

        print*, 'Vol is ',vol
        print*, 'depsum is ', depsum
        shine = 1 - depsum/vol
        print*, 'shinethrough is ',beam,ecomp, shine
        dep = dep/(1.-shine)
        J_nbeam = J_nbeam/(1.-shine)



        !print*, 'Beam number ',beam
        !do i=1,ncon
        !   print*, i, J_f(i), J_nbeam(i)
        !end do
        !print*, ' '

        J_nb = J_nb + J_nbeam


     end do

  end do

  rat = (J_nb(ncon-1) - J_nb(ncon-2))/(psiv(ncon-1) - psiv(ncon-2))
  J_nb(ncon) = J_nb(ncon-1) + rat*(psiv(ncon) - psiv(ncon-1))

  !print*, J_nb

  ! do i=1,nr
  !    if (ixout(i,nsym) .ne. 1) cycle
  !    psi = umax - u(i,nsym)
  !    do j=1,ncon
  !       if (psi .ge. psiv(j)) then

  !          rat = (psi - psiv(j))/(psiv(j-1) - psiv(j))
  !          exit
  !       end if
  !    end do


  !       rr = r(i)
  !       zz = z(nsym)
  !       bth = bp(rr,zz)

  !       fsi = fprof(psi,2)
  !       bphi = fsi/rr
  !       B_tot = sqrt( bth*bth+bphi*bphi )

  !       J_tot(i) =  (J_nb(j) + rat*(J_nb(j-1) - J_nb(j)) ) * B_tot


  !   end do



  !  end if
  write(nw,*) 'exiting loop2'
  write(nw,*) 'Completed nb2 contribution'

  avte =0.
  avne =0.
  count = 0
  do i=1,nr
     do j=1,nz
        if (ixout(i,j) .gt. 0) then
           count = count +1
           psi = umax - u(i,j)
           avte = avte + tempe(psi,0)
           avne = avne + dense(psi,0)

        end if
     end do
  end do

  avte = avte/count
  avne = avne/count

  print*, 'Average Temp', avte
  print*, 'Average Density',avne


end subroutine nbicur2









function dep2(con, volp, lambda, kap, kapp, del, delp, id,beam, ecomp)
  !calculates deposition profile from beams module
  !performs two integrals, over R for a fixed Z, and then integrates that
  !for all Z in the beam path

  use param
  implicit none

  double precision :: volp, lambda, kap, kapp, del, delp, dep2
  integer :: con, id, sgn, beam, ecomp

  double precision :: rr, rho, Zu, Rl, Ru, r_beam, rpm

  double precision :: wi, zti, zxi, rti, rxi

  double precision :: zterm, zsum, rterm, rsum, expterm, gauss

  double precision :: D0, D1, attenuation

  integer :: i, j, dpass


  !rho for each flux surface (outboard midplane for each flux sur
  rr = maxval(rpts(con,:))
  rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)


  r_beam = sig_r * 2

  Zu = min(rho*kap, r_beam)

  wi = pi/6
  zterm = 0.
  zsum = 0.


  !Z integration
  do i=1,6

     !Gaussian terms
     zti = cos((2*i -1)*pi/(2*6))
     zxi = (Zu/2) * (1+zti)

     !Calc R+ or R-
     if (id .eq. 1) then
        rpm = rcen + del + sqrt(rho**2 - (zxi/kap)**2)
        sgn = 1
     else if (id .eq. -1) then

        rpm = rcen + del - sqrt(rho**2 - (zxi/kap)**2)
        sgn = -1
     end if

     !Checks to see if R+- is below the tangency radius and if so cycles
     if (rpm .lt. (R_t(beam) - sqrt(r_beam**2-zxi**2) ) ) cycle


     !Terms in h(p) that involve Z
     zterm = (1 + kapp*zxi**2/(rho*kap**3))/sqrt(rho**2 - zxi**2/kap**2) + sgn*delp/rho

     !Multiplied by sqrt(1-ti^2)
     zterm = zterm * sqrt(1 - zti**2)


     !Limits for R integration
     Rl = R_t(beam) - sqrt(r_beam**2 - zxi**2)
     Ru = min(rpm, R_t(beam) + sqrt(r_beam**2 - zxi**2))
     !print*, 'Rl, Ru', Rl, Ru

     !R integration
     rsum = 0.
     rterm = 0.
     do j=1,6

        !Gaussian terms for R
        rti = cos( (2*j - 1) * pi / (2*6) )
        rxi = (Ru-Rl)*rti/2  + (Ru+Rl)/2

        !print*, i,j, rti,rxi

        !Checks to see if R+- is within integral
        !if (rpm .lt. rxi ) cycle

        !Calc attenuation and checks for double pass
        D0 = attenuation(rpm,rxi,0,beam, ecomp)
        dpass = 0
        if (rxi .ge. R_t(beam)) then
           dpass = 1

           D1 = attenuation(rpm,rxi,1,beam,ecomp)
        end if

        !Total atten
        expterm = exp(-D0) + dpass*exp(-(D0+2*D1))

        !if(con .eq. 25) then
        !   print*, con, expterm
        !end if


        !R terms in h(p)
        rterm = rpm/sqrt(rpm*2 - rxi**2) * expterm

!!! Check this again if numbers are weird
        !Gaussian for beam shape
        gauss = exp( - (zxi**2 + (R_t(beam) - rxi)**2)/sig_r**2) / (pi*sig_r**2 * (1 - exp(-r_beam**2/sig_r**2)) )



        !print*, D0, D1
        !R term
        rterm = rterm * wi * gauss * sqrt(1 - rti**2)

        rsum = rsum + rterm
     end do

     !multiply by limits (b-a/2)
     rsum = rsum * (Ru-Rl)/2.

     !Multiplies R and Z contribution of integral and adds to Z total
     zsum = zsum + zterm*rsum

  end do

  zsum = zsum * wi * 2 *Zu/2.
  !print*, rho, zsum

  dep2 = 2 * rho * vol * zsum / (volp*lambda)
  !if (beam .eq. 1)  print*, con, dep2, zsum
  !if (con .eq. 49 .or. con .eq. 12) then
  !   print*, 'Volp ',volp,rho, zsum, lambda, con
  !end if


end function dep2




function attenuation(rpm, rb, id, beam, ecomp)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, attenuation
  integer :: id, i, con, simfac, beam, ecomp

  double precision :: Ru, D, dense, invlam

  double precision :: rmin, rmax, r_int, rat, dr_nb

  D = 0.



  if (id .eq. 0) then

     Ru = rcen + (amin**2-Z_beam**2/elon**2)**0.5


     dr_nb = (Ru-rpm)/100

     do i= 0,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rpm + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b(beam)/ecomp)



        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))


     end do

  else if (id .eq. 1) then


     dr_nb = (rpm-rb)/100

     !stop atten go to inf
     do i= 0,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rb + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b(beam)/ecomp)

        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))

     end do


     end if





     attenuation = D*dr_nb/3.
     if (con .eq.1) then
        print*, 'atten ',attenuation
     end if

end function attenuation


function atten2(rpm, rb, id, beam)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, atten2
  integer :: id, i, con, simfac, beam

  atten2 =0.


end function atten2



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






function elec_return(J_f, zeff, Z_b,uval)
  !Calculates the current after the electron return
  ! current is accounted for

  use param
  implicit none

  double precision :: J_f, zeff, epsi, elec_return, uval
  double precision :: G, sqrt, Jnb
  integer :: Z_b

  call argen(uval, epsi,0)

  G = ( 1.55+(0.85/zeff) ) * sqrt(epsi)  - ( 0.2 + (1.55/zeff) )*epsi



  Jnb = (1 - ( Z_b/zeff * (1 - G) ) )
subroutine nbicur2()
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
  double precision :: psi,te, ne, dpsi, shine, depsum
  double precision :: tempe, densi, dense, avte, avne
  double precision :: v_c, y_c, Z_hat, I_func, E_c, I_0

  integer :: con, rb, zb, Z_b

  double precision, dimension(ncon) ::  volps, lambdas
  double precision :: shift, elong, kaps, kapps, dels, delps
  double precision, dimension(ncon) :: dpdrs

  double precision :: rr, zz, Ru, rhomax

  double precision :: depp, depm, dep2

  integer :: k, ik,err, beam, ecomp
  double precision :: rat, del, delp, kap, kapp, dpdr, volp, lambda

  double precision :: rrange, zrange, beam_int, elec_return, uval
  double precision, dimension(ncon) :: n_f, h_tot, j_f, dep, J_nbeam, rhos, rnormnb
  double precision, dimension(nr) :: J_tot


  double precision :: fprof, bphi, bth, fsi, B_tot, bp
  !See no. of cells in beam path
  integer :: count


  write(nw,*) 'Calculating Neutral beam current'

  vol = 0.
  do i=1,nr

     rr = r(i)
     do j=1,nz
        if (ixout(i,j) .le. 0) cycle

        vol = vol+rr*dr*dz

     end do
  end do

  vol = vol*2*pi

  !print*, vol
  n_f=0.
  count=0
  dep = 0.
  Z_b = 1


  !V beam (E_b in keV so change to J)

  !Calculate dpsi/drho and dvol/drho for each flux surface
  call dpsidrho(dpdrs, rnormnb)

  call dVdrho(volps)

  write(nw,*) 'calculated dVdrho'



  write(nw,*) 'calculated lambdas'
  !write(nw,*) lambdas
  !  if (icont .gt. -3) then

  J_nb= 0.
  J_nbeam = 0.


  do beam=1,1

     do ecomp = 1,3

        v_beam = sqrt(2.*E_b(beam)*1000*eq/(2.*mp * ecomp))
        call lam(lambdas, beam, ecomp)
        ! Calculate shafranov shift and elongation for each flux surface
        do con = 2, ncon

           psi = psiv(con)
           rr = maxval(rpts(con,:))

           if (rr .lt. R_t(beam)) cycle
           lambda = lambdas(con)
           volp = volps(con)
           dpdr = dpdrs(con)

           kap = elong(con,0)
           kapp = elong(con,1) * dpdr
           del = shift(con,0)
           delp = shift(con,1) * dpdr

           rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)
           rhos(con) = rho

           depp =  dep2(con, volp, lambda, kap, kapp, del, delp, 1, beam, ecomp)


           xi_b = R_t(beam)/(rho+rcen)

           !Pitch angle is different for inner and outer sections
           depp = depp * xi_b

           if ( rho**2 .ge. (Z_beam/kap)**2) then
              depm =  dep2(con, volp, lambda, kap, kapp, del, delp, -1, beam, ecomp)


              !Inner pitch angle
              xi_b = R_t(beam)/(rcen- rho)

              depm = depm * xi_b
           else
              depm = 0.
           end if

           dep(con) = depp + depm


           !write(nw,*) 'dep found', con, depp, depm


           ne = dense(psi,0)
           te = tempe(psi,0)

           coolog=log(sqrt(ne*1.0d-6)/te)
           coolog=24.-coolog
           tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)



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


           E_c = 1.2 * Z_b**(4./3.)  * 2*mp * te /((2.5*mp)**(2./3.) * me**(1./3.))
           v_c = sqrt(2.*E_c*eq/(2.*mp))
           y_c = v_c/v_beam

           I_0 = P_beam(beam)*P_frac(ecomp) * 1000 /(E_b(beam)/ecomp)

           n_f(con) = I_0/eq * dep(con) / vol



           J_f(con) = eq * Z_b * n_f(con) * tau_s  * v_beam * I_func(y_c, Z_hat)

           uval = umax - psi

           J_nbeam(con) = elec_return(J_f(con), zeff, Z_b, uval) * bdl(con) / bsqav(con)

           write(nw,*) con, beam, ecomp
           print*, dep(con), n_f(con)

           !dense(psi,0)*1e-14*pi*amin/cur



           !if (con .ge. 10 .and. beam .eq.1) then
           !   print*, 'Flux surface ', con
           !   print*, J_nbeam(con), dep(con)
           !end if
!!!!! Not sure if this is correct !!!!!!!!
           !J_nbeam(con) = J_nbeam(con)/rho


           !if(con.eq.49 .or. con .eq. 12) then
           !   print*, 'Dep and N_f ',dep(con), n_f(con)
           !end if

           !print*, con , J_nb(con)/rho, J_nb(con)



        end do




        J_tot = 0.



        depsum = 0.
!!! Need to see what happens at i=ncon (core)
        do i=2,ncon-1
           dpsi = (psiv(i-1) - psiv(i+1))/2.

           depsum = depsum + dep(i)*volps(i)*dpsi
           !print*, 'Dep sum', dep(i), volps(i), dpsi, depsum
           !print*, 'Dep is ', i, dep(i), J_nbeam(i)
        end do

        print*, 'Vol is ',vol
        print*, 'depsum is ', depsum
        shine = 1 - depsum/vol
        print*, 'shinethrough is ',beam,ecomp, shine
        dep = dep/(1.-shine)
        J_nbeam = J_nbeam/(1.-shine)



        !print*, 'Beam number ',beam
        !do i=1,ncon
        !   print*, i, J_f(i), J_nbeam(i)
        !end do
        !print*, ' '

        J_nb = J_nb + J_nbeam


     end do

  end do

  rat = (J_nb(ncon-1) - J_nb(ncon-2))/(psiv(ncon-1) - psiv(ncon-2))
  J_nb(ncon) = J_nb(ncon-1) + rat*(psiv(ncon) - psiv(ncon-1))

  !print*, J_nb

  ! do i=1,nr
  !    if (ixout(i,nsym) .ne. 1) cycle
  !    psi = umax - u(i,nsym)
  !    do j=1,ncon
  !       if (psi .ge. psiv(j)) then

  !          rat = (psi - psiv(j))/(psiv(j-1) - psiv(j))
  !          exit
  !       end if
  !    end do


  !       rr = r(i)
  !       zz = z(nsym)
  !       bth = bp(rr,zz)

  !       fsi = fprof(psi,2)
  !       bphi = fsi/rr
  !       B_tot = sqrt( bth*bth+bphi*bphi )

  !       J_tot(i) =  (J_nb(j) + rat*(J_nb(j-1) - J_nb(j)) ) * B_tot


  !   end do



  !  end if
  write(nw,*) 'exiting loop2'
  write(nw,*) 'Completed nb2 contribution'

  avte =0.
  avne =0.
  count = 0
  do i=1,nr
     do j=1,nz
        if (ixout(i,j) .gt. 0) then
           count = count +1
           psi = umax - u(i,j)
           avte = avte + tempe(psi,0)
           avne = avne + dense(psi,0)

        end if
     end do
  end do

  avte = avte/count
  avne = avne/count

  print*, 'Average Temp', avte
  print*, 'Average Density',avne


end subroutine nbicur2









function dep2(con, volp, lambda, kap, kapp, del, delp, id,beam, ecomp)
  !calculates deposition profile from beams module
  !performs two integrals, over R for a fixed Z, and then integrates that
  !for all Z in the beam path

  use param
  implicit none

  double precision :: volp, lambda, kap, kapp, del, delp, dep2
  integer :: con, id, sgn, beam, ecomp

  double precision :: rr, rho, Zu, Rl, Ru, r_beam, rpm

  double precision :: wi, zti, zxi, rti, rxi

  double precision :: zterm, zsum, rterm, rsum, expterm, gauss

  double precision :: D0, D1, attenuation

  integer :: i, j, dpass


  !rho for each flux surface (outboard midplane for each flux sur
  rr = maxval(rpts(con,:))
  rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)


  r_beam = sig_r * 2

  Zu = min(rho*kap, r_beam)

  wi = pi/6
  zterm = 0.
  zsum = 0.


  !Z integration
  do i=1,6

     !Gaussian terms
     zti = cos((2*i -1)*pi/(2*6))
     zxi = (Zu/2) * (1+zti)

     !Calc R+ or R-
     if (id .eq. 1) then
        rpm = rcen + del + sqrt(rho**2 - (zxi/kap)**2)
        sgn = 1
     else if (id .eq. -1) then

        rpm = rcen + del - sqrt(rho**2 - (zxi/kap)**2)
        sgn = -1
     end if

     !Checks to see if R+- is below the tangency radius and if so cycles
     if (rpm .lt. (R_t(beam) - sqrt(r_beam**2-zxi**2) ) ) cycle


     !Terms in h(p) that involve Z
     zterm = (1 + kapp*zxi**2/(rho*kap**3))/sqrt(rho**2 - zxi**2/kap**2) + sgn*delp/rho

     !Multiplied by sqrt(1-ti^2)
     zterm = zterm * sqrt(1 - zti**2)


     !Limits for R integration
     Rl = R_t(beam) - sqrt(r_beam**2 - zxi**2)
     Ru = min(rpm, R_t(beam) + sqrt(r_beam**2 - zxi**2))
     !print*, 'Rl, Ru', Rl, Ru

     !R integration
     rsum = 0.
     rterm = 0.
     do j=1,6

        !Gaussian terms for R
        rti = cos( (2*j - 1) * pi / (2*6) )
        rxi = (Ru-Rl)*rti/2  + (Ru+Rl)/2

        !print*, i,j, rti,rxi

        !Checks to see if R+- is within integral
        !if (rpm .lt. rxi ) cycle

        !Calc attenuation and checks for double pass
        D0 = attenuation(rpm,rxi,0,beam, ecomp)
        dpass = 0
        if (rxi .ge. R_t(beam)) then
           dpass = 1

           D1 = attenuation(rpm,rxi,1,beam,ecomp)
        end if

        !Total atten
        expterm = exp(-D0) + dpass*exp(-(D0+2*D1))

        !if(con .eq. 25) then
        !   print*, con, expterm
        !end if


        !R terms in h(p)
        rterm = rpm/sqrt(rpm*2 - rxi**2) * expterm

!!! Check this again if numbers are weird
        !Gaussian for beam shape
        gauss = exp( - (zxi**2 + (R_t(beam) - rxi)**2)/sig_r**2) / (pi*sig_r**2 * (1 - exp(-r_beam**2/sig_r**2)) )



        !print*, D0, D1
        !R term
        rterm = rterm * wi * gauss * sqrt(1 - rti**2)

        rsum = rsum + rterm
     end do

     !multiply by limits (b-a/2)
     rsum = rsum * (Ru-Rl)/2.

     !Multiplies R and Z contribution of integral and adds to Z total
     zsum = zsum + zterm*rsum

  end do

  zsum = zsum * wi * 2 *Zu/2.
  !print*, rho, zsum

  dep2 = 2 * rho * vol * zsum / (volp*lambda)
  !if (beam .eq. 1)  print*, con, dep2, zsum
  !if (con .eq. 49 .or. con .eq. 12) then
  !   print*, 'Volp ',volp,rho, zsum, lambda, con
  !end if


end function dep2




function attenuation(rpm, rb, id, beam, ecomp)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, attenuation
  integer :: id, i, con, simfac, beam, ecomp

  double precision :: Ru, D, dense, invlam

  double precision :: rmin, rmax, r_int, rat, dr_nb

  D = 0.



  if (id .eq. 0) then

     Ru = rcen + (amin**2-Z_beam**2/elon**2)**0.5


     dr_nb = (Ru-rpm)/100

     do i= 0,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rpm + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b(beam)/ecomp)



        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))


     end do

  else if (id .eq. 1) then


     dr_nb = (rpm-rb)/100

     !stop atten go to inf
     do i= 0,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rb + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b(beam)/ecomp)

        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))

     end do


     end if





     attenuation = D*dr_nb/3.
     if (con .eq.1) then
        print*, 'atten ',attenuation
     end if

end function attenuation


function atten2(rpm, rb, id, beam)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, atten2
  integer :: id, i, con, simfac, beam

  atten2 =0.


end function atten2



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






function elec_return(J_f, zeff, Z_b,uval)
  !Calculates the current after the electron return
  ! current is accounted for

  use param
  implicit none

  double precision :: J_f, zeff, epsi, elec_return, uval
  double precision :: G, sqrt, Jnb
  integer :: Z_b

  call argen(uval, epsi,0)

  G = ( 1.55+(0.85/zeff) ) * sqrt(epsi)  - ( 0.2 + (1.55/zeff) )*epsi



  Jnb = (1 - ( Z_b/zeff * (1 - G) ) )
subroutine nbicur2()
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
  double precision :: psi,te, ne, dpsi, shine, depsum
  double precision :: tempe, densi, dense, avte, avne
  double precision :: v_c, y_c, Z_hat, I_func, E_c, I_0

  integer :: con, rb, zb, Z_b

  double precision, dimension(ncon) ::  volps, lambdas
  double precision :: shift, elong, kaps, kapps, dels, delps
  double precision, dimension(ncon) :: dpdrs

  double precision :: rr, zz, Ru, rhomax

  double precision :: depp, depm, dep2

  integer :: k, ik,err, beam, ecomp
  double precision :: rat, del, delp, kap, kapp, dpdr, volp, lambda

  double precision :: rrange, zrange, beam_int, elec_return, uval
  double precision, dimension(ncon) :: n_f, h_tot, j_f, dep, J_nbeam, rhos, rnormnb
  double precision, dimension(nr) :: J_tot


  double precision :: fprof, bphi, bth, fsi, B_tot, bp
  !See no. of cells in beam path
  integer :: count


  write(nw,*) 'Calculating Neutral beam current'

  vol = 0.
  do i=1,nr

     rr = r(i)
     do j=1,nz
        if (ixout(i,j) .le. 0) cycle

        vol = vol+rr*dr*dz

     end do
  end do

  vol = vol*2*pi

  !print*, vol
  n_f=0.
  count=0
  dep = 0.
  Z_b = 1


  !V beam (E_b in keV so change to J)

  !Calculate dpsi/drho and dvol/drho for each flux surface
  call dpsidrho(dpdrs, rnormnb)

  call dVdrho(volps)

  write(nw,*) 'calculated dVdrho'



  write(nw,*) 'calculated lambdas'
  !write(nw,*) lambdas
  !  if (icont .gt. -3) then

  J_nb= 0.
  J_nbeam = 0.


  do beam=1,1

     do ecomp = 1,3

        v_beam = sqrt(2.*E_b(beam)*1000*eq/(2.*mp * ecomp))
        call lam(lambdas, beam, ecomp)
        ! Calculate shafranov shift and elongation for each flux surface
        do con = 2, ncon

           psi = psiv(con)
           rr = maxval(rpts(con,:))

           if (rr .lt. R_t(beam)) cycle
           lambda = lambdas(con)
           volp = volps(con)
           dpdr = dpdrs(con)

           kap = elong(con,0)
           kapp = elong(con,1) * dpdr
           del = shift(con,0)
           delp = shift(con,1) * dpdr

           rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)
           rhos(con) = rho

           depp =  dep2(con, volp, lambda, kap, kapp, del, delp, 1, beam, ecomp)


           xi_b = R_t(beam)/(rho+rcen)

           !Pitch angle is different for inner and outer sections
           depp = depp * xi_b

           if ( rho**2 .ge. (Z_beam/kap)**2) then
              depm =  dep2(con, volp, lambda, kap, kapp, del, delp, -1, beam, ecomp)


              !Inner pitch angle
              xi_b = R_t(beam)/(rcen- rho)

              depm = depm * xi_b
           else
              depm = 0.
           end if

           dep(con) = depp + depm


           !write(nw,*) 'dep found', con, depp, depm


           ne = dense(psi,0)
           te = tempe(psi,0)

           coolog=log(sqrt(ne*1.0d-6)/te)
           coolog=24.-coolog
           tau_s = 6.27e8 * A_beam*te**1.5/(Z_b*ne*1.0d-6*coolog)



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


           E_c = 1.2 * Z_b**(4./3.)  * 2*mp * te /((2.5*mp)**(2./3.) * me**(1./3.))
           v_c = sqrt(2.*E_c*eq/(2.*mp))
           y_c = v_c/v_beam

           I_0 = P_beam(beam)*P_frac(ecomp) * 1000 /(E_b(beam)/ecomp)

           n_f(con) = I_0/eq * dep(con) / vol



           J_f(con) = eq * Z_b * n_f(con) * tau_s  * v_beam * I_func(y_c, Z_hat)

           uval = umax - psi

           J_nbeam(con) = elec_return(J_f(con), zeff, Z_b, uval) * bdl(con) / bsqav(con)

           write(nw,*) con, beam, ecomp
           print*, dep(con), n_f(con)

           !dense(psi,0)*1e-14*pi*amin/cur



           !if (con .ge. 10 .and. beam .eq.1) then
           !   print*, 'Flux surface ', con
           !   print*, J_nbeam(con), dep(con)
           !end if
!!!!! Not sure if this is correct !!!!!!!!
           !J_nbeam(con) = J_nbeam(con)/rho


           !if(con.eq.49 .or. con .eq. 12) then
           !   print*, 'Dep and N_f ',dep(con), n_f(con)
           !end if

           !print*, con , J_nb(con)/rho, J_nb(con)



        end do




        J_tot = 0.



        depsum = 0.
!!! Need to see what happens at i=ncon (core)
        do i=2,ncon-1
           dpsi = (psiv(i-1) - psiv(i+1))/2.

           depsum = depsum + dep(i)*volps(i)*dpsi
           !print*, 'Dep sum', dep(i), volps(i), dpsi, depsum
           !print*, 'Dep is ', i, dep(i), J_nbeam(i)
        end do

        print*, 'Vol is ',vol
        print*, 'depsum is ', depsum
        shine = 1 - depsum/vol
        print*, 'shinethrough is ',beam,ecomp, shine
        dep = dep/(1.-shine)
        J_nbeam = J_nbeam/(1.-shine)



        !print*, 'Beam number ',beam
        !do i=1,ncon
        !   print*, i, J_f(i), J_nbeam(i)
        !end do
        !print*, ' '

        J_nb = J_nb + J_nbeam


     end do

  end do

  rat = (J_nb(ncon-1) - J_nb(ncon-2))/(psiv(ncon-1) - psiv(ncon-2))
  J_nb(ncon) = J_nb(ncon-1) + rat*(psiv(ncon) - psiv(ncon-1))

  !print*, J_nb

  ! do i=1,nr
  !    if (ixout(i,nsym) .ne. 1) cycle
  !    psi = umax - u(i,nsym)
  !    do j=1,ncon
  !       if (psi .ge. psiv(j)) then

  !          rat = (psi - psiv(j))/(psiv(j-1) - psiv(j))
  !          exit
  !       end if
  !    end do


  !       rr = r(i)
  !       zz = z(nsym)
  !       bth = bp(rr,zz)

  !       fsi = fprof(psi,2)
  !       bphi = fsi/rr
  !       B_tot = sqrt( bth*bth+bphi*bphi )

  !       J_tot(i) =  (J_nb(j) + rat*(J_nb(j-1) - J_nb(j)) ) * B_tot


  !   end do



  !  end if
  write(nw,*) 'exiting loop2'
  write(nw,*) 'Completed nb2 contribution'

  avte =0.
  avne =0.
  count = 0
  do i=1,nr
     do j=1,nz
        if (ixout(i,j) .gt. 0) then
           count = count +1
           psi = umax - u(i,j)
           avte = avte + tempe(psi,0)
           avne = avne + dense(psi,0)

        end if
     end do
  end do

  avte = avte/count
  avne = avne/count

  print*, 'Average Temp', avte
  print*, 'Average Density',avne


end subroutine nbicur2









function dep2(con, volp, lambda, kap, kapp, del, delp, id,beam, ecomp)
  !calculates deposition profile from beams module
  !performs two integrals, over R for a fixed Z, and then integrates that
  !for all Z in the beam path

  use param
  implicit none

  double precision :: volp, lambda, kap, kapp, del, delp, dep2
  integer :: con, id, sgn, beam, ecomp

  double precision :: rr, rho, Zu, Rl, Ru, r_beam, rpm

  double precision :: wi, zti, zxi, rti, rxi

  double precision :: zterm, zsum, rterm, rsum, expterm, gauss

  double precision :: D0, D1, attenuation

  integer :: i, j, dpass


  !rho for each flux surface (outboard midplane for each flux sur
  rr = maxval(rpts(con,:))
  rho = sqrt( (rr - (rcen + del))**2 + (Z_beam/kap)**2)


  r_beam = sig_r * 2

  Zu = min(rho*kap, r_beam)

  wi = pi/6
  zterm = 0.
  zsum = 0.


  !Z integration
  do i=1,6

     !Gaussian terms
     zti = cos((2*i -1)*pi/(2*6))
     zxi = (Zu/2) * (1+zti)

     !Calc R+ or R-
     if (id .eq. 1) then
        rpm = rcen + del + sqrt(rho**2 - (zxi/kap)**2)
        sgn = 1
     else if (id .eq. -1) then

        rpm = rcen + del - sqrt(rho**2 - (zxi/kap)**2)
        sgn = -1
     end if

     !Checks to see if R+- is below the tangency radius and if so cycles
     if (rpm .lt. (R_t(beam) - sqrt(r_beam**2-zxi**2) ) ) cycle


     !Terms in h(p) that involve Z
     zterm = (1 + kapp*zxi**2/(rho*kap**3))/sqrt(rho**2 - zxi**2/kap**2) + sgn*delp/rho

     !Multiplied by sqrt(1-ti^2)
     zterm = zterm * sqrt(1 - zti**2)


     !Limits for R integration
     Rl = R_t(beam) - sqrt(r_beam**2 - zxi**2)
     Ru = min(rpm, R_t(beam) + sqrt(r_beam**2 - zxi**2))
     !print*, 'Rl, Ru', Rl, Ru

     !R integration
     rsum = 0.
     rterm = 0.
     do j=1,6

        !Gaussian terms for R
        rti = cos( (2*j - 1) * pi / (2*6) )
        rxi = (Ru-Rl)*rti/2  + (Ru+Rl)/2

        !print*, i,j, rti,rxi

        !Checks to see if R+- is within integral
        !if (rpm .lt. rxi ) cycle

        !Calc attenuation and checks for double pass
        D0 = attenuation(rpm,rxi,0,beam, ecomp)
        dpass = 0
        if (rxi .ge. R_t(beam)) then
           dpass = 1

           D1 = attenuation(rpm,rxi,1,beam,ecomp)
        end if

        !Total atten
        expterm = exp(-D0) + dpass*exp(-(D0+2*D1))

        !if(con .eq. 25) then
        !   print*, con, expterm
        !end if


        !R terms in h(p)
        rterm = rpm/sqrt(rpm*2 - rxi**2) * expterm

!!! Check this again if numbers are weird
        !Gaussian for beam shape
        gauss = exp( - (zxi**2 + (R_t(beam) - rxi)**2)/sig_r**2) / (pi*sig_r**2 * (1 - exp(-r_beam**2/sig_r**2)) )



        !print*, D0, D1
        !R term
        rterm = rterm * wi * gauss * sqrt(1 - rti**2)

        rsum = rsum + rterm
     end do

     !multiply by limits (b-a/2)
     rsum = rsum * (Ru-Rl)/2.

     !Multiplies R and Z contribution of integral and adds to Z total
     zsum = zsum + zterm*rsum

  end do

  zsum = zsum * wi * 2 *Zu/2.
  !print*, rho, zsum

  dep2 = 2 * rho * vol * zsum / (volp*lambda)
  !if (beam .eq. 1)  print*, con, dep2, zsum
  !if (con .eq. 49 .or. con .eq. 12) then
  !   print*, 'Volp ',volp,rho, zsum, lambda, con
  !end if


end function dep2




function attenuation(rpm, rb, id, beam, ecomp)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, attenuation
  integer :: id, i, con, simfac, beam, ecomp

  double precision :: Ru, D, dense, invlam

  double precision :: rmin, rmax, r_int, rat, dr_nb

  D = 0.



  if (id .eq. 0) then

     Ru = rcen + (amin**2-Z_beam**2/elon**2)**0.5


     dr_nb = (Ru-rpm)/100

     do i= 0,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rpm + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b(beam)/ecomp)



        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))


     end do

  else if (id .eq. 1) then


     dr_nb = (rpm-rb)/100

     !stop atten go to inf
     do i= 0,100

        if (i.eq.0 .or. i.eq.100) then
           simfac=1
        else if (mod(i,2) .eq. 0) then
           simfac=2
        else
           simfac=4
        end if

        r_int = rb + i*dr_nb

        if (r_int .ge. r0) then
           do con=1,ncon

              rmax = rpts(con, int(npts/2))

              if (rmax .gt. r_int) cycle

              rat = (r_int - rmax)/( rpts(con-1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con-1) - psiv(con))

              exit
           end do

        else
           do con=ncon,1,-1

              rmin = rpts(con,1)

              if (rmin .gt. r_int) cycle

              rat = (r_int - rmin)/( rpts(con+1, int(npts/2)) - rmax)
              psi = psiv(con) + rat*(psiv(con+1) - psiv(con))

              exit
           end do
        end if

        invlam = dense(psi,0)/(2.8e17*E_b(beam)/ecomp)

        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))

     end do


     end if





     attenuation = D*dr_nb/3.
     if (con .eq.1) then
        print*, 'atten ',attenuation
     end if

end function attenuation


function atten2(rpm, rb, id, beam)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, atten2
  integer :: id, i, con, simfac, beam

  atten2 =0.


end function atten2



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






function elec_return(J_f, zeff, Z_b,uval)
  !Calculates the current after the electron return
  ! current is accounted for

  use param
  implicit none

  double precision :: J_f, zeff, epsi, elec_return, uval
  double precision :: G, sqrt, Jnb
  integer :: Z_b

  call argen(uval, epsi,0)

  G = ( 1.55+(0.85/zeff) ) * sqrt(epsi)  - ( 0.2 + (1.55/zeff) )*epsi



  Jnb = (1 - ( Z_b/zeff * (1 - G) ) )



 Jnb = Jnb *J_f


  elec_return = Jnb

end function elec_return
