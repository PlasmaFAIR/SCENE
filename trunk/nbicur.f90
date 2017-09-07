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

  double precision :: rr, zz, Ru, rhomax

  double precision :: depp, depm, dep2

  integer :: k, ik,err, beam
  double precision :: rat, del, delp, kap, kapp, dpdr, volp, lambda

  double precision :: rrange, zrange, beam_int, elec_return, uval
  double precision, dimension(ncon) :: n_f, h_tot, j_f, dep, J_nbeam, rhos
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
  call dpsidrho(dpdrs)

  call dVdrho(volps)

  write(nw,*) 'calculated dVdrho'



  write(nw,*) 'calculated lambdas'
  !write(nw,*) lambdas
  !  if (icont .gt. -3) then

  J_nb= 0.
  J_nbeam = 0.

  do beam=1,2
     v_beam = sqrt(2.*E_b(beam)*1000*eq/(2.*mp))
     call lam(lambdas, beam)
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

        depp =  dep2(con, volp, lambda, kap, kapp, del, delp, 1, beam)

        if ( rho**2 .ge. (Z_beam/kap)**2) then
           depm =  dep2(con, volp, lambda, kap, kapp, del, delp, -1, beam)
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

        n_f(con) = I_0(beam)/eq * dep(con) / vol

        xi_b = R_t(beam)/(rho+rcen)


        J_f(con) = eq * Z_b * n_f(con) * tau_s * xi_b * v_beam * I_func(y_c, Z_hat)

        uval = umax - psi

        J_nbeam(con) = elec_return(J_f(con), zeff, Z_b, uval) * bdl(con) / bsqav(con)

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
        dpsi = (psiv(i-1) + psiv(i+1))/2.
        !print*, i, dep(i), dpsi, volps(i)
        depsum = depsum + dep(i)*volps(i)*dpsi

     end do

     !print*, 'depsum is ', depsum
     shine = 1 - depsum/vol

     dep = dep/(1.-shine)
     J_nbeam = J_nbeam/(1.-shine)

     !do i=1,ncon
     !   print*, (ncon-i)*1./ncon, dep(i), dep(i)/rhos(i)
     !end do

     !print*, 'Beam number ',beam
     !do i=1,ncon
     !   print*, i, J_f(i), J_nbeam(i)
     !end do
     !print*, ' '
     
     J_nb = J_nb + J_nbeam




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

end subroutine nbicur2









function dep2(con, volp, lambda, kap, kapp, del, delp, id,beam)
  !calculates deposition profile from beams module
  !performs two integrals, over R for a fixed Z, and then integrates that
  !for all Z in the beam path

  use param
  implicit none

  double precision :: volp, lambda, kap, kapp, del, delp, dep2
  integer :: con, id, sgn, beam

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
        D0 = attenuation(rpm,rxi,0,beam)
        dpass = 0
        if (rxi .ge. R_t(beam)) then
           dpass = 1

           D1 = attenuation(rpm,rxi,1,beam)
        end if
        
        !Total atten
        expterm = exp(-D0) + dpass*exp(-(D0+2*D1))

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
  !print*, rho, vol, volp, lambda
  
  dep2 = 2 * rho * vol * zsum / (volp*lambda)

  !if (con .eq. 49 .or. con .eq. 12) then
  !   print*, 'Volp ',volp,rho, zsum, lambda, con
  !end if
  

end function dep2

  


function attenuation(rpm, rb, id, beam)

  use param
  implicit none

  double precision :: rpm, rb, psi, ne, attenuation
  integer :: id, i, con, simfac, beam

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

        invlam = dense(psi,0)/(2.8e17*E_b(beam))



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

        invlam = dense(psi,0)/(2.8e17*E_b(beam))
    
        D = D + simfac*(r_int *invlam/sqrt(r_int**2 - rb**2))

     end do
     

     end if
     


   
     
     attenuation = D*dr_nb/3.
     if (con .eq.1) then
        print*, 'atten ',attenuation
     end if
       
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



