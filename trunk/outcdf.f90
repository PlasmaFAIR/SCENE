subroutine write_netcdf()
  ! File ot write SCENE output to NETCDF format
  ! Initially focused for writing gs2 output

  use param
  use netcdf
  implicit none

  logical :: debug
  character(len=lrunname+4) :: cdf_file
  character(len=4) :: file_suffix
  integer :: ncid
  integer, parameter ::  oned=1, twod=2
  integer :: dimids0d(1), dimidsrho(oned), dimids2d(twod), dimidsR(oned)

  integer :: con, i
  double precision :: rat

  !Variable ID for different inputs
  integer :: rhopsi_varid, rhopsi_dimid
  integer :: rr_varid, rr_dimid, zz_varid, zz_dimid
  integer :: npts_varid, npts_dimid

  ! 0D profiles
  integer :: zerod_varid, zerod_dimid
  integer :: bcentr_varid, amin_varid

  ! 0D profile units
  character (len=*), parameter :: zerod_unit=''
  character (len=*), parameter :: bcentr_unit='T',amin_unit='m'

  ! 1d profile variable ID
  integer :: psi_varid
  integer :: te_varid, ti_varid, ne_varid, ni_varid
  integer :: Lte_varid, Lti_varid, Lne_varid, Lni_varid
  integer :: dshift_varid, shat_varid, pk_varid, eps_varid
  integer :: gs2_beta_prime_varid, zeff_varid, vnui_varid, vnue_varid
  integer :: sfac_varid, jtot_varid, jext_varid, jext2_varid
  integer :: jbs_varid, jdi_varid, jps_varid, jnbi_varid
  integer :: jdotb_varid, jbsdotb_varid, jnbdotb_varid
  integer :: fsi_varid, ffp_varid, pp_varid, dpsidrho_varid
  integer :: j_tor_varid, jbs_tor_varid, jnb_tor_varid, jex_tor_varid

  ! TGLF  1d profile variable ID 
  integer :: tglf_rmin_varid, tglf_rmaj_varid, tglf_q_varid, tglf_q_prime_varid
  integer :: tglf_p_prime_varid, tglf_drmajdx_varid, tglf_kappa_varid, tglf_s_kappa_varid
  integer :: tglf_delta_varid, tglf_s_delta_varid
  integer :: tglf_zs_i_varid, tglf_mass_i_varid, tglf_rlns_i_varid, tglf_rlts_i_varid
  integer :: tglf_taus_i_varid, tglf_as_i_varid, tglf_vpar_i_varid, tglf_vpar_shear_i_varid
  integer :: tglf_zs_e_varid, tglf_mass_e_varid, tglf_rlns_e_varid, tglf_rlts_e_varid
  integer :: tglf_taus_e_varid, tglf_as_e_varid, tglf_vpar_e_varid, tglf_vpar_shear_e_varid
  integer :: tglf_b_unit_varid
  
  
  ! 1d profiles Units
  character (len=*), parameter :: psi_unit="Wb/m2"
  character (len=*), parameter :: rhopsi_unit="", units="units"
  character (len=*), parameter :: rr_unit='m', zz_unit='m'
  character (len=*), parameter :: te_unit="eV", ti_unit='eV', ne_unit='m^-3', ni_unit='m^-3'
  character (len=*), parameter :: Lte_unit="", Lti_unit='', Lne_unit='', Lni_unit=''
  character (len=*), parameter :: dshift_unit='', shat_unit='', pk_unit='', eps_unit=''
  character (len=*), parameter :: gs2_beta_prime_unit='', zeff_unit='',vnui_unit='', vnue_unit=''
  character (len=*), parameter :: sfac_unit='', jtot_unit='kA/m2' ,jext_unit='kA/m2', jext2_unit='kA/m2'
  character (len=*), parameter :: jbs_unit='kA/m2', jdi_unit='kA/m2', jps_unit='kA/m2', jnbi_unit='kA/m2'
  character (len=*), parameter :: jdotb_unit='A/m2', jbsdotb_unit='A/m2', jnbdotb_unit='A/m2'
  character (len=*), parameter :: j_tor_unit='A/m2', jbs_tor_unit='A/m2', jnb_tor_unit='A/m2', jex_tor_unit='A/m2'
  character (len=*), parameter :: fsi_unit='', ffp_unit='', pp_unit='', dpsidrho_unit='Wb/m3'


  ! TGLF  1d profile variable ID 
  character (len=*), parameter :: tglf_rmin_unit='', tglf_rmaj_unit='', tglf_q_unit='', tglf_q_prime_unit=''
  character (len=*), parameter :: tglf_p_prime_unit='', tglf_drmajdx_unit='', tglf_kappa_unit='', tglf_s_kappa_unit=''
  character (len=*), parameter :: tglf_delta_unit='', tglf_s_delta_unit=''
  character (len=*), parameter :: tglf_zs_i_unit='', tglf_mass_i_unit='', tglf_rlns_i_unit='', tglf_rlts_i_unit=''
  character (len=*), parameter :: tglf_taus_i_unit='', tglf_as_i_unit='', tglf_vpar_i_unit='', tglf_vpar_shear_i_unit=''
  character (len=*), parameter :: tglf_zs_e_unit='', tglf_mass_e_unit='', tglf_rlns_e_unit='', tglf_rlts_e_unit=''
  character (len=*), parameter :: tglf_taus_e_unit='', tglf_as_e_unit='', tglf_vpar_e_unit='', tglf_vpar_shear_e_unit=''
  character (len=*), parameter :: tglf_b_unit_unit=''

  
  ! 1d Output data
  double precision :: psi, tempe, tempi, dense, densi, press, shift, fprof
  double precision, dimension(ncon) :: rhopsi, te, ti, ne, ni
  double precision, dimension(ncon) :: Lte, Lti, Lne, Lni
  double precision, dimension(ncon) :: dshift, shat, pk, eps
  double precision, dimension(ncon) :: gs2_beta_prime, zeff, vnui, vnue
  double precision, dimension(ncon) :: jdotb, jbsdotb, jnbdotb
  double precision, dimension(ncon) :: jbs_tor, jnb_tor, j_tor, jex_tor
  double precision, dimension(ncon) :: fsi, ffp, pp, bsmean
  
  double precision :: ne19, zni19, coolog, bcentr,colli,colle


  double precision :: tglf_zs_i, tglf_zs_e, tglf_mass_i, tglf_mass_e
  
  !TGLF 1D output data
  double precision, dimension(ncon) :: tglf_rmin, tglf_rmaj, tglf_q, tglf_q_prime
  double precision, dimension(ncon) :: tglf_p_prime, tglf_drmajdx, tglf_kappa, tglf_s_kappa
  double precision, dimension(ncon) :: tglf_delta, tglf_s_delta
  double precision, dimension(ncon) :: tglf_rlns_i, tglf_rlts_i
  double precision, dimension(ncon) :: tglf_taus_i, tglf_as_i, tglf_vpar_i, tglf_vpar_shear_i
  double precision, dimension(ncon) :: tglf_rlns_e, tglf_rlts_e
  double precision, dimension(ncon) :: tglf_taus_e, tglf_as_e, tglf_vpar_e, tglf_vpar_shear_e, tglf_b_unit


  ! 2d  variable ID
  integer :: rpts_varid, zpts_varid, bppts_varid

  ! 2d units
  character (len=*), parameter :: npts_unit=''
  character (len=*), parameter :: rpts_unit='m', zpts_unit='m', bppts_unit='T'

  
  double precision, dimension(ncon) :: dPsidrhos, rhos, dPsiNdrho
  integer, dimension(npts) :: pts
  double precision :: shafr, dshafr, rmaj, rmin
  double precision :: kappa, dkappa, p_prime, dPsidr
  double precision :: elong, tglf_shear, delta, ddelta, triang

  
  debug = .True.

  file_suffix = '.cdf'

  cdf_file = runname(1:lrunname)//file_suffix
  if (debug) print*, 'Saving to ',cdf_file

  ! Create (Overwrite) netcdf file
  call check(nf90_create(cdf_file, NF90_CLOBBER,ncid))
  write(nw,*) 'Create NetCDF Output'

  bcentr = mu0*rodi/(2.*pi*rcen)

  call dpsidrho(dPsidrhos,rhos)


  
  !Generate input data
  do con=1, ncon
     psi = psiv(con)
     rhopsi(con) = sqrt(psi/psiv(1))
     !Kinetic profiles
     te(con) = tempe(psi,0)
     ti(con) = tempi(psi,1,0)
     ne(con) = dense(psi,0)
     ni(con) = densi(psi,1,0)

     !Kinetic gradients
     Lte(con) = -tempe(psi,1)*umax/te(con)
     Lti(con) = -tempi(psi,1,1)*umax/ti(con)
     Lne(con) = -dense(psi,1)*umax/ne(con)
     Lni(con) = -densi(psi,1,1)*umax/ni(con)

     !Other GS2 inputs
     dshift(con) =shift(con,1)*umax/amin
     eps(con) = epsv(con)
     pk(con) = 2*amin/(rcen*sfac(con))
     shat(con) = qp(con)*psi/sfac(con)
     gs2_beta_prime(con) = press(psi,1)*8.0d-7*pi*umax/bcentr**2

     fsi(con) = fprof(psi,2)
     ffp(con) = -fprof(psi,1)
     pp(con) = -press(psi,1)
     bsmean(con) = bsqav(con)
     jdotb(con) = fsi(con)*pp(con)/bsmean(con) + ffp(con)/(fsi(con)*mu0)
     jbsdotb(con) = bsj(con)/sqrt(bsmean(con))
     jnbdotb(con) = J_nb(con)

     jnb_tor(con) = jnbdotb(con) * bphiav(con)
     jbs_tor(con) = jbsdotb(con)*bphiav(con)
     j_tor(con) = jdotb(con) * bphiav(con)
     
     !Zeff profile
     zeff(con) = zm
     if (imp.eq.1) then
        ne19 = ne(con)*1.e-19
        if (ne19.gt.0.) then
           zeff(con)=0.
           do i=1,nimp+1
              zni19=densi(psi,i,0)*1.0d-19
              zeff(con)=zeff(con)+(zni19*iz(i)**2)/ne19
           end do
        end if
     end if

     !Collisionality

     !Ions
     coolog=log(sqrt(ni(con)*1.0d-6)/ti(con))
     coolog=24.-coolog
     
     !From Wesson 2.15 (added by bhavin 21/03/18)
     ! colli = 6.6d17*zmas(1)**0.5 * (ti(con)/1000)**1.5/(ni(con)*iz(1)**4*coolog)
     
     colli = sqrt(2.0)*pi*ni(con) * iz(1)**4 * eq**4 * coolog &
          / ( sqrt(zmai*mp) * (ti(con)*bk)**1.5  * (4*pi*eps0)**2 )

     vnui(con) =colli * amin / sqrt(2*te(con)*bk/(zmai*mp)) 

     !Collisionality for electrons
     coolog=log(sqrt(ne(con)*1.0d-6)/te(con))
     coolog=24.-coolog

     !colle = 3.*(2.*pi)**1.5* eps0**2 * me**0.5 * (tempe(psi,0)*bk)**1.5 &
     !     / ( dense(psi,0) * eq**4 * coolog)
     
     colle = sqrt(2.0)*pi*ne(con) * eq**4 * coolog &
     	/ (sqrt(me) * (te(con)*bk)**1.5 * (4*pi*eps0)**2 )
     
     vnue(con) = colle * amin / sqrt(2*te(con)*bk/(zmai*mp))


     ! TGLF parameters
     ! Shafranov shift/elongation and derivative in r
     !Turn dpsi/drho into dPsi/dr
     dPsidr = dPsidrhos(con)/amin

     shafr = shift(con,0)
     dshafr = shift(con,1)*dPsidr
     kappa = elong(con,0)
     dkappa = elong(con,1)*dPsidr
     delta = triang(con,0)
     ddelta = triang(con,1)*dPsidr

     !Major radius 
     rmaj = rcen + shafr
     rmin = maxval(rpts(con,:)) - rmaj

     if (con .eq. ncon) then
        rat = (tglf_rmin(ncon-1) -tglf_rmin(ncon-2))/(psiv(ncon-1) - psiv(ncon-2))
        rmin = tglf_rmin(ncon-1) + rat*psiv(ncon-1)
     end if
     
     p_prime = press(psi, 1)*dpsidr

     
     tglf_rmin(con) = rmin/amin
     tglf_rmaj(con) = rmaj/amin
     tglf_q(con) = sfac(con)
     tglf_shear = rmin * qp(con) * dpsidr / tglf_q(con)
     tglf_q_prime(con) = tglf_q(con)**2 * tglf_shear / tglf_rmin(con)**2


     tglf_drmajdx(con) = dshafr

     tglf_kappa(con) = kappa
     tglf_s_kappa(con) = rmin*dkappa/kappa

     tglf_delta(con) = delta
     tglf_s_delta(con) = rmin*ddelta/sqrt(1-delta**2)

     
     tglf_b_unit(con) = tglf_q(con)*dpsidr/(rmin*2.0*pi)

     tglf_p_prime(con) = tglf_q(con) * amin**2 * p_prime * 2 * mu0 / (rmin*tglf_b_unit(con)**2)

  end do

  jex_tor = j_tor - jnb_tor - jbs_tor

  tglf_zs_i = 1.0
  tglf_zs_e = -1.0
  tglf_mass_i = 1.0
  tglf_mass_e = me/(zmas(1)*mp)
  tglf_vpar_i = 0.0
  tglf_vpar_e = 0.0
  tglf_vpar_shear_i = 0.0
  tglf_vpar_shear_e = 0.0

  !Covert GS2s dT/dPsiN to dT/dr
  ! dT/dr = dT/dPsiN * dPsi/drho / (a*umax)
  dPsiNdrho = dPsidrhos/umax
  
  tglf_rlns_i = Lni*dPsiNdrho
  tglf_rlns_e = Lne*dPsiNdrho
  tglf_rlts_i = Lti*dPsiNdrho
  tglf_rlts_e = Lte*dPsiNdrho

  tglf_as_i = ni/ne
  tglf_taus_i = ti/te

  tglf_as_e = 1.0
  tglf_taus_e = 1.0
  
  tglf_q_prime(ncon) = 0.0
  tglf_p_prime(ncon) = 0.0

  
  !Define Radial dimension
  call check( nf90_def_dim(ncid, "rho_psi", ncon, rhopsi_dimid) )
  !Define Co-ordinate variable and its units
  call check( nf90_def_var(ncid, "rho_psi", NF90_REAL, rhopsi_dimid, rhopsi_varid))
  call check( nf90_put_att(ncid, rhopsi_varid, units, rhopsi_unit))

  call check( nf90_def_dim(ncid, "npts", npts, npts_dimid))
  call check( nf90_def_var(ncid, "npts", NF90_REAL, npts_dimid, npts_varid))
  call check( nf90_put_att(ncid, npts_varid, units, npts_unit))
  
  call check( nf90_def_dim(ncid, "R", nr, rr_dimid))
  call check( nf90_def_var(ncid, "R", NF90_REAL, rr_dimid, rr_varid))
  call check( nf90_put_att(ncid, rr_varid, units, rr_unit))
  
  
  !Define "units" for 0D data
  call check( nf90_def_dim(ncid, "ZeroD", 1, zerod_dimid) )
  call check( nf90_def_var(ncid, "ZeroD", NF90_INT, zerod_dimid, zerod_varid))
  call check( nf90_put_att(ncid, zerod_varid, units, zerod_unit))
  if (debug) print*, 'Assigned variable names'
  
  
  dimids0d = (/ zerod_dimid /)
  dimidsrho = (/rhopsi_dimid/)
  dimidsR = (/rr_varid/)

  dimids2d = (/rhopsi_varid, npts_varid/)

  if (debug) print*, 'Assigned dimensions'
  !Assign 0d Variables:
  call check( nf90_def_var(ncid, "bcentr", NF90_REAL, dimids0d, bcentr_varid) ) 
  call check( nf90_def_var(ncid, "amin", NF90_REAL, dimids0d, amin_varid) )

  if (debug) print*, 'Assigned 0D variables'
  
  !Assign 1d Variables:
  call check( nf90_def_var(ncid, "Psi", NF90_REAL, dimidsrho, psi_varid) ) 
  call check( nf90_def_var(ncid, "Te", NF90_REAL, dimidsrho, te_varid) )
  call check( nf90_def_var(ncid, "Ti", NF90_REAL, dimidsrho, ti_varid) )
  call check( nf90_def_var(ncid, "Ne", NF90_REAL, dimidsrho, ne_varid) )
  call check( nf90_def_var(ncid, "Ni", NF90_REAL, dimidsrho, ni_varid) )
  call check( nf90_def_var(ncid, "LTe", NF90_REAL, dimidsrho, Lte_varid) )
  call check( nf90_def_var(ncid, "LTi", NF90_REAL, dimidsrho, Lti_varid) )
  call check( nf90_def_var(ncid, "Lne", NF90_REAL, dimidsrho, Lne_varid) )
  call check( nf90_def_var(ncid, "Lni", NF90_REAL, dimidsrho, Lni_varid) )
  call check( nf90_def_var(ncid, "shift", NF90_REAL, dimidsrho, dshift_varid) )
  call check( nf90_def_var(ncid, "eps", NF90_REAL, dimidsrho, eps_varid) )
  call check( nf90_def_var(ncid, "pk", NF90_REAL, dimidsrho, pk_varid) )
  call check( nf90_def_var(ncid, "shat", NF90_REAL, dimidsrho, shat_varid) )
  call check( nf90_def_var(ncid, "gs2_beta_prime", NF90_REAL, dimidsrho, gs2_beta_prime_varid) )
  call check( nf90_def_var(ncid, "zeff", NF90_REAL, dimidsrho, zeff_varid) )
  call check( nf90_def_var(ncid, "vnui", NF90_REAL, dimidsrho, vnui_varid) )
  call check( nf90_def_var(ncid, "vnue", NF90_REAL, dimidsrho, vnue_varid) )
  call check( nf90_def_var(ncid, "Sfac", NF90_REAL, dimidsrho, sfac_varid) )
  call check( nf90_def_var(ncid, "Jtot.B", NF90_REAL, dimidsrho, jdotb_varid) )
  call check( nf90_def_var(ncid, "Jbs.B", NF90_REAL, dimidsrho, jbsdotb_varid) )
  call check( nf90_def_var(ncid, "Jnb.B", NF90_REAL, dimidsrho, jnbdotb_varid) )
  call check( nf90_def_var(ncid, "J_tor", NF90_REAL, dimidsrho, j_tor_varid) )
  call check( nf90_def_var(ncid, "Jbs_tor", NF90_REAL, dimidsrho, jbs_tor_varid) )
  call check( nf90_def_var(ncid, "Jnb_tor", NF90_REAL, dimidsrho, jnb_tor_varid) )
  call check( nf90_def_var(ncid, "Jex_tor", NF90_REAL, dimidsrho, jex_tor_varid) )
  call check( nf90_def_var(ncid, "Fpsi", NF90_REAL, dimidsrho, fsi_varid) )
  call check( nf90_def_var(ncid, "FFP", NF90_REAL, dimidsrho, ffp_varid) )
  call check( nf90_def_var(ncid, "Pp", NF90_REAL, dimidsrho, pp_varid) )
  call check( nf90_def_var(ncid, "dPsidrho", NF90_REAL, dimidsrho, dpsidrho_varid) )

  call check( nf90_def_var(ncid, "JTOT", NF90_REAL, dimidsR, jtot_varid) )
  call check( nf90_def_var(ncid, "JEXT", NF90_REAL, dimidsR, jext_varid) )
  call check( nf90_def_var(ncid, "JEXT2", NF90_REAL, dimidsR, jext2_varid) )
  call check( nf90_def_var(ncid, "JBS", NF90_REAL, dimidsR, jbs_varid) )
  call check( nf90_def_var(ncid, "JDIA", NF90_REAL, dimidsR, jdi_varid) )
  call check( nf90_def_var(ncid, "JPS", NF90_REAL, dimidsR, jps_varid) )
  call check( nf90_def_var(ncid, "JNBI", NF90_REAL, dimidsR, jnbi_varid) )



  !TGLF 1d variables assigment
  call check( nf90_def_var(ncid, "TGLF_ZS_i", NF90_REAL, dimids0d, tglf_zs_i_varid) )
  call check( nf90_def_var(ncid, "TGLF_MASS_i", NF90_REAL, dimids0d, tglf_mass_i_varid) )
  call check( nf90_def_var(ncid, "TGLF_ZS_e", NF90_REAL, dimids0d, tglf_zs_e_varid) )
  call check( nf90_def_var(ncid, "TGLF_MASS_e", NF90_REAL, dimids0d, tglf_mass_e_varid) )

  
  call check( nf90_def_var(ncid, "TGLF_B_UNIT", NF90_REAL, dimidsrho, tglf_b_unit_varid) )
   

  call check( nf90_def_var(ncid, "TGLF_RLNS_i", NF90_REAL, dimidsrho, tglf_rlns_i_varid) )
  call check( nf90_def_var(ncid, "TGLF_RLTS_i", NF90_REAL, dimidsrho, tglf_rlts_i_varid) )
  call check( nf90_def_var(ncid, "TGLF_TAUS_i", NF90_REAL, dimidsrho, tglf_taus_i_varid) )
  call check( nf90_def_var(ncid, "TGLF_AS_i", NF90_REAL, dimidsrho, tglf_as_i_varid) )
  call check( nf90_def_var(ncid, "TGLF_VPAR_i", NF90_REAL, dimidsrho, tglf_vpar_i_varid) )
  call check( nf90_def_var(ncid, "TGLF_VPAR_SHEAR_i", NF90_REAL, dimidsrho, tglf_vpar_shear_i_varid) )


  
  call check( nf90_def_var(ncid, "TGLF_RLNS_e", NF90_REAL, dimidsrho, tglf_rlns_e_varid) )
  call check( nf90_def_var(ncid, "TGLF_RLTS_e", NF90_REAL, dimidsrho, tglf_rlts_e_varid) )
  call check( nf90_def_var(ncid, "TGLF_TAUS_e", NF90_REAL, dimidsrho, tglf_taus_e_varid) )
  call check( nf90_def_var(ncid, "TGLF_AS_e", NF90_REAL, dimidsrho, tglf_as_e_varid) )
  call check( nf90_def_var(ncid, "TGLF_VPAR_e", NF90_REAL, dimidsrho, tglf_vpar_e_varid) )
  call check( nf90_def_var(ncid, "TGLF_VPAR_SHEAR_e", NF90_REAL, dimidsrho, tglf_vpar_shear_e_varid) )

  call check( nf90_def_var(ncid, "TGLF_RMIN", NF90_REAL, dimidsrho, tglf_rmin_varid) )
  call check( nf90_def_var(ncid, "TGLF_RMAJ", NF90_REAL, dimidsrho, tglf_rmaj_varid) )
  call check( nf90_def_var(ncid, "TGLF_Q", NF90_REAL, dimidsrho, tglf_q_varid) )
  call check( nf90_def_var(ncid, "TGLF_Q_PRIME", NF90_REAL, dimidsrho, tglf_q_prime_varid) )
  call check( nf90_def_var(ncid, "TGLF_P_PRIME", NF90_REAL, dimidsrho, tglf_p_prime_varid) )
  call check( nf90_def_var(ncid, "TGLF_DRMAJDX", NF90_REAL, dimidsrho, tglf_drmajdx_varid) )
  call check( nf90_def_var(ncid, "TGLF_KAPPA", NF90_REAL, dimidsrho, tglf_kappa_varid) )
  call check( nf90_def_var(ncid, "TGLF_S_KAPPA", NF90_REAL, dimidsrho, tglf_s_kappa_varid) )
  call check( nf90_def_var(ncid, "TGLF_DELTA", NF90_REAL, dimidsrho, tglf_delta_varid) )
  call check( nf90_def_var(ncid, "TGLF_S_DELTA", NF90_REAL, dimidsrho, tglf_s_delta_varid) )

  if (debug) print*, 'Defined TGLF variables'
  
  call check( nf90_def_var(ncid, "rpts", NF90_REAL, dimids2d, rpts_varid) )

  call check( nf90_def_var(ncid, "zpts", NF90_REAL, dimids2d, zpts_varid) )

  call check( nf90_def_var(ncid, "bppts", NF90_REAL, dimids2d, bppts_varid) )

    
  if (debug) print*, 'Defined variables'

  !Assign Units

  !Assign 0d units
  call check( nf90_put_att(ncid,bcentr_varid, units, bcentr_unit)) 
  call check( nf90_put_att(ncid,amin_varid, units, amin_unit)) 
  
  !Assign 1d units
  call check( nf90_put_att(ncid,psi_varid, units, psi_unit)) 
  call check( nf90_put_att(ncid,te_varid, units, te_unit))
  call check( nf90_put_att(ncid,ti_varid, units, ti_unit))
  call check( nf90_put_att(ncid,ne_varid, units, ne_unit))
  call check( nf90_put_att(ncid,ni_varid, units, ni_unit))
  call check( nf90_put_att(ncid,Lte_varid, units, Lte_unit))
  call check( nf90_put_att(ncid,Lti_varid, units, Lti_unit))
  call check( nf90_put_att(ncid,Lne_varid, units, Lne_unit))
  call check( nf90_put_att(ncid,Lni_varid, units, Lni_unit))
  call check( nf90_put_att(ncid,dshift_varid, units, dshift_unit))
  call check( nf90_put_att(ncid,eps_varid, units, eps_unit))
  call check( nf90_put_att(ncid,pk_varid, units, pk_unit))
  call check( nf90_put_att(ncid,shat_varid, units, shat_unit))
  call check( nf90_put_att(ncid,gs2_beta_prime_varid, units, gs2_beta_prime_unit))
  call check( nf90_put_att(ncid,zeff_varid, units, zeff_unit))
  call check( nf90_put_att(ncid,vnui_varid, units, vnui_unit))
  call check( nf90_put_att(ncid,vnue_varid, units, vnue_unit))  
  call check( nf90_put_att(ncid,sfac_varid, units, sfac_unit))
  call check( nf90_put_att(ncid,jdotb_varid, units, jdotb_unit))
  call check( nf90_put_att(ncid,jbsdotb_varid, units, jbsdotb_unit))
  call check( nf90_put_att(ncid,jnbdotb_varid, units, jnbdotb_unit))
  call check( nf90_put_att(ncid,j_tor_varid, units, j_tor_unit))
  call check( nf90_put_att(ncid,jbs_tor_varid, units, jbs_tor_unit))
  call check( nf90_put_att(ncid,jnb_tor_varid, units, jnb_tor_unit))
  call check( nf90_put_att(ncid,jex_tor_varid, units, jex_tor_unit))
  
  call check( nf90_put_att(ncid,jtot_varid, units, jtot_unit))
  call check( nf90_put_att(ncid,jext_varid, units, jext_unit))
  call check( nf90_put_att(ncid,jext2_varid, units, jext2_unit))
  call check( nf90_put_att(ncid,jbs_varid, units, jbs_unit))
  call check( nf90_put_att(ncid,jdi_varid, units, jdi_unit))
  call check( nf90_put_att(ncid,jps_varid, units, jps_unit))
  call check( nf90_put_att(ncid,jnbi_varid, units, jnbi_unit))
  call check( nf90_put_att(ncid,fsi_varid, units, fsi_unit))
  call check( nf90_put_att(ncid,ffp_varid, units, ffp_unit))
  call check( nf90_put_att(ncid,pp_varid, units, pp_unit))
  call check( nf90_put_att(ncid,dpsidrho_varid, units, dpsidrho_unit))

  !Put tglf parameters units
  call check( nf90_put_att(ncid, tglf_zs_i_varid, units, tglf_zs_i_unit) )
  call check( nf90_put_att(ncid, tglf_zs_e_varid, units, tglf_zs_e_unit) )
  call check( nf90_put_att(ncid, tglf_mass_i_varid, units, tglf_mass_i_unit) )
  call check( nf90_put_att(ncid, tglf_mass_e_varid, units, tglf_mass_e_unit) )

  call check( nf90_put_att(ncid, tglf_rlns_i_varid, units, tglf_rlns_i_unit) )
  call check( nf90_put_att(ncid, tglf_rlts_i_varid, units, tglf_rlts_i_unit) )
  call check( nf90_put_att(ncid, tglf_taus_i_varid, units, tglf_taus_i_unit) )
  call check( nf90_put_att(ncid, tglf_as_i_varid, units, tglf_as_i_unit) )
  call check( nf90_put_att(ncid, tglf_vpar_i_varid, units, tglf_vpar_i_unit) )
  call check( nf90_put_att(ncid, tglf_vpar_shear_i_varid, units, tglf_vpar_shear_i_unit) )

  call check( nf90_put_att(ncid, tglf_rlns_e_varid, units, tglf_rlns_e_unit) )
  call check( nf90_put_att(ncid, tglf_rlts_e_varid, units, tglf_rlts_e_unit) )
  call check( nf90_put_att(ncid, tglf_taus_e_varid, units, tglf_taus_e_unit) )
  call check( nf90_put_att(ncid, tglf_as_e_varid, units, tglf_as_e_unit) )
  call check( nf90_put_att(ncid, tglf_vpar_e_varid, units, tglf_vpar_e_unit) )
  call check( nf90_put_att(ncid, tglf_vpar_shear_e_varid, units, tglf_vpar_shear_e_unit) )


  call check( nf90_put_att(ncid, tglf_rmaj_varid, units, tglf_rmaj_unit) )
  call check( nf90_put_att(ncid, tglf_rmin_varid, units, tglf_rmin_unit) )
  call check( nf90_put_att(ncid, tglf_q_varid, units, tglf_q_unit) )
  call check( nf90_put_att(ncid, tglf_q_prime_varid, units, tglf_q_prime_unit) )
  call check( nf90_put_att(ncid, tglf_p_prime_varid, units, tglf_p_prime_unit) )
  call check( nf90_put_att(ncid, tglf_drmajdx_varid, units, tglf_drmajdx_unit) )
  call check( nf90_put_att(ncid, tglf_kappa_varid, units, tglf_kappa_unit) )
  call check( nf90_put_att(ncid, tglf_s_kappa_varid, units, tglf_s_kappa_unit))
  call check( nf90_put_att(ncid, tglf_delta_varid, units, tglf_delta_unit) )
  call check( nf90_put_att(ncid, tglf_s_delta_varid, units, tglf_s_delta_unit))
  call check( nf90_put_att(ncid, tglf_b_unit_varid, units, tglf_b_unit_unit))

  call check( nf90_put_att(ncid, rpts_varid, units, rpts_unit))
  call check( nf90_put_att(ncid, zpts_varid, units, zpts_unit))
  call check( nf90_put_att(ncid, bppts_varid, units, bppts_unit))
  
  if (debug) print*, 'Assigned variable attributes'
  
  call check(nf90_enddef(ncid))

  !Put in value for 0D 'co-ordinate'
  call check(nf90_put_var(ncid, zerod_dimid, 1))
  
  !Put Radial coordinates in
  call check(nf90_put_var(ncid, rhopsi_dimid, rhopsi))

  do i=1, npts
     pts(i) = i
  end do
  !Put poloidal points
  call check(nf90_put_var(ncid, npts_dimid, pts))
  
  !Put in 0d variable data
  call check(nf90_put_var(ncid, bcentr_varid, bcentr))
  call check(nf90_put_var(ncid, amin_varid, amin))
  call check(nf90_put_var(ncid, tglf_zs_i_varid,tglf_zs_i) )
  call check(nf90_put_var(ncid, tglf_zs_e_varid, tglf_zs_e) )
  call check(nf90_put_var(ncid, tglf_mass_i_varid, tglf_mass_i) )
  call check(nf90_put_var(ncid, tglf_mass_e_varid, tglf_mass_e) )

  
  !Put in 1d variable data
  call check(nf90_put_var(ncid, psi_varid, psiv))
  call check(nf90_put_var(ncid, te_varid, te))
  call check(nf90_put_var(ncid, ti_varid, ti)) 
  call check(nf90_put_var(ncid, ne_varid, ne))
  call check(nf90_put_var(ncid, ni_varid, ni))
  call check(nf90_put_var(ncid, Lte_varid, Lte))
  call check(nf90_put_var(ncid, Lti_varid, Lti))
  call check(nf90_put_var(ncid, Lne_varid, Lne))
  call check(nf90_put_var(ncid, Lni_varid, Lni))
  call check(nf90_put_var(ncid, dshift_varid, dshift))
  call check(nf90_put_var(ncid, eps_varid, eps))
  call check(nf90_put_var(ncid, pk_varid, pk))
  call check(nf90_put_var(ncid, shat_varid, shat))
  call check(nf90_put_var(ncid, gs2_beta_prime_varid, gs2_beta_prime))
  call check(nf90_put_var(ncid, zeff_varid, zeff))
  call check(nf90_put_var(ncid, vnui_varid, vnui)) 
  call check(nf90_put_var(ncid, vnue_varid, vnue))
  call check(nf90_put_var(ncid, sfac_varid, sfac)) 
  call check(nf90_put_var(ncid, jdotb_varid, jdotb)) 
  call check(nf90_put_var(ncid, jbsdotb_varid, jbsdotb)) 
  call check(nf90_put_var(ncid, jnbdotb_varid, jnbdotb))
  call check(nf90_put_var(ncid, j_tor_varid, j_tor)) 
  call check(nf90_put_var(ncid, jbs_tor_varid, jbs_tor)) 
  call check(nf90_put_var(ncid, jnb_tor_varid, jnb_tor))
  call check(nf90_put_var(ncid, jex_tor_varid, jex_tor))
  call check(nf90_put_var(ncid, fsi_varid, fsi)) 
  call check(nf90_put_var(ncid, ffp_varid, ffp)) 
  call check(nf90_put_var(ncid, pp_varid, pp))
  call check(nf90_put_var(ncid, dpsidrho_varid, dPsidrhos)) 

  call check(nf90_put_var(ncid, jtot_varid, gradj(:,nsym)/1000) )
  call check(nf90_put_var(ncid, jext_varid, exph(:,nsym)/1000) )
  call check(nf90_put_var(ncid, jext2_varid, exph2(:,nsym)/1000)) 
  call check(nf90_put_var(ncid, jbs_varid, bsph(:,nsym)/1000) )
  call check(nf90_put_var(ncid, jdi_varid, diph(:,nsym)/1000) )
  call check(nf90_put_var(ncid, jps_varid, psph(:,nsym)/1000) )
  call check(nf90_put_var(ncid, jnbi_varid, nbph(:,nsym)/1000) )

  print*, 'Putting TGLF attrs'
  !Put 1D TGLF variables

  
  call check( nf90_put_var(ncid, tglf_rlns_i_varid, tglf_rlns_i) )
  call check( nf90_put_var(ncid, tglf_rlts_i_varid, tglf_rlts_i) )
  call check( nf90_put_var(ncid, tglf_taus_i_varid, tglf_taus_i) )
  call check( nf90_put_var(ncid, tglf_as_i_varid, tglf_as_i) )
  call check( nf90_put_var(ncid, tglf_vpar_i_varid, tglf_vpar_i) )
  call check( nf90_put_var(ncid, tglf_vpar_shear_i_varid, tglf_vpar_shear_i) )

  call check( nf90_put_var(ncid, tglf_rlns_e_varid, tglf_rlns_e) )
  call check( nf90_put_var(ncid, tglf_rlts_e_varid, tglf_rlts_e) )
  call check( nf90_put_var(ncid, tglf_taus_e_varid, tglf_taus_e) )
  call check( nf90_put_var(ncid, tglf_as_e_varid, tglf_as_e) )
  call check( nf90_put_var(ncid, tglf_vpar_e_varid, tglf_vpar_e) )
  call check( nf90_put_var(ncid, tglf_vpar_shear_e_varid, tglf_vpar_shear_e) )

  print*, 'Put TGLF species values'
  
  call check( nf90_put_var(ncid, tglf_rmaj_varid, tglf_rmaj) )
  call check( nf90_put_var(ncid, tglf_rmin_varid, tglf_rmin) )
  call check( nf90_put_var(ncid, tglf_q_varid, tglf_q) )
  call check( nf90_put_var(ncid, tglf_q_prime_varid, tglf_q_prime) )
  call check( nf90_put_var(ncid, tglf_p_prime_varid, tglf_p_prime) )
  call check( nf90_put_var(ncid, tglf_drmajdx_varid, tglf_drmajdx) )
  call check( nf90_put_var(ncid, tglf_kappa_varid, tglf_kappa) )
  call check( nf90_put_var(ncid, tglf_s_kappa_varid, tglf_s_kappa))
  call check( nf90_put_var(ncid, tglf_delta_varid, tglf_delta) )
  call check( nf90_put_var(ncid, tglf_s_delta_varid, tglf_s_delta))
  call check( nf90_put_var(ncid, tglf_b_unit_varid, tglf_b_unit))

  if (debug) print*, 'Put variables values in'

  ! 2d variables
  call check( nf90_put_var(ncid, rpts_varid, rpts ))
  call check( nf90_put_var(ncid, zpts_varid, zpts ))
  call check( nf90_put_var(ncid, bppts_varid, bppts ))

  
  
  call check(nf90_close(ncid))
  print *, "*** SUCCESS writing NETCDF file "

end subroutine write_netcdf


subroutine check(status)

  use netcdf
  integer, intent ( in) :: status
  
  if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop 2
  end if
end subroutine check
       
