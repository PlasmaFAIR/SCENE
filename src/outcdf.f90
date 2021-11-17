module netcdf_interface
  implicit none
contains

  !> Write SCENE output to netCDF file
  subroutine write_netcdf()
    use, intrinsic :: iso_fortran_env, only : real64
    use param
    use neasyf, only: neasyf_open, neasyf_close, neasyf_dim, neasyf_write
    use profiles_mod, only: tempe, tempi, dense, densi, press, shift, fprof, dpsidrho, dVdrho, rhotor, elong, triang
    implicit none

    logical :: debug
    character(len=lrunname + 4) :: cdf_file
    character(len=4) :: file_suffix
    integer :: ncid
    integer :: dimidsrho(1), dimids2d(2), dimidsR(1)

    integer :: con, i
    real(real64) :: rat

    !Variable ID for different inputs
    integer :: rhopsi_dimid
    integer :: rr_dimid
    integer :: npts_dimid

    ! 1d Output data
    real(real64) :: psi
    real(real64), dimension(ncon) :: rhopsi, te, ti, ne, ni, nHe, tHe
    real(real64), dimension(ncon) :: Lte, Lti, Lne, Lni
    real(real64), dimension(ncon) :: dshift, shat, pk, eps
    real(real64), dimension(ncon) :: gs2_beta_prime, zeff, vnui, vnue
    real(real64), dimension(ncon) :: jdotb, jbsdotb, jnbdotb
    real(real64), dimension(ncon) :: jbs_tor, jnb_tor, j_tor, jex_tor
    real(real64), dimension(ncon) :: fsi, ffp, pp, bsmean

    real(real64) :: ne19, zni19, coolog, bcentr, colli, colle

    real(real64) :: tglf_zs_i, tglf_zs_e, tglf_mass_i, tglf_mass_e

    !TGLF 1D output data
    real(real64), dimension(ncon) :: tglf_rmin, tglf_rmaj, tglf_q, tglf_q_prime
    real(real64), dimension(ncon) :: tglf_p_prime, tglf_drmajdx, tglf_kappa, tglf_s_kappa
    real(real64), dimension(ncon) :: tglf_delta, tglf_s_delta
    real(real64), dimension(ncon) :: tglf_rlns_i, tglf_rlts_i
    real(real64), dimension(ncon) :: tglf_taus_i, tglf_as_i, tglf_vpar_i, tglf_vpar_shear_i
    real(real64), dimension(ncon) :: tglf_rlns_e, tglf_rlts_e
    real(real64), dimension(ncon) :: tglf_taus_e, tglf_as_e, tglf_vpar_e, tglf_vpar_shear_e, tglf_b_unit

    real(real64), dimension(ncon) :: dPsidrhos, rhos, dPsiNdrho
    real(real64) :: shafr, dshafr, rmaj, rmin
    real(real64) :: kappa, dkappa, p_prime, dPsidr
    real(real64) :: tglf_shear, delta, ddelta

    debug = .False.

    file_suffix = '.cdf'

    cdf_file = runname(1:lrunname)//file_suffix
    if (debug) print *, 'Saving to ', cdf_file

    ncid = neasyf_open(cdf_file, "w")

    if (debug) write (nw, *) 'Create NetCDF Output'

    bcentr = mu0 * rodi / (2.*pi * rcen)

    call dpsidrho(dPsidrhos, rhos)

    !Generate input data
    do con = 1, ncon
      psi = psiv(con)
      rhopsi(con) = sqrt(psi / psiv(1))
      !Kinetic profiles
      te(con) = tempe(psi, 0)
      ti(con) = tempi(psi, 1, 0)
      ne(con) = dense(psi, 0)
      ni(con) = densi(psi, 1, 0)

      if (imp .eq. 1) then
        nHe(con) = densi(psi, 2, 0)
        tHe(con) = tempi(psi, 2, 0)
      end if

      !Kinetic gradients
      Lte(con) = -tempe(psi, 1) * umax / te(con)
      Lti(con) = -tempi(psi, 1, 1) * umax / ti(con)
      Lne(con) = -dense(psi, 1) * umax / ne(con)
      Lni(con) = -densi(psi, 1, 1) * umax / ni(con)

      !Other GS2 inputs
      dshift(con) = shift(con, 1) * umax / amin
      eps(con) = epsv(con)
      pk(con) = 2 * amin / (rcen * sfac(con))
      shat(con) = qp(con) * psi / sfac(con)
      gs2_beta_prime(con) = press(psi, 1) * 8.0d-7 * pi * umax / bcentr**2

      fsi(con) = fprof(psi, 2)
      ffp(con) = -fprof(psi, 1)
      pp(con) = -press(psi, 1)
      bsmean(con) = bsqav(con)
      jdotb(con) = fsi(con) * pp(con) / bsmean(con) + ffp(con) / (fsi(con) * mu0)
      jbsdotb(con) = bsj(con) / sqrt(bsmean(con))
      jnbdotb(con) = J_nb(con)
      bphiav(con) = fsi(con) * rsqinv(con) / rinv(con)
      jnb_tor(con) = jnbdotb(con) * bphiav(con)
      jbs_tor(con) = jbsdotb(con) * bphiav(con)
      j_tor(con) = jdotb(con) * bphiav(con)
      jex_tor(con) = (ffp(con) / (fsi(con) * mu0) - jbsdotb(con) + fsi(con) * pp(con) / bsmean(con)) * bphiav(con)

      !Zeff profile
      zeff(con) = zm
      if (imp .eq. 1) then
        ne19 = ne(con) * 1.e-19
        if (ne19 .gt. 0.) then
          zeff(con) = 0.
          do i = 1, nimp + 1
            zni19 = densi(psi, i, 0) * 1.0d-19
            zeff(con) = zeff(con) + (zni19 * iz(i)**2) / ne19
          end do
        end if
      end if

      !Collisionality

      !Ions
      coolog = log(sqrt(ni(con) * 1.0d-6) / ti(con))
      coolog = 24.-coolog

      !From Wesson 2.15 (added by bhavin 21/03/18)
      ! colli = 6.6d17*zmas(1)**0.5 * (ti(con)/1000)**1.5/(ni(con)*iz(1)**4*coolog)

      colli = sqrt(2.0) * pi * ni(con) * iz(1)**4 * eq**4 * coolog &
              / (sqrt(zmai * mp) * (ti(con) * bk)**1.5 * (4 * pi * eps0)**2)

      vnui(con) = colli * amin / sqrt(2 * te(con) * bk / (zmai * mp))

      !Collisionality for electrons
      coolog = log(sqrt(ne(con) * 1.0d-6) / te(con))
      coolog = 24.-coolog

      !colle = 3.*(2.*pi)**1.5* eps0**2 * me**0.5 * (tempe(psi,0)*bk)**1.5 &
      !     / ( dense(psi,0) * eq**4 * coolog)

      colle = sqrt(2.0) * pi * ne(con) * eq**4 * coolog &
              / (sqrt(me) * (te(con) * bk)**1.5 * (4 * pi * eps0)**2)

      vnue(con) = colle * amin / sqrt(2 * te(con) * bk / (zmai * mp))

      ! TGLF parameters
      ! Shafranov shift/elongation and derivative in r
      !Turn dpsi/drho into dPsi/dr
      dPsidr = dPsidrhos(con) / amin

      shafr = shift(con, 0)
      dshafr = shift(con, 1) * dPsidr
      kappa = elong(con, 0)
      dkappa = elong(con, 1) * dPsidr
      delta = triang(con, 0)
      ddelta = triang(con, 1) * dPsidr

      !Major radius
      rmaj = rcen + shafr
      rmin = maxval(rpts(con, :)) - rmaj

      if (con .eq. ncon) then
        rat = (tglf_rmin(ncon - 1) - tglf_rmin(ncon - 2)) / (psiv(ncon - 1) - psiv(ncon - 2))
        rmin = tglf_rmin(ncon - 1) + rat * psiv(ncon - 1)
      end if

      p_prime = press(psi, 1) * dpsidr

      tglf_rmin(con) = rmin / amin
      tglf_rmaj(con) = rmaj / amin
      tglf_q(con) = sfac(con)
      tglf_shear = rmin * qp(con) * dpsidr / tglf_q(con)
      tglf_q_prime(con) = tglf_q(con)**2 * tglf_shear / tglf_rmin(con)**2

      tglf_drmajdx(con) = dshafr

      tglf_kappa(con) = kappa
      tglf_s_kappa(con) = rmin * dkappa / kappa

      tglf_delta(con) = delta
      tglf_s_delta(con) = rmin * ddelta / sqrt(1 - delta**2)

      tglf_b_unit(con) = tglf_q(con) * dpsidr / (rmin * 2.0 * pi)

      tglf_p_prime(con) = tglf_q(con) * amin**2 * p_prime * 2 * mu0 / (rmin * tglf_b_unit(con)**2)

    end do

    tglf_zs_i = 1.0
    tglf_zs_e = -1.0
    tglf_mass_i = 1.0
    tglf_mass_e = me / (zmas(1) * mp)
    tglf_vpar_i = 0.0
    tglf_vpar_e = 0.0
    tglf_vpar_shear_i = 0.0
    tglf_vpar_shear_e = 0.0

    !Covert GS2s dT/dPsiN to dT/dr
    ! dT/dr = dT/dPsiN * dPsi/drho / (a*umax)
    dPsiNdrho = dPsidrhos / umax

    tglf_rlns_i = Lni * dPsiNdrho
    tglf_rlns_e = Lne * dPsiNdrho
    tglf_rlts_i = Lti * dPsiNdrho
    tglf_rlts_e = Lte * dPsiNdrho

    tglf_as_i = ni / ne
    tglf_taus_i = ti / te

    tglf_as_e = 1.0
    tglf_taus_e = 1.0

    tglf_q_prime(ncon) = 0.0
    tglf_p_prime(ncon) = 0.0

    !Define Radial dimension
    call neasyf_dim(ncid, "rho_psi", values=rhopsi, dimid=rhopsi_dimid, &
                    description="Radial coordinate")
    call neasyf_dim(ncid, "npts", dim_size=npts, dimid=npts_dimid, &
                    description="Poloidal points")
    call neasyf_dim(ncid, "R", values=R, dimid=rr_dimid, &
                    units="m", description="Major radius")

    if (debug) print *, 'Assigned variable names'

    dimidsrho = [rhopsi_dimid]
    dimidsR = [rr_dimid]
    dimids2d = [rhopsi_dimid, npts_dimid]

    if (debug) print *, 'Assigned dimensions'
    !Assign 0d Variables:
    call neasyf_write(ncid, "bcentr", bcentr, units="T")
    call neasyf_write(ncid, "amin", amin, units="m")

    if (debug) print *, 'Assigned 0D variables'

    !Assign 1d Variables:
    call neasyf_write(ncid, "Psi", psiv, dim_ids=dimidsrho, units="Wb/m2")
    call neasyf_write(ncid, "Te", te, dim_ids=dimidsrho, units="eV")
    call neasyf_write(ncid, "Ti", ti, dim_ids=dimidsrho, units="eV")
    call neasyf_write(ncid, "Ne", ne, dim_ids=dimidsrho, units="m^-3")
    call neasyf_write(ncid, "Ni", ni, dim_ids=dimidsrho, units="m^-3")
    if (imp .eq. 1) then
      call neasyf_write(ncid, "N_He", nHe, dim_ids=dimidsrho, units="m^-3")
      call neasyf_write(ncid, "T_He", tHe, dim_ids=dimidsrho, units="keV")
    end if

    call neasyf_write(ncid, "LTe", Lte, dim_ids=dimidsrho)
    call neasyf_write(ncid, "LTi", Lti, dim_ids=dimidsrho)
    call neasyf_write(ncid, "Lne", Lne, dim_ids=dimidsrho)
    call neasyf_write(ncid, "Lni", Lni, dim_ids=dimidsrho)
    call neasyf_write(ncid, "shift", dshift, dim_ids=dimidsrho)
    call neasyf_write(ncid, "eps", eps, dim_ids=dimidsrho)
    call neasyf_write(ncid, "pk", pk, dim_ids=dimidsrho)
    call neasyf_write(ncid, "shat", shat, dim_ids=dimidsrho)
    call neasyf_write(ncid, "gs2_beta_prime", gs2_beta_prime, dim_ids=dimidsrho, units="")
    call neasyf_write(ncid, "zeff", zeff, dim_ids=dimidsrho)
    call neasyf_write(ncid, "vnui", vnui, dim_ids=dimidsrho)
    call neasyf_write(ncid, "vnue", vnue, dim_ids=dimidsrho)
    call neasyf_write(ncid, "Sfac", sfac, dim_ids=dimidsrho)
    call neasyf_write(ncid, "Jtot.B", jdotb, dim_ids=dimidsrho, units='A/m2')
    call neasyf_write(ncid, "Jbs.B", jbsdotb, dim_ids=dimidsrho, units='A/m2')
    call neasyf_write(ncid, "Jnb.B", jnbdotb, dim_ids=dimidsrho, units='A/m2')
    call neasyf_write(ncid, "J_tor", j_tor, dim_ids=dimidsrho, units='A/m2')
    call neasyf_write(ncid, "Jbs_tor", jbs_tor, dim_ids=dimidsrho, units='A/m2')
    call neasyf_write(ncid, "Jnb_tor", jnb_tor, dim_ids=dimidsrho, units='A/m2')
    call neasyf_write(ncid, "Jex_tor", jex_tor, dim_ids=dimidsrho, units='A/m2')
    call neasyf_write(ncid, "Fpsi", fsi, dim_ids=dimidsrho)
    call neasyf_write(ncid, "FFP", ffp, dim_ids=dimidsrho)
    call neasyf_write(ncid, "Pp", pp, dim_ids=dimidsrho)
    call neasyf_write(ncid, "dPsidrho", dPsidrhos, dim_ids=dimidsrho, units='Wb/m3')

    call neasyf_write(ncid, "JTOT", gradj(:, nsym) / 1000, dim_ids=dimidsR, units="kA/m2")
    call neasyf_write(ncid, "JEXT", exph(:, nsym) / 1000, dim_ids=dimidsR, units="kA/m2")
    call neasyf_write(ncid, "JEXT2", exph2(:, nsym) / 1000, dim_ids=dimidsR, units="kA/m2")
    call neasyf_write(ncid, "JBS", bsph(:, nsym) / 1000, dim_ids=dimidsR, units="kA/m2")
    call neasyf_write(ncid, "JDIA", diph(:, nsym) / 1000, dim_ids=dimidsR, units="kA/m2")
    call neasyf_write(ncid, "JPS", psph(:, nsym) / 1000, dim_ids=dimidsR, units="kA/m2")
    call neasyf_write(ncid, "JNBI", nbph(:, nsym) / 1000, dim_ids=dimidsR, units="kA/m2")

    !TGLF 1d variables assigment
    call neasyf_write(ncid, "TGLF_ZS_i", tglf_zs_i)
    call neasyf_write(ncid, "TGLF_MASS_i", tglf_mass_i)
    call neasyf_write(ncid, "TGLF_ZS_e", tglf_zs_e)
    call neasyf_write(ncid, "TGLF_MASS_e", tglf_mass_e)

    call neasyf_write(ncid, "TGLF_B_UNIT", tglf_b_unit, dim_ids=dimidsrho)

    call neasyf_write(ncid, "TGLF_RLNS_i", tglf_rlns_i, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_RLTS_i", tglf_rlts_i, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_TAUS_i", tglf_taus_i, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_AS_i", tglf_as_i, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_VPAR_i", tglf_vpar_i, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_VPAR_SHEAR_i", tglf_vpar_shear_i, dim_ids=dimidsrho)

    call neasyf_write(ncid, "TGLF_RLNS_e", tglf_rlns_e, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_RLTS_e", tglf_rlts_e, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_TAUS_e", tglf_taus_e, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_AS_e", tglf_as_e, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_VPAR_e", tglf_vpar_e, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_VPAR_SHEAR_e", tglf_vpar_shear_e, dim_ids=dimidsrho)

    call neasyf_write(ncid, "TGLF_RMIN", tglf_rmin, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_RMAJ", tglf_rmaj, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_Q", tglf_q, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_Q_PRIME", tglf_q_prime, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_P_PRIME", tglf_p_prime, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_DRMAJDX", tglf_drmajdx, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_KAPPA", tglf_kappa, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_S_KAPPA", tglf_s_kappa, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_DELTA", tglf_delta, dim_ids=dimidsrho)
    call neasyf_write(ncid, "TGLF_S_DELTA", tglf_s_delta, dim_ids=dimidsrho)

    if (debug) print *, 'Defined TGLF variables'

    call neasyf_write(ncid, "rpts", rpts, dim_ids=dimids2d, units='m')
    call neasyf_write(ncid, "zpts", zpts, dim_ids=dimids2d, units='m')
    call neasyf_write(ncid, "bppts", bppts, dim_ids=dimids2d, units='T')

    call neasyf_close(ncid)
    if (debug) print *, "*** SUCCESS writing NETCDF file "

  end subroutine write_netcdf

end module netcdf_interface
