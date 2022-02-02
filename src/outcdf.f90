module netcdf_interface
  implicit none
contains

  !> Write SCENE output to netCDF file
  subroutine write_netcdf(header)
    use, intrinsic :: iso_fortran_env, only: real64
    use param
    use neasyf, only: neasyf_open, neasyf_close, neasyf_dim, neasyf_write, neasyf_metadata
    use profiles_mod, only: tempe, tempi, dense, densi, press, shift, fprof, dpsidrho, dVdrho, rhotor, elong, triang
    use git_version, only: get_git_version
    use header_mod, only: header_type
    implicit none

    !> Header with run UUID and creation timestamp
    type(header_type), intent(in) :: header

    logical :: debug
    character(len=:), allocatable :: cdf_file
    character(len=*), parameter :: file_suffix = ".cdf"
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

    allocate(character(len=len_trim(runname) + len_trim(file_suffix))::cdf_file)
    cdf_file = trim(runname) // file_suffix
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

    call neasyf_metadata(ncid, &
                         software_name="SCENE", &
                         software_version=get_git_version(), &
                         title=title, &
                         file_id=header%run_uuid, &
                         date_created=header%date_time)

    !Define Radial dimension
    call neasyf_dim(ncid, "rho_psi", values=rhopsi, dimid=rhopsi_dimid, &
                    long_name="Radial coordinate")
    call neasyf_dim(ncid, "npts", dim_size=npts, dimid=npts_dimid, &
                    long_name="Poloidal points")
    call neasyf_dim(ncid, "R", values=R, dimid=rr_dimid, &
                    units="m", long_name="Major radius")

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

    call neasyf_write(ncid, "rpts", rpts, dim_ids=dimids2d, units='m', &
         long_name="Flux surface R points")
    call neasyf_write(ncid, "zpts", zpts, dim_ids=dimids2d, units='m', &
         long_name="Flux surface Z points")
    call neasyf_write(ncid, "bppts", bppts, dim_ids=dimids2d, units='T')

    call write_input_parameters(ncid)
    call write_0D_quantities(ncid)

    call neasyf_close(ncid)
    if (debug) print *, "*** SUCCESS writing NETCDF file "

  end subroutine write_netcdf

  !> Write the input parameters to a new group in given file
  subroutine write_input_parameters(file_id)
    use param, only: rcen, tokeps, ibdry, nfm, elon, tri, quad, kval, step, &
                     ipr, igr, fpow, powj, fpow3, fpow4, af3, af4, fpow1, &
                     fpow2, af1, af0, af2, psic, ipswtch, ppow, p0, pa, &
                     pped, ppo1, ape, pedg, nipow, ni0, nia, naa, nbb, niped, &
                     niedg, ane, npo1, ate, tpo1, tpoe, te0, tea, teped, teedg, &
                     tpoi, ti0, tia, tiped, tiedg, xitb, cur, paux, bpol, pfac, &
                     betan, dil, rodi, imat, ifast, scl, omega, frac, zm, zmai, &
                     icontour, imp, itot, fast, neo, nco, ncon, npts, nouter, &
                     ninner, npass, errit, errffd, icont, nbi, nimp
    use balpar, only: nbal, ibal, nturns, nchi0, nchi, chi0val, lamges
    use neasyf, only: neasyf_write, neasyf_error
    use netcdf, only: nf90_def_grp

    !> NetCDF ID of parent file
    integer, intent(in) :: file_id

    ! NetCDF ID of the parameter group
    integer :: group_id

    call neasyf_error(nf90_def_grp(file_id, "inputs", group_id), &
                      ncid=file_id, message="creating inputs group")

    call neasyf_write(group_id, "af0", af0, &
         long_name="Coefficients used in ffp and J profile")
    call neasyf_write(group_id, "af1", af1, &
         long_name="Coefficients used in ffp and J profile")
    call neasyf_write(group_id, "af2", af2, &
         long_name="Coefficients used in ffp and J profile")
    call neasyf_write(group_id, "af3", af3, &
         long_name="Coefficients used in ffp and J profile")
    call neasyf_write(group_id, "af4", af4, &
         long_name="Coefficients used in ffp and J profile")
    call neasyf_write(group_id, "ane", ane, &
         long_name="Flatten core density for ips==6 or 21")
    call neasyf_write(group_id, "ape", ape, &
         long_name="Flatten core pressure if ips==21")
    call neasyf_write(group_id, "ate", ate, &
         long_name="Flatten core temperature for ips==6")
    call neasyf_write(group_id, "betn", betan, &
         long_name="Target value of betan (bpol re-scaled if betan>0)")
    call neasyf_write(group_id, "bpol", bpol, &
         long_name="Approximate starting beta-poloidal for p0==1")
    call neasyf_write(group_id, "chi0", chi0val, &
         long_name="")
    call neasyf_write(group_id, "cur", cur, &
         long_name="Total toroidal plasma current", units="MA")
    call neasyf_write(group_id, "dil", dil, &
         long_name="Ion dilution factor for imp==0")
    call neasyf_write(group_id, "elon", elon, &
         long_name="Elongation of plasma boundary")
    call neasyf_write(group_id, "eps", tokeps, &
         long_name="Inverse aspect ratio of boundary")
    call neasyf_write(group_id, "erfd", errffd, &
         long_name="Target error between input and output currents")
    call neasyf_write(group_id, "erit", errit, &
         long_name="Target error for Grad-Shafranov solver")
    call neasyf_write(group_id, "fast", fast, &
         long_name="Calculate beta due to fast alpha particles")
    call neasyf_write(group_id, "fpo1", fpow1, &
         long_name="Initial ff' profile (increase to peak near axis)")
    call neasyf_write(group_id, "fpo2", fpow2, &
         long_name="Initial ff' profile (increase to peak near edge)")
    call neasyf_write(group_id, "fpo3", fpow3, &
         long_name="Initial ff' profile (increase to peak near axis)")
    call neasyf_write(group_id, "fpo4", fpow4, &
         long_name="Initial ff' profile (increase to peak near edge)")
    call neasyf_write(group_id, "fpow", fpow, &
         long_name="Initial ff' profile (on-axis peaking)")
    call neasyf_write(group_id, "frac", frac, &
         long_name="Control rate of convergence of Grad-Shafranov solver")
    call neasyf_write(group_id, "ibal", ibal, &
         long_name="The flux surface analysed for ballooning stability for nbal==1")
    call neasyf_write(group_id, "ibdr", ibdry, &
         long_name="Boundary shape: 0 is default (Sykes) shape; 1 is more general shape")
    call neasyf_write(group_id, "icon", icont, &
         long_name="If 1, use a 'warm' start for the equilibrium solver")
    call neasyf_write(group_id, "ictr", icontour, &
         long_name="Controls contouring algorithm; 0: linear interpolation, &
         &1: linear interpolation + smoothing, 2: bi-cubic spline interpolation")
    call neasyf_write(group_id, "ifas", ifast, &
         long_name="If 1, bootstrap current includes alpha-particle contribution")
    call neasyf_write(group_id, "igr", igr, &
         long_name="Controls graphics output (requires GHOST)")
    call neasyf_write(group_id, "imat", imat, &
         long_name="Bootstrap model: 0: Hirshman, 1: Hirshman-Sigmar")
    call neasyf_write(group_id, "imp", imp, &
         long_name="Enable impurities")
    call neasyf_write(group_id, "ipr", ipr, &
         long_name="Controls printout to screen")
    call neasyf_write(group_id, "ipsw", ipswtch, &
         long_name="Profile form")
    call neasyf_write(group_id, "itot", itot, &
         long_name="If 1, use ff-prime and p-prime to derive equilibrium, else use specified driven current profile")
    call neasyf_write(group_id, "kval", kval, &
         long_name="If ibdr==1, kval->1 provides DND; kval=0 limiter")
    call neasyf_write(group_id, "lam", lamges, &
         long_name="")
    call neasyf_write(group_id, "naa", naa, &
         long_name="")
    call neasyf_write(group_id, "nbal", nbal, &
         long_name="")
    call neasyf_write(group_id, "nbb", nbb, &
         long_name="")
    call neasyf_write(group_id, "nbi", nbi, &
         long_name="")
    call neasyf_write(group_id, "nchi", nchi, &
         long_name="")
    call neasyf_write(group_id, "nci0", nchi0, &
         long_name="")
    call neasyf_write(group_id, "nco", nco, &
         long_name="Collisionality model")
    call neasyf_write(group_id, "ncon", ncon, &
         long_name="Number of flux surfaces evaluated")
    call neasyf_write(group_id, "nedg", niedg, &
         long_name="Pedestal density gradient for imp==1")
    call neasyf_write(group_id, "neo", neo, &
         long_name="Enable neoclassical current")
    call neasyf_write(group_id, "nfm", nfm, &
         long_name="Number of Fourier modes for ibdr==2")
    call neasyf_write(group_id, "ni0", ni0, &
         long_name="Central main ion density for imp==1")
    call neasyf_write(group_id, "nia", nia, &
         long_name="Edge main ion density fraction for imp==1")
    call neasyf_write(group_id, "ninn", ninner, &
         long_name="Number of inner iterations in Grad-Shafranov solver")
    call neasyf_write(group_id, "nimp", nimp, &
         long_name="Number of impurity species")
    call neasyf_write(group_id, "nout", nouter, &
         long_name="Number of iterations in Grad-Shafranov solver")
    call neasyf_write(group_id, "npas", npass, &
         long_name="Maximum number of equilibrium iterations")
    call neasyf_write(group_id, "nped", niped, &
         long_name="Pedestal main ion density fraction for imp==1")
    call neasyf_write(group_id, "npo1", npo1, &
         long_name="")
    call neasyf_write(group_id, "npow", nipow, &
         long_name="Main ion density profile for imp==1")
    call neasyf_write(group_id, "npts", npts, &
         long_name="Number of poloidal points on flux surface")
    call neasyf_write(group_id, "ntrn", nturns, &
         long_name="Number of periods in 2*pi along field line for ballooning calculation &
         &(nturns in positive and negative directions)")
    call neasyf_write(group_id, "omeg", omega, &
         long_name="Accelerate equilibrium convergence (proportional to psi)")
    call neasyf_write(group_id, "p0", p0, &
         long_name="Extra pressure to help get beta up")
    call neasyf_write(group_id, "pa", pa, &
         long_name="mu0* (SI units) edge pressure for imp==0")
    call neasyf_write(group_id, "paux", paux, &
         long_name="Auxiliary heating power for taue", units="MW")
    call neasyf_write(group_id, "pedg", pedg, &
         long_name="Pedestal pressure gradient for imp==0")
    call neasyf_write(group_id, "pfac", pfac, &
         long_name="")
    call neasyf_write(group_id, "powj", powj, &
         long_name="Exponent of driven current profile")
    call neasyf_write(group_id, "pped", pped, &
         long_name="mu0* pedestal pressure (approx) (SI units) for imp==0")
    call neasyf_write(group_id, "ppo1", ppo1, &
         long_name="")
    call neasyf_write(group_id, "ppow", ppow, &
         long_name="Electron pressure profile for imp==0")
    call neasyf_write(group_id, "psic", psic, &
         long_name="Extent in psi where bootstrap current is filled in")
    call neasyf_write(group_id, "quad", quad, &
         long_name="Quadracity (squareness) for ibdr==1")
    call neasyf_write(group_id, "rcen", rcen, &
         long_name="Geometric axis of plasma boundary", units="m")
    call neasyf_write(group_id, "rodi", rodi, &
         long_name="Rod current (gives vacuum B-field", units="MA")
    call neasyf_write(group_id, "scl", scl, &
         long_name="Guess at current scaling")
    call neasyf_write(group_id, "step", step, &
         long_name="Mesh size", units="m")
    call neasyf_write(group_id, "te0", te0, &
         long_name="Central electron temperature", units="eV")
    call neasyf_write(group_id, "tea", tea, &
         long_name="Edge electron temperature", units="eV")
    call neasyf_write(group_id, "tedg", teedg, &
         long_name="Pedestal electron temperature gradient")
    call neasyf_write(group_id, "tepd", teped, &
         long_name="Pedestal electron temperature", units="eV")
    call neasyf_write(group_id, "ti0", ti0, &
         long_name="Central ion temperature", units="eV")
    call neasyf_write(group_id, "tia", tia, &
         long_name="Edge ion temperature", units="eV")
    call neasyf_write(group_id, "tidg", tiedg, &
         long_name="Pedestal ion temperature gradient")
    call neasyf_write(group_id, "tipd", tiped, &
         long_name="Pedestal ion temperature", units="eV")
    call neasyf_write(group_id, "tpo1", tpo1, &
         long_name="Used to flatten core temperature")
    call neasyf_write(group_id, "tpoe", tpoe, &
         long_name="Electron temperature profile")
    call neasyf_write(group_id, "tpoi", tpoi, &
         long_name="Ion temperature profile")
    call neasyf_write(group_id, "tri", tri, &
         long_name="Triangularity of plasma boundary")
    call neasyf_write(group_id, "xitb", xitb, &
         long_name="Position of transport barrier in normalised flux, from the edge")
    call neasyf_write(group_id, "zm", zm, &
         long_name="Charge of main species ion (<ZEFF)")
    call neasyf_write(group_id, "zmai", zmai, &
         long_name="Main ion species mass (c.f. H)")

    call write_impurity_input_parameters(group_id)

  end subroutine write_input_parameters

  subroutine write_impurity_input_parameters(input_group_id)
    use param, only: iz, zmas, ztpow, zt0, zta, ztped, ztedg, znpow, zn0, zna, znped, znedg, nimp
    use neasyf, only: neasyf_write, neasyf_error
    use netcdf, only: nf90_def_grp
    !> NetCDF ID of the input parameter group
    integer, intent(in) :: input_group_id

    ! NetCDF ID of the parameter group
    integer :: impurity_group_id
    character(len=11) :: impurity_group_name
    integer :: impurity

    do impurity = 2, nimp + 1
      write(impurity_group_name, '("impurity_", i2.2)') impurity - 1
      call neasyf_error(nf90_def_grp(input_group_id, impurity_group_name, impurity_group_id), &
                        ncid=input_group_id, message="creating inputs group")

      call neasyf_write(impurity_group_id, "Z", iz(impurity), &
                        long_name="Impurity charge")
      call neasyf_write(impurity_group_id, "M", zmas(impurity), &
                        long_name="Impurity mass", units="Mp")
      call neasyf_write(impurity_group_id, "at", ztpow(impurity), &
                        long_name="Impurity temperature profile")
      call neasyf_write(impurity_group_id, "T0", zt0(impurity), &
                        long_name="Central impurity temperature")
      call neasyf_write(impurity_group_id, "Ta", zta(impurity), &
                        long_name="Edge impurity temperature")
      call neasyf_write(impurity_group_id, "Tped", ztped(impurity), &
                        long_name="Pedestal impurity temperature")
      call neasyf_write(impurity_group_id, "Tedg", ztedg(impurity), &
                        long_name="Pedestal impurity temperature gradient")
      call neasyf_write(impurity_group_id, "an", znpow(impurity), &
                        long_name="Impurity density profile")
      call neasyf_write(impurity_group_id, "n0", zn0(impurity), &
                        long_name="Central impurity density")
      call neasyf_write(impurity_group_id, "na", zna(impurity), &
                        long_name="Edge impurity density")
      call neasyf_write(impurity_group_id, "nped", znped(impurity), &
                        long_name="Pedestal impurity density")
      call neasyf_write(impurity_group_id, "nedg", znedg(impurity), &
                        long_name="Pedestal impurity density gradient")
    end do

  end subroutine write_impurity_input_parameters

  !> Write some 0D quantities
  !>
  !> These are the equivalent to the variables written in [[getdata_mod:popcon]]
  subroutine write_0D_quantities(file_id)
    use param, only : avt, avti, avel, nebar, avio, zm, zeffav, area, beta, betap, betexp, &
         bmfus, bmshine, conft, cur, betlim, hipb98y1, hipb98y2, ncon, mu0, negw, paux, &
         pfus, petty, pi, r0, rcen, rli2, rli3, rodi, taue, tauh, totbs, totdi, totex, totex2, &
         vol, sfac, totnb, totps
    use profiles_mod, only : dense, densi, fprof
    use neasyf, only: neasyf_write

    !> NetCDF ID of parent file
    integer, intent(in) :: file_id

    call neasyf_write(file_id, 'volume_average_Te', avt, units='eV', &
         long_name="Volume average electron temperature")
    call neasyf_write(file_id, 'volume_average_Ti', avti, units='eV', &
         long_name="Volume average electron temperature")

    call neasyf_write(file_id, 'volume_average_Ne', avel, units='m^-3', &
         long_name="Volume average electron density")
    call neasyf_write(file_id, 'central_Ne', dense(0.0d0, 0), units='m^-3', &
         long_name="Central electron density")
    call neasyf_write(file_id, 'line_average_Ne', nebar*1.0d19, units='m^-3', &
         long_name="Line average electron density")
    call neasyf_write(file_id, 'N_gw', negw, &
         long_name="Line average electron density normalised to Greenwald limit")

    call neasyf_write(file_id, 'volume_average_Ni', avio*1.0d19, units='m^-3', &
         long_name="Volume average ion density")
    call neasyf_write(file_id, 'central_Ni', densi(0.0d0, 1, 0), units='m^-3', &
         long_name="Central ion density")

    call neasyf_write(file_id, 'Z_i', zm, long_name="Charge on main ions")
    call neasyf_write(file_id, 'Zeff', zeffav, &
         long_name="Volume average effective charge")

    call neasyf_write(file_id, 'B tor (mag)', fprof(0.0d0, 2)/r0, units='T')
    call neasyf_write(file_id, 'B tor (geo)', fprof(0.0d0, 2)/rcen, units='T')
    call neasyf_write(file_id, 'Vac B tor (geo)', mu0*rodi/(2.*pi*rcen), units='T')

    call neasyf_write(file_id, 'Tor. tot cur', cur/1000000., units='MA', &
         long_name="Total toroidal current")
    call neasyf_write(file_id, 'Tor. bs cur', totbs/1000000., units='MA', &
         long_name="Bootstrap toroidal current")
    call neasyf_write(file_id, 'Tor. nb cur', totnb/1000000., units='MA', &
         long_name="Neutral beam toroidal current")
    call neasyf_write(file_id, 'Tor. ps cur', totps/1000000., units='MA', &
         long_name="Pfirsch-Schluter toroidal current")
    call neasyf_write(file_id, 'Tor. dia cur', totdi/1000000., units='MA', &
         long_name="Diamagnetic toroidal current")
    call neasyf_write(file_id, 'Tor. ext cur', (totex+totex2)/1000000., units='MA', &
         long_name="Externally driven toroidal current")

    call neasyf_write(file_id, 'Rod Current', rodi/1000000., units='MA')

    call neasyf_write(file_id, 'beta', beta, units='%', &
         long_name="Plasma beta")
    call neasyf_write(file_id, 'beta_norm', 3.5*betexp/betlim, &
         long_name="Normalised beta")
    call neasyf_write(file_id, 'beta_poloidal', betap, &
         long_name="Poloidal beta")

    call neasyf_write(file_id, 'q0', sfac(ncon), &
         long_name="Safety factor at magnetic axis")
    call neasyf_write(file_id, 'qa', sfac(1), &
         long_name="Safety factor at edge")
    call neasyf_write(file_id, 'qmin', minval(sfac,1), &
         long_name="Minimum value of safety factor")

    call neasyf_write(file_id, 'Pfus',pfus*1.0d-6*5, units='MW', &
         long_name="Alpha heating power")
    call neasyf_write(file_id, 'Beam fus',bmfus, units='MW', &
         long_name="Beam fusion power")
    call neasyf_write(file_id, 'Aux Pow',paux, units='MW', &
         long_name="Auxiallary heating power")
    call neasyf_write(file_id, 'NBP eff.', bmfus/(paux))
    call neasyf_write(file_id, 'NBCD eff.', totnb/(1.0d6*paux))
    call neasyf_write(file_id, 'Shine through', bmshine, units='MW')
    call neasyf_write(file_id, 'Q', ((pfus*1.0d-6*5)+bmfus)/paux)

    call neasyf_write(file_id, 'H_IPB98(y1)', hipb98y1, &
         long_name="Confinement time relative to scaling laws")
    call neasyf_write(file_id, 'H_IPB98(y2)', hipb98y2, &
         long_name="Confinement time relative to scaling laws")
    call neasyf_write(file_id, 'H_Petty08',   petty)
    call neasyf_write(file_id, 'Tau_e', taue*1.0e3, units='ms', &
         long_name="Energy confinement time")
    call neasyf_write(file_id, 'Tau_h', tauh, &
         long_name="Alpha confinement time")

    call neasyf_write(file_id, 'Area', area, units='m^2', &
         long_name="Plasma cross-sectional area")
    call neasyf_write(file_id, 'Volume', vol, units='m^3', &
         long_name="Plasma volume")
    call neasyf_write(file_id, 'li(3)', rli3, &
         long_name="Internal inductance")
    call neasyf_write(file_id, 'li(2)', rli2, &
         long_name="Internal inductance")
    call neasyf_write(file_id, 'Energy', conft/1.0e6, units='MJ')
  end subroutine write_0D_quantities

end module netcdf_interface
