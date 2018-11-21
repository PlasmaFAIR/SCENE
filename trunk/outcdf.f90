subroutine write_netcdf()
  ! File ot write SCENE output to NETCDF format
  ! Initially focused for writing gs2 output

  use param
  use netcdf
  implicit none

  logical :: debug
  character(len=12)  :: cdf_file
  character(len=4) :: file_suffix
  integer :: ncid
  integer, parameter :: ndims=1
  integer :: dimids(ndims)

  integer :: con, i


  !Variable ID for different inputs
  integer :: rhopsi_varid, rhopsi_dimid

  integer :: psi_varid
  integer :: te_varid, ti_varid, ne_varid, ni_varid
  integer :: Lte_varid, Lti_varid, Lne_varid, Lni_varid
  integer :: dshift_varid, shat_varid, pk_varid, eps_varid
  integer :: beta_varid, zeff_varid, vnui_varid, vnue_varid
  
  !Units
  character (len=*), parameter :: psi_unit="Wb/m2"
  character (len=*), parameter :: rhopsi_unit="", units="units"
  character (len=*), parameter :: te_unit="keV", ti_unit='keV', ne_unit='m^-3', ni_unit='m^-3'
  character (len=*), parameter :: Lte_unit="", Lti_unit='', Lne_unit='', Lni_unit=''
  character (len=*), parameter :: dshift_unit='', shat_unit='', pk_unit='', eps_unit=''
  character (len=*), parameter :: beta_unit='', zeff_unit='',vnui_unit='', vnue_unit=''

  
  !Output data
  double precision :: psi, tempe, tempi, dense, densi, press, shift
  double precision, dimension(ncon) :: rhopsi, te, ti, ne, ni
  double precision, dimension(ncon) :: Lte, Lti, Lne, Lni
  double precision, dimension(ncon) :: dshift, shat, pk, eps
  double precision, dimension(ncon) :: gs2_beta_prime, zeff, vnui, vnue


  double precision :: ne19, zni19, coolog, bcentr,colli,colle

  debug = .True.

  file_suffix = '.cdf'
  cdf_file = runname(1:lrunname)//file_suffix
  if (debug) print*, 'Saving to ',cdf_file, runname

  
  ! Create (Overwrite) netcdf file
  call check(nf90_create(cdf_file, NF90_CLOBBER,ncid))
  write(nw,*) 'Create NetCDF Output'

  bcentr = mu0*rodi/(2.*pi*rcen)
  
  !Generate input data
  do con=1,ncon
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
     Lne(con) = -dense(psi,0)*umax/ne(con)
     Lni(con) = -densi(psi,1,0)*umax/ni(con)

     !Other GS2 inputs
     dshift(con) =shift(con,1)*umax/amin
     eps(con) = epsv(con)
     pk(con) = 2*amin/(rcen*sfac(con))
     shat(con) = qp(con)*psi/sfac(con)
     gs2_beta_prime(con) = press(psi,1)*8.0d-7*pi*umax/bcentr**2

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
     coolog=log(sqrt(ni(con)*1.0d-6)/ti(con))
     coolog=24.-coolog
     
     !From Wesson 2.15 (added by bhavin 21/03/18)
     colli = 6.6d17*zmas(1)**0.5 * (ti(con)/1000)**1.5/(ni(con)*iz(1)**4*coolog)
     vnui(con) =amin/ (sqrt(2*te(con)*bk/mp) *colli)

     !vss = sqrt(2.0)*pi*ni * iz(i)**4 * eq**4 * coolog &
     !	/ ( sqrt(zmas(i)*mp) * (ti*eq/bk)**1.5  * (4*pi*eps0)**2 )

     !Collisionality for electrons
     coolog=log(sqrt(ne(con)*1.0d-6)/te(con))
     coolog=24.-coolog
     !vss = sqrt(2.0)*pi*dense(psi,0) * eq**4 * coolog &
     !	/ (sqrt(me) * (tempe(psi,0)/1000.)**1.5 * (4*pi*eps0)**2 )
     
     colle = 3.*(2.*pi)**1.5* eps0**2 * me**0.5 * (tempe(psi,0)*bk)**1.5 &
          / ( dense(psi,0) * eq**4 * coolog)
     vnue(con) = amin/ (sqrt(2*te(con)*bk/mp) *colle)

  
  end do

  
  !Define Radial dimension
  call check( nf90_def_dim(ncid, "rho_psi", ncon, rhopsi_dimid) )
  if (debug) print*, 'Assigned variable name'

  !Define Co-ordinate variable and its units
  call check( nf90_def_var(ncid, "rho_psi", NF90_REAL, rhopsi_dimid, rhopsi_varid))
  call check( nf90_put_att(ncid, rhopsi_varid, units, rhopsi_unit))
  dimids = (/rhopsi_dimid/)

  !Assign Variables: Electron Temp
  call check( nf90_def_var(ncid, "Psi", NF90_REAL, dimids, psi_varid) ) 
  call check( nf90_def_var(ncid, "Te", NF90_REAL, dimids, te_varid) )
  call check( nf90_def_var(ncid, "Ti", NF90_REAL, dimids, ti_varid) )
  call check( nf90_def_var(ncid, "Ne", NF90_REAL, dimids, ne_varid) )
  call check( nf90_def_var(ncid, "Ni", NF90_REAL, dimids, ni_varid) )
  call check( nf90_def_var(ncid, "LTe", NF90_REAL, dimids, Lte_varid) )
  call check( nf90_def_var(ncid, "LTi", NF90_REAL, dimids, Lti_varid) )
  call check( nf90_def_var(ncid, "Lne", NF90_REAL, dimids, Lne_varid) )
  call check( nf90_def_var(ncid, "Lni", NF90_REAL, dimids, Lni_varid) )
  call check( nf90_def_var(ncid, "shift", NF90_REAL, dimids, dshift_varid) )
  call check( nf90_def_var(ncid, "eps", NF90_REAL, dimids, eps_varid) )
  call check( nf90_def_var(ncid, "pk", NF90_REAL, dimids, pk_varid) )
  call check( nf90_def_var(ncid, "shat", NF90_REAL, dimids, shat_varid) )
  call check( nf90_def_var(ncid, "gs2_beta_prime", NF90_REAL, dimids, beta_varid) )
  call check( nf90_def_var(ncid, "zeff", NF90_REAL, dimids, zeff_varid) )
  call check( nf90_def_var(ncid, "vnui", NF90_REAL, dimids, vnui_varid) )
  call check( nf90_def_var(ncid, "vnue", NF90_REAL, dimids, vnue_varid) )
  
  if (debug) print*, 'Defined variables'

  !Assign Units
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
  call check( nf90_put_att(ncid,beta_varid, units, beta_unit))
  call check( nf90_put_att(ncid,zeff_varid, units, zeff_unit))
  call check( nf90_put_att(ncid,vnui_varid, units, vnui_unit))
  call check( nf90_put_att(ncid,vnue_varid, units, vnue_unit))  
  

  if (debug) print*, 'Assigned variable attributes'
  
  call check(nf90_enddef(ncid))

  !Put Radial coordinates in
  call check(nf90_put_var(ncid, rhopsi_dimid, rhopsi))

  !Put in variable data
  call check(nf90_put_var(ncid, psi_varid, psi))
  call check(nf90_put_var(ncid, te_varid, te))
  call check(nf90_put_var(ncid, ti_varid, ti)) 
  call check(nf90_put_var(ncid, ne_varid, ne))
  call check(nf90_put_var(ncid, ni_varid, ni))
  call check(nf90_put_var(ncid, ne_varid, Lte))
  call check(nf90_put_var(ncid, Lte_varid, Lti))
  call check(nf90_put_var(ncid, Lti_varid, Lne))
  call check(nf90_put_var(ncid, Lni_varid, Lni))
  call check(nf90_put_var(ncid, dshift_varid, dshift))
  call check(nf90_put_var(ncid, eps_varid, eps))
  call check(nf90_put_var(ncid, pk_varid, pk))
  call check(nf90_put_var(ncid, shat_varid, shat))
  call check(nf90_put_var(ncid, beta_varid, beta))
  call check(nf90_put_var(ncid, zeff_varid, zeff))
  call check(nf90_put_var(ncid, vnui_varid, vnui)) 
  call check(nf90_put_var(ncid, vnue_varid, vnue)) 

  if (debug) print*, 'Put variables values in'
  

  
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
       
