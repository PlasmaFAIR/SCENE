! SCENE: Tokamak magnetic equilibrium solver
! Copyright (C) 2021 Howard Wilson, Bhavin Patel, David Dickinson, Peter Hill
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!> SPDX-License-Identifier: GPL-3.0-only
program scene
  use, intrinsic :: iso_fortran_env, only : output_unit
  use bootstap, only : fastbs
  use equilibrium, only : equil, xarea
  use ff_mod, only : ffdgen
  use flux_average, only: flxav
  use flux_plot, only : initialise_graphs, make_all_graphs
  use flux_surface, only : flxsur
  use getdata_mod, only : getdata
  use hirsig_mod, only : erfun
  use netcdf_interface, only : write_netcdf
  use nbi_mod, only : nbicur
  use output_mod, only : output
  use param
  use usrcal_mod, only : usrcal
  use setnt_mod, only : setnt
  use scene_errors, only : error_msg
  use scene_init, only : init
  use toroidal_current_mod, only : torcur
  use header_mod, only : header_type
  implicit none
  integer i,j,icur,niter
  double precision errcur,bpold,bnold,bpnew
  logical :: debug
  type(header_type) :: header

  header = header_type()

  debug = .false.
  nbet=0
  ! Perform initialisations
  write(output_unit, *) 'Starting SCENE'
  write(output_unit, '(a)') header%to_string(comment_character=" ")
  call error_msg('SCENE has begun', -1)

  call init
  ! Initialise error function array for use in bootstrap routine, hirsig
  ! In call bootstrap routine hirsig, don;'t forget you need to load up
  !  error function first...and don't call it errf!
  !  Junk that needs to be taken account of before bootstrap, etc
  !  Calculate error function array
  call erfun

!!$C  Initialise ICUR parameter...if ICUR=-1 trapped particle
!!$C  effects are not included in the plasma conductivity.
!!$      ICUR=0
  if (igr.gt.0) then
    call initialise_graphs()
  end if

  beta_convergence: do
    ! Call SCENE equilibrium routines
    equilibrium_convergence: do niter = 1, npass

      if (ipr.eq.0) write(nw,*)' on the ',niter,' iteration'
      if (debug) write(6,*) 'Printing poloidal/toroidal graphs'
      call equil(niter)


      if (debug) write(6,*)' done equil'
      ! Convert equilibrium to flux surface coordinates
      call flxsur
      if (debug) write(6,*)' done flxsur'
      ! Perform a number of useful averages over flux surfaces
      call flxav
      ! Calculate alpha-particle bootstrap current
      call fastbs
      if (debug) write(6,*)' done flxav'

      j_nb = 0.

      !Include NBI in current calc
      if (nbi.ge.1) call nbicur()

      if (debug) print*, 'done nbicur'
      !  Calculate currents (and read in externally applied current
      !  profile if itot=0)
      call torcur(icur)

      !do i=1,nr
      !   print*, i, nbph(i,nsym)
      !end do

      !  Calculate total current contributions....
      call xarea(bsph,totbs)
      call xarea(absph,alfbs)
      call xarea(psph,totps)
      call xarea(diph,totdi)
      call xarea(exph,totex)
      call xarea(nbph,totnb)
      call xarea(exph2,totex2)
      call xarea(gradj,totgs)

      ! Deriving equilibrium from ff-prime and p-prime profiles
      if (itot /= 0) exit equilibrium_convergence

      !  First calculate scaling factor on external current to give
      !  required total current (this will be the loop voltage in the
      ! case of a neoclassical ohmic profile)

      ! Flag on whether to include NBI in FFDGEN calc
      if (nbi .eq. 2) then
        vloop=(totgs-totbs-totps-totdi-totex2-totnb)/(totex)
      else
        vloop=(totgs-totbs-totps-totdi-totex2)/(totex)
      end if

      if (debug) print*, ' '
      if (debug) print*, 'Vloop = ', vloop
      if (debug) print*, ' '
      if (abs(neo).eq.1) then

        if (nbi .eq. 2) then
          vnobs=(totgs-totps-totdi-totex2-totnb)/(totex)
        else
          vnobs=(totgs-totps-totdi-totex2)/(totex)
        end if

        call xarea(spit,spiti)
        vspit=(totgs-totps-totdi-totex2)/spiti
      end if
      !  Label exph as total driven current...
      do i=1,nr
        do j=1,nz
          exph(i,j)=exph(i,j)*vloop+exph2(i,j)
        end do
      end do
      totex=totex*vloop+totex2

      !if (debug) print*, totgs, totbs, totps, totdi, totex2, totnb, totex

      !  calculate error in ff'
      call ffdgen(1,icur,errcur)
      if (debug) write(6,*)' done ffdgen'

      if (errcur <= errffd) then
        ! We're converged!
        exit equilibrium_convergence
      end if

      !  need to iterate further on ff'
      !  up-date ff'
      call ffdgen(0,icur,errcur)
      !  up-date n,T mesh if appropriate
      !if (debug) write(6,*)' in setnt'
      if (igr.eq.5) call setnt
      !if (debug) write(6,*)' out setnt'
      icont=-3
    end do equilibrium_convergence

    ! We reached npass iterations without convergence
    if (errcur > errffd) then
      write(nw,*)'warning***ffd has not converged to specified'
      write(nw,*)'accuracy. accuracy=',errcur,' requested=',errffd
    end if

    !  equilibrium converged to specified accuracy...
    !  generate output
    call output

    ! If we're not trying to converge in betan, we're done
    if (betan <= 0.) then
      exit beta_convergence
    end if

    nbet=nbet+1
    if (debug) write(6,*)' betan=',3.5*betexp/betlim

    if (abs(betan-3.5*betexp/betlim)/betan <= 0.01) then
      ! Converged in betan
      exit beta_convergence
    end if

    if (nbet.eq.1) then
      bpold=bpol
      bnold=betlim/(3.5*betexp)
      bpol=bpol*betan*betlim/(3.5*betexp)
      if (debug) write(6,*)' nbet=',nbet,' betan=',3.5*betexp/betlim, &
           ' target=',betan
    else
      bpnew=(betan*(bpold-bpol)+bnold*bpol-bpold*3.5*betexp/betlim)  &
           /(bnold-3.5*betexp/betlim)
      bpold=bpol
      bnold=3.5*betexp/betlim
      bpol=bpnew
      if (debug) write(6,*)' nbet=',nbet,' betan=',3.5*betexp/betlim, &
           ' target=',betan
    end if

  end do beta_convergence

  if (ipswtch.eq.3) call setnt
  if (debug) write(6,*)' done setnt'

  !Doesn't include NBI in current calc
  !if (nbi.eq.1) call nbicur

 !Write data to files
  call getdata

  !  Plot out some useful figs
  if (igr.gt.0) then
    call make_all_graphs(debug)
  end if

  !  Plot stability plots if igr set to 3
  !      if (igr.eq.3) call stab
  !      if (igr.eq.3) call epsplot
  !  Call user interface routine
  call usrcal
  if (debug) print*, 'Calling NETCDF writer'
  call write_netcdf()

  call error_msg('Clean exit', 0)

end program scene
