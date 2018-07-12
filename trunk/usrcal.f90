      subroutine usrcal
!     *****************
!
!  This is the user interface with SCENE.  Call to user subroutines can be
!  made from here.....the equilibrium is now converged to the required
!  accuracy.
!
      use param
      use balpar
      implicit none
      integer icur,i,n1,n2,j
      double precision extapp,extapp2,rat,scale,eps
      double precision fprof
      double precision psi1,psi2,eps1,eps2,temp0,temp1,temp2,tempi,psi
      double precision tempe,dense,densi
      double precision tlen(ncon),epsn(ncon),psinorm(ncon)
      double precision arr0(ncon),arr1(ncon),arr2(ncon),arr3(ncon),arr4(ncon)
      double precision pow(3,ncon),dum,voltst,pprof(3,ncon),ptprof(ncon)
      character(len=8) text
!
!  Write temperature and density data for radiation calc
      open(unit=29,file=runname(1:lrunname)//'.radat', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.radat'
         stop
      endif
      write(29,*)' Electrons: psi, Te, ne'
      do i=1,ncon
        psi=psiv(i)
        write(29,*)psi/umax,tempe(psi,0),dense(psi,0)
      end do
      do j=1,nimp+1
        write(29,*)' Ion species, charge=',iz(j)
        write(29,*)' Psi, Ti, ni'
        do i=1,ncon
          psi=psiv(i)
          write(29,*)psi/umax,tempi(psi,j,0),densi(psi,j,0)
        end do
      end do
      close(29)
!!$! Read in power densities from Martin O'Mullane's code...
!!$      open(unit=29,file='pdenw.dat', &
!!$           status='unknown',iostat=ios)
!!$      if(ios.ne.0) then
!!$         write(6,*) 'problem opening pdenw.dat'
!!$         stop
!!$      endif
!!$      read(29,*)text
!!$      write(6,*)text
!!$      read(29,*)text
!!$      write(6,*)text
!!$      read(29,*)text
!!$      write(6,*)text
!!$      do i=1,ncon
!!$        read(29,*)psi,dum,dum,pow(1,i),pow(2,i),pow(3,i)
!!$      end do
!!$      close(29)
!!$!  Integrate over volume (skip axis point)
!!$      powsum=0.
!!$      do j=1,3
!!$         voltst=0.
!!$         powtot(j)=0.
!!$         do i=1,ncon-2
!!$           voltst=voltst+0.5*(rnorm(i)+rnorm(i+1))*(psiv(i)-psiv(i+1))
!!$           powtot(j)=powtot(j)+0.5*(rnorm(i)*pow(j,i)+rnorm(i+1)*pow(j,i+1))*  &
!!$                                    (psiv(i)-psiv(i+1))
!!$         end do
!!$         powtot(j)=2.*pi*powtot(j)
!!$         powsum=powsum+powtot(j)
!!$         voltst=voltst*2.*pi
!!$         write(6,*)' Impurity=',j,'  voltst=',voltst,' power=',powtot(j)
!!$      end do
!!$      write(6,*)' Total radiated power=',powsum,' MWm**-3'
      do j=1,3
        pprof(j,ncon-1)=0.
        pprof(j,ncon)=0.
        do i=ncon-2,1,-1
           pprof(j,i)=pprof(j,i+1)+0.5*(rnorm(i)*pow(j,i)+rnorm(i+1)*pow(j,i+1))*  &
                                    (psiv(i)-psiv(i+1))*2.*pi
        end do
      end do
      do i=1,ncon
        ptprof(i)=0.
        do j=1,3
          ptprof(i)=ptprof(i)+pprof(j,i)
        end do
!        write(6,*)' i=',i,' ptprof=',ptprof(i),' pprof=',pprof(3,i)
      end do
!      call tstplt(ncon,psiv,ptprof,0.0d0,0.0d0)
      icur=1
      call extj(tokeps,umax,extapp,extapp2,icur)
      rat=extapp*vloop*sqrt(bsqav(1))/bsj(1)
      !write(6,*)' Ratio of edge Ohmic to bootstrap current=',rat
      !write(6,*)' Ohmic edge current, J||/B=',extapp*vloop,extapp2
      !write(6,*)' Bootstrapc edge current, J||/B=',bsj(1)/sqrt(bsqav(1))
      !write(6,*)' sqrt(bsqav)=',sqrt(bsqav(1))
      call gs2
!      call astradat
!!$      call bandat
!!$      write(6,*)' done bandat'
      call tokameq
!      call bandit
!      call curplt
!!$      write(6,*)' done bandit'
!      call flxorb
      call helena

      call geqdsk
!!$      call lynton
!!$      write(6,*)' done lynton'
!!$!      call fastbs
!!$      write(6,*)' done fastbs'
! Data for Preinhalter
      open(unit=73,file=runname(1:lrunname)//'.ebweq', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.ebweq'
         stop
      endif
      write(73,*)'NR, NZ'
      write(73,*)nr,nz
      write(73,*)'R(m), Z(m), psi(R,Z) Tm**2'
      do i=1,nr
        do j=1,nz
          write(73,*)r(i),z(j),u(i,j)
        end do
      end do
      write(73,*)'psi (Tm**2), f(psi)=R*Bphi (Tm); psi=0 is edge'
      do i=1,ncon
       write(73,*)umax-psiv(i),fprof(psiv(i),2)
      end do
      close(73)
      call baleq
!      call chease_dat
      call mercier
!      call gfile
      call elite_data
      call elite2_data
!      call caxout
!      call iaea_plots
      if (nbal.ge.1) then
        allocate(eq1(2*npts*nturns),eq2(2*npts*nturns),eqd1(2*npts*nturns))
        allocate(chi(2*npts*nturns),fb(2*npts*nturns))
        allocate(lambda(ncon))
      end if
      if (nbal.eq.1) call balloon(ibal)
      if (nbal.gt.1) then
        write(6,*)' input min and max flux surfaces (1=edge)'
        read(5,*)n1,n2
        open(92,file='balloon')
!        call mercier
        do i=n1,n2
!          write(6,*)' in mervcier'
!          write(6,*)' out mercier'
          call balloon(i)
        end do
        close(92)
      end if
!      do i=1,ncon-1
!        call mercier(i)
!      end do
! plot out temperature gradient length scale...
      write(6,*)' doing loop'
      do i=2,ncon-1
        psi=psiv(i)
        psi1=psiv(i-1)
        eps1=epsv(i-1)
        psi2=psiv(i+1)
        eps2=epsv(i+1)
!        write(6,*)' eps1=',eps1,' eps2=',eps2,' temp1=',temp1,' temp2=',temp2,' temp0=',temp0
        temp1=tempi(psi1,1,0)
        temp2=tempi(psi2,1,0)
        temp0=tempi(psi,1,0)
        tlen(i)=-temp0*(eps2-eps1)/(temp2-temp1)
        epsn(i)=epsv(i)/tokeps
      end do
      do i=1,ncon
        psi=psiv(i)
        rat=psi/umax
        psinorm(i)=rat
        arr0(i)=vloop*af0*(1.-rat)**powj
        arr1(i)=((1.-rat)**(fpow1))*((psi/umax)**fpow2)
        scale=((fpow1/(fpow1+fpow2))**(fpow1))*((fpow2/(fpow2+fpow1))**fpow2)
        arr1(i)=af1*vloop*arr1(i)/scale
        arr2(i)=af2*vloop*((((1.-psi/umax)/0.9)**(powj/10))     &
             *(((psi/umax)/0.1)**(0.111*powj/10.)))
        arr3(i)=((1.-rat)**(fpow3))*((psi/umax)**fpow4)
         scale=((fpow3/(fpow3+fpow4))**(fpow3))*((fpow4/(fpow4+fpow3))**fpow4)
        arr3(i)=af3*vloop*arr3(i)/scale
        eps=epsv(i)
        call extj(eps,psi,extapp,extapp2,icur)
        arr4(i)=extapp2
      end do
!      call tstplt(ncon,psinorm,sfac,0.2d0,0.0d0)
!      call tstplt5(ncon,psinorm,arr0,arr1,arr2,arr3,arr4,sfac)
   end subroutine usrcal
