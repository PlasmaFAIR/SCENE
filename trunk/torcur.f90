      subroutine torcur(icur)
!     ***********************
!
!  This subroutine calculates the toroidal current profiles of the
!  input, bootstrap, diamagnetic and pfirsch-schluter components.
!
      use param
      implicit none
!
      integer i,j,k,ik,icur,jcur
      double precision rr,zz,psi,bth,pd,fsi,ffd,fd,bphi,bsq,bmod
      double precision alintp,press,fprof,bp
      double precision rat,bsmean,bstrap,eps,bra,absj
      double precision extapp,extapp2,spitc,extapp3
!
      if (icont.gt.-3) then
        allocate( bsph(nr,nz),psph(nr,nz),diph(nr,nz),exph(nr,nz),   &
                exph2(nr,nz),gradj(nr,nz),spit(nr,nz),absph(nr,nz) )
      end if
!
      do i=1,nr
        do j=1,nz
          if (ixout(i,j).le.0) then
!  zero current outside plasma...
            bsph(i,j)=0.
            psph(i,j)=0.
            diph(i,j)=0.
            exph(i,j)=0.
            gradj(i,j)=0.
            exph2(i,j)=0.
          else
!  calculate currents inside plasma...
            rr=r(i)
            zz=z(j)
            psi=umax-u(i,j)
            bth=bp(rr,zz)
            pd=press(psi,1)
            fsi=fprof(psi,2)
            ffd=fprof(psi,1)
            fd=ffd/fsi
            bphi=fsi/rr
            bsq=bphi*bphi+bth*bth
            bmod=sqrt(bsq)
!  extrapolate between flux surfaces:
	    ik=0
	    if (psi.gt.umax) then
	      write(nw,*)'error***psi value greater than maximum'
	      write(nw,*)'problem in torcur'
	      stop
	    end if
	    do 5 k=1,ncon
	      if (psi.ge.psiv(k)) goto 5
	      ik=ik+1
 5          continue
	    if (ik.eq.0) then
	      ik=1
	      rat=0.
	    else if (ik.lt.ncon) then
	      rat=(psi-psiv(ik))/(psiv(ik+1)-psiv(ik))
	    else
	      write(nw,*)'error***cannot interpolate a psi value'
	      write(nw,*)'problem in torcur'
	      stop
	    end if
	    bsmean=bsqav(ik)+rat*(bsqav(ik+1)-bsqav(ik))
	    bstrap=bsj(ik)+rat*(bsj(ik+1)-bsj(ik))
	    absj=ajdotb(ik)+rat*(ajdotb(ik+1)-ajdotb(ik))
	    eps=epsv(ik)+rat*(epsv(ik+1)-epsv(ik))
!   bootstrap...
	    bsph(i,j)=bstrap*bphi/sqrt(bsmean)
! and alphs contribution
	    absph(i,j)=absj*bphi/bsmean
	    bra=1.-bsq/bsmean
!   pfirsch-schluter...
	    psph(i,j)=-fsi*fsi*pd*bra/(rr*bsq)
!   diamagnetic...
	    diph(i,j)=-rr*bth*bth*pd/bsq
!   grad--shafranov current (should converge to total j, as specified)
!  (note, if itot=1 total specified current = gradj -> no convergence
!  required).
	    gradj(i,j)=-rr*pd-ffd/(rr*mu0)
!  externally applied current
	    if (itot.eq.0) then
              icur=1
              if (neo.lt.0) icur=-1
!  external current profile calculated in extj, returned in extapp
	      call extj(eps,psi,extapp,extapp2,icur)
	      if (abs(neo).eq.1) then
		jcur=-1
		call extj(eps,psi,spitc,extapp3,jcur)
	      end if
	      exph(i,j)=extapp*fsi/rr
              exph2(i,j)=extapp2*fsi/rr
	      if (abs(neo).eq.1) spit(i,j)=spitc*fsi/rr
	    else
!  total current profile given in gradj
	      exph(i,j)=gradj(i,j)-psph(i,j)-diph(i,j)-bsph(i,j)
              exph2(i,j)=0.
	    end if
	  end if
        end do
      end do
  end subroutine torcur
!
!**********************************************************************
