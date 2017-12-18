      subroutine fastbs
!     *****************
!
!  Calculates the bootstrap current arising from alpha-particles
!  Calculation by C T Hsu et al PFB 4 (1992) 4023
!
      use param
      use coldat
      implicit none
!
      double precision tempe,tempi,dense,densi,fprof,press
      double precision fcp,fcd,psi
! Alpha particle properties...
      double precision amass, alen,zal,v0,vc
      double precision te,ti,ted,tid,tau,ne19,ni19,nid,ned,ne,ni
      double precision ablog,sigv,srce,coolog,tslowsgn,sigd,srcd
      double precision zeff,falnum,tslow
      double precision f1d,f2d,ss,arg
      double precision sgn,pdal,fl0,fl1,fl2
      double precision x1d,x2d,denom,fmua,fsi,pq
      double precision femu,bot
      integer l,nh,k
!
      nh=229
      open(unit=nh,file=runname(1:lrunname)//'.lynton', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.lynton'
         stop
      endif
      !write(6,*)' writing fast ptcle info for lynton'
      write(nh,*)' psi_norm, vc, S*taus/(4*pi); f_alpha=(S*taus/(4*pi))/(v**3+vc**3)H(v0-v)'
!  alpha mass (relative to proton mass)
      amass=4.
!  alpha energy (in mev)
      alen=3.5
!  alpha charge
      zal=2.
!  alpha birth speed
      v0=1.383e7*sqrt(alen/amass)
      do k=1,ncon-1
        psi=psiv(k)
        ss=psi/umax
        fcp=1.-ftrap(k)
        fcd=1.-ftrapd(k)
!  electron and ion temperature (in ev)
        ti=tempi(psi,1,0)
        te=tempe(psi,0)
        tid=tempi(psi,1,1)
        ted=tempe(psi,1)
        tau=te/ti
!  alpha critical speed
        vc=3.979d4*sqrt(te)/amass**(1./3.)
!  and densities (ni=nd+nt=2*nd)
        ne19=1.0d-19*dense(psi,0)
        ni19=1.0d-19*densi(psi,1,0)
        ned=1.0d-19*dense(psi,1)
        nid=1.0d-19*densi(psi,1,1)
!  dt fusion cross-section
        ablog=(abs(log(1.45d-5*ti)))**2.25
        sigv=9.07d-22*exp(-0.476*ablog)
        srce=0.25d+38*ni19*ni19*sigv
!  coulomb logarithm (take same for ions and electrons)
        coolog=24.-log(1.0d6*sqrt(10.*ne19)/te)
!  alpha slowing down time
        tslow=6.316d-5*amass*te**1.5/(zal**2*ne19*coolog)
! alpha disbn fun = (s tau_s/4 pi)/(v^3+v_c^3) (v<v_0)
        falnum=srce*tslow/(4.*pi)
        write(nh,10)ss,vc,falnum
!  zeff
       	zeff=zm
        if (imp.eq.1) then
          ne=ne19*1.0d19
	  if (ne.gt.0.) then
	    zeff=0.
	    do l=1,nimp+1
	      ni=densi(psi,l,0)
	      zeff=zeff+(ni*iz(l)**2)/ne
            end do
          end if
        end if
        pq=3.*zeff/5.
!  derivative of cross section
        sgn=1.
        if (ti.lt.1./1.45d-5) sgn=-1.
        sigd=-sgn*(tid/ti)*1.071*(abs(log(1.45e-5*ti)))**1.25
        srcd=2.*nid/ni19+sigd
! Note T' terms from f1d are taken into f2d...f1d and f2d are the coefficients
! in dP_alpha/dpsi of the two velocity integrals
        f1d=srcd-ned/ne19
        f2d=1.5*ted/te
!  alpha pressure gradient (wrt psi)...returned in pdal
!  velocity-space integrals for trapped alpha fraction
!  returned in fl1 and fl2
        call pdgen(f1d,f2d,falnum,vc,v0,amass,pdal,fl0,fl1,fl2)
!  fast particle pressure
        fastp(k)=4.*pi*amass*mp*falnum*fl1*v0**2
        x1d=-v0*(srcd+1.5*ted/te-ned/ne19)
        x2d=-1.5*v0*ted/te
        denom=fl1*x1d+fl2*x2d
        fmua=fl1*x1d*fl1/(pq*fl1+4.*fcp*fl0)
        fmua=fmua-fl2*fl2*x2d/(pq*fl2+4.*fcp*fl1)
        fmua=1.-fcd+pq*(fcd-fcp)*fmua/denom
        fsi=fprof(psi,2)
!  Add in electron screening factor:
        femu=(mu1(k)*(mu3(k)+2.*(nuee(k)+13./(8.*nuei(k))))- &
             mu2(k)*(mu2(k)-3./(2.*nuei(k))))                      
        bot=((mu1(k)+nuei(k))*                               &
                (mu3(k)+2.*(nuee(k)+13.*nuei(k)/8.))    &
              -(mu2(k)-1.5*nuei(k))**2)
        femu=femu/bot
        ajdotb(k)=-fmua*fsi*pdal*(1.-(zal/zeff)*(1.-femu))
! If ifast=1 add fast particle bootstrap to thermal
        if (ifast.eq.1) bsj(k)=bsj(k)+ajdotb(k)/sqrt(bsqav(k))
      end do
! On axis bootstrap current is neglected...
      ajdotb(ncon)=0.
 10   format(3e14.6)
      return
  end subroutine fastbs
!
!**********************************************************************
      subroutine pdgen(f1d,f2d,falnum,vc,v0,amass,pdal,fl0,fl1,fl2)
!     *************************************************************
!
!   Calculates the psi-derivative of the fusion alpha-particle pressure (PDAL)
!   and velocity space integrals required for effective trapped alpha
!   fraction (FL1 and FL2)
!
      implicit none
      double precision rat,v0,vc,vmin,dv,v,pdal
      double precision f1,f2,f3,f4,fl0,fl1,fl2
      double precision f1d,f2d,falnum,amass
      integer nv,i

      rat=v0/vc
      vmin=v0
      if (vc.lt.vmin) vmin=vc
      dv=vmin/500
      nv=v0/dv
      dv=v0/(nv-1)
      f1=0.0d0
      f2=0.0d0
      f3=0.0d0
      f4=0.0d0
      do i=1,nv-1
        v=(i-1)*dv+dv/2.
        f1=f1+dv*v**4/(v**3+vc**3)
        f2=f2+dv*v**7/(v**3+vc**3)**2
        f3=f3+dv*v**4/(v**3+vc**3)**2
        f4=f4+dv*v**4
      end do
      fl0=f4*rat**3/(3.*v0**5)
      fl1=f1/(3.*v0**2)
      fl2=f3*v0/(3.*rat**3)
      pdal=7.0e-27*(f1d*f1+f2d*f2)*falnum*amass
      return
   end subroutine pdgen







