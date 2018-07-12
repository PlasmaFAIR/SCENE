      subroutine mercier
!     **********************
!
!  Calculates various flux surface averages, including the Mercier
!  coefficient, on flux surface label nf: 1<nf<ncon
!
      use param
      implicit none
      integer i,j,k,nf
      double precision qd,rr,bth,psi,fsi,ant,dm,bsq,pd,dres,hfun
      double precision fprof,press,pp,neoterm(ncon),hegterm(ncon),alfs
      double precision av1,av2,av3,av4,av5,av6,bvacu,bpol2,    &
                       glassterm(ncon),volavp,bdel
      double precision marr1(npts),marr2(npts),marr3(npts),marr4(npts)   &
                       ,marr5(npts),marr6(npts)

      do nf=1,ncon
      qd=qp(nf)
      psi=psiv(nf)
      fsi=fprof(psi,2)
      pd=mu0*press(psi,1)
      pp=press(psi,0)
      do i=1,npts
!        write(6,*)' i=',i
        rr=rpts(nf,i)
        bth=bppts(nf,i)
        bsq=bth**2+(fsi/rr)**2
        marr1(i)=fsi/(rr**2*bth**3)
        marr2(i)=marr1(i)*bsq
        marr3(i)=marr1(i)/bsq
        marr4(i)=(pd+0.5*dbsqdpsi(nf,i))/(bsq*bth)
        marr5(i)=2.*pi*nup(nf,i)/(bsq*circumf(nf))
        marr6(i)=bsq/bth
      end do
      call flxint(marr1,nf,ant)
      av1=ant/(2.*pi)
      call flxint(marr2,nf,ant)
      av2=ant/(2.*pi)
      call flxint(marr3,nf,ant)
      av3=ant/(2.*pi)
      call flxint(marr4,nf,ant)
      av4=ant/(2.*pi)
      call flxint(marr5,nf,ant)
      av5=ant/(2.*pi)
      call flxint(marr6,nf,ant)
      av6=ant/(2.*pi)
!      write(6,*)' av1=',av1,' av2-',av2,' av3=',av3
!      write(6,*)' av4=',av4,' av5=',av5
! Mercier interchange parameter
      dm=(pd/(fsi*qd**2))*(fsi*qd*av1-fsi*pd*av1**2+       &
         av2*(fsi*pd*av3+2.*av4-fsi*av5))
      !write(6,*)' psi=',psi/psiv(1),' Dm=',dm
! GGJ H-function
      hfun=pd*(av1-av2/bsqav(nf))/qd
! resistive interchange parameter, DR
      dres=dm-hfun+hfun**2
      write(716,*)' psi=',psi/psiv(1),' DM=',dm,' DR=',dres
      bvacu=mu0*rodi/(2.*pi*rcen)
      volavp=betexp/(2.*mu0*100./(bvacu*bvacu))
      bpol2=2.*mu0*pp*(circumf(1)/(mu0*cur))**2
!      write(6,*)' bpol2=',bpol2
      neoterm(nf)=3.2*bsj(nf)*sqrt(bsqav(nf))*sfac(nf)*bpol2/(qp(nf)*fsi*pp)
      glassterm(nf)=6.*dres
! Evaluate modification due to Hegna (PoP 6 (1999) 3980)
      alfs=0.5+sqrt(0.25-dm)
      hegterm(nf)=glassterm(nf)/(alfs-hfun)
      write(82,*)' psiN=',psi/psiv(1),' q=',sfac(nf),' bsterm=',neoterm(nf), &
              ' glass term=',glassterm(nf),' ratio=',neoterm(nf)/glassterm(nf)
      write(82,*)' psiN=',psi/psiv(1),' q=',sfac(nf),' bsterm=',neoterm(nf),&
                ' hegna term=',hegterm(nf),' ratio=',neoterm(nf)/hegterm(nf)
      write(82,*)' Lbs=',mu0*bsj(nf)*sqrt(bsqav(nf))/(fsi*pd)
      write(82,*)' DR=',dres,' dm=',dm
      if ((nf.gt.1).and.(nf.lt.ncon)) write(82,*)' qp=',qp(nf),' qp diff=',(sfac(nf+1)-sfac(nf-1))/(psiv(nf+1)-psiv(nf-1))
      write(82,*)' Lp/psi0=',mu0*pp/(pd*umax),' bpol2=',bpol2,' Lq/psi0=',sfac(nf)/(qp(nf)*umax)
      bdel=-av2*2.*pi*rnorm(nf)*mu0*bsj(nf)*sqrt(bsqav(nf))/    &
            (4.*fsi*qd*av6*pi**2)
        !write(6,*)' f90 boot/old boot=',neoterm(nf)/bdel
      end do
!      call tstplt2(ncon,sfac,neoterm,hegterm)
      return
  end subroutine mercier
