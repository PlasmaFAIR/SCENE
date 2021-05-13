module idball_output
  implicit none
contains
      subroutine baleq
!     ****************
!
!  outputs data file for ballooning stability code, idball
!
!----------------------------------------------------------------------
!
      use equilibrium, only : bp
      use param
      implicit none
      double precision ptab(ncon),fptab(ncon),ftab(ncon)
      double precision press,fprof
      double precision bvac0,bvacg,bet0,betg,b0,b0g
      double precision psimin,psimax,dpsi,psi,psih
      double precision fsi,pt
      double precision rat,rr,zz,bth
      double precision aidbal,alpha,dqdpsi,sh,qq
      integer i,j,k,krev,jr,n,jp
!
      integer mgf,mgx,mgy,nopts,ifail
      logical dpflag
      double precision xgrid(nr),ygrid(nz/2+1)
      double precision psig(nr*(nz/2+1)),qf(ncon),ff(ncon),dpdpsi(ncon)
      double precision r0gd,bt0gd
      double precision psival(nr*(nz/2+1)),rgr(nr*(nz/2+1)),zgr(nr*(nz/2+1))
!
     integer, allocatable:: LCELL(:,:), LNEXT(:)
     double precision:: XMIN, YMIN, DX, DY, RMAX
     double precision, allocatable:: RW(:), A(:,:)
     double precision, external:: CS2VAL
     integer:: NCC, NRR, NWW
     !
     logical :: debug
     
     debug = .false.
     
      open(66,file='alfsh.dat')
      pt=press(0.0d0,0)
      bvac0=mu0*rodi/(2.*pi*r0)
      bvacg=mu0*rodi/(2.*pi*rcen)
      bet0=2.*mu0*pt/bvac0**2
      betg=2.*mu0*pt/bvacg**2
      psimax=umax
      psimin=0.
      dpsi=(psimax-psimin)/(ncon-1)
      do k=1,ncon
    krev=ncon-k+1
    psi=psimax-(k-1)*dpsi
        psih=psimax-psi
    fsi=fprof(psi,2)
    ftab(krev)=fsi*1.0e6
    pt=press(psi,1)
    ptab(krev)=pt
        if (k.lt.ncon) then
          jr=0
          if (k.gt.1) then
            do 55 j=nr,1,-1
              if (r(j).lt.r0) goto 55
              if (psih.gt.u(j,nsym)) jr=j
 55         continue
            if (jr.eq.0) then
              write(6,*)'error in baleq***could not find flux surface'
              write(6,*)' psih=',psih,' umax=',umax
              stop
            end if
            rat=(psih-u(jr,nsym))/(u(jr-1,nsym)-u(jr,nsym))
            rr=r(jr)+rat*(r(jr-1)-r(jr))
          else
            rr=r2
          end if
          zz=z(nsym)
          bth=bp(rr,zz)
          aidbal=2.*(rr-r0)**2/bth
          alpha=-mu0*aidbal*pt
!  evaluate dqdpsi
          jp=1
          do 67 j=1,ncon
            if (psi.gt.psiv(j)) goto 67
            jp=j
 67       continue
          if (jp.eq.ncon) jp=jp-1
          rat=(psi-psiv(jp))/(psiv(jp+1)-psiv(jp))
          qq=sfac(jp)+rat*(sfac(jp+1)-sfac(jp))
          dqdpsi=qp(jp)+rat*(qp(jp+1)-qp(jp))
          sh=(rr-r0)*r0*bth*dqdpsi/qq
          write(66,15)psi/psimax,qq,alpha,sh
 15       format(' psi=',e13.5,' q=',f8.4,' alpha=',e13.5,' sh=',e13.5)
        end if
    n=1
    do 88 j=1,ncon
      if (psi.gt.psiv(j)) goto 88
      n=j
 88     continue
    if (n.eq.ncon) n=ncon-1
    rat=(psi-psiv(n))/(psiv(n+1)-psiv(n))
    qq=sfac(n)+rat*(sfac(n+1)-sfac(n))
    fptab(krev)=qq
      end do
      b0=mu0*rodi/(2.*pi*rcen)
      b0g=b0*1.0e4
      z(nsym)=0.
!!$      open(40,file='baleq.dat')
!!$!      write(40,'(a)') 'hagevl: no of x y points, coord of magax: '
!!$      write(40,*)nr,nz/2+1,ncon
!!$      write(40,*)'-------------------------------------'
!!$      write(40,64)(r(i)*100.,i=1,nr)
!!$      write(40,*)'-------------------------------------'
!!$      write(40,64)(z(i)*100.,i=nsym,nz)
!!$      write(40,*)'-------------------------------------'
!!$      write(40,64)((-u(i,j),i=1,nr),j=nsym,nz)
!!$      write(40,*)'-------------------------------------'
!!$      write(40,66)(fptab(i),ftab(i),i=1,ncon)
!!$      write(40,*)'-------------------------------------'
!!$      psimin=-umax
!!$      write(40,68)psimin,b0g,rcen*100.
!!$      write(40,*)'-------------------------------------'
!!$      write(40,64)(ptab(i),i=1,ncon)
!!$ 64   format(1p5e13.5)
!!$ 66   format(1p2e13.5)
!!$ 68   format(1p3e13.5)
!  paste in Appel's interpolator routine....
      mgx=nr
      mgy=nz/2+1
      mgf=ncon
      do i=1,mgx
        xgrid(i)=r(i)*100
      end do
      k=0
      do j=nsym,nz
        ygrid(j-nsym+1)=z(j)*100.
        k=k+1
      end do
      do i=1,mgx
        do j=nsym,nz
          k=j-nsym+1
          psig(mgy*(i-1)+k)=-u(i,j)
        end do
      end do
      do k=1,mgf
        qf(k)=fptab(k)
        ff(k)=ftab(k)
        dpdpsi(k)=ptab(k)
      end do
      bt0gd=b0g
      r0gd=rcen*100.
      dpflag=.true.
!
      k=0
      do 100 i=1,mgx
         do 90 j=1,mgy
            if (psig(mgy*(i-1)+j).lt.0.d0) then
               k=k+1
               psival(k)=psig(mgy*(i-1)+j)
               rgr(k)=xgrid(i)
               zgr(k)=ygrid(j)
            endif
 90      continue
 100  continue
!
      if (debug) print *,'make interpolant',mgx*mgy,k
      nopts=k
      ifail=0
      ! Setup for interpolation
      NCC=MIN(17,nopts-1)
      NWW=MIN(30,nopts-1)
      NRR=ceiling( sqrt(dble(nopts)/dble(3.0)) )
      allocate(LCELL(NRR,NRR))
      allocate(LNEXT(nopts))
      allocate(RW(nopts))
      allocate(A(9,nopts))
      call CSHEP2(nopts,rgr,zgr,psival,NCC,NWW,NRR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,ifail)
      if(ifail.NE.0) then
          write(6,*) 'CSHEP2 ifail not zero!'
          stop
      end if
      ! End of interpolation setup
!      call e01saf(nopts,rgr,zgr,psival,triang,grads,ifail)
!
      !print *,'reconstruct psi'
      do 200 i=1,mgx
         do 190 j=1,mgy
            ifail=1
!            call e01sbf(nopts,rgr,zgr,psival,triang,grads,    &
!                xgrid(i),ygrid(j),psig(mgy*(i-1)+j),ifail)
            psig(mgy*(i-1)+j)=CS2VAL(xgrid(i),ygrid(j),        &
                  nopts,rgr,zgr,psival,NRR,LCELL,LNEXT,        &
                  XMIN,YMIN,DX,DY,RMAX,RW,A)
 190     continue
 200  continue

      if (debug) print *,'write out'
!      open(unit=40,file='interp.baleq',status='new')
      open(unit=40,file='interp.baleq')
         write(40,*) mgx,mgy,mgf
!         write(40,*)' R***'
         write(40,*) (real(xgrid(i)),i=1,mgx)
!         write(40,*)' Z***'
         write(40,*) (real(ygrid(i)),i=1,mgy)
!         write(40,*)' psi***'
         write(40,*) ((real(psig(mgy*(i-1)+j)),i=1,mgx),j=1,mgy)
!         write(40,*)' q***'
         write(40,*) (real(qf(i)),real(ff(i)),i=1,mgf)
!         write(40,*)' psimin, etc***'
         write(40,*) real(psimin),real(bt0gd),real(r0gd)
!         write(40,*)' pprime'
      if (dpflag) then
         write(40,*) (dpdpsi(i),i=1,mgf)
      endif
      close (40)


      return
  end subroutine baleq
end module idball_output
