      subroutine bandat
!     *****************
!
! writes out flux surface data file for use in bandit
!
      use param
      implicit none
      double precision fax,baxis,dpsi,psi,rat,rr,zz,bth,fsi,bphi,bt,bmod
      double precision psi1,psi2,rn1,ri1,rn2,ri2,vint,aint
      double precision drl,dzl,dl
      double precision psiban(ncon)
      double precision fprof,bp
      double precision rban(npts,ncon),zban(npts,ncon),btban(npts,ncon)
      double precision bban(npts,ncon),sban(npts,ncon)
      double precision circ(ncon),qban(ncon),vban(ncon),aban(ncon)
      integer iban(ncon)
      integer nban,k,ik,i,k1,k2
      real xp(npts),zp(npts),xmin,xmax,zmin,zmax
!      dimension btban(400,50),bban(400,50)
!      dimension sban(400,50),vban(50),aban(50),circ(50),qban(50)
!      dimension xp(400),y1(400),y2(400),y3(400)
!
      do i=1,npts
        rpts(ncon,i)=r0
        zpts(ncon,i)=0.
      end do
      psi=0.
      fax=fprof(psi,2)
      baxis=fax/r0
!  spacing for bandit flux surfaces (same as nfreya surfaces)
      dpsi=umax/ncon
!  no bandit flux surfaces
      nban=ncon
      do k=1,nban
    psiban(k)=k*dpsi
!  store coords of bandit flux surfaces in rban,zban
        ik=1
        do 10 i=1,ncon
          if (psiv(i).lt.psiban(k)) goto 10
          ik=i
 10     continue
        if (ik.eq.ncon) ik=ik-1
        rat=(psiban(k)-psiv(ik))/(psiv(ik+1)-psiv(ik))
        iban(k)=npts
        do i=1,npts
          rban(i,k)=rpts(ik,npts+1-i)+                               &
                    rat*(rpts(ik+1,npts+1-i)-rpts(ik,npts+1-i))
          zban(i,k)=zpts(ik,npts+1-i)+                               &
                   rat*(zpts(ik+1,npts+1-i)-zpts(ik,npts+1-i))
        end do
      end do
      do 30 k=1,nban
    psi=psiban(k)
    fsi=fprof(psi,2)
    do i=1,npts
      rr=rban(i,k)
      zz=zban(i,k)
      btban(i,k)=fsi/rr
      bth=bp(rr,zz)
      bphi=fsi/rr
      bban(i,k)=sqrt(bth*bth+bphi*bphi)
        end do
!  distance along field line and circumference of flux surface
    do i=1,npts
      if (i.eq.1) then
        sban(i,k)=0.
        circ(k)=0.
      else
        drl=rban(i,k)-rban(i-1,k)
        dzl=zban(i,k)-zban(i-1,k)
        dl=sqrt(drl**2+dzl**2)
        circ(k)=circ(k)+dl
        bmod=0.5*(bban(i,k)+bban(i-1,k))
        bt=0.5*(btban(i,k)+btban(i-1,k))
        bth=sqrt(bmod**2-bt**2)
        sban(i,k)=sban(i-1,k)+dl*bmod/bth
      end if
        end do
!  calculate areas and volumes between flux surfaces
    psi1=psiban(k)
    if (k.eq.1) then
      psi2=0.
    else
      psi2=psiban(k-1)
    end if
    k1=1
    k2=1
    do 35 ik=1,ncon
      if (psi1.gt.psiv(ik)) goto 35
      k1=ik
 35     continue
    do 36 ik=1,ncon
      if (psi2.gt.psiv(ik)) goto 36
      k2=ik
 36     continue
    if (k1.eq.ncon) k1=k1-1
    if (k2.eq.ncon) k2=k2-1
    rat=(psi1-psiv(k1))/(psiv(k1+1)-psiv(k1))
    rn1=rnorm(k1)+rat*(rnorm(k1+1)-rnorm(k1))
    ri1=rinv(k1)+rat*(rinv(k1+1)-rinv(k1))
    qban(k)=sfac(k1)+rat*(sfac(k1+1)-sfac(k1))
    rat=(psi2-psiv(k2))/(psiv(k2+1)-psiv(k2))
    rn2=rnorm(k2)+rat*(rnorm(k2+1)-rnorm(k2))
    ri2=rinv(k2)+rat*(rinv(k2+1)-rinv(k2))
    dpsi=abs(psi2-psi1)
    vint=0.5*(rn1+rn2)
    aint=0.5*(ri2*rn2+ri1*rn1)
    vban(k)=2.*pi*vint*dpsi
    aban(k)=aint*dpsi
 30   continue
      open(50,file='flxinf.dat')
      write(50,*)'no flux surfaces'
      write(50,1030)nban
      write(50,*)'current, r0(mag axis) b0(mag axis)'
      write(50,1010)cur,r0,baxis
      write(50,*)'no. points on flux surface'
      write(50,1040)(iban(k),k=1,nban)
      write(50,*)'distance along field line'
      write(50,1020)((sban(i,k),i=1,iban(k)),k=1,nban)
      write(50,*)'r coords'
      write(50,1020)((rban(i,k),i=1,iban(k)),k=1,nban)
      write(50,*)'z coords'
      write(50,1020)((zban(i,k),i=1,iban(k)),k=1,nban)
      write(50,*)'bphi'
      write(50,1020)((btban(i,k),i=1,iban(k)),k=1,nban)
      write(50,*)'b'
      write(50,1020)((bban(i,k),i=1,iban(k)),k=1,nban)
      write(50,*)'psi'
      write(50,1020)(psiban(k),k=1,nban)
      write(50,*)'area'
      write(50,1020)(aban(k),k=1,nban)
      write(50,*)'volume'
      write(50,1020)(vban(k),k=1,nban)
      write(50,*)'safety factor'
      write(50,1020)(qban(k),k=1,nban)
      write(50,*)'circumference of flux surface'
      write(50,1020)(circ(k),k=1,nban)
 1010 format(1p3e13.5)
 1020 format(1p5e13.5)
 1030 format(i5)
 1040 format(8i5)
      close(50)
!!$      call pspace(0.05,0.95,0.05,0.95)
!!$      zmin=z(1)
!!$      zmax=z(nz)
!!$      xmin=r(1)
!!$      xmax=r(nr)
!!$      call map(xmin,xmax,zmin,zmax)
!!$      do k=1,ncon
!!$        do i=1,npts
!!$          xp(i)=rban(i,k)
!!$          zp(i)=zban(i,k)
!!$        end do
!!$        call ptjoin(xp,zp,1,npts,0)
!!$      end do
!!$      call frame
      return
      end subroutine bandat
!
