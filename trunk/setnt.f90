module setnt_mod
  implicit none
contains
      subroutine setnt
!     *****************
!
!   plots the flux contours around which variables are integrated.
!   ....and the ecrh footprint
!
      use param
      implicit none
      integer :: i,jgot,j,ishot,ntdat
      double precision :: pltr(1000),pltt(1000),pltn(1000)
      double precision :: psi_new(ncon),ne_new(ncon),te_new(ncon),p_new(ncon)
      double precision :: rr1,rr2,ne1,te1,ne2,te2,rat,rr
      double precision :: timsht


!  read in temperature and density profiles from data file
     open(20)
     read(20,*)ishot,timsht,ntdat
     if (ntdat.gt.1000) then
       write(6,*)' ERROR**** need to increase dimensions of'
       write(6,*)' pltr,pltt and pltn in setnt and graphs'
       stop
     end if
     do i=1,ntdat
       read(20,*)pltr(i),pltn(i),pltt(i)
     end do
     close(20)
! Automatic matching
! Update mesh of ne, Te, and ensure ipswtch=3
!     ipswtch=3
     deallocate(psi_m,ne_m,te_m,p_m)
     nterp=ncon
     allocate(psi_m(nterp),ne_m(nterp),te_m(nterp),p_m(nterp))
     do i=1,ncon
       psi_new(i)=psiv(i)/umax
       if (psi_new(i).ge.1.) then
!  edge values
         psi_new(i)=1.
         ne_new(i)=0.001
         te_new(i)=0.001
       else if (psi_new(i).le.0.) then
          psi_new(i)=0.
          rr=r0
          jgot=1
!  Find TS mesh point:
          do 29 j=1,ntdat
            if (pltr(j).gt.rr) goto 29
            jgot=j
 29       continue
          if (jgot.eq.ntdat) jgot=jgot-1
          rat=(rr-pltr(jgot))/(pltr(jgot+1)-pltr(jgot))
          te_new(i)=pltt(jgot)+rat*(pltt(jgot+1)-pltt(jgot))
          ne_new(i)=pltn(jgot)+rat*(pltn(jgot+1)-pltn(jgot))
       else
! find R-value of
         rr1=rpts(i,1)
         rr2=rpts(i,npts/2+1)
         jgot=1
          do 31 j=1,ntdat
            if (pltr(j).gt.rr1) goto 31
            jgot=j
 31       continue
          if (jgot.eq.ntdat) jgot=jgot-1
          rat=(rr1-pltr(jgot))/(pltr(jgot+1)-pltr(jgot))
          te1=pltt(jgot)+rat*(pltt(jgot+1)-pltt(jgot))
          ne1=pltn(jgot)+rat*(pltn(jgot+1)-pltn(jgot))
          jgot=1
          do 32 j=1,ntdat
            if (pltr(j).gt.rr2) goto 32
            jgot=j
 32       continue
          if (jgot.eq.ntdat) jgot=jgot-1
          rat=(rr2-pltr(jgot))/(pltr(jgot+1)-pltr(jgot))
          te2=pltt(jgot)+rat*(pltt(jgot+1)-pltt(jgot))
          ne2=pltn(jgot)+rat*(pltn(jgot+1)-pltn(jgot))
          te_new(i)=0.5*(te1+te2)
          ne_new(i)=0.5*(ne1+ne2)
       end if
     end do
     do i=1,nterp
       p_new(i)=bk*te_new(i)*ne_new(i)
     end do
     do i=1,nterp
       te_m(i)=te_new(i)
       p_m(i)=p_new(i)
       ne_m(i)=ne_new(i)
!!$       write(6,*)' i=',i
!!$       psi_m(i)=psi_new(i)
!!$       te_m(i)=te_new(i)
!!$       p_m(i)=0.
!!$       nsmooth=0
!!$       ilo=i-nsmooth
!!$       if (ilo.lt.1) ilo=1
!!$       iup=i+(i-ilo)
!!$       if (iup.gt.nterp) then
!!$         iup=nterp
!!$         ilo=i-(iup-i)
!!$       end if
!!$       do j=ilo,iup
!!$!         p_m(i)=p_m(i)+p_new(j)*pscl/(scl*bpol*umax)
!!$         p_m(i)=p_m(i)+p_new(j)
!!$       end do
!!$       if (iup.gt.ilo) p_m(i)=p_m(i)/(iup-ilo+1)
!!$       ne_m(i)=p_m(i)/(2.*bk*te_m(i))
!!$       write(6,*)' psi=',psi_m(i),' te=',te_m(i),' ne=',ne_m(i)
     end do
!
!
      return
end subroutine setnt
end module setnt_mod
