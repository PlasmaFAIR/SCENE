      subroutine elite_data
!     *********************
!
!  This subroutine provides data to the ELITE edge stability code in the
!  required format.  Note, ELITE requires data to be in CGS units!
!
      use param
      implicit none
      character(len=8) dums(10)
      character(len=12) ctitle
      double precision arr(ncon)
      double precision psi
      double precision press,fprof,dense,tempe,tempi
      integer ndsk,i,j,ncon1,ncstrt
!
      ncstrt=1
      ncon1=ncon-ncstrt+1
      ctitle='scene.dskbal'
      ndsk=20
      write(6,*)' opening file...'
      open(ndsk,file=ctitle)
!  A dummy string (must start 'iend....')
      do i=1,10
        dums(i)='noend   '
      end do
      do i=1,2
         write(ndsk,10) (dums(j),j=1,10)
      end do
      write(ndsk,10) ' &end   ',(dums(j),j=1,9)
      write(ndsk,10) (dums(j),j=1,10)
!  number of flux surfaces, and no. points on flux surfaces
      write(ndsk,15)ncon1,npts
!  dummy text
      dums(1)='psi'
      write(ndsk,10) (dums(i),i=1,10)
!  psi values...
      do i=1,ncon1
        arr(i)=1.0d8*psiv(ncon-i+1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='pprime'
      write(ndsk,10) (dums(j),j=1,10)
!  p-prime: note p variable in eqbm part of ELITE is actually p/(4*pi)!
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=1.0d-7*press(psi,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='f(psi)'
      write(ndsk,10) (dums(j),j=1,10)
!  f(psi)
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=1.0d6*fprof(psi,2)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='ffprime'
      write(ndsk,10) (dums(j),j=1,10)
!  f*f-prime(psi)
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=1.0d4*fprof(psi,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='chipsi'
      write(ndsk,10) (dums(j),j=1,10)
!  chipsi, corrects for non-constant psi-mesh
      do i=2,ncon1
        if (i.lt.ncon) then
          arr(i)=0.5*(psiv(ncon-i)-psiv(ncon-i+2))*(ncon-1.)*1.0d8
        else
!  Use quadratic approx for end point....
          arr(ncon)=2.*(0.75*psiv(1)-psiv(2)+0.25*psiv(3))*(ncon-1.)*1.0d8
        end if
      end do
      arr(1)=arr(2)
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='q'
      write(ndsk,10) (dums(j),j=1,10)
!  q(psi)
      do i=1,ncon1
        arr(i)=sfac(ncon-i+1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='R'
      write(ndsk,10) (dums(j),j=1,10)
!  R-pts of the flux surfaces
      do j=npts/2+1,1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=r0*100.
          else
            arr(i)=rpts(ncon-i+1,j)*100.
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
      do j=npts-1,npts/2+1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=r0*100.
          else
            arr(i)=rpts(ncon-i+1,j)*100.
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
!  dummy text
      dums(1)='Z'
      write(ndsk,10) (dums(j),j=1,10)
!  Z-pts of the flux surfaces
      do j=npts/2+1,1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=0.0d0
          else
            arr(i)=zpts(ncon-i+1,j)*100.
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
      do j=npts-1,npts/2+1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=0.0d0
          else
            arr(i)=zpts(ncon-i+1,j)*100.
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
!  dummy text
      dums(1)='ne'
      write(ndsk,10) (dums(j),j=1,10)
!  ne:
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=1.0d-6*dense(psi,0)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='nep'
      write(ndsk,10) (dums(j),j=1,10)
!  ne-prime:
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=1.0d-14*dense(psi,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='Te'
      write(ndsk,10) (dums(j),j=1,10)
!  Te:
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=tempe(psi,0)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='Tep'
      write(ndsk,10) (dums(j),j=1,10)
!  Te-prime
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=1.0d-8*tempe(psi,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dums(1)='Ti'
      write(ndsk,10) (dums(j),j=1,10)
!  Ti:
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=tempi(psi,1,0)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!      dums(1)='Tip'
      write(ndsk,10) (dums(j),j=1,10)
!  Ti-prime
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=1.0d-8*tempi(psi,1,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      close(ndsk)
      write(6,*)' file closed...'
 10   format(10a8)
 15   format(2i5)
 200  format(5g20.12)
!      do i=1,npts
!        write(801,*)rpts(1,i),zpts(1,i)
!      end do
   end subroutine elite_data
!
!----------------------------------------------------------------
      subroutine elite2_data
!     *********************
!
!  This subroutine provides data to the ELITE edge stability code in the
!  required format. New ELITE format 11/11/02, shape='eqbm'
!
      use param
      implicit none
      character(len=8) dum
      character(len=12) ctitle
      double precision arr(ncon)
      double precision psi
      double precision press,fprof,dense,tempe,tempi
      double precision x1,x2,x3,p1,p2,p3,aa,bb
      integer ndsk,i,j,ncon1,ncstrt
!
      ncstrt=1
      ncon1=ncon-ncstrt+1
      ndsk=21
      open(unit=ndsk,file=runname(1:lrunname)//'.elite', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.elite'
         stop
      endif
      write(ndsk,*)'Input data for ELITE, shape =eqbm format'
!  number of flux surfaces, and no. points on flux surfaces
      write(ndsk,15)ncon1,npts
!  psi values...
      dum='psi'
      write(ndsk,10) dum
      do i=1,ncon1
        arr(i)=psiv(ncon-i+1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dum='pprime'
      write(ndsk,10) dum
!  p-prime: note p variable in eqbm part of ELITE is actually p/(4*pi)!
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=mu0*press(psi,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
! Need derivative of p-prime
!  dummy text
      dum='ppprime'
      write(ndsk,10) dum
      do i=1,ncon1
        if (i.eq.1) then
          x1=psiv(ncon-2)
          x2=psiv(ncon-1)
          x3=psiv(ncon)
        else if (i.eq.ncon1) then
          x1=psiv(ncon-ncon1+1)
          x2=psiv(ncon-ncon1+2)
          x3=psiv(ncon-ncon1+3)
        else
          x1=psiv(ncon-i)
          x2=psiv(ncon-i+1)
          x3=psiv(ncon-i+2)
        end if
        p1=mu0*press(x1,1)
        p2=mu0*press(x2,1)
        p3=mu0*press(x3,1)
        aa=((p1-p2)/(x1-x2)-(p2-p3)/(x2-x3))/(x1-x3)
        bb=(p1-p2)/(x1-x2)-aa*(x1+x2)
        if (i.eq.1) then
          arr(i)=2.*aa*x3+bb
        else if (i.eq.ncon1) then
          arr(i)=2.*aa*x1+bb
        else
          arr(i)=2.*aa*x2+bb
        end if
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dum='f(psi)'
      write(ndsk,10)dum
!  f(psi)
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=fprof(psi,2)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dum='ffprime'
      write(ndsk,10) dum
!  f*f-prime(psi)
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=fprof(psi,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
! Need derivative of p-prime
!  dummy text
      dum='ffpprime'
      write(ndsk,10) dum
      do i=1,ncon1
        if (i.eq.1) then
          x1=psiv(ncon-2)
          x2=psiv(ncon-1)
          x3=psiv(ncon)
        else if (i.eq.ncon1) then
          x1=psiv(ncon-ncon1+1)
          x2=psiv(ncon-ncon1+2)
          x3=psiv(ncon-ncon1+3)
        else
          x1=psiv(ncon-i)
          x2=psiv(ncon-i+1)
          x3=psiv(ncon-i+2)
        end if
        p1=fprof(x1,1)
        p2=fprof(x2,1)
        p3=fprof(x3,1)
        aa=((p1-p2)/(x1-x2)-(p2-p3)/(x2-x3))/(x1-x3)
        bb=(p1-p2)/(x1-x2)-aa*(x1+x2)
        if (i.eq.1) then
          arr(i)=2.*aa*x3+bb
        else if (i.eq.ncon1) then
          arr(i)=2.*aa*x1+bb
        else
          arr(i)=2.*aa*x2+bb
        end if
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dum='q'
      write(ndsk,10) dum
!  q(psi)
      do i=1,ncon1
        arr(i)=sfac(ncon-i+1)
!        write(6,*)' i=',i,' ncon=',ncon,' q=',arr(i)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  dummy text
      dum='R'
      write(ndsk,10) dum
!  R-pts of the flux surfaces
      do j=npts/2+1,1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=r0
          else
            arr(i)=rpts(ncon-i+1,j)
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
      do j=npts-1,npts/2+1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=r0
          else
            arr(i)=rpts(ncon-i+1,j)
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
!  dummy text
      dum='Z'
      write(ndsk,10) dum
!  Z-pts of the flux surfaces
      do j=npts/2+1,1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=0.0d0
          else
            arr(i)=zpts(ncon-i+1,j)
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
      do j=npts-1,npts/2+1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=0.0d0
          else
            arr(i)=zpts(ncon-i+1,j)
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
!  dummy text
      dum='Bp'
      write(ndsk,10) dum
!  R-pts of the flux surfaces
      do j=npts/2+1,1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=0.
          else
            arr(i)=bppts(ncon-i+1,j)
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
      do j=npts-1,npts/2+1,-1
        do i=1,ncon1
          if (i.eq.1) then
            arr(i)=0.
          else
            arr(i)=bppts(ncon-i+1,j)
          end if
        end do
        do i=1,ncon1/5
          write(ndsk,200) arr((i*5)-4:i*5)
        end do
        write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      end do
!  ne:
      dum='ne'
      write(ndsk,10) dum
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=dense(psi,0)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  ne-prime:
      dum='neprime'
      write(ndsk,10) dum
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=dense(psi,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  Te:
      dum='Te'
      write(ndsk,10) dum
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=tempe(psi,0)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  Te-prime
      dum='Teprime'
      write(ndsk,10) dum
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=tempe(psi,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  Ti:
      dum='Ti'
      write(ndsk,10) dum
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=tempi(psi,1,0)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
!  Ti-prime
      dum='Tiprime'
      write(ndsk,10) dum
      do i=1,ncon1
        psi=psiv(ncon-i+1)
        arr(i)=tempi(psi,1,1)
      end do
      do i=1,ncon1/5
        write(ndsk,200) arr((i*5)-4:i*5)
      end do
      write(ndsk,200) arr((ncon1/5)*5+1:ncon1)
      close(ndsk)
      write(6,*)' file closed...'
 10   format(a8)
 15   format(2i5)
 200  format(5g20.12)
!      do i=1,npts
!        write(801,*)rpts(1,i),zpts(1,i)
!      end do
   end subroutine elite2_data
