      subroutine chease_dat
!     *********************
!
!  Subroutine to write out data to Chease equilibrium code
!  H R Wilson 3/7/02
!
      use param
      implicit none
      double precision rbound(npts),zbound(npts)
      double precision fcsm(ncon),rppf(ncon),rfun(ncon)
      double precision psi,fprof,press,rz0,pedge0,b0exp,ffpval,psip,psim
      integer i,j,nsttp
!
      b0exp=mu0*rodi/(2.*pi*r0)
      open(unit=48,file=runname(1:lrunname)//'.chease', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.chease'
         stop
      endif
!  Write aspect ratio
      write(48,996)tokeps
!  Vertical position of magnetic axis (normalised)
      rz0=0.
      write(48,996)rz0
!  Edge pressure (normalised)
      pedge0=0.
      write(48,996)pedge0
!  nsttp tells CHEASE we are writing the ff' profile
      nsttp=1
!  Load up boundary....
      do i=1,npts
        rbound(i)=rpts(1,i)/r0
        zbound(i)=zpts(1,i)/r0
      end do
      write(48,998)npts
      write(48,999)(rbound(i),zbound(i),i=1,npts)
      write(48,998)ncon,nsttp
!  Load up s p-prime and ff-prime arrays
      do i=1,ncon
        j=ncon-i+1
        fcsm(j)=psiv(i)/umax
        psi=psiv(i)
        if (fcsm(j).lt.-1.0e-4) then
          write(6,*)' ERROR in chease output, psi outside range'
          write(6,*)' k=',i,' psiv=',psiv(i),' umax=',umax
          stop
        end if
        if (fcsm(j).gt.1.0001) then
          write(6,*)' ERROR in chease output, psi outside range'
          write(6,*)' k=',i,' psiv=',psiv(i),' umax=',umax
          stop
        end if
        if (fcsm(j).gt.1.) fcsm(j)=1.
        if (fcsm(j).lt.0.)fcsm(j)=0.
        fcsm(j)=sqrt(fcsm(j))
        rppf(j)=press(psi,1)*mu0*r0**2/b0exp
        rfun(j)=fprof(psi,1)/b0exp
      end do
      write(48,996)(fcsm(i),i=1,ncon)
      write(48,996)(rppf(i),i=1,ncon)
      write(48,996)(rfun(i),i=1,ncon)
      write(48,*)
      write(48,*)
      write(48,*)b0exp,'      B0EXP'
      write(48,*)r0,'      R0EXP'
      write(48,*)'psi_n, pressure (Nm**-2), f(psi)=R*B_phi (Tm)'
      do i=1,ncon
        psi=psiv(i)
        write(48,909)psi/umax,press(psi,0),fprof(psi,2)
      end do
      write(48,*)'psi_n, p-prime, f-prime'
      do i=1,ncon
        psi=psiv(i)
        write(48,909)psi/umax,press(psi,1),fprof(psi,3)
      end do
 996  format(e18.8)
 998  format(i5)
 999  format(2e18.8)
 909  format(3e18.8)
      close(48)
      write(84,*)' ipass=',ipass,' ipswtch=',ipswtch
      do i=2,ncon-1
        psi=psiv(i)
        psip=psiv(i+1)
        psim=psiv(i-1)
        ffpval=0.5*(fprof(psip,2)**2-fprof(psim,2)**2)/(psip-psim)
        write(84,*)' psi=',psi/umax,' ffp=',fprof(psi,1),' ffpval=',ffpval
        ffpval=fprof(psi,2)*(fprof(psip,2)-fprof(psim,2))/(psip-psim)
        write(84,*)' ffpval2=',ffpval
      end do
      return
  end subroutine chease_dat
