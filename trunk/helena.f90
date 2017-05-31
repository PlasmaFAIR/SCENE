      subroutine helena
!     *****************
!
!  Code provides output to reconstruct input data for HELENA equilibrium code
!  11/4/01
!
      use param
      implicit none
      integer nh,ishape,mharm,imesh,ias,igam,ipai,npr1,npl1,nchi
      integer i,n1,n2,n3,n4,nh1
      double precision quadh,rhalf,zup,arg
      double precision press,fprof,amix,rat
      double precision psi,dpsi,pp0,ffp0,bv0,bb,gampsi,xv(ncon),yv(ncon),yv2(ncon)
!
! Evaluate quadracity
      zup=elon*rcen*tokeps
      rhalf=-1.
      do i=npts/2+1,1,-1
        if (rhalf.lt.0.) then
          if (zpts(1,i).gt.zup/2.) then
            rat=(zup/2.-zpts(1,i))/(zpts(1,i+1)-zpts(1,i))
            rhalf=rpts(1,i)+rat*(rpts(1,i+1)-rpts(1,i))
            write(6,*)' rat=',rat,' rhalf=',rhalf,' rptsi=',rpts(1,i),' rptsi-1=',rpts(1,i-1)
          end if
        end if
      end do
      rhalf=rhalf-rcen
      if ((rhalf.lt.0.).or.(rhalf.gt.(rcen*tokeps))) then
        write(6,*)' Error in calculating quadh in helena'
        stop
      end if
      arg=rhalf/(tokeps*rcen)
      quadh=(2./sqrt(3.))*(acos(arg)-pi/6.-tri/2.)
      write(6,*)' in helena'
      psi=0.
      pp0=press(psi,1)
      ffp0=fprof(psi,1)
      bb=mu0*rcen**2*pp0/ffp0
!      bb=2.*tokeps*bb/(1.+bb)
      nh=110
      open(unit=nh,file=runname(1:lrunname)//'.helena', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.helena'
         stop
      endif
      nh1=120
      open(unit=nh1,file=runname(1:lrunname)//'.profs', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.profs'
         stop
      endif
      write(nh1,*)' normalised flux, dp/dpsi, fdfdpsi'
      write(nh1,*)ncon
      do i=ncon,1,-1
        psi=psiv(i)
        write(nh1,55)psi/umax,press(psi,1),fprof(psi,1)
      end do
 55   format(3e14.6)
      close(nh1)
!      open(nh,file='helena.dat')
      ishape=1
      mharm=128
      imesh=1
!      imesh=0
      ias=0
      write(nh,*)' Equilibrium data for HELENA'
      write(nh,*)' PPF SCENE'
      write(nh,*)' '
      write(nh,*)'&SHAPE'
      write(nh,10)elon,tri,quadh,mharm,ishape,imesh,ias
 10   format('ELLIP= ',F6.4,', TRIA=',F6.4,', QUAD=',F6.4,', MHARM=',I3,   &
             ', ISHAPE=',I3,', IMESH=',I2,', IAS=',I1,',')
      write(nh,*)'&END'
      write(nh,*)'&PROFILE'
      igam=7
!      igam=4
      ipai=7
      write(nh,20)igam,ipai,ncon
 20   format('IGAM=',I4,', IPAI=',I4,', NPTS=',I5,',')
      dpsi=umax/(ncon-1)
      do i=1,ncon
        psi=(i-1)*dpsi
        xv(i)=sqrt(psi/umax)
        yv(i)=fprof(psi,1)
        yv2(i)=press(psi,1)
        gampsi=(mu0*rcen**2*press(psi,1)+fprof(psi,1))/(mu0*rcen**2*pp0+ffp0)
        if (i.lt.10) then
          write(nh,30)i,press(psi,1)/pp0,i,fprof(psi,1)/ffp0
!          write(nh,30)i,press(psi,1)/pp0,i,gampsi
        else if (i.lt.100) then
          write(nh,40)i,press(psi,1)/pp0,i,fprof(psi,1)/ffp0
!          write(nh,40)i,press(psi,1)/pp0,i,gampsi
        else
          write(nh,50)i,press(psi,1)/pp0,i,fprof(psi,1)/ffp0
!          write(nh,50)i,press(psi,1)/pp0,i,gampsi
        end if
      end do
 30   format('DPR(',I1,')=',E14.5,', DF2(',I1,')=',E14.5,',')
 40   format('DPR(',I2,')=',E14.5,', DF2(',I2,')=',E14.5,',')
 50   format('DPR(',I3,')=',E14.5,', DF2(',I3,')=',E14.5,',')
      write(nh,*)'&END'
      write(nh,*)' &PHYS'
      bv0=mu0*rodi/(2.*pi*rcen)
      write(nh,60)tokeps,bb,mu0*cur/(tokeps*rcen*bv0)
 60   format(' EPS=',F8.4,', B=',E14.5,', XIAB=',E14.5)
      write(nh,*)'&END'
      write(nh,*)'&NUM'
      nchi=1+2**8
      n1=51
      n2=(2**6)+1
      n3=101
      n4=500
      amix=0.9
      write(nh,70)n1,n2,n3,nchi,nchi,n4,amix
 70   format(' NR=',i4,', NP=',i4,', NRMAP=',i4,', NPMAP=',i4,   &
             ', NCHI=',i4,', NITER=',i4,' AMIX=',f7.4)
      write(nh,*)'&END'
      npr1=0
      npl1=1
      write(nh,*)'&PRI'
      write(nh,*)' NPR1=   ',npr1
      write(nh,*)'&END'
      write(nh,*)'&PLOT'
      write(nh,*)' NPL1=   ',npl1
      write(nh,*)'&END'
      write(nh,*)'&BALL'
      write(nh,*)'&END'
      write(nh,*)'&MERC'
      write(nh,*)'&END'
!       call tstplt(ncon,xv,yv)
!      call tstplt(ncon,xv,yv2)
!
      write(nh,*)' '
      write(6,*)' out helena'
   end subroutine helena
