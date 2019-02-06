      subroutine flxplt
!     *****************
!
!   plots the flux contours around which variables are integrated.
!   ....and the ecrh footprint
!
      use param
      implicit none
      character(len=8) ctim
      character(len=12) cdat
      character(len=35) string
      character(len=30) string2
      character(len=2) string1
      real xp(nr*nz),yp(nr*nz),clevls(50),ftr(4),ftz(4)
      real rdim(nr*nz),zdim(nr*nz)
      real rlo,rup,zlo,zup,ypos,zl,zu,rl,ru,zav
      real zr(nr),zz(nz),elong,spsi
      real sval,dum
      real xhalf,yhalf,xx,yx,dzsep,aa,bb
      double precision psi,qaxis
      integer i,nc,itst,lmax,l,ipos,iq,nh,nreadpts,nk,k,ndmr
      real, dimension(:), allocatable:: rbdr,zbdr
!!$!
!!$!  corners of ecrh footprint
!!$      if (iso.eq.1) then
!!$        ftr(1)=rf-drf/2.
!!$        ftr(2)=rf-drf/2.
!!$        ftr(3)=rf+drf/2.
!!$        ftr(4)=rf+drf/2.
!!$        ftz(1)=zf-dzf/2.
!!$        ftz(2)=zf+dzf/2.
!!$        ftz(3)=zf+dzf/2.
!!$        ftz(4)=zf-dzf/2.
!!$      end if
! Read in free boundary
      nh=39
!      open(unit=nh,file='hidelta04_rob4_2.dat', &
!!$      open(unit=nh,file='step15_fiesta.txt', &
!!$           status='unknown',iostat=ios)
!!$      if(ios.ne.0) then
!!$         write(6,*) 'problem opening rz8694.dat'
!!$         stop
!!$      endif
!!$      read(nh,*)nreadpts
!!$      if (nreadpts.gt.(nr*nz)) then
!!$        write(6,*)' error*****nead to define larger rdim and zdim arrays in flxplt'
!!$        write(6,*)' nr*nz=',nr*nz,' nreadpts=',nreadpts
!!$        stop
!!$      end if
!!$      do i=1,nreadpts
!!$        read(nh,*)rdim(i),zdim(i)
!!$!        write(6,*)' i=',i,' R=',rdim(i),' Z=',zdim(i)
!!$      end do
!*************************************************
      if (ibdry.eq.2) then
        open(unit=nh,file='bdy.txt',status='unknown',iostat=ios)
        if(ios.ne.0) then
          write(6,*) 'problem opening bdy.txt for boundary R,Z'
          stop
        endif
        read(nh,*)nreadpts
        allocate(rbdr(nreadpts),zbdr(nreadpts))
        do i=1,nreadpts
          read(nh,*)rbdr(i),zbdr(i)
        end do
        close(nh)
      end if
      qaxis=sfac(ncon)
      iq=int(qaxis)+1
!!$!
!  Write out all input parameters:
      call pspace(0.85,1.3,0.1,0.92)
      call map(0.,1.,0.,1.)
      call border
      call ctrmag(15)
      call undlin(1)
      call lincol(2)
      call plotcs(0.15,0.95,'Input Parameters')
      call lincol(0)
      call undlin(0)
      call ctrmag(12)
      call plotcs(0.02,0.9,' rcen=')
      sval=sngl(rcen)
      call typenf(sval,3)
      call typecs(' eps=')
      sval=sngl(tokeps)
      call typenf(sval,3)
      call typecs(' elon=')
      sval=elon
      call typenf(sval,3)
      call plotcs(0.02,0.87,' tri=')
      sval=tri
      call typenf(sval,3)
      call typecs(' step=')
      sval=sngl(step)
      call typenf(sval,3)
      call typecs(' powj=')
      sval=sngl(powj)
      call typenf(sval,3)
      call plotcs(0.02,0.84,' fpow=')
      sval=sngl(fpow)
      call typenf(sval,3)
      call typecs(' fpo1=')
      sval=sngl(fpow1)
      call typenf(sval,3)
      call typecs(' fpo2=')
      sval=sngl(fpow2)
      call typenf(sval,3)
      call plotcs(0.02,0.81,' af0=')
      sval=sngl(af0)
      call typenf(sval,3)
      call typecs(' af1=')
      sval=sngl(af1)
      call typenf(sval,3)
      call typecs(' af2=')
      sval=sngl(af2)
      call typenf(sval,3)
      call plotcs(0.02,0.78,' ipsw=')
      call typeni(ipswtch,2)
      call typecs(' ppow=')
      sval=sngl(ppow)
      call typenf(sval,3)
      call typecs(' p0=')
      sval=sngl(p0)
      call typenf(sval,3)
      call typecs(' pa=')
      sval=sngl(pa)
      call typenf(sval,3)
      sval=sngl(pped)
      call plotcs(0.02,0.75,' pped=')
      call typenf(sval,3)
      call typecs(' pedg=')
      sval=sngl(pedg)
      call typenf(sval,3)
      call plotcs(0.02,0.72,' npow=')
      sval=sngl(nipow)
      call typenf(sval,3)
      call typecs(' ni0=')
      sval=sngl(ni0*1.0e-19)
      call typenf(sval,3)
      call typecs(' nia=')
      sval=sngl(nia*1.0e-19)
      call typenf(sval,3)
      call plotcs(0.02,0.69,' nped=')
      sval=sngl(niped*1.0e-19)
      call typenf(sval,3)
      call typecs(' nedg=')
      sval=sngl(niedg)
      call typenf(sval,3)
      call plotcs(0.02,0.66,' tpoe=')
      sval=sngl(tpoe)
      call typenf(sval,3)
      call typecs(' te0=')
      sval=sngl(te0)
      call typenf(sval,3)
      call typecs(' tea=')
      sval=sngl(tea)
      call typenf(sval,3)
      call plotcs(0.02,0.63,' tepd=')
      sval=sngl(teped)
      call typenf(sval,3)
      call typecs(' tedg=')
      sval=sngl(teedg)
      call typenf(sval,3)
      call plotcs(0.02,0.60,' tpoi=')
      sval=sngl(tpoi)
      call typenf(sval,3)
      call typecs(' ti0=')
      sval=sngl(ti0)
      call typenf(sval,3)
      call typecs(' tia=')
      sval=sngl(tia)
      call typenf(sval,3)
      call plotcs(0.02,0.57,' tipd=')
      sval=sngl(tiped)
      call typenf(sval,3)
      call typecs(' tidg=')
      sval=sngl(tiedg)
      call typenf(sval,3)
      call plotcs(0.02,0.54,' cur=')
      sval=sngl(cur*1.0e-6)
      call typenf(sval,3)
      call typecs(' paux=')
      if (neo.eq.1) then
        sval=sngl(paux-vloop*(totex-totex2)*1.0d-6)
      else
        sval=sngl(paux)
      end if  
      call typenf(sval,3)
      call typecs(' bpol=')
      sval=sngl(bpol)
      call typenf(sval,4)
      call plotcs(0.02,0.51,' rodi=')
      sval=sngl(rodi*1.0d-6)
      call typenf(sval,3)
      call typecs(' imat=')
      call typeni(imat,3)
      call typecs(' scl=')
      sval=sngl(inscl)
      call typene(sval,3)
      call plotcs(0.02,0.48,' frac=')
      sval=sngl(frac)
      call typenf(sval,3)
      call typecs(' zm=')
      sval=sngl(zm)
      call typenf(sval,1)
      call typecs(' zmai=')
      sval=sngl(zmai)
      call typenf(sval,2)
      call typecs(' omeg=')
      sval=sngl(omega)
      call typenf(sval,3)
      call plotcs(0.02,0.45,' imp=')
      call typeni(imp,2)
      call typecs(' itot=')
      call typeni(itot,2)
      call typecs(' neo=')
      call typeni(neo,2)
      call typecs(' nco=')
      call typeni(nco,2)
      call typecs(' ncon=')
      call typeni(ncon,3)
      call plotcs(0.02,0.42,' npts=')
      call typeni(npts,5)
      call typecs(' nout=')
      call typeni(nouter,5)
      call typecs(' ninn=')
      call typeni(ninner,2)
      call typecs(' npas=')
      call typeni(npass,2)
      call plotcs(0.02,0.39,' erit=')
      sval=sngl(errit)
      call typene(sval,3)
      call typecs(' erfd=')
      sval=sngl(errffd)
      call typene(sval,3)
      call typecs(' icon=')
      call typeni(icont,2)
      call plotcs(0.02,0.36,' dil=')
      sval=sngl(dil)
      call typenf(sval,4)
      call typecs(' pfac=')
      sval=sngl(pfac)
      call typenf(sval,4)
      call typecs(' ate=')
      sval=sngl(ate)
      call typenf(sval,4)
      call typecs(' tpo1=')
      sval=sngl(tpo1)
      call typenf(sval,4)
      call plotcs(0.02,0.33,' kval=')
      sval=sngl(kval)
      call typenf(sval,6)
      call typecs(' ibdr=')
      call typeni(ibdry,2)
      call typecs(' quad=')
      sval=sngl(quad)
      call typenf(sval,5)
      call typecs(' psic=')
      sval=psic
      call typenf(sval,4)
      call plotcs(0.02,0.3,' af3=')
      sval=sngl(af3)
      call typenf(sval,4)
      call typecs(' af4=')
      sval=sngl(af4)
      call typenf(sval,4)
      call typecs(' fpo3=')
      sval=sngl(fpow3)
      call typenf(sval,4)
      call typecs(' fpo4=')
      sval=sngl(fpow4)
      call typenf(sval,4)
      call ctrmag(15)
      call lincol(2)
      call plotcs(0.02,0.26,' Impurity information')
      call lincol(0)
      call ctrmag(12)
      call plotcs(0.02,0.23,' Number of impurities=')
      call typeni(nimp,2)
      call ctrmag(10)
      if (nimp.gt.0) then
        ypos=0.2
        do i=2,nimp+1
          call plotcs(0.02,ypos,' Z=')
          call typeni(iz(i),2)
          call typecs(' M/Mp=')
          sval=sngl(zmas(i))
          call typenf(sval,2)
          ypos=ypos-0.025
          call plotcs(0.02,ypos,' at,T0,Ta,Tped,Tpedg=')
          sval=sngl(ztpow(i))
          call typenf(sval,2)
          sval=sngl(zt0(i))
          call typenf(sval,2)
          sval=sngl(zta(i))
          call typenf(sval,2)
          sval=sngl(ztped(i))
          call typenf(sval,2)
          sval=sngl(ztedg(i))
          call typenf(sval,2)
          ypos=ypos-0.025
          call plotcs(0.02,ypos,' an,n0,na,nped,npedg=')
          sval=sngl(znpow(i))
          call typenf(sval,2)
          sval=sngl(zn0(i)*1.0d-19)
          call typenf(sval,2)
          sval=sngl(zna(i)*1.0d-19)
          call typenf(sval,2)
          sval=sngl(znped(i)*1.0d-19)
          call typenf(sval,2)
          sval=sngl(znedg(i))
          call typenf(sval,2)
          ypos=ypos-0.025
        end do
      end if
!
!---------------------------------------------------------------------
!
!!$	  if (itst.eq.0) then
!!$	    call plotcs(0.03,ypos,' edge')
!!$	  end if
!!$	  itst=(ncon-nc)/1
!!$	  if (itst.eq.0) then
!!$	    call plotcs(0.03,ypos,' centre')
!!$	  end if
!!$	  call ctrmag(15)
!!$          call plotne(0.4,ypos,spsi,4)
      elong=sngl(2.*zs/(r2-r1))
      do i=1,nr
        zr(i)=sngl(r(i))
      end do
      do i=1,nz
        zz(i)=sngl(z(i))
      end do
      do 20 nc=ncon,1,-1
        psi=psiv(nc)
        spsi=sngl(psi)
!       qq=sfac(nc)
        itst=1*((nc-1)/1)-(nc-1)
        if (itst.eq.0) then
          zlo=0.16
          zup=0.75
          rup=0.81
          rlo=rup-(zup-zlo)/elong
          call pspace(rlo,rup,zlo,zup)
          call map(zr(1),zr(nr),zz(1),zz(nz))
          call border
          call axorig(zr(1),zz(1))
          call ctrmag(10)
          call axessi(0.,0.)
          call ctrmag(12)
          lmax=npts
          do  l=1,lmax
            xp(l)=sngl(rpts(nc,l))
            yp(l)=sngl(zpts(nc,l))
          end do
	  if (sfac(nc).gt.iq) then
	    iq=iq+1
	    if (iq.le.3) then 
              call thick(2)
              call lincol(2)
            end if
	  end if
!	  if (modb.eq.-1) then
!	    if (nc.ne.1) goto 15
!	  end if
          if (nc.eq.1) call thick(5)
!          if (iso.eq.1) call ptjoin(ftr,ftz,1,4,-1)
          call ptjoin(xp,yp,1,lmax,-1)
          call lincol(2)
          call broken(2,10,2,10)
!          call ptjoin(rbdr,zbdr,1,nreadpts,-1)
          call full
!!$          if (nc.eq.1) then
!!$!  Set up boundary (Preinhalter data)
!!$          r1=0.2338
!!$          r2=1.3269
!!$          xx=0.53
!!$          yx=1.11**2
!!$          xhalf=0.297
!!$          yhalf=0.588**2
!!$          aa=((xhalf-r1)/yhalf-(xx-r1)/yx)/(yhalf-yx)
!!$          bb=(xx-r1-aa*yx**2)/yx
!!$          nreadpts=20
!!$          dzsep=sqrt(yx)/(nreadpts-1)
!!$          k=0
!!$          do i=1,nreadpts
!!$            k=k+1
!!$            zdim(k)=(i-1)*dzsep
!!$            rdim(k)=r1+aa*zdim(k)**4+bb*zdim(k)**2
!!$            write(6,*)' k=',k,' R=',rdim(k),' Z=',zdim(k)
!!$          end do
!!$          write(6,*)'**************'
!!$          xhalf=1.12
!!$          yhalf=0.552**2
!!$          aa=((xhalf-r2)/yhalf-(xx-r2)/yx)/(yhalf-yx)
!!$          bb=(xx-r2-aa*yx**2)/yx
!!$          nreadpts=20
!!$          dzsep=sqrt(yx)/(nreadpts-1)
!!$          do i=2,nreadpts
!!$            k=k+1
!!$            zdim(k)=sqrt(yx)-(i-1)*dzsep
!!$            rdim(k)=r2+aa*zdim(k)**4+bb*zdim(k)**2
!!$            write(6,*)' k=',k,' R=',rdim(k),' Z=',zdim(k)
!!$          end do
!!$          write(6,*)' *************'
!!$          nk=k
!!$          do i=1,nk-1
!!$            k=k+1
!!$            rdim(k)=rdim(nk-i)
!!$            zdim(k)=-zdim(nk-i)
!!$            write(6,*)' k=',k,' R=',rdim(k),' Z=',zdim(k)
!!$          end do
!!$          nk=k
!!$          call lincol(2)
!!$          call broken(10,5,10,5)
!!$          call ptjoin(rdim,zdim,1,nk,-1)
!!$          call full
!!$          end if
 15       call thick(1)
          call lincol(0)
!!$          call pspace(0.85,1.2,0.1,0.87)
!!$          call map(0.,1.,0.,1.)
!!$          call border
!!$          call ctrmag(15)
!!$          ipos=nc/1
!!$          ypos=1.-(ipos+4)*0.04
!!$	  call undlin(1)
!!$          call plotcs(0.15,0.95,'poloidal flux surfaces')
!!$!	  call plotcs(0.15,0.95,'safety factor')
!!$          call undlin(0)
!!$          call ctrmag(12)
!!$          call pcscen(0.5,0.88,'t m')
!!$	  call supfix
!!$          call typecs('2')
!!$	  call normal
!!$	  itst=nc/1
!!$	  if (itst.eq.0) then
!!$	    call plotcs(0.03,ypos,' edge')
!!$	  end if
!!$	  itst=(ncon-nc)/1
!!$	  if (itst.eq.0) then
!!$	    call plotcs(0.03,ypos,' centre')
!!$	  end if
!!$	  call ctrmag(15)
!!$          call plotne(0.4,ypos,spsi,4)
!	  if (nc.lt.ncon+1) call plotnf(0.4,ypos,qq,2)
	  zl=zlo-0.06
	  zu=zlo-0.03
	  call pspace(0.75,0.81,zl,zu)
	  call map(0.,1.,0.,1.)
	  call pcscen(0.5,0.5,' R (m)')
	  rl=rlo-0.055
	  ru=rl+0.03
	  zav=(zlo+zup)/2.
	  zl=zav-0.2
	  zu=zav+0.2
	  call pspace(rl,ru,zl,zu)
	  call map(0.,1.,0.,1.)
	  call ctrori(90.)
	  call pcscen(0.5,0.5,' Z (m)')
	  call ctrori(0.)
	end if
 20   continue
      call pspace(0.15,1.3,0.1,0.92)
      call border
!  title
      call pspace(0.15,0.85,0.8,0.92)
      call border
      call map(0.,1.,0.,1.)
      call ctrmag(26)
      call thick(2)
      call plotcs(0.35,0.8,title)
      call thick(1)
      call enqdat(cdat)
      call enqtim(ctim)
      call ctrmag(16)
      call plotcs(0.01,0.5,' SCENE 6.1 ')
      call typecs('queries: H.R. Wilson, Culham ext 6902')
      write(string,fmt='('' Run time:'',a8,'' on '',a12)') ctim,cdat
!      call pcsend(0.98,0.6,string)
      call plotcs(0.01,0.2,string)
!!$      call thick(2)
!!$      call plotcs(0.25,0.25,title)
!!$      call thick(1)
!      if (abs(modb).eq.1) then
!!  mod b contours
!	zlo=0.16
!	zup=0.75
!	rup=0.81
!	rlo=rup-(zup-zlo)/elong
!	call pspace(rlo,rup,zlo,zup)
!	call map(zr(1),zr(nr),zz(1),zz(nz))
!	do 50 i=1,ncon
!	   clevls(i)=bmin+(i-1.)*(bmax-bmin)/(ncon-1.)
! 50     continue
!	call contra(bmoda,1,nr,400,1,nz,400,clevls,1,10)
!      end if
      call frame

      elong=sngl(2.*zs/(r2-r1))
      do i=1,nr
        zr(i)=sngl(r(i))
      end do
      do i=1,nz
        zz(i)=sngl(z(i))
      end do
      do 147 nc=ncon,1,-1
        psi=psiv(nc)
        spsi=sngl(psi)
!       qq=sfac(nc)
        itst=1*((nc-1)/1)-(nc-1)
        if (itst.eq.0) then
          zlo=0.05
          zup=0.95
          rup=0.81
          rlo=rup-(zup-zlo)/elong
          call pspace(rlo,rup,zlo,zup)
          call map(zr(1),zr(nr),zz(1),zz(nz))
          call border
          call ctrori(90.)
          call axorig(zr(1),zz(1))
          call ctrmag(10)
          call axessi(0.,0.)
          call ctrori(0.)
          call ctrmag(12)
          lmax=npts
          do  l=1,lmax
            xp(l)=sngl(rpts(nc,l))
            yp(l)=sngl(zpts(nc,l))
          end do
	  if (sfac(nc).gt.iq) then
	    iq=iq+1
	    if (iq.le.3) then 
              call thick(2)
              call lincol(2)
            end if
	  end if
!	  if (modb.eq.-1) then
!	    if (nc.ne.1) goto 15
!	  end if
          if (nc.eq.1) call thick(5)
!          if (iso.eq.1) call ptjoin(ftr,ftz,1,4,-1)
          call ptjoin(xp,yp,1,lmax,-1)
          call lincol(2)
          call broken(2,10,2,10)
!          call ptjoin(rbdr,zbdr,1,nreadpts,-1)
          call full
!!$          if (nc.eq.1) then
!!$!  Set up boundary (Preinhalter data)
!!$          r1=0.2338
!!$          r2=1.3269
!!$          xx=0.53
!!$          yx=1.11**2
!!$          xhalf=0.297
!!$          yhalf=0.588**2
!!$          aa=((xhalf-r1)/yhalf-(xx-r1)/yx)/(yhalf-yx)
!!$          bb=(xx-r1-aa*yx**2)/yx
!!$          nreadpts=20
!!$          dzsep=sqrt(yx)/(nreadpts-1)
!!$          k=0
!!$          do i=1,nreadpts
!!$            k=k+1
!!$            zdim(k)=(i-1)*dzsep
!!$            rdim(k)=r1+aa*zdim(k)**4+bb*zdim(k)**2
!!$            write(6,*)' k=',k,' R=',rdim(k),' Z=',zdim(k)
!!$          end do
!!$          write(6,*)'**************'
!!$          xhalf=1.12
!!$          yhalf=0.552**2
!!$          aa=((xhalf-r2)/yhalf-(xx-r2)/yx)/(yhalf-yx)
!!$          bb=(xx-r2-aa*yx**2)/yx
!!$          nreadpts=20
!!$          dzsep=sqrt(yx)/(nreadpts-1)
!!$          do i=2,nreadpts
!!$            k=k+1
!!$            zdim(k)=sqrt(yx)-(i-1)*dzsep
!!$            rdim(k)=r2+aa*zdim(k)**4+bb*zdim(k)**2
!!$            write(6,*)' k=',k,' R=',rdim(k),' Z=',zdim(k)
!!$          end do
!!$          write(6,*)' *************'
!!$          nk=k
!!$          do i=1,nk-1
!!$            k=k+1
!!$            rdim(k)=rdim(nk-i)
!!$            zdim(k)=-zdim(nk-i)
!!$            write(6,*)' k=',k,' R=',rdim(k),' Z=',zdim(k)
!!$          end do
!!$          nk=k
!!$          call lincol(2)
!!$          call broken(10,5,10,5)
!!$          call ptjoin(rdim,zdim,1,nk,-1)
!!$          call full
!!$          end if
          call thick(1)
          call lincol(0)
!!$          call pspace(0.85,1.2,0.1,0.87)
!!$          call map(0.,1.,0.,1.)
!!$          call border
!!$          call ctrmag(15)
!!$          ipos=nc/1
!!$          ypos=1.-(ipos+4)*0.04
!!$	  call undlin(1)
!!$          call plotcs(0.15,0.95,'poloidal flux surfaces')
!!$!	  call plotcs(0.15,0.95,'safety factor')
!!$          call undlin(0)
!!$          call ctrmag(12)
!!$          call pcscen(0.5,0.88,'t m')
!!$	  call supfix
!!$          call typecs('2')
!!$	  call normal
!!$	  itst=nc/1
!!$	  if (itst.eq.0) then
!!$	    call plotcs(0.03,ypos,' edge')
!!$	  end if
!!$	  itst=(ncon-nc)/1
!!$	  if (itst.eq.0) then
!!$	    call plotcs(0.03,ypos,' centre')
!!$	  end if
!!$	  call ctrmag(15)
!!$          call plotne(0.4,ypos,spsi,4)
!	  if (nc.lt.ncon+1) call plotnf(0.4,ypos,qq,2)
	  zl=zlo-0.06
	  zu=zlo-0.03
	  call pspace(0.75,0.81,zl,zu)
	  call map(0.,1.,0.,1.)
	  call pcscen(0.5,0.5,' R (m)')
	  rl=rlo-0.055
	  ru=rl+0.03
	  zav=(zlo+zup)/2.
	  zl=zav-0.2
	  zu=zav+0.2
	  call pspace(rl,ru,zl,zu)
	  call map(0.,1.,0.,1.)
	  call ctrori(90.)
	  call pcscen(0.5,0.5,' Z (m)')
	  call ctrori(0.)
	end if
 147   continue
!      call pspace(0.15,1.3,0.1,0.92)
!      call border
      call frame

   end subroutine flxplt
!
!**********************************************************************
!
      subroutine graphs
!     *****************
        !
!    Plots various profiles and outputs run parameters.
!
      use param
      implicit none
      character(len=8) ctim
      character(len=9) cdat
      character(len=30) string
      character(len=34) txt
      integer nrob
      integer i,j,l,lmid,ii,ig,jj,ix,lp,icur
      double precision press,dense,densi,fprof,tempe,tempi,bp
      double precision rstar,zstar,tric,elonc,rhalf,arg,quadc
      double precision psi,rat,bt0,bv0,pcen,ne0
      double precision curgs,rin,rout,rr
      double precision fcop,rc,rescop,pden,fsi
      double precision pdedg,ffdedg,extapp,extapp2
      double precision bth,bphi,bsq,paxis,baxis,bigs,dl
      integer ishot,ntdat
      real seplo,sepup,plo,pup,ffplo,ffpup,pplo,ppup
      real seps(ncon),sffp(ncon),sp(ncon),spp(ncon)
      real sxlo,sxup,sylo,syup,sy1lo,sy1up,sprup
      real sx(nr+5),sx1(nr+5),sy(nr+5),sy1(nr+5),sy2(nr+5),sy3(nr+5),spr(nr+5)
      real sxx(1000),syy(1000)
      real sxpos,sxposup,sypos
      real svar,szero,rlo,rup
      real dilo,diup,bslo,bsup,pslo,psup,tolo,toup
      double precision timsht
      double precision jedge,jaxis,pd,fd,fval
      double precision psilo,psiup,plod,pupd,flo,fup,rlod,rupd
      double precision q95,bthlo,bthup,betlo,betup,rv
      double precision dpsi,ffp
      double precision psi_new(ncon),ne_new(ncon),te_new(ncon)
      double precision rr1,rr2,ne1,te1,ne2,te2
      double precision tgot,ngot,psilst,qgot,rpos,zz
      integer jgot,nh
      real pltr(1000),pltt(1000),pltn(1000),pltter(1000),pltner(1000), &
           pltp(1000)
      real tmax,rnmax,te,ne,rmin,rmax,dum,jrat,pmax
      real alf95
      integer isk,na
!
!
!      call tester
!--------------------------------------------------------------------------
      na=122
      open(unit=na,file=runname(1:lrunname)//'.cairns', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.cairns'
         stop
      endif
     write(na,*)' R(m), ne (m**-3), B (T), J_total, J_aux (kA m**-2), T (keV)'
     do 50 i=1,nr
       if (ixout(i,nsym).le.0) goto 50
       psi=umax-u(i,nsym)
       rr=r(i)
       zz=0.
       sy(i)=sngl(gradj(i,nsym)/1000.)
       sy1(i)=sngl(exph(i,nsym)/1000.)
       sy2(i)=sngl(sqrt((fprof(psi,2)/rr)**2+(bp(rr,zz))**2))
       write(na,14)rr,dense(psi,0),sy2(i),sy(i),sy1(i),tempe(psi,0)/1000.
 50  continue
 14  format(6e14.6)
     close(na)
     write(6,*)' in graphs'
     rin=rpts(1,1)
     rout=rpts(1,npts/2+1)
!  density/temperature profile
      lmid=0
      j=0
      do 10 i=1,nr
        if (ixout(i,nsym).le.0) goto 10
        if (r(i).gt.r0) then
!  insert axis value
          if (j.eq.0) then
            lmid=lmid+1
            j=1
            sx(lmid+1)=sngl(r0)
            psi=0.
            sy(lmid+1)=sngl(dense(psi,0))
            sy1(lmid+1)=sngl(tempe(psi,0)/1000.)
            sy2(lmid+1)=sngl(tempi(psi,1,0)/1000.)
            spr(lmid+1)=sngl(press(psi,0))
            write(6,*)' Central pressure=',press(psi,0),' M**-2'
            if (ipswtch.eq.-1) then
              sy(lmid+1)=sngl(press(psi,0))
            end if
          end if
        end if
        lmid=lmid+1
        psi=umax-u(i,nsym)
        sx(lmid+1)=sngl(r(i))
        sy(lmid+1)=sngl(dense(psi,0))
        sy1(lmid+1)=sngl(tempe(psi,0)/1000.)
        sy2(lmid+1)=sngl(tempi(psi,1,0)/1000.)
        spr(lmid+1)=sngl(press(psi,0))
        if (ipswtch.eq.-1) then
          sy(lmid+1)=sngl(press(psi,0))
        end if
 10   continue
      lmid=lmid+2
      sx(1)=sngl(rin)
      sx(lmid)=sngl(rout)
      psi=psiv(1)
      sy(1)=sngl(dense(psi,0))
      sy1(1)=sngl(tempe(psi,0)/1000.)
      sy2(1)=sngl(tempi(psi,1,0)/1000.)
      spr(1)=sngl(press(psi,0))
      if (ipswtch.eq.-1) then
        sy(1)=sngl(press(psi,0))
      end if
      sy(lmid)=sngl(dense(psi,0))
      sy1(lmid)=sngl(tempe(psi,0)/1000.)
      sy2(lmid)=sngl(tempi(psi,1,0)/1000.)
      spr(lmid)=sngl(press(psi,0))
      write(73,*)' R(m)       p(Nm**-2)'
      sprup=0.
      do i=1,lmid
        write(73,*)sx(i),spr(i)
        if (spr(i).gt.sprup) sprup=spr(i)
      end do
      nh=112
      open(unit=nh,file=runname(1:lrunname)//'.rob', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.rob'
         stop
      endif
      write(nh,*)' number of data points=',lmid
      write(nh,*)' R(m), ne(m**-3), Te(keV)'
      do i=1,lmid
        write(nh,11)sx(i),sy(i),sy1(i)
      end do
 11   format(3e14.6)
      if (ipswtch.eq.-1) then
        sy(lmid)=sngl(press(psi,0))
      end if
      call pspace(0.11,0.45,0.47,0.81)
      sxlo=sx(1)
      sxup=sx(lmid)
      call range(sy,lmid,sylo,syup)
      syup=sylo+1.3*(syup-sylo)
      sylo=0.
      call map(sxlo,sxup,sylo,syup)
      sxpos=sxlo+0.5*(sxup-sxlo)
      sxposup=sxlo+0.75*(sxup-sxlo)
      sypos=0.95*syup
      call ctrmag(12)
      if (ipswtch.eq.-1) then
        call positn(sxpos,sypos)
        call join(sxposup,sypos)
        call plotcs(sxposup,sypos,' pressure')
        call axorig(sxlo,0.)
        call ctrmag(6)
        call yaxisi(0.)
        call ctrmag(14)
        call curveo(sx,sy,1,lmid)
        call border
        call axorig(sxlo,0.)
        call ctrmag(6)
        call axessi(0.,0.)
        call ctrmag(14)
      else
        call lincol(0)
        call positn(sxpos,sypos)
        call broken(10,5,10,5)
        call join(sxposup,sypos)
        call full
        call plotcs(sxposup,sypos,' e temp.')
        sypos=0.9*syup
        call lincol(4)
        call positn(sxpos,sypos)
        call broken(2,10,2,10)
        call join(sxposup,sypos)
        call full
        call plotcs(sxposup,sypos,' i temp.')
        print*, 'Sxpos is: ',sxpos
        print*, 'Sypos is: ',sypos
        print*, 'Syposup is: ',sxposup
        sypos=0.85*syup
        call lincol(0)
        call positn(sxpos,sypos)
        call join(sxposup,sypos)
        call plotcs(sxposup,sypos,' e dens.')
        sypos=0.81*syup
        call lincol(3)
        call positn(sxpos,sypos)
        call broken(10,5,2,5)
        call join(sxposup,sypos)
        call full
        call plotcs(sxposup,sypos,' press')
        call lincol(0)
        call axorig(sxlo,0.)
        call ctrmag(6)
        call yaxisi(0.)
        call ctrmag(14)
        call curveo(sx,sy,1,lmid)
        syup=sngl(1.3*te0/1000.)
        sy1up=sngl(1.3*ti0/1000.)
        if (sy1up.gt.syup) syup=sy1up
        call map(sxlo,sxup,0.,syup)
        call lincol(2)
        call broken(10,5,10,5)
        call curveo(sx,sy1,1,lmid)
        call lincol(0)
        call broken(2,10,2,10)
        call curveo(sx,sy2,1,lmid)
        call lincol(0)
        call full
        call axorig(sxup,0.)
        call border
        call ctrmag(6)
        call xaxisi(0.)
        call ctrmag(14)
        call annotp(0,1)
        call ctrmag(6)
        call yaxisi(0.)
        call ctrmag(14)
        call annotp(0,0)
        call axorig(0.,0.)
!  pressure
        syup=1.3*sprup
        call map(sxlo,sxup,0.,syup)
        call lincol(3)
        call broken(10,5,2,5)
        call curveo(sx,spr,1,lmid)
        call lincol(0)
        call full
      end if
      call pspace(0.05,0.07,0.47,0.81)
      call map(0.,1.,0.,1.)
      call ctrori(90.)
      sypos=0.5
      sxpos=0.4
      if (ipswtch.eq.-1) then
        call pcscen(sxpos,sypos,'pressure')
      else
        call pcscen(sxpos,sypos,'e dens. m')
        call supfix
        call typecs('-3')
        call normal
        call pspace(0.465,0.535,0.47,0.81)
        call map(0.,1.,0.,1.)
        sxpos=0.5
        call pcscen(sxpos,sypos,' temp.  kev')
      end if
      call ctrori(0.)
      call pspace(0.38,0.42,0.42,0.44)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call pcscen(sxpos,sypos,'R(m)')
      psi=0.5*umax
      dpsi=0.001*umax
      ffp=fprof(psi,2)*(fprof(psi+dpsi,2)-fprof(psi-dpsi,2))/(2.*dpsi)
      write(6,*)'CHECK****ffp-num=',ffp,' ffp=',fprof(psi,1)
      write(6,*)' ffp on axis=',fprof(0.0d0,1)
      ffp=(press(psi+dpsi,0)-press(psi-dpsi,0))/(2.*dpsi)
      write(6,*)'CHECK***p-prim-num=',ffp,' p-prome=',press(psi,1)
!-------------------------------------------------------------------------
!  Collisionality plot
      call pspace(0.57,0.91,0.47,0.81)
      call border
      sxlo=sx(1)
      sxup=sx(lmid)
      lmid=0
      j=0
      do 15 i=1,nr
        psi=umax-u(i,nsym)
        if (psi.gt.psiv(2)) goto 15
        if (r(i).gt.r0) then
!  insert axis value
          if (j.eq.0) then
            lmid=lmid+1
            j=1
            sx(lmid)=sngl(r0)
            psi=0.
            sy(lmid+1)=sngl(tnue(ncon))
          end if
        end if
        lmid=lmid+1
        ig=1
        do 5 ii=1,ncon
          if (psiv(ii).lt.psi) goto 5
          ig=ii
 5      continue
        if (ig.eq.ncon) ig=ig-1
        rat=(psi-psiv(ig))/(psiv(ig+1)-psiv(ig))
        sy(lmid+1)=sngl(tnue(ig)+rat*(tnue(ig+1)-tnue(ig)))
 15   continue
      lmid=lmid+2
      sy(1)=sngl(tnue(1))
      sy(lmid)=sngl(tnue(1))
      sxlo=sx(1)
      sxup=sx(lmid)
      call range(sy,lmid,sylo,syup)
      if (syup.gt.5.) syup=5.
      if (sylo.gt.0.) sylo=0.
      call map(sxlo,sxup,sylo,syup)
      call axorig(sxlo,sylo)
      call ctrmag(8)
      call axessi(0.,0.)
      call ctrmag(14)
      call ptjoin(sx,sy,1,lmid,1)
      call pspace(0.51,0.535,0.47,0.81)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call ctrori(90.)
      call pcscen(sxpos,sypos,' ')
      call ctrfnt(2)
      call typenc(110)
      call suffix
      call typenc(250)
      call ctrfnt(0)
      call typecs('e')
      call normal
      call ctrori(0.)
      call pspace(0.85,0.9,0.42,0.44)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call pcscen(sxpos,sypos,' R (m)')
!-------------------------------------------------------------------------
!  Plot of toroidal field
      lmid=0
      j=0
      do 30 i=1,nr
        if (ixout(i,nsym).le.0) goto 30
        if (r(i).gt.r0) then
!  insert axis value
          if (j.eq.0) then
            lmid=lmid+1
            j=1
            sx(lmid)=sngl(r0)
            psi=0.
            sy(lmid)=sngl(fprof(psi,2)/r(i))
            psi=umax
            sy1(lmid)=sngl(fprof(psi,2)/r(i))
          end if
        end if
        lmid=lmid+1
        psi=umax-u(i,nsym)
        sx(lmid)=sngl(r(i))
        sy(lmid)=sngl(fprof(psi,2)/r(i))
        psi=umax
        sy1(lmid)=sngl(fprof(psi,2)/r(i))
 30   continue
      sxlo=sx(1)
      sxup=sx(lmid)
      call pspace(0.11,0.45,0.08,0.42)
      call range(sy,lmid,sylo,syup)
      call range(sy1,lmid,sy1lo,sy1up)
      if (sy1lo.lt.sylo) sylo=sy1lo
      if (sy1up.gt.syup) syup=sy1up
      call map(sxlo,sxup,sylo,syup)
      call axorig(sxlo,sylo)
      call ctrmag(6)
      call axessi(0.,0.)
      call ctrmag(14)
      call curveo(sx,sy,1,lmid)
      call broken(10,5,10,5)
      call curveo(sx,sy1,1,lmid)
      call full
      call border
      call pspace(0.05,0.07,0.15,0.39)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.9
      call ctrori(90.)
      call pcscen(sxpos,sypos,'B')
      call suffix
      call ctrfnt(2)
      call typenc(102)
      call ctrfnt(0)
      call normal
      call typecs(' (T)')
      call ctrori(0.)
      call pspace(0.38,0.42,0.03,0.05)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call pcscen(sxpos,sypos,' R (m)')
!-------------------------------------------------------------------------
!   safety factor plot
      call pspace(0.57,0.91,0.08,0.42)
      call border
      lmid=0
      j=0
      do 25 i=1,nr
        psi=umax-u(i,nsym)
        if (psi.gt.psiv(1)) goto 25
        if (r(i).gt.r0) then
!  insert axis value
          if (j.eq.0) then
            lmid=lmid+1
            j=1
            sx(lmid+1)=sngl(r0)
            psi=0.
            sy(lmid+1)=sngl(sfac(ncon))
          end if
        end if
        lmid=lmid+1
        ig=1
        do 24 ii=1,ncon-1
          if (psiv(ii).lt.psi) goto 24
          ig=ii
 24      continue
        if (ig.eq.ncon-1) ig=ig-1
        rat=(psi-psiv(ig))/(psiv(ig+1)-psiv(ig))
        sy(lmid+1)=sngl(sfac(ig)+rat*(sfac(ig+1)-sfac(ig)))
        sx(lmid+1)=sngl(r(i))
 25   continue
      lmid=lmid+2
      sy(1)=sngl(sfac(1))
      sy(lmid)=sngl(sfac(1))
      sx(1)=sngl(rin)
      sx(lmid)=sngl(rout)
      sxlo=sx(1)
      sxup=sx(lmid)
      sxlo=sx(1)
      sxup=sx(lmid)
      call range(sy,lmid,sylo,syup)
      if (sylo.gt.0.) sylo=0.
      syup=1.1*syup
      call map(sxlo,sxup,sylo,syup)
      call axorig(sxlo,sylo)
      call ctrmag(6)
      call axessi(0.,0.)
      call ctrmag(14)
      call curveo(sx,sy,1,lmid)
      call pspace(0.51,0.535,0.15,0.41)
      call map(0.,1.,0.,1.)
      sypos=0.2
      sxpos=0.4
      call ctrori(90.)
      call plotcs(sxpos,sxpos,'q')
      call ctrori(0.)
      call pspace(0.85,0.9,0.03,0.05)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call pcscen(sxpos,sypos,' R (m)')
!------------------------------------------------------------------------
!   run details
      call pspace(0.95,1.4,0.07,0.95)
      call map(0.,1.,0.,1.)
      call border
!  toroidal and vacuum magnetic field on axis
      l=0
      do 80 i=1,nr
	if (r(i).gt.rcen) goto 80
	l=i
 80   continue
      rat=(rcen-r(l))/(r(l+1)-r(l))
      psi=umax-(u(l,nsym)+rat*(u(l+1,nsym)-u(l,nsym)))
      bt0=fprof(psi,2)/rcen
      bv0=mu0*rodi/(2.*pi*rcen)
      pcen=press(0.0d0,0)
      ne0=dense(0.0d0,0)
      ni0=densi(0.0d0,1,0)
      call ctrmag(12)
      call undlin(1)
      call pcscen(0.5,0.95,'Parameters')
      call undlin(0)
      call ctrmag(12)
      call plotcs(0.03,0.9,' Central e temp.=')
      svar=sngl(te0/1000.)
      call typenf(svar,4)
      call typecs(' keV')
      call plotcs(0.03,0.87,' Cent. main i temp.=')
      svar=sngl(ti0/1000.)
      call typenf(svar,4)
      call typecs(' keV')
      call plotcs(0.03,0.84,' Vol. av. e temp.=')
      svar=(avt/1000.)
      call typenf(svar,4)
      call typecs(' keV')
      call plotcs(0.03,0.81,' Cent.dens.(e,i)=')
      svar=sngl(ne0)
      call typene(svar,3)
      call typecs(', ')
      svar=sngl(ni0)
      call typene(svar,3)
      call typecs(' m')
      call supfix
      call typecs('-3')
      call normal
      call plotcs(0.03,0.78,' vol.av.el. density=')
      svar=sngl(avel)
      call typene(svar,3)
      call typecs(' m')
      call supfix
      call typecs('-3')
      call normal
      call plotcs(0.03,0.75,' line av.el. density=')
      svar=sngl(nebar*1.0d19)
      call typene(svar,3)
      call typecs(' m')
      call supfix
      call typecs('-3')
      call normal
      if ((ipswtch.eq.2).or.(ipswtch.eq.0)) then
        call plotcs(0.03,0.72,' Pedestal width, T=')
!!$        do i=nr,1,-1
!!$          if (ixout(i,nsym).eq.1) then
!!$          write(6,*)' dpsi/dR=',(u(i,nsym)-u(i-1,nsym))/(r(i)-r(i-1)),   &
!!$                    ' r*bth=',r2*bp(r2,0.0d0)
!!$          psi=umax-u(i,nsym)
!!$          write(6,*)' r=',r(i),' temp=',tempe(psi,0)/1000.,' psi=',psi
!!$          end if
!!$        end do
        svar=sngl(umax/(r2*bp(r2,0.0d0)*teedg))
        call typenf(svar,3)
        call typecs(', ')
        svar=sngl(teped/1000.)
        call typenf(svar,3)
        call typecs(' keV')
      end if
      call plotcs(0.03,0.69,' Curr. index=')
      if (itot.eq.0) then
	if (neo.eq.1) then
	  call typecs('neoclassical')
	else if (neo.eq.-1) then
	  call typecs('spitzer')
	else
	  call typecs('user-defined')
	end if
      else
	call typecs("ff' profile specified")
      end if
      call plotcs(0.03,0.66,' Z=')
      svar=sngl(zm)
      call typenf(svar,3)
      call typecs(' Zeff=')
      svar=sngl(zeffav)
      call typenf(svar,3)
      call plotcs(0.03,0.63,' Toroidal B (geo axis)=')
      svar=sngl(bt0)
      call typenf(svar,4)
      call typecs(' T')
! Toroidal magnetic field at magnetic axis
      psi=0.
      bt0=fprof(psi,2)/r0
      call plotcs(0.03,0.6,' Toroidal B (mag axis)=')
      svar=sngl(bt0)
      call typenf(svar,4)
      call typecs(' T')
      call plotcs(0.03,0.57,' Vacuum B (geo axis)=')
      svar=sngl(bv0)
      call typenf(svar,4)
      call typecs(' T')
      call plotcs(0.03,0.54,' Tor. b/s current=')
      svar=sngl(totbs/1000000.)
      call typenf(svar,4)
      call typecs(' MA')
      if (ifast.eq.1) then
       call typecs(" (inc alpha's)")
      else
       call typecs(" (exc alpha's)")
      end if
      call plotcs(0.03,0.51,' Toroidal current=')
      svar=sngl(cur/1000000.)
      call typenf(svar,4)
      call typecs(' MA')
!  total beta...
      call plotcs(0.03,0.48,' ')
      call ctrfnt(2)
      call typenc(98)
      call ctrfnt(0)
      call typecs('=')
      svar=sngl(beta)
      call typenf(svar,3)
      call ctrfnt(2)
      call typenc(37)
      call typecs('   ')
      call typenc(98)
      call ctrfnt(0)
      call suffix
      call typecs('N')
      call normal
      call typecs('=')
!  normalised beta
      svar=sngl(3.5*betexp/betlim)
      call typenf(svar,3)
      call ctrfnt(2)
      call typecs('   ')
      call typenc(98)
      call ctrfnt(0)
      call suffix
      call typecs('exp')
      call normal
      call typecs('=')
      svar=sngl(betexp)
      call typenf(svar,3)
      call ctrfnt(2)
      call typenc(37)
      call ctrfnt(0)
!  beta current...
      call plotcs(0.03,0.45,' ')
      call ctrfnt(2)
      call typenc(98)
      call ctrfnt(0)
      call suffix
      call typecs('I')
      call normal
      call typecs('=')
      svar=sngl(betai)
      call typenf(svar,4)
!  beta poloidal...
      call typecs('  ')
      call ctrfnt(2)
      call typenc(98)
      call ctrfnt(0)
      call suffix
      call typecs('p')
      call normal
      call typecs('=')
      svar=sngl(betap)
      call typenf(svar,4)
!  beta on axis
      paxis=press(0.0d0,0)
      baxis=fprof(0.0d0,2)/r0
      svar=sngl(200.*mu0*paxis/baxis**2)
      call ctrfnt(2)
      call typecs('   ')
      call typenc(98)
      call ctrfnt(0)
      call suffix
      call typecs('axis')
      call normal
      call typecs('=')
      call typenf(svar,2)
      call ctrfnt(2)
      call typenc(37)
      call ctrfnt(0)
! poloidal flux...
      call plotcs(0.03,0.42,' ')
      call ctrfnt(2)
      call typenc(118)
      call ctrfnt(0)
      call suffix
      call typecs('0')
      call normal
      call typecs('=')
      svar=0.
      call typene(svar,3)
      call typecs('  ')
      call ctrfnt(2)
      call typenc(118)
      call ctrfnt(0)
      call suffix
      call typecs('a')
      call normal
      call typecs('=')
      svar=sngl(umax)
      call typene(svar,3)
      call typecs(' Tm')
      call supfix
      call typecs('-2')
      call normal
      call plotcs(0.03,0.39,' q')
      call suffix
      call typecs('0')
      call normal
      call typecs('=')
      svar=sngl(sfac(ncon))
      call typenf(svar,3)
      call typecs('  ')


      if (fast.eq.1) then
         call ctrfnt(2)
         call typenc(98)
         call ctrfnt(0)
         call suffix
         call typecs('fast')
         call normal
         call typecs('=')
         svar=fastb
         call typenf(svar,3)
      end if
         

      call plotcs(0.03,0.36,' q')
      call suffix
      call typecs('a')
      call normal
!      iqw=ncon-int(0.95*ncon)+1
!      if (iqw.ne.1) then
!	rat=0.95*ncon-int(0.95*ncon)
!	spl=sfac(iqw)+rat*(sfac(iqw-1)-sfac(iqw))
!      else
!	spl=sfac(1)
!      end if
      svar=sngl(sfac(1))
      call typenf(svar,3)
      call plotcs(0.03,0.33,' R')
      call suffix
      call typecs('0')
      call normal
      call typecs('=')
      svar=sngl(r0)
      call typenf(svar,3)
      call typecs('m')
!  impurity info
      call ctrmag(16)
      call plotcs(0.03,0.3,'impurity info')
      call ctrmag(12)
      call plotcs(0.03,0.27,'no impurity species=')
      call typeni(nimp,2)
! Ratio of edge to central current density
      psi=psiv(1)
      fval=fprof(psi,2)
      pd=press(psi,1)
      fd=fprof(psi,1)/fval
      jedge=-(fval*pd/sqrt(bsqav(1))+fd*sqrt(bsqav(1))/mu0)
!      write(6,*)' pd=',pd,' fd=',fd,' fval=',fval,' jedge=',jedge
!!$      psi=psiv(ncon)
!!$      fval=fprof(psi,2)
!!$      pd=press(psi,1)
!!$      fd=fprof(psi,1)/fval
!!$      jaxis=-(fval*pd/sqrt(bsqav(ncon))+fd*sqrt(bsqav(ncon))/mu0)
      jrat=sngl(jedge*area/cur)
      call plotcs(0.03,0.24,'jedge/j-avge=')
      call typene(jrat,4)
      isk=0
      do 90 i=1,ncon
        if (isk.gt.0) goto 90
        psi=psiv(i)/umax
        if (psi.lt.0.95) isk=i
 90   continue
      psilo=psiv(isk)
      psiup=psiv(isk-1)
      rat=(0.95*umax-psilo)/(psiup-psilo)
      write(6,*)' psilo=',psilo/umax,' psiup=',psiup/umax,' rat=',rat
      psi=0.95*umax
      fval=fprof(psi,2)
      pd=press(psi,1)
      fd=fprof(psi,1)/fval
      bsq=bsqav(isk)+rat*(bsqav(isk-1)-bsqav(isk))
      q95=sfac(isk)+rat*(sfac(isk-1)-sfac(isk))
      jedge=-(fval*pd/sqrt(bsq)+fd*sqrt(bsq)/mu0)
      jrat=sngl(jedge*area/cur)
      call plotcs(0.03,0.21,'j95/j-avge=')
      call typene(jrat,4)
      isk=0
      do 91 i=nr,1,-1
        if (isk.gt.0) goto 91
        psi=1.-u(i,nsym)/umax
        if (psi.lt.0.95) isk=i
 91   continue
      psilo=umax-u(isk,nsym)
      psiup=umax-u(isk+1,nsym)
      if (psiup.gt.umax) then
        psiup=umax
        rupd=r2
        write(6,*)' r2=',r2
      else
        rupd=r(isk+1)
      end if
      psi=0.95*umax
      plod=bk*dense(psilo,0)*tempe(psilo,0)
      pupd=bk*dense(psiup,0)*tempe(psiup,0)
      flo=fprof(psilo,2)
      fup=fprof(psiup,2)
      rlod=r(isk)
      bthlo=bp(rlod,0.0d0)
      bthup=bp(rupd,0.0d0)
      betup=2.*mu0*pupd/bv0**2
      betlo=2.*mu0*plod/bv0**2
      psilo=psilo/umax
      psiup=psiup/umax
      rat=(0.95-psilo)/(psiup-psilo)
      rv=rlod+rat*(rupd-rlod)
      write(6,*)' psilo=',psilo,' psiup=',psiup,' rat=',rat
      write(6,*)' betup=',betup,' betlo=',betlo,' rupd=',rupd,' rlod=',rlod
      write(6,*)' plod=',plod,' pupd=',pupd,' bthup=',bthup,' bthlo=',bthlo
      write(6,*)' fup=',fup,' flo=',flo
      write(6,*)' rv=',rv,' q95=',q95
      write(6,*)' dp/dR=',(pupd-plod)/(rupd-rlod),   &
                2.*mu0*rupd*bthup*press(psi,1)
      write(6,*)' bv0=',bv0,' q95=',q95,' r0=',rcen
      alf95=sngl(-rcen*q95**2*(betup-betlo)/(rupd-rlod))
      call plotcs(0.03,0.18,'alpha95=')
      call typene(alf95,4)
      call plotcs(0.03,0.15,'p peaking=')
      svar=sngl(ppeak)
      call typenf(svar,4)
!!$      if (nimp.gt.0) then
!!$      call plotcs(0.03,0.27,' z')
!!$      call plotcs(0.18,0.27,' n')
!!$      call suffix
!!$      call typecs('0')
!!$      call normal
!!$      call typecs('/n')
!!$      call suffix
!!$      call typecs('i0')
!!$      call normal
!!$      call plotcs(0.38,0.27,' t')
!!$      call suffix
!!$      call typecs('0')
!!$      call normal
!!$      call typecs('(ev)')
!!$      call plotcs(0.57,0.27,' ')
!!$      call ctrfnt(2)
!!$      call typenc(97)
!!$      call suffix
!!$      call typecs('t')
!!$      call normal
!!$      call plotcs(0.74,0.27,' ')
!!$      call typenc(97)
!!$      call ctrfnt(0)
!!$      call suffix
!!$      call typecs('n')
!!$      call normal
!!$      call plotcs(0.9,0.27,'m')
!!$      ds=0.03
!!$      do 50 i=1,nimp
!!$      yy=0.27-i*ds
!!$      call plotcs(0.01,yy,' ')
!!$      call typeni(iz(i+1),2)
!!$      call typecs('  ')
!!$      call typenf(zn(i+1),4)
!!$      call typecs(' ')
!!$      call typenf(zt(i+1),1)
!!$      call typecs(' ')
!!$      call typenf(zat(i+1),2)
!!$      call typecs(' ')
!!$      call typenf(zan(i+1),2)
!!$      call typecs(' ')
!!$      call typenf(zmas(i+1),2)
!!$ 50   continue
!!$      end if
      call plotcs(0.03,0.07,' tauhe*/taue=')
      svar=sngl(tauh)
      call typenf(svar,4)
      if (nco.ge.1) then
	call plotcs(0.03,0.03,'collisional model')
      else
	call plotcs(0.03,0.03,'collisionless model')
      end if
!---------------------------------------------------------------------
! job  title
      call ctrmag(26)
      call pspace(0.03,0.95,0.83,0.95)
      call border
      call map(0.,1.,0.,1.)
      call thick(2)
      call plotcs(0.35,0.8,title)
      call thick(1)
!!$      call enqdat(cdat)
!!$      call enqtim(ctim)
!!$      call ctrmag(16)
!!$      call plotcs(0.01,0.5,'SCENE 6.0 ')
!!$      call typecs('queries: H.R. wilson, Culham ext 3902')
!!$      write(string,fmt='(''Run time:'',a8,'' on '',a9)') ctim,cdat
!!$!      call pcsend(0.98,0.6,string)
!!$      call plotcs(0.01,0.2,string)
      call pspace(0.03,1.41,0.03,0.95)
      call border
      call frame
!
!-------------------------------------------------------------------
!
! page 3 plots
!
!--------------------------------------------------------------------
!  Set normalisations for plots...
      dilo=0.
      diup=0.
      bslo=0.
      bsup=0.
      pslo=0.
      psup=0.
      tolo=0.
      toup=0.
      jj=nsym
      do i=1,nr
	if (sngl(diph(i,jj)).gt.diup) diup=sngl(diph(i,jj))
	if (sngl(diph(i,jj)).lt.dilo) dilo=sngl(diph(i,jj))
	if (sngl(psph(i,jj)).gt.psup) psup=sngl(psph(i,jj))
	if (sngl(psph(i,jj)).lt.pslo) pslo=sngl(psph(i,jj))
	if (sngl(bsph(i,jj)).gt.bsup) bsup=sngl(bsph(i,jj))
	if (sngl(bsph(i,jj)).lt.bslo) bslo=sngl(bsph(i,jj))
	if (sngl(exph(i,jj)).gt.toup) toup=sngl(exph(i,jj))
	if (sngl(exph(i,jj)).lt.tolo) tolo=sngl(exph(i,jj))
	if (sngl(exph2(i,jj)).gt.toup) toup=sngl(exph2(i,jj))
	if (sngl(exph2(i,jj)).lt.tolo) tolo=sngl(exph2(i,jj))
	if (sngl(gradj(i,jj)).gt.toup) toup=sngl(gradj(i,jj))
	if (sngl(gradj(i,jj)).lt.tolo) tolo=sngl(gradj(i,jj))
      end do
      tolo=tolo/1000.
      toup=toup/1000.
      toup=tolo+1.3*(toup-tolo)
      bslo=bslo/1000.
      bsup=bsup*1.5/1000.
      pslo=pslo/1000.
      psup=psup/1000.
      dilo=dilo/1000.
      diup=diup/1000.
      rlo=sngl(rcen*(1.-tokeps))
      rup=sngl(rcen*(1.+tokeps))
      call ctrmag(14)
!  toroidal current plot
!      call filnam('curplot.grd')
!      call thick(2)
      call pspace(0.56,0.9,0.47,0.81)
      call border
      call map(rlo,rup,tolo,toup)
      sxpos=rlo+0.4*(rup-rlo)
      sxposup=rlo+0.65*(rup-rlo)
      sypos=tolo+0.95*(toup-tolo)
      call positn(sxpos,sypos)
      call broken(10,5,2,5)
      call join(sxposup,sypos)
      call full
      call plotcs(sxposup,sypos,' External')
      sypos=tolo+0.89*(toup-tolo)
      call broken(2,7,2,7)
      call positn(sxpos,sypos)
      call join(sxposup,sypos)
      call full
      call plotcs(sxposup,sypos,' Aux J')
      sypos=tolo+0.83*(toup-tolo)
      call positn(sxpos,sypos)
      call join(sxposup,sypos)
      call plotcs(sxposup,sypos,' Total J')
      call axorig(rlo,tolo)
      call ctrmag(6)
      szero=0.
      call axessi(szero,szero)
      call ctrmag(14)
!!$      do i=1,nr
!!$	sy(i)=sngl(gradj(i,jj)/1000.)
!!$	sy1(i)=sngl(exph(i,jj)/1000.)
!!$	sx(i)=sngl(r(i))
!!$      end do
      write(73,*)' R,          Jphi(MAm**-2)'
      do  i=1,nr
        sx1(i)=sngl(r(i))
        sy1(i)=sngl(gradj(i,jj)/1000.)
        sy2(i)=sngl(exph(i,jj)/1000.)
        sy3(i)=sngl(exph2(i,jj)/1000.)
        if (ixout(i,jj).le.0) sy1(i)=0.
        if (ixout(i,jj).le.0) sy2(i)=0.
        if (ixout(i,jj).le.0) sy3(i)=0.
      end do
!  add in edge values
      lp=0
      pdedg=press(psiv(1),1)
      ffdedg=fprof(psiv(1),1)
      do i=1,nr
        sx(i+lp)=sx1(i)
        sy(i+lp)=sy1(i)
        rr=dble(sx1(i))
        if (rr.gt.rin) then
          if (lp.eq.0) then
            lp=1
            sx(i+lp)=sx1(i)
            sy(i+lp)=sy1(i)
            sx(i)=sngl(rin)
	    sy(i)=(-rin*pdedg-ffdedg/(rin*mu0))/1000.
          end if
        end if
        if (rr.gt.rout) then
          if (lp.eq.1) then
            lp=2
            sx(i+lp)=sx1(i)
            sy(i+lp)=sy1(i)
            sx(i+1)=sngl(rout)
            sy(i+1)=(-rout*pdedg-ffdedg/(rout*mu0))/1000.
          end if
        end if
      end do
      lp=0
      icur=1
      if (neo.lt.0) icur=-1
      call extj(epsv(1),psiv(1),extapp,extapp2,icur)
      fsi=fprof(psiv(1),2)
! total driven current
      do i=1,nr
        sx(i+lp)=sx1(i)
        sy1(i+lp)=sy2(i)
        rr=dble(sx1(i))
        if (rr.gt.rin) then
          if (lp.eq.0) then
            lp=1
            sx(i+lp)=sx1(i)
            sy1(i+lp)=sy2(i)
            sx(i)=sngl(rin)
	    sy1(i)=sngl((vloop*extapp+extapp2)*fsi/(1000.*rin))
            write(6,*)' R=',sx(i),' exJ=',sy1(i)
            write(6,*)' extapp=',extapp,' fsi=',fsi,' extapp2=',extapp2
          end if
        end if
        if (rr.gt.rout) then
          if (lp.eq.1) then
            lp=2
            sx(i+lp)=sx1(i)
            sy1(i+lp)=sy2(i)
            sx(i+1)=sngl(rout)
	    sy1(i+1)=sngl((vloop*extapp+extapp2)*fsi/(1000.*rout))
          end if
        end if
      end do
! auxilliary current drive
      lp=0
      do i=1,nr
        sx(i+lp)=sx1(i)
        sy2(i+lp)=sy3(i)
        rr=dble(sx1(i))
        if (rr.gt.rin) then
          if (lp.eq.0) then
            lp=1
            sx(i+lp)=sx1(i)
            sy2(i+lp)=sy3(i)
            sx(i)=sngl(rin)
	    sy2(i)=sngl(extapp2*fsi/(1000.*rin))
          end if
        end if
        if (rr.gt.rout) then
          if (lp.eq.1) then
            lp=2
            sx(i+lp)=sx1(i)
            sy2(i+lp)=sy3(i)
            sx(i+1)=sngl(rout)
	    sy2(i+1)=sngl(extapp2*fsi/(1000.*rout))
          end if
        end if
      end do
      call xarea(gradj,curgs)
      call lincol(2)
      call broken(2,7,2,7)
!  Auxiliary
      call curveo(sx,sy2,1,nr+2)
      call lincol(3)
      call broken(10,5,2,5)
!  External
      write(nh,*)' No. points in Driven J-profile=',nr+2
      write(nh,13)jcore/1.0d6,(jtotal-jcore)/1.0d6
      write(nh,*)' R(m), J_aux_phi, J_tot_phi (kA m**-2)'
      do i=1,nr+2
        write(nh,12)sx(i),sy1(i),sy(i)
      end do
 12   format(3e14.6)
 13   format('Core current driven in MA (20%psi)=',e14.6,' edge CD=',e14.6)
      call curveo(sx,sy1,1,nr+2)
      call lincol(0)
      call full
!  Total
      do i=1,nr+2
        write(73,*)sx(i),sy(i)*1.0e-3
      end do
      close(73)
      call curveo(sx,sy,1,nr+2)
!      open(24,file='JCTF-NBI-ECCD.nbi')
!      read(24,56)nrob
! 56   format(' No. points in Driven J-profile=       ',i6)
!      write(6,*)' nrob=',nrob
!      read(24,*)txt
!      read(24,*)txt
!      do i=1,nrob
!        read(24,*)sxx(i),syy(i)
!        write(6,*)' sx=',sxx(i),' sy=',syy(i)
!      end do
!      close(24)
!      call lincol(4)
!      call curveo(sxx,syy,1,nrob)
      call lincol(0)
      call pspace(0.5,0.525,0.55,0.81)
      call map(0.,1.,0.,1.)
      sypos=0.2
      sxpos=0.4
      call ctrori(90.)
      call plotcs(sxpos,sypos,'J')
      call suffix
      call ctrfnt(2)
      call typenc(102)
      call ctrfnt(0)
      call normal
      call typecs(' (kA m')
      call supfix
      call typecs('-2')
      call normal
      call typecs(')')
      call ctrori(0.)
      call pspace(0.85,0.9,0.425,0.44)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call pcscen(sxpos,sypos,' R(m)') 
!      call frame
!      call filnam('graphs.grd')
!  plot of bootstrap current
      if (neo.ge.0) then
	call pspace(0.1,0.44,0.08,0.41)
	call map(rlo,rup,bslo,bsup)
	call axorig(rlo,bslo)
	call ctrmag(6)
	call axessi(0.,0.)
	call ctrmag(14)
	lp=0
	do  i=1,nr
	 if (r(i).lt.r0) then
           sx1(i)=sngl(r(i))
	   sy1(i)=sngl(bsph(i,jj)/1000.)
           if (ixout(i,jj).le.0) sy(i)=0.
	 else
	   if (lp.eq.1) then
	     lp=lp+1
	     sy1(i)=0.
	     sx1(i)=r0
	   else
	     sy1(i)=sngl(bsph(i-1,jj)/1000.)
	     sx1(i)=sngl(r(i-1))
             if (ixout(i-1,jj).le.0) sy1(i)=0.
	   end if
	 end if
        end do
	sx1(nr+1)=sngl(r(nr))
	sy1(nr+1)=0.
!  add in edge values
        lp=0
        do i=1,nr+1
          sx(i+lp)=sx1(i)
          sy(i+lp)=sy1(i)
          rr=dble(sx1(i))
          if (rr.gt.rin) then
            if (lp.eq.0) then
              lp=1
              sx(i+lp)=sx1(i)
              sy(i+lp)=sy1(i)
              sx(i)=sngl(rin)
              sy(i)=sngl(bsj(1)*fsi/(1000.*rin*sqrt(bsqav(1))))
            end if
          end if
          if (rr.gt.rout) then
            if (lp.eq.1) then
              lp=2
              sx(i+lp)=sx1(i)
              sy(i+lp)=sy1(i)
              sx(i+1)=sngl(rout)
              sy(i+1)=sngl(bsj(1)*fsi/(1000.*rout*sqrt(bsqav(1))))
            end if
          end if
        end do
	call curveo(sx,sy,1,nr+3)
	call border
	call pspace(0.05,0.07,0.15,0.39)
	call map(0.,1.,0.,1.)
	sypos=0.5
	sxpos=0.5
	call ctrori(90.)
	call pcscen(sxpos,sypos,'J')
	call suffix
	call typecs('b')
	call normal
	call typecs(' (kA m')
	call supfix
	call typecs('-2')
	call normal
	call typecs(')')
	call ctrori(0.)
	call pspace(0.38,0.42,0.04,0.05)
	call map(0.,1.,0.,1.)
	sypos=0.5
	sxpos=0.5
	call pcscen(sxpos,sypos,' R (m)')
      end if
!  Pfirsch--Schluter current plots
      call pspace(0.1,0.44,0.47,0.81)
      call map(rlo,rup,pslo,psup)
      call axorig(rlo,pslo)
      call ctrmag(6)
      call axessi(0.,0.)
      call ctrmag(14)
      do i=1,nr
       sy1(i)=sngl(psph(i,jj)/1000.)
       sx1(i)=sngl(r(i))
       if (ixout(i,jj).le.0) sy(i)=0.
      end do
!  add in edge values
      lp=0
      do i=1,nr
        sx(i+lp)=sx1(i)
        sy(i+lp)=sy1(i)
        rr=dble(sx1(i))
        if (rr.gt.rin) then
          if (lp.eq.0) then
            lp=1
            sx(i+lp)=sx1(i)
            sy(i+lp)=sy1(i)
            sx(i)=sngl(rin)
            bphi=fsi/rin
            bth=bp(rin,0.0d0)
            bsq=bphi*bphi+bth*bth
	    sy(i)=-fsi*fsi*pdedg*(1.-bsq/bsqav(1))/(1000.*rin*bsq)
          end if
        end if
        if (rr.gt.rout) then
          if (lp.eq.1) then
            lp=2
            sx(i+lp)=sx1(i)
            sy(i+lp)=sy1(i)
            sx(i+1)=sngl(rout)
            bphi=fsi/rout
            bth=bp(rout,0.0d0)
            bsq=bphi*bphi+bth*bth
	    sy(i+1)=-fsi*fsi*pdedg*(1.-bsq/bsqav(1))/(1000.*rout*bsq)
          end if
        end if
      end do
      call curveo(sx,sy,1,nr+2)
      call border
      call pspace(0.05,0.07,0.55,0.81)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call ctrori(90.)
      call pcscen(sxpos,sypos,'J')
      call suffix
      call typecs('ps')
      call normal
      call typecs(' (kA m')
      call supfix
      call typecs('-2')
      call normal
      call typecs(')')
      call ctrori(0.)
      call pspace(0.38,0.42,0.425,0.44)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call pcscen(sxpos,sypos,' R (m)')
!  diamagnetic current plots
      call pspace(0.56,0.9,0.08,0.41)
      call map(rlo,rup,dilo,diup)
      call axorig(rlo,dilo)
      call ctrmag(6)
      call axessi(0.,0.)
      call ctrmag(14)
      do i=1,nr
       sy1(i)=sngl(diph(i,jj)/1000.)
       sx1(i)=sngl(r(i))
       if (ixout(i,jj).le.0) sy(i)=0.
      end do
!  add in edge values
      lp=0
      do i=1,nr
        sx(i+lp)=sx1(i)
        sy(i+lp)=sy1(i)
        rr=dble(sx1(i))
        if (rr.gt.rin) then
          if (lp.eq.0) then
            lp=1
            sx(i+lp)=sx1(i)
            sy(i+lp)=sy1(i)
            sx(i)=sngl(rin)
            bphi=fsi/rin
            bth=bp(rin,0.0d0)
            bsq=bphi*bphi+bth*bth
	    sy(i)=-rin*bth*bth*pdedg/(1000.*bsq)
          end if
        end if
        if (rr.gt.rout) then
          if (lp.eq.1) then
            lp=2
            sx(i+lp)=sx1(i)
            sy(i+lp)=sy1(i)
            sx(i+1)=sngl(rout)
            bphi=fsi/rout
            bth=bp(rout,0.0d0)
            bsq=bphi*bphi+bth*bth
	    sy(i+1)=-rout*bth*bth*pdedg/(1000.*bsq)
          end if
        end if
      end do
      call curveo(sx,sy,1,nr+2)
      call border
      call pspace(0.5,0.52,0.15,0.39)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call ctrori(90.)
      call pcscen(sxpos,sypos,'J')
      call suffix
      call typecs('dia')
      call normal
      call typecs(' (kA m')
      call supfix
      call typecs('-2')
      call normal
      call typecs(')')
      call ctrori(0.)
      call pspace(0.8,0.85,0.03,0.05)
      call map(0.,1.,0.,1.)
      sypos=0.5
      sxpos=0.5
      call pcscen(sxpos,sypos,' R (m)')
!   Run details
      call pspace(0.95,1.4,0.07,0.87)
      call map(0.,1.,0.,1.)
      call border
      call ctrmag(23)
      call undlin(1)
      call pcscen(0.5,0.95,'Parameters')
      call undlin(0)
      call ctrmag(14)
      svar=sngl(curgs/1000000.)
      call plotcs(0.03,0.9,' Toroidal current=')
      call typenf(svar,4)
      call typecs(' MA')
      call ctrmag(14)
      call plotcs(0.03,0.86,' Tor. b/s current=')
      svar=sngl(totbs/1000000.)
      call typenf(svar,4)
      call typecs(' MA')
      if (ifast.eq.1) then
       call typecs(" (inc alpha's)")
      else
       call typecs(" (exc alpha's)")
      end if
      call plotcs(0.03,0.82,' Tor. p/s current=')
      svar=sngl(totps/1000000.)
      call typenf(svar,4)
      call typecs(' MA')
      call plotcs(0.03,0.78,' Tor. dia current=')
      svar=sngl(totdi/1000000.)
      call typenf(svar,4)
      call typecs(' MA')
      call plotcs(0.03,0.74,' Tor. Ohmic current=')
      if (neo.eq.1) then
        svar=sngl((totex-totex2)/1000000.)
      else
        svar=0.
      end if
      call typenf(svar,4)
      call typecs(' MA')
      call plotcs(0.03,0.7,' Tor. aux. current=')
      if (neo.eq.1) then
        svar=sngl(totex2/1000000.)
      else
        svar=sngl(totex/1000000.)
      end if
      call typenf(svar,4)
      call typecs(' MA')
      call ctrmag(14)
      call plotcs(0.03,0.65,' alpha,aux pow=')
      svar=sngl(pfus*1.0d-6)
      call typenf(svar,2)
      call typecs(',')
      svar=sngl(paux)
      call typenf(svar,2)
      call typecs(' MW')
      call plotcs(0.03,0.6,' l')
      call suffix
      call typecs('i')
      call normal
      call typecs('(3),(2)')
      call typecs('=')
      svar=sngl(rli3)
      call typenf(svar,4)
      call typecs(', ')
      svar=sngl(rli2)
      call typenf(svar,4)
      call plotcs(0.03,0.55,' area, vol=')
      svar=sngl(area)
      call typenf(svar,4)
      call typecs('m')
      call supfix
      call typecs('2')
      call normal
      svar=sngl(vol)
      call typenf(svar,4)
      call typecs('m')
      call supfix
      call typecs('3')
      call normal
!  Plasma surface area
      bigs=0.
      do i=1,npts-1
        dl=sqrt((rpts(1,i)-rpts(1,i+1))**2+(zpts(1,i)-zpts(1,i+1))**2)
        bigs=bigs+pi*(rpts(1,i)+rpts(1,i+1))*dl
      end do
      call plotcs(0.03,0.5,' 1.5')
      call ctrfnt(2)
      call typenc(194)
      call ctrfnt(0)
      call typecs('p dV =')
      svar=sngl(conft/1.0e6)
      call typenf(svar,4)
      call typecs(' MJ')
      if (abs(neo).eq.1) then
        if (ibv.eq.0) then
          call plotcs(0.03,0.45,' Loop volts=')
        else
          call plotcs(0.03,0.45,' Bv-dot=')
        end if
        svar=sngl(vloop)
	call typenf(svar,4)
        if (ibv.eq.0) then
	  call typecs(' V')
        else
	  call typecs(' Ts')
          call supfix
          call typecs('-1')
          call normal
        end if
!	call plotcs(0.03,0.4,' loop v (j')
!	call suffix
!	call typecs('bs')
!	call normal
!	call typecs('=0) = ')
!	call typenf(vnobs,4)
!	call typecs(' v')
	if (neo.eq.1) then
	  call plotcs(0.03,0.4,' Loop v (Spitz.)=')
          svar=sngl(vspit)
	  call typenf(svar,4)
	  call typecs(' V')
	end if
      end if
      call plotcs(0.03,0.35,' H')
      call suffix
      call typecs('IPB98(y,1)')
      call normal
      call typecs('=')
      svar=sngl(hipb98y1)
      call typenf(svar,3)
      call typecs(', H')
      call suffix
      call typecs('IPB98(y,2)')
      call normal
      svar=sngl(hipb98y2)
      call typenf(svar,3)
      call plotcs(0.03,0.31,' n/n')
      call suffix
      call typecs('GW')
      call normal
      call typecs('=')
      svar=sngl(negw)
      call typenf(svar,3)
      call ctrfnt(2)
      call typecs(', t')
      call ctrfnt(0)
      call suffix
      call typecs('E')
      call normal
      call typecs('=')
      svar=sngl(taue*1.0d3)
      call typenf(svar,3)
      call typecs('ms')
! centre column power density
      fcop=0.6
      rc=rcen*(1.-tokeps)
      rescop=2.5d-8
      pden=rescop*elon*tokeps*rcen/(fcop*pi*rc**2)**2
      pden=rodi*rodi*pden*1.0e-6
      call plotcs(0.03,0.28,'C C power density=')
      svar=sngl(pden)
      call typenf(svar,3)
      call typecs('MWm')
      call supfix
      call typecs('2')
      call normal
      zstar=0
      do i=1,npts/2
        if (zstar.lt.zpts(1,i)) then
          zstar=zpts(1,i)
          rstar=rpts(1,i)
        end if
      end do
      elonc=zstar/(tokeps*rcen)
      tric=(1.-rstar/rcen)/tokeps
! Evaluate quadracity
      rhalf=-1.
      do i=npts/2+1,1,-1
        if (rhalf.lt.0.) then
          if (zpts(1,i).gt.zstar/2.) then
            rat=(zstar/2.-zpts(1,i))/(zpts(1,i+1)-zpts(1,i))
            rhalf=rpts(1,i)+rat*(rpts(1,i+1)-rpts(1,i))
          end if
        end if
      end do
      rhalf=rhalf-rcen
      if ((rhalf.lt.0.).or.(rhalf.gt.(rcen*tokeps))) then
        write(6,*)' Error in calculating quadc in graphs'
        stop
      end if
      arg=rhalf/(tokeps*rcen)
      quadc=(2./sqrt(3.))*(acos(arg)-pi/6.-tri/2.)
      call plotcs(0.03,0.25,' Elong=')
      svar=sngl(elonc)
      call typenf(svar,3)
      call typecs(' Bishop k=')
      svar=sngl(kval)
      call typenf(svar,3)
      if (abs(tric).gt.1.) then
         write(6,*)'WARNING****Could not calculate triangularity'
         call plotcs(0.03,0.22,' Tri=*****')
      else
        call plotcs(0.03,0.22,' Tri=')
        svar=sngl(asin(tric))
      end if
      call typenf(svar,3)
      call typecs(' quad=')
      svar=quadc
      call typenf(svar,4)
      call plotcs(0.03,0.19,' Aspect ratio=')
      svar=sngl(1./tokeps)
      call typenf(svar,3)
      call plotcs(0.03,0.16,' Geometric axis=')
      svar=sngl(rcen)
      call typenf(svar,3)
      call typecs('m')
      call plotcs(0.03,0.13,' Rod I=')
      svar=sngl(rodi/1000000.)
      call typenf(svar,3)
      call typecs('MA')
      call plotcs(0.03,0.1,' alpha b/strap=')
      svar=sngl(alfbs/1000000.)
      call typenf(svar,4)
      call typecs(' MA')
      call plotcs(0.03,0.07,' alpha press. fraction=')
      svar=sngl(palpha)
      call typenf(svar,4)
      call plotcs(0.03,0.04,' plasma surface area=')
      svar=sngl(bigs)
      call typenf(svar,4)
      call typecs('m')
      call supfix
      call typecs('2')
      call normal
      call plotcs(0.06,0.01,'Tritium/year =')
      svar=sngl(pfus*0.28238964*1.0d-6)
      call typenf(svar,2)
      call typecs('kg/year')
      call ctrmag(26)
      call pspace(0.08,0.9,0.83,0.87)
      call border
      call map(0.,1.,0.,1.)
      call thick(2)
      call plotcs(0.35,0.25,title)
      call thick(1)
      call pspace(0.03,1.41,0.03,0.89)
      call border
      call frame
!      call filnam('flxplot.grd')
!
!------------------------------------------------------------------
!
      if (igr.eq.5) then
!  read in temperature and density profiles from data file
!  and plot versus radius
      call filnam('graphs.grd')
 20   open(20,file='mast_8694.dat')
      read(20,*)ishot,timsht,ntdat
      write(6,*)' ishot=',ishot,' timsht=',timsht,ntdat
!!$      if (ntdat.gt.nr) then
!!$        write(6,*)'WARNING***only plotting 1st ',nr,' data points'
!!$        write(6,*)' Increase dimensions of pltr,pltt and pltn for more'
!!$        ntdat=nr
!!$      end if
!!$      do i=1,ntdat
!!$        read(20,*)pltr(i),pltn(i),pltt(i)
!!$        pltn(i)=pltn(i)*1.0d19
!!$        pltt(i)=pltt(i)*1000.
!!$        pltp(i)=(1.+ti0/te0)*bk*pltt(i)*pltn(i)
!!$!        write(6,*)' i=',i,' R=',pltr(i),' T=',pltt(i),' ne=',pltn(i)
!!$      end do
      do i=1,ntdat
        read(20,*)pltr(i),pltn(i),pltt(i)
        write(6,*)' R=',pltr(i),' n=',pltn(i),' T=',pltt(i)
        pltn(i)=pltn(i)*1.0d19
        pltt(i)=pltt(i)*1000.
      end do
      close(20)
!!$    open(20,file='mast_4985Te.dat')
!!$      do i=1,ntdat
!!$        read(20,*)pltr(i),pltt(i)
!!$        write(6,*)' R=',pltr(i),' T=',pltt(i)
!!$        pltt(i)=pltt(i)*1000.
!!$        pltp(i)=(1.+ti0/te0)*bk*pltt(i)*pltn(i)
!!$      end do
      do  i=1,ntdat-1
        pltter(i)=0.
        pltner(i)=0.
        if (i.eq.1) then
          tmax=pltt(i)+pltter(i)
          rnmax=pltn(i)+pltner(i)
        else
          if (pltr(i).lt.(tokeps*rcen+rcen)) then
            te=pltt(i)+pltter(i)
            if (te.gt.tmax) tmax=te
            ne=pltn(i)+pltner(i)
            if (ne.gt.rnmax) rnmax=ne
          end if
        end if
      end do
      rmax=r2*1.1
      rnmax=1.1*rnmax
      tmax=1.1*tmax
      pmax=2.*bk*tmax*rnmax
      rmin=0.
      lmid=0
      j=0
      tgot=-1.
      ngot=-1.
      rpos=rcen*(1.+tokeps)-0.028
      do 40 i=1,nr
        if (ixout(i,nsym).le.0) goto 40
        if (r(i).gt.r0) then
!  insert axis value
          if (j.eq.0) then
            lmid=lmid+1
            j=1
            sx(lmid+1)=sngl(r0)
            psi=0.
            sy(lmid+1)=sngl(dense(psi,0))
            sy1(lmid+1)=sngl(tempe(psi,0))
            sy2(lmid+1)=sngl(press(psi,0))
          end if
        end if
        lmid=lmid+1
        psi=umax-u(i,nsym)
        sx(lmid+1)=sngl(r(i))
        sy(lmid+1)=sngl(dense(psi,0))
        sy1(lmid+1)=sngl(tempe(psi,0))
        sy2(lmid+1)=sngl(press(psi,0))
        if (r(i).gt.rpos) then
          if (tgot.lt.0.) then
            rat=(rpos-r(i-1))/(r(i)-r(i-1))
            psilst=umax-u(i-1,nsym)
            tgot=tempe(psilst,0)+rat*(tempe(psi,0)-tempe(psilst,0))
            ngot=dense(psilst,0)+rat*(dense(psi,0)-dense(psilst,0))
            psilst=psilst+rat*(psi-psilst)
          end if
        end if
 40   continue
      do 41 i=1,ncon
        if (psilst.gt.psiv(i)) goto 41
        ii=i
 41   continue
      if (ii.ge.ncon) ii=ncon-1
      rat=(psilst-psiv(ii))/(psiv(ii+1)-psiv(ii))
      qgot=sfac(ii)+rat*(sfac(ii+1)-sfac(ii))
      write(6,*)'-------------------------------------'
      write(6,*)' rpos=',rpos,' psi=',psilst/umax
      write(6,*)' Te=',tgot,' ne19=',ngot*1.0d-19,' q=',qgot
      write(6,*)'-------------------------------------'
      lmid=lmid+2
      sx(1)=sngl(rin)
      sx(lmid)=sngl(rout)
      psi=psiv(1)
      sy(1)=sngl(dense(psi,0))
      sy1(1)=sngl(tempe(psi,0))
      sy2(1)=sngl(press(psi,0))
      sy(lmid)=sngl(dense(psi,0))
      sy1(lmid)=sngl(tempe(psi,0))
      sy2(lmid)=sngl(press(psi,0))
      call ctrmag(15)
      call pspace(0.1,0.5,0.2,0.7)
      call border
      call map(rmin,rmax,0.,tmax)
      call axorig(0.,0.)
      call ctrmag(12)
      call axessi(0.,0.)
      call errbar(pltr,pltt,pltter,pltter,1,ntdat,235,1)
      call lincol(3)
      call curveo(sx,sy1,1,lmid)
      call lincol(0)
      call pspace(0.1,0.5,0.7,0.8)
      call border
      call map(0.,1.,0.,1.)
      call pcscen(0.5,0.5,'TEMPERATURE PROFILE')
      call pspace(0.55,0.95,0.2,0.7)
      call border
      call map(rmin,rmax,0.,rnmax)
      call axorig(0.,0.)
      call axessi(0.,0.)
      call errbar(pltr,pltn,pltner,pltner,1,ntdat,235,1)
      call lincol(3)
      call curveo(sx,sy,1,lmid)
      call lincol(0)
      call pspace(0.55,0.95,0.7,0.8)
      call border
      call map(0.,1.,0.,1.)
      call pcscen(0.5,0.5,'DENSITY PROFILE')
!!$      call pspace(1.0,1.4,0.2,0.7)
!!$      call border
!!$      call map(rmin,rmax,0.,pmax)
!!$      call axorig(0.,0.)
!!$      call axessi(0.,0.)
!!$      call errbar(pltr,pltp,pltner,pltner,1,ntdat,235,1)
!!$      call lincol(3)
!!$      call curveo(sx,sy2,1,lmid)
!!$      call lincol(0)
!!$      call pspace(1.0,1.4,0.7,0.8)
!!$      call border
!!$      call map(0.,1.,0.,1.)
!!$      call pcscen(0.5,0.5,'PRESSURE PROFILE')
      call frame
      close(20)
!!$      do i=1,ncon
!!$        sx(i)=psiv(i)/umax
!!$        sy(i)=-press(psiv(i),1)
!!$        if (i.eq.1) then
!!$          pmax=sy(1)
!!$        else
!!$          if (pmax.lt.sy(i)) pmax=sy(i)
!!$        end if
!!$      end do
!!$      call pspace(0.1,0.5,0.2,0.7)
!!$      call border
!!$      call map(0.,1.,0.,pmax)
!!$      call axorig(0.,0.)
!!$      call ctrmag(12)
!!$      call axessi(0.,0.)
!!$      call lincol(3)
!!$      call curveo(sx,sy,1,ncon)
!!$      call lincol(0)
!!$      call pspace(0.1,0.5,0.7,0.8)
!!$      call border
!!$      call map(0.,1.,0.,1.)
!!$      call pcscen(0.5,0.5,'P-prime')
      do i=1,ncon
        sx(i)=psiv(i)/umax
        sy(i)=dense(psiv(i),0)
        sy2(i)=ne_m(i)
        if (i.eq.1) then
          pmax=sy(1)
        else
          if (pmax.lt.sy(i)) pmax=sy(i)
        end if
      end do
!!$      call pspace(0.1,0.5,0.2,0.7)
!!$      call border
!!$      call map(0.,1.,0.,pmax)
!!$      call axorig(0.,0.)
!!$      call ctrmag(12)
!!$      call axessi(0.,0.)
!!$      call errbar(sx,sy2,pltter,pltter,1,ncon,235,1)
!!$      call lincol(3)
!!$      call curveo(sx,sy,1,ncon)
!!$      call lincol(0)
!!$      call pspace(0.1,0.5,0.7,0.8)
!!$      call border
!!$      call map(0.,1.,0.,1.)
!!$      call pcscen(0.5,0.5,'n(psi)')
!!$      do i=1,ncon
!!$        sx(i)=psiv(i)/umax
!!$        sy(i)=bk*dense(psiv(i),0)*tempe(psiv(i),0)
!!$        sy2(i)=p_m(i)
!!$        if (i.eq.1) then
!!$          pmax=sy(1)
!!$        else
!!$          if (pmax.lt.sy(i)) pmax=sy(i)
!!$        end if
!!$      end do
!!$      call pspace(0.55,0.95,0.2,0.7)
!!$      call border
!!$      call map(0.,1.,0.,pmax)
!!$      call axorig(0.,0.)
!!$      call ctrmag(12)
!!$      call axessi(0.,0.)
!!$      call errbar(sx,sy2,pltter,pltter,1,ncon,235,1)
!!$      call lincol(3)
!!$      call curveo(sx,sy,1,ncon)
!!$      call lincol(0)
!!$      call pspace(0.55,0.95,0.7,0.8)
!!$      call border
!!$      call map(0.,1.,0.,1.)
!!$      call pcscen(0.5,0.5,'P(psi)')
!!$      do i=1,ncon
!!$        sx(i)=sngl(psiv(i)/umax)
!!$        sy(i)=sngl(tempe(psiv(i),0))
!!$        sy2(i)=sngl(te_m(i))
!!$        if (i.eq.1) then
!!$          pmax=sy(1)
!!$        else
!!$          if (pmax.lt.sy(i)) pmax=sy(i)
!!$        end if
!!$      end do
!!$      call pspace(1.05,1.4,0.2,0.7)
!!$      call border
!!$      call map(0.,1.,0.,pmax)
!!$      call axorig(0.,0.)
!!$      call ctrmag(12)
!!$      call axessi(0.,0.)
!!$      call errbar(sx,sy2,pltter,pltter,1,ncon,235,1)
!!$      call lincol(3)
!!$      call curveo(sx,sy,1,ncon)
!!$      call lincol(0)
!!$      call pspace(1.05,1.4,0.7,0.8)
!!$      call border
!!$      call map(0.,1.,0.,1.)
!!$      call pcscen(0.5,0.5,'T(psi)')
!!$      call frame
      end if
      j=0
      do i=1,nr
        if (ixout(i,nsym).gt.0) then
        j=j+1
        sx(j)=sngl(r(i))
        psi=umax-u(i,nsym)
        fsi=fprof(psi,2)
        bphi=fsi/r(i)
        bth=bp(r(i),0.0d0)
        bsq=bphi**2+bth**2
        sy(j)=sngl(1.76d+11*sqrt(bsq))
        sy2(j)=2.*sy(j)
        if (j.eq.1) then
          pmax=sy(1)
        else
          if (pmax.lt.sy(j)) pmax=sy(j)
        end if
        sy1(j)=sngl(56.4*sqrt(dense(psi,0)))
          if (pmax.lt.sy1(j)) pmax=sy1(j)
        end if
      end do
      pmax=4.*pmax
      call pspace(0.2,0.8,0.2,0.8)
      call border
      call map(0.,sx(j),0.,pmax)
      call axorig(0.,0.)
      call ctrmag(12)
      call axessi(0.,0.)
      call broken(10,5,10,5)
      call lincol(3)
      call curveo(sx,sy1,1,j)
      call lincol(0)
      call full
      call curveo(sx,sy,1,j)
      call curveo(sx,sy2,1,j) 
      do i=1,j
        sy2(i)=3.*sy(i)
        sy(i)=4.*sy(i)
      end do
      call curveo(sx,sy,1,j)
      call curveo(sx,sy2,1,j) 
      call pspace(0.2,0.8,0.8,0.9)
      call border
      call map(0.,1.,0.,1.)
      call pcscen(0.5,0.5,'Wce')
      call frame
!
!
!---------------------------------------------------------
!








      end subroutine graphs










