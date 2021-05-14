      subroutine output
!     *****************
!
!  routine to calculate any outputs required
!
      use param
      implicit none
      double precision dense,densi
      double precision pdiam,dlo,dup,psi,dhlo,dhup,dimp
      double precision arad,circar,cirkap,htpow,plasi,bvaxis,risomas
      double precision gwden,aa,bvacu,nhebar
      integer i
!
      !  calculate beta values...
      call betas
!  total auxiliary heating power (in MW):
      if (neo.eq.1) then
        paux=paux+1.0e-6*(totex-totex2)*vloop
      end if
      nebar=0.
      nhebar=0.
      pdiam=0.
!  line averaged density
      do 10 i=2,nr
        if (ixout(i,nsym).le.0) goto 10
        psi=umax-u(i,nsym)
        if (ixout(i-1,nsym).eq.1) then
          dlo=dup
          dhlo=dhup
        else
          dlo=0.
          dhlo=0.
        end if
        dup=dense(psi,0)
        nebar=nebar+0.5*(dlo+dup)*dr
        if (nimp .ge. 1) then
           dhup=densi(psi,2,0)
           nhebar=nhebar+0.5*(dhlo+dhup)*dr
        end if
        if (nimp .ge. 2) then
           dimp=densi(psi,3,0)
        end if        
!        write(6,*)' R=',r(i),' nz/ne=',dimp/dup        
        pdiam=pdiam+dr
 10   continue
      !write(6,*)' nHE/ne=',nhebar/nebar
      nebar=1.0d-19*nebar/pdiam
!  Confinement time relative to scaling laws
      arad=tokeps*rcen
      circar=pi*arad**2
      cirkap=area/circar
!  Heating power
      htpow=pfus*1.0e-6+paux
      plasi=cur*1.0e-6
      bvaxis=mu0*rodi/(2.*pi*rcen)
      risomas=2.5
!  H factors for y1 and y2 scaling laws
      !write(6,*)' plasi=',plasi,' bvaxis=',bvaxis,' nebar=',nebar
      !write(6,*)' htpow=',htpow,' rcen=',rcen,' cirkap=',cirkap
      !write(6,*)' tokeps=',tokeps,' risisomas=',risomas
      hipb98y1=1./(0.0503*(plasi**0.91)*(bvaxis**0.15)*(nebar**0.44)*  &
      (htpow**(-0.65))*(rcen**2.05)*(cirkap**0.72)*(tokeps**0.57)      &
      *(risomas**0.13)*htpow*1.0d6/conft)
      !write(6,*)' tauipb=',conft/(htpow*1.0d6)/hipb98y1
      hipb98y2=1./(0.0562*(plasi**0.93)*(bvaxis**0.15)*(nebar**0.41)*  &
      (htpow**(-0.69))*(rcen**1.97)*(cirkap**0.78)*(tokeps**0.58)      &
      *(risomas**0.19)*htpow*1.0d6/conft)    
!  Energy and helium confinement times
      taue=conft/(htpow*1.0d6)
      !write(6,*)' taue=',taue,' hipb98y1=',hipb98y1
      tauh=3.5*1.0d6*1.602d-19*confa/(pfus*taue)

      ! Petty 2008 scaling law
      petty = taue/(0.052*(plasi**0.75)*(bvaxis**0.3)*(nebar**0.32)* &
           (htpow**(-0.47))*(rcen**2.09)*(cirkap**0.88)* &
           (tokeps**0.84)*(risomas**0.0))
! Density normalised to Greenwald
      gwden=plasi/(pi*arad*arad)
      negw=nebar/(10.*gwden)
!  calculate beta limit
      aa=tokeps*rcen
      bvacu=mu0*rodi/(2.*pi*rcen)
      betlim=3.5*cur*1.e-6/(aa*bvacu)
  end subroutine output
!
!**********************************************************************
!
      subroutine betas
!     ****************
!
!   Calculates beta value as a percentage.
!   Also calculates different definitions of poloidal beta.
!
      use param
      implicit none
!
      double precision bptot,ptota,ptotv,btotv
      double precision arg,pden,zni, efus, xsec
      double precision rr,zz,uu,p,psi,dalf
      double precision press,densi,dense,tempe,tempi,fprof,bp
      double precision te,ti,ne,dion,rat
      double precision bth,bphi,btot,bvacu,psicut
      double precision pimp3
      double precision, dimension(:), allocatable:: pimp3v
      integer i,j,ip,l,k
      logical :: debug

      debug = .false.
!
!
!  Calculate volume/area integrals of pressure and B-field
      bptot=0.
      ptota=0.
      ptotv=0.
      btotv=0.
      vol=0.
      avel=0.
      avio=0.
      avt=0.
      avti=0.
      area=0.
      pfus=0.
      confa=0.
      if (imp.eq.1) zeffav=0.
      ip=0
      jcore=0.
      jtotal=0.
      psicut=0.2*umax
      palpha=0.
      pimp3=0.
      do i=1,nr
        do 10 j=1,nz
          if (ixout(i,j).le.0) goto 10
          rr=r(i)
          zz=z(j)
          uu=u(i,j)
          psi=umax-uu
          if (psi.le.0.) then
            k=ncon
          else if (psi.ge.umax) then
            k=1
          else
            k=1
            do 15 l=1,ncon
              if (psi.gt.psiv(l)) goto 15
              k=l
 15          continue
            if (k.eq.ncon) k=ncon-1
            rat=(psi-psiv(k))/(psiv(k+1)-psiv(k))
            if ((rat.lt.0.).or.(rat.gt.1.)) then
              write(6,*)' ERROR in BETAS***rat=',rat
              write(6,*)' k=',k,' psi=',psi
              write(6,*)' psiv(k)=',psiv(k),' psiv(k+1)=',psiv(k+1)
              stop
            end if
          end if
          palpha=palpha+(fastp(k)+rat*(fastp(k+1)-fastp(k)))*rr*dr*dz
!  calculate core auxiliary current drive, and total current drive to check
          jtotal=jtotal+exph(i,j)*dr*dz
          if (psi.lt.psicut) jcore=jcore+exph(i,j)*dr*dz
!  calculate pressure
          p=press(psi,0)
          ptota=ptota+p*dr*dz
          ptotv=ptotv+p*rr*dr*dz
          if (nimp.ge.2) then
             pimp3=pimp3+bk*densi(psi,3,0)*tempi(psi,3,0)*rr*dr*dz
          end if
          
          if (nimp .ge. 1) then
             !   calculate helium density
             dalf=densi(psi,2,0)
          else
            dalf=0.5*(1.-dil)*dense(psi,0)
          end if
          confa=confa+dalf*rr*dr*dz*2.*pi
          bth=bp(rr,zz)
          bphi=fprof(psi,2)/rr
          btot=bphi*bphi+bth*bth
          bptot=bptot+bth*bth*rr*dr*dz
          btotv=btotv+btot*rr*dr*dz
          te=tempe(psi,0)
          avt=avt+te*rr*dr*dz

          if (te.gt.1.d-10) then
            ti=tempi(psi,1,0)
            ne=dense(psi,0)
            dion=1.0d-19*densi(psi,1,0)
!   include inconsistent dilution factor if no impurities...
            if (imp.eq.0) dion=dil*ne*1.0d-19/zm
            arg=-0.476*(abs(log(1.45d-5*ti)))**2.25
          else
            ne=0.
            dion=0.
            ti=0.
            arg=1.
         end if
          avti=avti+ti*rr*dr*dz
          !Alternate Way to calc fusion power
          !cross section taken from Wesson
          
          !xsec = 1.1e-24*(ti/1000)**2
          !efus = 3.5e6*eq
          !print*, 'Fusion Power Calc !!!'
          !print*, ti/1000, dion*1.0e19, efus
          !pden = xsec*0.25*(dion*1.0e19)**2*efus
          pden=1.27d4*dion**2*exp(arg)
          pfus=pfus+2.*pi*pden*rr*dr*dz
          avel=avel+ne*rr*dr*dz
          avio=avio+dion*rr*dr*dz
          if (imp.eq.1) then
            if (ne.gt.0.) then
              do 5 l=1,nimp+1
                zni=densi(psi,l,0)
                zeffav=zeffav+(zni*rr*dr*dz*iz(l)**2)/ne
 5            continue
            end if
          end if
          vol=vol+rr*dr*dz
          area=area+dr*dz
 10     continue
      end do


      if (fast.eq.1) then
         !write(6,*)' pimp3/palpha=',pimp3/palpha
         allocate (pimp3v(ncon))
         pimp3v=0.
         do k=1,ncon
            psi=psiv(k)
            pimp3v(k)=bk*densi(psi,3,0)*tempi(psi,3,0)
         end do
!         call tstplt2(ncon,psiv,pimp3v,fastp)
      end if


!Fast beta calculation
      if (fast.eq.1) then
         !Write(6,*)'--------------------------------------------'
         !Write(6,*)'Calculating fast beta'
         !Write(6,*)'--------------------------------------------'
!         if (nimp.lt.2) then
!            Write(6,*)'  Error calculating fast beta'
!            Write(6,*)' Input fast alpha as second impurity, check number of impurities'
!            Stop
!         else
!            do i=1,nr
!               do 30 j=1,nz
!                  if (ixout(i,j).le.0) goto 30
!                  rr=r(i)
!                  zz=z(j)
!                  uu=u(i,j)
!                  psi=umax-uu
!                  if (psi.le.0) then
!                     k=ncon
!                  else if (psi.ge.umax) then
!                     k=1
!                  else
!                     k=1
!                     do 35 l=1,ncon
!                        if (psi.gt.psiv(l)) goto 35
!                        k=l
!35                      continue
!                        if (k.eq.ncon) k=ncon-1
!                        rat=(psi-psiv(k))/(psiv(k+1)-psiv(k))
!                        if ((rat.lt.0.).or.(rat.gt.1.)) then
!                           Write(6,*)' *****ERROR IN FAST BETA***** rat=',rat
!                           Write(6,*)' k=',k,' psi=',psi
!                           Write(6,*)' psiv(k)=',psiv(k),' psiv(k+1)=',psiv(k+1)
!                           Stop
!                        end if
!                  end if
!                  pfast=pressi(psi,3,0)
!                  pfastv=pfastv+pfast*rr*dr*dz
!30               continue
!            end do
            fastb=2.*mu0*pimp3*100/btotv
!            Write(6,*)'********************************************'
!            Write(6,*)' fastb=',fastb
!            Write(6,*)'********************************************'
!         end if
         if (debug) write(6,*)' Fast calculation finished'
      end if

!  pressure peaking factor
      psi=0.
      ppeak=press(psi,0)*vol/ptotv
!  Fast alpha fraction
      palpha=palpha/ptotv
!  Calculate beta values
      betai=8.*pi*ptota/(mu0*cur*cur)
      betap=2.*mu0*ptotv/bptot
      conft=1.5*2.*pi*ptotv
      beta=2.*mu0*ptotv*100./btotv
      bvacu=mu0*rodi/(2.*pi*rcen)
      betexp=2.*mu0*ptotv*100./(vol*bvacu*bvacu)
      avt=avt/vol
      avti=avti/vol
      avel=avel/vol
      avio=avio/vol
      if (imp.eq.0) then
        zeffav=zm
      else
        zeffav=zeffav/vol
      end if
!  internal inductance
      rli2=bptot*circumf(1)**2/(mu0*mu0*cur*cur*vol)
      rli3=4.*pi*bptot/(mu0*mu0*cur*cur*rcen)
!  plasma volume
      vol=2.*pi*vol
  end subroutine betas
