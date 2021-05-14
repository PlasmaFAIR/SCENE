      subroutine hirsh(psi,rj0,x,bsq,k,bstrap)
!     **********************************************
!
!  Bootstrap current as calculated by Hirshman (Phys.Fluids 31
!  (1988) 3150). Ion current of main ion species only included.
!  Impurities included through their contribution to the electron
!  density (and thereby the electron contribution to the bootstrap
!  current)
!
      use param
      use profiles_mod, only : tempe, tempi, dense, densi
      implicit none

      double precision ted,te,tid,ti,ne,tau,zni,zeff,zb
      double precision pe,pip,ped,pid,ned,znid,dox
      double precision rl310,rl320,alfi0
      double precision x,rj0,psi,bsq,bstrap
      double precision a13,b13,c13,a23,b23,c23
      double precision f1,f2
      double precision rl31,rl32,fa,alfi
      double precision a1e,a1i,a2e,a2i
      integer l,k
!
      ted=tempe(psi,1)
      te=tempe(psi,0)
      tid=tempi(psi,1,1)
      ti=tempi(psi,1,0)
      ne=dense(psi,0)
      zeff=zm
      if (imp.eq.1) then
        if (ne.gt.0.) then
          zeff=0.
          do l=1,nimp+1
            zni=densi(psi,l,0)
            zeff=zeff+(zni*iz(l)**2)/ne
          end do
        else
          nw=10
          write(nw,*)'error*** problem in hirsh, ne=0'
          write(nw,*)'cannot evaluate zeff'
          stop
        end if
      end if
      zb=zeff
      tau=te/ti
      pe=ne*bk*te
! calculate electron and ion contributions to pressure gradient drives...
      ned=dense(psi,1)
      ped=bk*(ne*ted+ned*te)
      zni=densi(psi,1,0)
      znid=densi(psi,1,1)
      pip=zni*bk*ti
      pid=bk*(zni*tid+znid*ti)
      dox=1.414*zb+zb*zb+x*(0.754+2.657*zb+2.*zb*zb)+x*x*    &
          (0.348+1.243*zb+zb*zb)
      rl310=rj0*x*(0.754+2.21*zb+zb*zb+x*(0.348+1.243*zb+zb*zb))  &
           /dox
      rl320=-rj0*x*(0.884+2.074*zb)/dox
      alfi0=-1.172/(1.+0.462*x)
      if (nco.eq.1) then
!  collisionality coeficients from hinton,hazeltine
        a13=0.027*zm*zm-0.211*zm+1.204
        b13=0.14*zm*zm-0.87*zm+1.8
        c13=0.097*zm*zm-0.67*zm+1.643
        a23=0.01*zm*zm-0.08*zm+0.64
        b23=0.088*zm*zm-0.535*zm+1.057
        c23=0.06*zm*zm-0.41*zm+0.97
        f1=(1.+a13*sqrt(tnue(k))+b13*tnue(k))*(1.+c13*cnue(k))
        f2=(1.+a23*sqrt(tnue(k))+b23*tnue(k))*(1.+c23*cnue(k))
        rl31=rl310/f1
        rl32=(rl320+2.5*rl310)/f2-2.5*rl31
        fa=(1.+cnui(k)**2)*(1.+cnue(k)**2)
        alfi=((alfi0+0.35*sqrt(tnui(k)))/(1.+0.7*sqrt(tnui(k)))   &
             +2.1*cnui(k)**2)/fa
      else
!  collisionless model
        rl31=rl310
        rl32=rl320
        alfi=alfi0
      end if
      a1e=-ped/pe
      a1i=-pid/pip
      a2e=-ted/te
      a2i=-tid/ti
!      a2i=0.
!      a2e=0.
      bstrap=(rl31*(a1e+(a1i+alfi*a2i)/(tau*zm))+rl32*a2e)/sqrt(bsq)
      write(6,*)' psin=',psi/umax,' l31=',rl31/rj0,' lb=',bstrap*sqrt(bsq)/(rj0*a1e)
      write(6,*)' l310=',rl310,' l31=',rl31,' f1=',f1,' a13=',a13
      write(6,*)' b13=',b13,' c13=',c13,' zm=',zm,' tnue=',tnue(k),' cnue=',cnue(k)
    end  subroutine hirsh
