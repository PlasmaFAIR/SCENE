      subroutine tester
!     *****************
!
!    Plots various profiles and outputs run parameters.
!
      use param
      implicit none
      integer i,j
      double precision press,dense,densi,fprof,tempe,tempi,bp
      double precision psi,psi1,psi2,psi3,pv,pp,pd,ppp,pdd
!
      do i=2,ncon-1
        psi=psiv(i)
        psi1=psiv(i-1)
        psi2=psiv(i+1)
        pp=(press(psi1,0)-press(psi2,0))/(psi1-psi2)
        ppp=(press(psi1,1)-press(psi2,1))/(psi1-psi2)
        pd=press(psi,1)
        pdd=press(psi,2)
        write(26,*)' ******i=',i,' p=',press(psi,0),' ipswtch=',ipswtch
        write(26,*)' pp=',pp,' pd=',pd
        write(26,*)' ppp=',ppp,' pdd=',pdd
      end do
      return
      end
        
