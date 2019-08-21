!===================full range (not expanded)
!geq
      subroutine feqM(rh,tt,uh,vh,wh,feq)
      use var_inc
      implicit none
      real, dimension(0:npop-1) :: feq,eux,euy,euz,csqq
      real rh,tt,uh,vh,wh

      eux(:) = cix(:) - uh
      euy(:) = ciy(:) - vh
      euz(:) = ciz(:) - wh
      csqq(:) = (eux(:)*eux(:) + euy(:)*euy(:) + euz(:)*euz(:))/(2.0*R*tt)
      feq(:) = tp(:)*rh*exp( -csqq(:) )*(T_ref/tt)**1.5   !! 1.5 = D/2
      return
      end

!gS hS
      subroutine Shakhov(rh,tt,uh,vh,wh,qqx,qqy,qqz,g_Sh,h_Sh)
      use var_inc
      implicit none
      real, dimension(0:npop-1) :: feq,g_Sh,h_Sh,cq,cn,fac
      real rh,tt,uh,vh,wh,uv,qqx,qqy,qqz

      call feqM(rh,tt,uh,vh,wh,feq)

      cq(:)=((cix(:)-uh)*qqx+(ciy(:)-vh)*qqy+(ciz(:)-wh)*qqz)/(5.0*rh*R*R*tt*tt)
      cn(:)=((cix(:)-uh)*(cix(:)-uh)+(ciy(:)-vh)*(ciy(:)-vh)+(ciz(:)-wh)*(ciz(:)-wh))/(R*tt) - 5.0
      fac(:) = (1.-Pr)*cq(:)
      g_Sh(:)=feq(:)*(1.0 + fac(:)*(cn(:)+cL))
      h_Sh(:)=(cLK+fac(:)*(cn(:)*cLK+cL*(cLK+2)) )*R*TT*feq(:)

      return
      end

!----------------------------------------------------
      SUBROUTINE setVelocity
      use mpi
      use var_inc
      implicit none
      integer k
      real vr, vs, vt, w0, wr, ws, wt,tpsum

! Use in Shan et al (2006) page 438, correct the missing 2 in t^2 before
! sqrt(15)
! E^27_{3,7}, top sign
! I used the first (top) sign choice 
      c0=sqrt(15.0)
      vr=(15.0+c0)/2.0
      vs=6.0-c0
      vt=9.0+2.0*c0
      vr=sqrt(vr)
      vs=sqrt(vs)
      vt=sqrt(vt)
      w0 = (720.0+8.0*c0)/2205.0
      wr = (270.0-46.0*c0)/15435.0
      ws = (162.0+41.0*c0)/6174.0
      wt = (783.0-202.0*c0)/24696.0

      c0=sqrt(R*T_ref)
      vr = vr*c0
      vs = vs*c0
      vt = vt*c0
      Vmax=sqrt(3.0)*vt

      cix=(/0.,vr,-vr,0., 0.,0., 0.,vs, vs,-vs,-vs,vs, vs,-vs,-vs,0., 0., 0., 0., &
           vt, vt, vt, vt,-vt,-vt,-vt,-vt/)
      ciy=(/0.,0., 0.,vr,-vr,0., 0.,vs,-vs, vs,-vs,0., 0., 0., 0.,vs, vs,-vs,-vs, &
           vt, vt,-vt,-vt, vt, vt,-vt,-vt/)
      ciz=(/0.,0., 0., 0,  0,vr,-vr,0., 0., 0., 0.,vs,-vs, vs,-vs,vs,-vs, vs,-vs, &
           vt,-vt, vt,-vt, vt,-vt, vt,-vt/)
      ipopp=(/0,2,  1, 4,  3, 6, 5,10,  9,  8,  7,14, 13, 12, 11, 18, 17, 16, 15, &
           26, 25, 24, 23, 22, 21, 20, 19/)
      tp =(/w0,wr,wr, wr, wr,wr, wr,ws,ws, ws, ws, ws, ws, ws, ws,ws, ws, ws, ws, &
           wt, wt, wt, wt, wt, wt, wt, wt/)


      if(myid.eq.0)then
      tpsum = 0.0
      do k=0,npop-1
      tpsum = tpsum + tp(k)
      end do
      write(*,*)'tpsum=',tpsum
      end if
! different from ZLG's half GH code, tp was no multiplied by c0
! tp is only used in feq, one single place
      do k = 0,npop-1
      tp(k) = tp(k)*exp((cix(k)*cix(k)+ciy(k)*ciy(k)+ciz(k)*ciz(k))/(2.0*c0*c0))
      end do

      RETURN
      END SUBROUTINE setVelocity

