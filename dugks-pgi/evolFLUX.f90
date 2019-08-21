      SUBROUTINE EVOL
      use mpi
      use var_inc
      implicit none 
      real gdiff,hdiff
      real time_start2,time_end2

      real factor,facy,facz
      real rh,uh,vh,wh,xp0,yp0,zp0
      real,dimension(0:npop-1) :: gx,gy,gz,grad
      real,dimension(0:npop-1) :: f9,fixl,fbl,fbxp,fil,force1,feq
!
      integer ix,iy,iz,k 
      real qx, qy, qz, cux, cuy, cuz, tau, mu, w, T_factor, u_local, Qgh 
      real Frho, Fux, Fuy, Fuz, Fe

      
! ----  x direction ----
      do iz = 1,lz
      do iy = 1,ly

      call InterpX(iy,iz)

      do ix = 1,lx
      cx = dt/dx(ix)
      g(:,ix,iy,iz)=(4.0*gp(:,ix,iy,iz)-g(:,ix,iy,iz))/3.0 &
                  +cix(:)*cx*(xg_face(:,ix)-xg_face(:,ix+1))
      h(:,ix,iy,iz)=(4.0*hp(:,ix,iy,iz)-h(:,ix,iy,iz))/3.0 &
                  +cix(:)*cx*(xh_face(:,ix)-xh_face(:,ix+1))
      end do
! -- also accumulate conservative quantities 
      do ix = 1,lx
      Frho=0.0
      Fux=0.0
      Fuy=0.0
      Fuz=0.0
      Fe=0.0 
        do k=0,npop-1
        gdiff = xg_face(k,ix)-xg_face(k,ix+1)
        hdiff = xh_face(k,ix)-xh_face(k,ix+1)
        Frho=Frho+cix(k)*gdiff
        Fux=Fux+cix(k)*cix(k)*gdiff
        Fuy=Fuy+cix(k)*ciy(k)*gdiff
        Fuz=Fuz+cix(k)*ciz(k)*gdiff
        Fe=Fe+cix(k)*((cix(k)*cix(k)+ciy(k)*ciy(k)+ciz(k)*ciz(k))*gdiff + hdiff)
        end do
      rho(ix,iy,iz)=rho(ix,iy,iz)+Frho*cx
      Jx(ix,iy,iz)=Jx(ix,iy,iz)+Fux*cx
      Jy(ix,iy,iz)=Jy(ix,iy,iz)+Fuy*cx
      Jz(ix,iy,iz)=Jz(ix,iy,iz)+Fuz*cx
      E(ix,iy,iz)=E(ix,iy,iz)+Fe*cx
      end do

      end do
      end do

! ----  y direction ----
      do iz = 1,lz
      do ix = 1,lx

      call InterpY(ix,iz)

      do iy = 1,ly
      cy = dt/dy(iy)
      g(:,ix,iy,iz)=g(:,ix,iy,iz)+ciy(:)*cy*(yg_face(:,iy)-yg_face(:,iy+1))
      h(:,ix,iy,iz)=h(:,ix,iy,iz)+ciy(:)*cy*(yh_face(:,iy)-yh_face(:,iy+1))
      end do
! -- accumulate conservative quantities
      do iy = 1,ly
      Frho=0.0
      Fux=0.0
      Fuy=0.0
      Fuz=0.0
      Fe=0.0
        do k=0,npop-1
        gdiff = yg_face(k,iy)-yg_face(k,iy+1)
        hdiff = yh_face(k,iy)-yh_face(k,iy+1)
        Frho=Frho+ciy(k)*gdiff
        Fux=Fux+ciy(k)*cix(k)*gdiff
        Fuy=Fuy+ciy(k)*ciy(k)*gdiff
        Fuz=Fuz+ciy(k)*ciz(k)*gdiff
        Fe=Fe+ciy(k)*((cix(k)*cix(k)+ciy(k)*ciy(k)+ciz(k)*ciz(k))*gdiff +hdiff)
        end do
      rho(ix,iy,iz)=rho(ix,iy,iz)+Frho*cy
      Jx(ix,iy,iz)=Jx(ix,iy,iz)+Fux*cy
      Jy(ix,iy,iz)=Jy(ix,iy,iz)+Fuy*cy
      Jz(ix,iy,iz)=Jz(ix,iy,iz)+Fuz*cy
      E(ix,iy,iz)=E(ix,iy,iz)+Fe*cy
      end do

      end do
      end do

! ----  z direction ----
      do iy = 1,ly
      do ix = 1,lx

      call InterpZ(ix,iy)

      do iz = 1,lz
      cz = dt/dz(iz)
      g(:,ix,iy,iz)=g(:,ix,iy,iz)+ciz(:)*cz*(zg_face(:,iz)-zg_face(:,iz+1))
      h(:,ix,iy,iz)=h(:,ix,iy,iz)+ciz(:)*cz*(zh_face(:,iz)-zh_face(:,iz+1))
      end do
! -- accumulate conservative quantities
      do iz = 1,lz
      Frho=0.0
      Fux=0.0
      Fuy=0.0
      Fuz=0.0
      Fe=0.0
        do k=0,npop-1
        gdiff=zg_face(k,iz)-zg_face(k,iz+1)
        hdiff=zh_face(k,iz)-zh_face(k,iz+1)
        Frho=Frho+ciz(k)*gdiff
        Fux=Fux+ciz(k)*cix(k)*gdiff
        Fuy=Fuy+ciz(k)*ciy(k)*gdiff
        Fuz=Fuz+ciz(k)*ciz(k)*gdiff
        Fe=Fe+ciz(k)*((cix(k)*cix(k)+ciy(k)*ciy(k)+ciz(k)*ciz(k))*gdiff + hdiff)
        end do
      rho(ix,iy,iz)=rho(ix,iy,iz)+Frho*cz
      Jx(ix,iy,iz)=Jx(ix,iy,iz)+Fux*cz
      Jy(ix,iy,iz)=Jy(ix,iy,iz)+Fuy*cz
      Jz(ix,iy,iz)=Jz(ix,iy,iz)+Fuz*cz
      E(ix,iy,iz)=E(ix,iy,iz)+Fe*cz

! now update hydro variables using conservative variables
      ux(ix,iy,iz)=Jx(ix,iy,iz)/rho(ix,iy,iz)
      uy(ix,iy,iz)=Jy(ix,iy,iz)/rho(ix,iy,iz)
      uz(ix,iy,iz)=Jz(ix,iy,iz)/rho(ix,iy,iz)
      Te(ix,iy,iz)=(E(ix,iy,iz)-rho(ix,iy,iz)*( ux(ix,iy,iz)*ux(ix,iy,iz)+uy(ix,iy,iz)*uy(ix,iy,iz) &
                    +uz(ix,iy,iz)*uz(ix,iy,iz)  ))/((D+cLK)*rho(ix,iy,iz)*R)
      end do

      end do
      end do

      return 
      END SUBROUTINE EVOL

!-----------------------------------------------------------
      SUBROUTINE InterpX(iy,iz)
      use var_inc
      implicit none
      real fL, fR, delta, x, y, z, ux_face, uy_face, uz_face
      real rho_face, T_face, E_face, qx, qy, qz, Qgh
      real T_factor, tau, wb, mu, cux, cuy, cuz
      integer  iy,iz,ix,k
      real,dimension(0:npop-1) :: g_Sh,h_Sh
      real sgn

! x-direction
      call SlopeX(iy,iz)

      do ix=1,lx+1
        do k=0,npop-1
        delta=0.5*(1.0+sgn(cix(k))) 
        x=0.5*cix(k)*dt
        y=0.5*ciy(k)*dt
        z=0.5*ciz(k)*dt
        fL=gp(k,ix-1,iy,iz)+sgx(k,ix-1)*(0.5*dx(ix-1)-x)-sgy(k,ix-1)*y-sgz(k,ix-1)*z
        fR=gp(k,ix,iy,iz)-sgx(k,ix)*(0.5*dx(ix)+x)-sgy(k,ix)*y-sgz(k,ix)*z
        xg_face(k,ix)=delta*fL+(1.0-delta)*fR

        fL=hp(k,ix-1,iy,iz)+shx(k,ix-1)*(0.5*dx(ix-1)-x)-shy(k,ix-1)*y-shz(k,ix-1)*z
        fR=hp(k,ix,iy,iz)-shx(k,ix)*(0.5*dx(ix)+x)-shy(k,ix)*y-shz(k,ix)*z
        xh_face(k,ix)=delta*fL+(1.0-delta)*fR
        end do
      end do

!  the original f at interface
!  W at interface:
     do ix=1,lx+1
       ux_face=0.0
       uy_face=0.0
       uz_face=0.0
       rho_face=0.0
       E_face=0.0
        do k=0,npop-1
        rho_face=rho_face+xg_face(k,ix)
        ux_face=ux_face+cix(k)*xg_face(k,ix)
        uy_face=uy_face+ciy(k)*xg_face(k,ix)
        uz_face=uz_face+ciz(k)*xg_face(k,ix)
        E_face = E_face + (cix(k)*cix(k)+ciy(k)*ciy(k)+ciz(k)*ciz(k))*xg_face(k,ix)+xh_face(k,ix)
        end do
        ux_face=ux_face/rho_face
        uy_face=uy_face/rho_face
        uz_face=uz_face/rho_face
        T_face=(E_face-rho_face*(ux_face*ux_face+uy_face*uy_face+uz_face*uz_face)) &
                   /((D+cLK)*rho_face*R)
! - get relaxation time
      T_factor=exp(omega*log(T_face/T_ref))
      mu=mu_ref*T_factor 
      mu=mu_ref*T_factor 
      tau=mu/(rho_face*R*T_face) 
      wb=0.5*dt/(2*tau+0.5*dt) 
!- heat flux
      qx=0.0
      qy=0.0
      qz=0.0
           do k=0,npop-1
           cux=cix(k)-ux_face 
           cuy=ciy(k)-uy_face 
           cuz=ciz(k)-uz_face 
           Qgh=(cux*cux+cuy*cuy+cuz*cuz)*xg_face(k,ix)+xh_face(k,ix)
           qx=qx+cux*Qgh 
           qy=qy+cuy*Qgh 
           qz=qz+cuz*Qgh 
           end do
      qx=(tau*0.5*qx)/(tau+0.25*dt*Pr) ! correction from \bar{f}
      qy=(tau*0.5*qy)/(tau+0.25*dt*Pr) ! correction from \bar{f}
      qz=(tau*0.5*qz)/(tau+0.25*dt*Pr) ! correction from \bar{f}

      call Shakhov(rho_face,T_face,ux_face,uy_face,uz_face,qx,qy,qz,g_Sh,h_Sh)
      xg_face(:,ix)=xg_face(:,ix)-wb*(xg_face(:,ix)-g_Sh(:))
      xh_face(:,ix)=xh_face(:,ix)-wb*(xh_face(:,ix)-h_Sh(:))
      end do

      return
      end subroutine InterpX
!-----------------------------------------------------------
      SUBROUTINE InterpY(ix,iz)
      use var_inc
      implicit none
      real fL, fR, delta, x, y, z, ux_face, uy_face, uz_face
      real rho_face, T_face, E_face, qx, qy, qz, Qgh
      real T_factor, tau, wb, mu, cux, cuy, cuz
      integer  iy,iz,ix,k
      real,dimension(0:npop-1) :: g_Sh,h_Sh
      real sgn

! x-direction
      call SlopeY(ix,iz)

      do iy=1,ly+1
        do k=0,npop-1
        delta=0.5*(1.0+sgn(ciy(k))) 
        x=0.5*cix(k)*dt
        y=0.5*ciy(k)*dt
        z=0.5*ciz(k)*dt
        fL=gp(k,ix,iy-1,iz)+sgy(k,iy-1)*(0.5*dy(iy-1)-y)-sgx(k,iy-1)*x-sgz(k,iy-1)*z
        fR=gp(k,ix,iy,iz)-sgy(k,iy)*(0.5*dy(iy)+y)-sgx(k,iy)*y-sgz(k,iy)*z
        yg_face(k,iy)=delta*fL+(1.0-delta)*fR

        fL=hp(k,ix,iy-1,iz)+shy(k,iy-1)*(0.5*dy(iy-1)-y)-shx(k,iy-1)*x-shz(k,iy-1)*z
        fR=hp(k,ix,iy,iz)-shy(k,iy)*(0.5*dy(iy)+y)-shx(k,iy)*x-shz(k,iy)*z
        yh_face(k,iy)=delta*fL+(1.0-delta)*fR
        end do
      end do

!  the original f at interface
!  W at interface:
     do iy=1,ly+1
       ux_face=0.0
       uy_face=0.0
       uz_face=0.0
       rho_face=0.0
       E_face=0.0
        do k=0,npop-1
        rho_face=rho_face+yg_face(k,iy)
        ux_face=ux_face+cix(k)*yg_face(k,iy)
        uy_face=uy_face+ciy(k)*yg_face(k,iy)
        uz_face=uz_face+ciz(k)*yg_face(k,iy)
        E_face = E_face + (cix(k)*cix(k)+ciy(k)*ciy(k)+ciz(k)*ciz(k))*yg_face(k,iy)+yh_face(k,iy)
        end do
        ux_face=ux_face/rho_face
        uy_face=uy_face/rho_face
        uz_face=uz_face/rho_face
        T_face=(E_face-rho_face*(ux_face*ux_face+uy_face*uy_face+uz_face*uz_face)) &
                   /((D+cLK)*rho_face*R)
! - get relaxation time
      T_factor=exp(omega*log(T_face/T_ref))
      mu=mu_ref*T_factor 
      mu=mu_ref*T_factor 
      tau=mu/(rho_face*R*T_face) 
      wb=0.5*dt/(2*tau+0.5*dt) 
!- heat flux
      qx=0.0
      qy=0.0
      qz=0.0
           do k=0,npop-1
           cux=cix(k)-ux_face 
           cuy=ciy(k)-uy_face 
           cuz=ciz(k)-uz_face 
           Qgh=(cux*cux+cuy*cuy+cuz*cuz)*yg_face(k,iy)+yh_face(k,iy)
           qx=qx+cux*Qgh 
           qy=qy+cuy*Qgh 
           qz=qz+cuz*Qgh 
           end do
      qx=(tau*0.5*qx)/(tau+0.25*dt*Pr) ! correction from \bar{f}
      qy=(tau*0.5*qy)/(tau+0.25*dt*Pr) ! correction from \bar{f}
      qz=(tau*0.5*qz)/(tau+0.25*dt*Pr) ! correction from \bar{f}

      call Shakhov(rho_face,T_face,ux_face,uy_face,uz_face,qx,qy,qz,g_Sh,h_Sh)
      yg_face(:,iy)=yg_face(:,iy)-wb*(yg_face(:,iy)-g_Sh(:))
      yh_face(:,iy)=yh_face(:,iy)-wb*(yh_face(:,iy)-h_Sh(:))
      end do

      return
      end subroutine InterpY
!-----------------------------------------------------------
      SUBROUTINE InterpZ(ix,iy)
      use var_inc
      implicit none
      real fL, fR, delta, x, y, z, ux_face, uy_face, uz_face
      real rho_face, T_face, E_face, qx, qy, qz, Qgh
      real T_factor, tau, wb, mu, cux, cuy, cuz
      integer  iy,iz,ix,k
      real,dimension(0:npop-1) :: g_Sh,h_Sh
      real sgn

! x-direction
      call SlopeZ(ix,iy)

      do iz=1,lz+1
        do k=0,npop-1
        delta=0.5*(1.0+sgn(ciz(k))) 
        x=0.5*cix(k)*dt
        y=0.5*ciy(k)*dt
        z=0.5*ciz(k)*dt
        fL=gp(k,ix,iy,iz-1)+sgz(k,iz-1)*(0.5*dz(iz-1)-z)-sgy(k,iz-1)*y-sgx(k,iz-1)*x
        fR=gp(k,ix,iy,iz)-sgz(k,iz)*(0.5*dz(iz)+z)-sgy(k,iz)*y-sgx(k,iz)*x
        zg_face(k,iz)=delta*fL+(1.0-delta)*fR

        fL=hp(k,ix,iy,iz-1)+shz(k,iz-1)*(0.5*dz(iz-1)-z)-shy(k,iz-1)*y-shx(k,iz-1)*x
        fR=hp(k,ix,iy,iz)-shz(k,iz)*(0.5*dz(iz)+z)-shy(k,iz)*y-shx(k,iz)*x
        zh_face(k,iz)=delta*fL+(1.0-delta)*fR
        end do
      end do

!  the original f at interface
!  W at interface:
     do iz=1,lz+1
       ux_face=0.0
       uy_face=0.0
       uz_face=0.0
       rho_face=0.0
       E_face=0.0
        do k=0,npop-1
        rho_face=rho_face+zg_face(k,iz)
        ux_face=ux_face+cix(k)*zg_face(k,iz)
        uy_face=uy_face+ciy(k)*zg_face(k,iz)
        uz_face=uz_face+ciz(k)*zg_face(k,iz)
        E_face = E_face + (cix(k)*cix(k)+ciy(k)*ciy(k)+ciz(k)*ciz(k))*zg_face(k,iz)+zh_face(k,iz)
        end do
        ux_face=ux_face/rho_face
        uy_face=uy_face/rho_face
        uz_face=uz_face/rho_face
        T_face=(E_face-rho_face*(ux_face*ux_face+uy_face*uy_face+uz_face*uz_face)) &
                   /((D+cLK)*rho_face*R)
! - get relaxation time
      T_factor=exp(omega*log(T_face/T_ref))
      mu=mu_ref*T_factor 
      mu=mu_ref*T_factor 
      tau=mu/(rho_face*R*T_face) 
      wb=0.5*dt/(2*tau+0.5*dt) 
!- heat flux
      qx=0.0
      qy=0.0
      qz=0.0
           do k=0,npop-1
           cux=cix(k)-ux_face 
           cuy=ciy(k)-uy_face 
           cuz=ciz(k)-uz_face 
           Qgh=(cux*cux+cuy*cuy+cuz*cuz)*zg_face(k,iz)+zh_face(k,iz)
           qx=qx+cux*Qgh 
           qy=qy+cuy*Qgh 
           qz=qz+cuz*Qgh 
           end do
      qx=(tau*0.5*qx)/(tau+0.25*dt*Pr) ! correction from \bar{f}
      qy=(tau*0.5*qy)/(tau+0.25*dt*Pr) ! correction from \bar{f}
      qz=(tau*0.5*qz)/(tau+0.25*dt*Pr) ! correction from \bar{f}

      call Shakhov(rho_face,T_face,ux_face,uy_face,uz_face,qx,qy,qz,g_Sh,h_Sh)
      zg_face(:,iz)=zg_face(:,iz)-wb*(zg_face(:,iz)-g_Sh(:))
      zh_face(:,iz)=zh_face(:,iz)-wb*(zh_face(:,iz)-h_Sh(:))
      end do

      return
      end subroutine InterpZ

!-----------------------------------------------------------
      SUBROUTINE SlopeX(iy,iz)
      use var_inc
      implicit none
      integer  iy,iz,ix,k
      real sL,sR
      real sgn

!   df/dx
      do ix=0,lx+1
       do k=0,npop-1
       sL=(gp(k,ix,iy,iz)-gp(k,ix-1,iy,iz))/(xc(ix)-xc(ix-1)) 
       sR=(gp(k,ix+1,iy,iz)-gp(k,ix,iy,iz))/(xc(ix+1)-xc(ix))
          if(I_Limiter) then
           sgx(k,ix)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           sgx(k,ix) = 0.5*(SL+SR)
          end if

       sL=(hp(k,ix,iy,iz)-hp(k,ix-1,iy,iz))/(xc(ix)-xc(ix-1)) 
       sR=(hp(k,ix+1,iy,iz)-hp(k,ix,iy,iz))/(xc(ix+1)-xc(ix))
          if(I_Limiter) then 
          shx(k,ix)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else 
          shx(k,ix) = 0.5*(SL+SR)
          end if
       enddo
      end do

!  df/dy
      do ix=0,lx+1
       do k=0,npop-1
       sL=(gp(k,ix,iy,iz)-gp(k,ix,iy-1,iz))/(yc(iy)-yc(iy-1))
       sR=(gp(k,ix,iy+1,iz)-gp(k,ix,iy,iz))/(yc(iy+1)-yc(iy))
          if(I_Limiter) then 
           sgy(k,ix)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           sgy(k,ix) = 0.5*(SL+SR)
          end if

       sL=(hp(k,ix,iy,iz)-hp(k,ix,iy-1,iz))/(yc(iy)-yc(iy-1))
       sR=(hp(k,ix,iy+1,iz)-hp(k,ix,iy,iz))/(yc(iy+1)-yc(iy))
          if(I_Limiter) then 
           shy(k,ix)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           shy(k,ix) = 0.5*(SL+SR)
          end if
       enddo
      end do

!  df/dz
      do ix=0,lx+1
       do k=0,npop-1
       sL=(gp(k,ix,iy,iz)-gp(k,ix,iy,iz-1))/(zc(iz)-zc(iz-1))
       sR=(gp(k,ix,iy,iz+1)-gp(k,ix,iy,iz))/(zc(iz+1)-zc(iz))
          if(I_Limiter) then 
           sgz(k,ix)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           sgz(k,ix) = 0.5*(SL+SR)
          end if

       sL=(hp(k,ix,iy,iz)-hp(k,ix,iy,iz-1))/(zc(iz)-zc(iz-1))
       sR=(hp(k,ix,iy,iz+1)-hp(k,ix,iy,iz))/(zc(iz+1)-zc(iz))
          if(I_Limiter) then 
           shz(k,ix)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           shz(k,ix) = 0.5*(SL+SR)
          end if
       enddo
      end do

! Boundary treatment: 
! this part is no longer needed due to 2-layer expansion
!      sgx(:,0)=0.0
!      shx(:,0)=0.0
!      sgy(:,0)=sgy(:,1)
!      shy(:,0)=shy(:,1)
!      sgz(:,0)=sgz(:,1)
!      shz(:,0)=shz(:,1)

!      sgx(:,lx+1)=0.0
!      shx(:,lx+1)=0.0
!      sgy(:,lx+1)=sgy(:,lx)
!      shy(:,lx+1)=shy(:,lx)
!      sgz(:,lx+1)=sgz(:,lx)
!      shz(:,lx+1)=shz(:,lx)

      return
      end subroutine SlopeX

!-----------------------------------------------------------
      SUBROUTINE SlopeY(ix,iz)
      use var_inc
      implicit none
      integer  iy,iz,ix,k
      real sL,sR
      real sgn

!   df/dx
      do iy=0,ly+1
       do k=0,npop-1
       sL=(gp(k,ix,iy,iz)-gp(k,ix-1,iy,iz))/(xc(ix)-xc(ix-1)) 
       sR=(gp(k,ix+1,iy,iz)-gp(k,ix,iy,iz))/(xc(ix+1)-xc(ix))
          if(I_Limiter) then 
           sgx(k,iy)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           sgx(k,iy) = 0.5*(SL+SR)
          end if

       sL=(hp(k,ix,iy,iz)-hp(k,ix-1,iy,iz))/(xc(ix)-xc(ix-1)) 
       sR=(hp(k,ix+1,iy,iz)-hp(k,ix,iy,iz))/(xc(ix+1)-xc(ix))
          if(I_Limiter) then 
           shx(k,iy)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           shx(k,iy) = 0.5*(SL+SR)
          end if
       enddo
      end do

!  df/dy
      do iy=0,ly+1
       do k=0,npop-1
       sL=(gp(k,ix,iy,iz)-gp(k,ix,iy-1,iz))/(yc(iy)-yc(iy-1))
       sR=(gp(k,ix,iy+1,iz)-gp(k,ix,iy,iz))/(yc(iy+1)-yc(iy))
          if(I_Limiter) then 
           sgy(k,iy)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           sgy(k,iy) = 0.5*(SL+SR)
          end if

       sL=(hp(k,ix,iy,iz)-hp(k,ix,iy-1,iz))/(yc(iy)-yc(iy-1))
       sR=(hp(k,ix,iy+1,iz)-hp(k,ix,iy,iz))/(yc(iy+1)-yc(iy))
          if(I_Limiter) then 
           shy(k,iy)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           shy(k,iy) = 0.5*(SL+SR)
          end if
       enddo
      end do

!  df/dz
      do iy=0,ly+1
       do k=0,npop-1
       sL=(gp(k,ix,iy,iz)-gp(k,ix,iy,iz-1))/(zc(iz)-zc(iz-1))
       sR=(gp(k,ix,iy,iz+1)-gp(k,ix,iy,iz))/(zc(iz+1)-zc(iz))
          if(I_Limiter) then 
           sgz(k,iy)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           sgz(k,iy) = 0.5*(SL+SR)
          end if

       sL=(hp(k,ix,iy,iz)-hp(k,ix,iy,iz-1))/(zc(iz)-zc(iz-1))
       sR=(hp(k,ix,iy,iz+1)-hp(k,ix,iy,iz))/(zc(iz+1)-zc(iz))
          if(I_Limiter) then 
           shz(k,iy)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           shz(k,iy) = 0.5*(SL+SR)
          end if
       enddo
      end do

! Boundary treatment
!      if(indy.eq.0)then
!      sgx(:,0)=sgx(:,1)
!      shx(:,0)=shx(:,1)
!      sgy(:,0)=0.0
!      shy(:,0)=0.0
!      sgz(:,0)=sgz(:,1)
!      shz(:,0)=shz(:,1)
!      end if

!      if(indy.eq.(nprocY-1))then
!      sgx(:,ly+1)=sgx(:,ly)
!      shx(:,ly+1)=shx(:,ly)
!      sgy(:,ly+1)=0.0
!      shy(:,ly+1)=0.0
!      sgz(:,ly+1)=sgz(:,ly)
!      shz(:,ly+1)=shz(:,ly)
!      end if

      return
      end subroutine SlopeY

!-----------------------------------------------------------
      SUBROUTINE SlopeZ(ix,iy)
      use var_inc
      implicit none
      integer  iy,iz,ix,k
      real sL,sR
      real sgn

!   df/dx
      do iz=0,lz+1
       do k=0,npop-1
       sL=(gp(k,ix,iy,iz)-gp(k,ix-1,iy,iz))/(xc(ix)-xc(ix-1)) 
       sR=(gp(k,ix+1,iy,iz)-gp(k,ix,iy,iz))/(xc(ix+1)-xc(ix))
          if(I_Limiter) then 
           sgx(k,iz)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           sgx(k,iz) = 0.5*(SL+SR)
          end if

       sL=(hp(k,ix,iy,iz)-hp(k,ix-1,iy,iz))/(xc(ix)-xc(ix-1)) 
       sR=(hp(k,ix+1,iy,iz)-hp(k,ix,iy,iz))/(xc(ix+1)-xc(ix))
          if(I_Limiter) then 
           shx(k,iz)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           shx(k,iz) = 0.5*(SL+SR)
          end if 
       enddo
      end do

!  df/dy
      do iz=0,lz+1
       do k=0,npop-1
       sL=(gp(k,ix,iy,iz)-gp(k,ix,iy-1,iz))/(yc(iy)-yc(iy-1))
       sR=(gp(k,ix,iy+1,iz)-gp(k,ix,iy,iz))/(yc(iy+1)-yc(iy))
          if(I_Limiter) then 
           sgy(k,iz)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           sgy(k,iz) = 0.5*(SL+SR)
          end if

       sL=(hp(k,ix,iy,iz)-hp(k,ix,iy-1,iz))/(yc(iy)-yc(iy-1))
       sR=(hp(k,ix,iy+1,iz)-hp(k,ix,iy,iz))/(yc(iy+1)-yc(iy))
          if(I_Limiter) then 
           shy(k,iz)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           shy(k,iz) = 0.5*(SL+SR)
          end if
       enddo
      end do

!  df/dz
      do iz=0,lz+1
       do k=0,npop-1
       sL=(gp(k,ix,iy,iz)-gp(k,ix,iy,iz-1))/(zc(iz)-zc(iz-1))
       sR=(gp(k,ix,iy,iz+1)-gp(k,ix,iy,iz))/(zc(iz+1)-zc(iz))
          if(I_Limiter) then 
           sgz(k,iz)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           sgz(k,iz) = 0.5*(SL+SR)
          end if

       sL=(hp(k,ix,iy,iz)-hp(k,ix,iy,iz-1))/(zc(iz)-zc(iz-1))
       sR=(hp(k,ix,iy,iz+1)-hp(k,ix,iy,iz))/(zc(iz+1)-zc(iz))
          if(I_Limiter) then 
           shz(k,iz)=(sgn(sR)+sgn(sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+1.0e-10)
          else
           shz(k,iz) = 0.5*(SL+SR)
          end if
       enddo
      end do

! Boundary treatment
!      if(indz.eq.0)then
!      sgx(:,0)=sgx(:,1)
!      shx(:,0)=shx(:,1)
!      sgy(:,0)=sgy(:,1)
!      shy(:,0)=shy(:,1)
!      sgz(:,0)=0.0
!      shz(:,0)=0.0
!      end if

!      if(indz.eq.(nprocZ-1))then
!      sgx(:,lz+1)=sgx(:,lz)
!      shx(:,lz+1)=shx(:,lz)
!      sgy(:,lz+1)=sgy(:,lz)
!      shy(:,lz+1)=shy(:,lz)
!      sgz(:,lz+1)=0.0
!      shz(:,lz+1)=0.0 
!      end if

      return
      end subroutine SlopeZ

