      subroutine interface
      use mpi
      use var_inc
      implicit none 

      real factor1, factor2, qqq, Umax, Umax1, T_factor,mu
      real qx,qy,qz,cux,cuy,cuz,tau
      integer  ix,iy,iz,ip,k
      real,dimension(0:npop-1) :: g_Sh,h_Sh,force1
      real,allocatable,dimension(:,:,:,:):: tmpfL, tmpfR, tmphL, tmphR
      real,allocatable,dimension(:,:,:,:):: tmpfU, tmpfD, tmphU, tmphD
      real rh,uh,vh,wh,tt,bw2,Qgh

! ========================================dt hdt=0.5dt
      Umax1=MAXVAL(te(1:lx,1:ly,1:lz))
      Umax1=Vmax+sqrt(gam*R*Umax1)
      CALL MPI_REDUCE(Umax1,Umax,1,MPI_REAL8,MPI_MAX,0, &
                                   mpi_comm_world,ierr)
      dt=CFL*dx(1)/Umax
      CALL MPI_BCAST (dt,1,MPI_REAL8,0,mpi_comm_world,ierr)
      hdt = 0.5*dt
!=================================================
 
      do iz =1,lz
      do iy =1,ly
      do ix = 1,lx
      !rho ux uy uz T
      rh = rho(ix,iy,iz)
      uh = ux(ix,iy,iz)
      vh = uy(ix,iy,iz)
      wh = uz(ix,iy,iz)
      tt = te(ix,iy,iz)
      !mu(T)  tau(rho,T)
      T_factor=exp(omega*log(tt/T_ref))
      mu=mu_ref*T_factor
      tau=mu/(rh*R*tt)
      !(3h)/(2tau+dt)
      bw2=3.*hdt/(2.*tau+dt)
!===========================================      
! Compute heat flux qx,qy,qz
      qx = 0.0
      qy = 0.0
      qz = 0.0
      do ip=0,npop-1
      cux = cix(ip)-uh
      cuy = ciy(ip)-vh
      cuz = ciz(ip)-wh
      Qgh=(cux*cux+cuy*cuy+cuz*cuz)*g(ip,ix,iy,iz)+h(ip,ix,iy,iz) 
      qx = qx + cux*Qgh
      qy = qy + cuy*Qgh
      qz = qz + cuz*Qgh
      end do
      qx=(tau*0.5*qx)/(tau+0.5*dt*Pr)
      qy=(tau*0.5*qy)/(tau+0.5*dt*Pr)
      qz=(tau*0.5*qz)/(tau+0.5*dt*Pr)
!==============================================
!gbar+ hbar+  at cell center   tn
!gS(rho,T,u,q)   hS(rho,T,u,q)
      call Shakhov(rh,tt,uh,vh,wh,qx,qy,qz,g_Sh,h_Sh)

      gp(:,ix,iy,iz) = g(:,ix,iy,iz)-bw2*(g(:,ix,iy,iz)-g_sh(:))
      hp(:,ix,iy,iz) = h(:,ix,iy,iz)-bw2*(h(:,ix,iy,iz)-h_sh(:))
      enddo
      enddo
      enddo
!===================================================
!  DATA COMMUNICATION
      allocate(tmpfL(0:npop-1,lx,ly,2))
      allocate(tmpfR(0:npop-1,lx,ly,2))
      allocate(tmpfU(0:npop-1,lx,2,-1:lz+2))
      allocate(tmpfD(0:npop-1,lx,2,-1:lz+2))
      allocate(tmphL(0:npop-1,lx,ly,2))
      allocate(tmphR(0:npop-1,lx,ly,2))
      allocate(tmphU(0:npop-1,lx,2,-1:lz+2))
      allocate(tmphD(0:npop-1,lx,2,-1:lz+2))


      call exchng5z(gp(:,1:lx,1:ly,1:2),tmpfR,gp(:,1:lx,1:ly,lz-1:lz),tmpfL)
      call exchng5z(hp(:,1:lx,1:ly,1:2),tmphR,hp(:,1:lx,1:ly,lz-1:lz),tmphL)

      gp(:,1:lx,1:ly,-1:0) = tmpfL(:,:,:,1:2)
      gp(:,1:lx,1:ly,lz+1:lz+2) = tmpfR(:,:,:,1:2)
      hp(:,1:lx,1:ly,-1:0) = tmphL(:,:,:,1:2)
      hp(:,1:lx,1:ly,lz+1:lz+2) = tmphR(:,:,:,1:2)

      call exchng5y(gp(:,1:lx,1:2,-1:lz+2),tmpfU,gp(:,1:lx,ly-1:ly,-1:lz+2),tmpfD)
      call exchng5y(hp(:,1:lx,1:2,-1:lz+2),tmphU,hp(:,1:lx,ly-1:ly,-1:lz+2),tmphD)

      gp(:,1:lx,-1:0,:) = tmpfD(:,:,:,:)
      gp(:,1:lx,ly+1:ly+2,:) = tmpfU(:,:,:,:)
      hp(:,1:lx,-1:0,:) = tmphD(:,:,:,:)
      hp(:,1:lx,ly+1:ly+2,:) = tmphU(:,:,:,:)

! then local extension in x
! If periodic in x direction
! periodic
      gp(:,-1:0,:,:) = gp(:,lx-1:lx,:,:)
      gp(:,lx+1:lx+2,:,:) = gp(:,1:2,:,:)
      hp(:,-1:0,:,:) = hp(:,lx-1:lx,:,:)
      hp(:,lx+1:lx+2,:,:) = hp(:,1:2,:,:)
!
      deallocate(tmpfL)
      deallocate(tmpfR)
      deallocate(tmpfU)
      deallocate(tmpfD)
      deallocate(tmphL)
      deallocate(tmphR)
      deallocate(tmphU)
      deallocate(tmphD)                                
      end subroutine interface
!===================================================================
!===========================================================================

      subroutine exchng5z(tmp1,tmp2,tmp3,tmp4)
      use mpi
      use var_inc
      implicit none

!                     1          2       3        4
!                     5          6       7        8
!      call exchng5(f(:,:,:,1),tmpfR,f(:,:,:,lz),tmpfL)


      integer ileny, ilenz
      real, dimension(0:npop-1,lx,ly,2)    :: tmp1, tmp2, tmp3, tmp4


      integer status_array(MPI_STATUS_SIZE,4), req(4)

      ilenz = npop*lx*ly*2

      call MPI_IRECV(tmp2,ilenz,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp4,ilenz,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp1,ilenz,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp3,ilenz,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      end subroutine exchng5z
!===========================================================================
!===================================================================

      subroutine exchng5y(tmp5l,tmp6,tmp7l,tmp8)
      use mpi
      use var_inc
      implicit none

!                     5          6       7        8
!      call exchng5(f(:,:,1,:),tmpfU,f(:,:,ly,:),tmpfD)


      integer ileny, ilenz
      real, dimension(0:npop-1,lx,2,-1:lz+2):: tmp5l, tmp6, tmp7l, tmp8

      integer status_array(MPI_STATUS_SIZE,4), req(4)

      ileny = npop*lx*(lz+4)*2

      call MPI_IRECV(tmp6,ileny,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp8,ileny,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp5l,ileny,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp7l,ileny,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

      end subroutine exchng5y
!===========================================================================

