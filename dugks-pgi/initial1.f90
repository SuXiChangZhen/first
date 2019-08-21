!==================================================================
      SUBROUTINE initpop!only geq and heq 
      use var_inc
      use mpi
      implicit none

      real u9,v9,w9,rho9,t9,Umax1,Umax
      real,dimension(0:npop-1) :: feq
      integer ip,ix,iy,iz
!=============density
      rho = 1.0
!============temperature
      Te = T_ref
!===========momentum
      Jx = rho*ux
      Jy = rho*uy
      Jz = rho*uz
!==========total energy rho*(u^2+(3+K)RT)
      E = rho*(ux*ux+uy*uy+uz*uz+(D+cLK)*R*te)

!====================================================
      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx

      rho9 = rho(ix,iy,iz)
      u9 = ux(ix,iy,iz)
      v9 = uy(ix,iy,iz)
      w9 = uz(ix,iy,iz)
      t9 = te(ix,iy,iz)

      call feqM(rho9,t9,u9,v9,w9,feq)
!g~
      g(:,ix,iy,iz)=feq(:)
!h~
      h(:,ix,iy,iz)=feq(:)*cLK*R*Te(ix,iy,iz)

      enddo
      enddo
      enddo
!======================================================
!/////////////////////////////////////////dt
      Umax1=MAXVAL(te(1:lx,1:ly,1:lz))
      Umax1=Vmax+sqrt(gam*R*Umax1)
      CALL MPI_REDUCE(Umax1,Umax,1,MPI_REAL8,MPI_MAX,0, &
                                   mpi_comm_world,ierr)
      dt=CFL*dx(1)/Umax
      CALL MPI_BCAST (dt,1,MPI_REAL8,0,mpi_comm_world,ierr)
      if(myid.eq.0) write(*,*)'dt = ', dt
!////////////////////////////////////////////////////////////
      END SUBROUTINE initpop
!==================================================================
! NOTE: here the irand() is replaced by random_seed() and
! random_number() function to make the subroutine suitable for the intel
! compiler on yellowstone
      SUBROUTINE initrand(inseed)
      use mpi
      use var_inc
      use curand
      implicit none

      integer :: rng 
      integer i, k, inseed
      integer, dimension(nproc):: iseed
      integer, dimension(1):: seed, oldseed, myseed
      double precision:: temprand
      k = 1

      call random_seed

      if(myid == 0)then
        seed = inseed
        write(*,*)'inseed',size(seed)
        call random_seed(PUT = seed(1:k))
        write(*,*)'rand_1'
        do i = 1,nproc
          call random_seed(get = oldseed(1:k))
          write(*,*)'seeddone'
          call random_number(temprand)
          iseed(i) = oldseed(1)
        end do
      end if

      call MPI_BCAST(iseed,nproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      myseed(1) = iseed(myid + 1)
      call random_seed(put = myseed(1:k))

      END SUBROUTINE initrand
!==================================================================
      SUBROUTINE initvel
      use mpi   
      use var_inc
      implicit none

      real, dimension(lx+2,lly,lz) :: tmp
      real ek, e_t, k9, kpeak, Ek9
      integer ik

      call gaussian(vx)
      call gaussian(vy)
      call gaussian(vz)

      tmp = ( kx*vx + ky*vy + kz*vz ) / k2
      vx  = vx - kx*tmp
      vy  = vy - ky*tmp
      vz  = vz - kz*tmp

      call symmetrize(vx)
      call symmetrize(vy)
      call symmetrize(vz)

      where(ik2 > nek)
        vx = 0.0
        vy = 0.0
        vz = 0.0
      end where

      tmp = vx*vx + vy*vy + vz*vz
      if(indy.eq.0) then
         tmp(:,1,:) = 0.5 * tmp(:,1,:)
         tmp(:,2,:) = 0.5 * tmp(:,2,:)
      endif

      kpeak = 4.0
      do ik = 1,nek
!  Form 1
        k9 = (real(ik)/kpeak)**2
        Ek9 = 0.011*(real(ik))**4*exp(-2.0*k9)
!  Form 2
!       EK9 = 0.0d0
!       if(ik.ge.3 .and. ik.le.8)then
!       k9 = real(ik)
!       Ek9 = 1.1474d-2*k9**4*exp(-0.14*k9*k9)
!       end if

        ek = 0.0
        e_t = 0.0
        ek = sum(tmp, mask=(ik2 == ik))
        !write(*,*)'tmp as shown',size(tmp)
        !write(*,*)'Pre_mpiar',ek
        call MPI_ALLREDUCE(ek,e_t,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
        !write(*,*)'post_mpiar'
        e_t = sqrt(Ek9/e_t)

        where(ik2 == ik)
          vx = vx * e_t
          vy = vy * e_t
          vz = vz * e_t
        end where
      end do

      write(*,*)'3crx'
      call mpifft3DCR(vx)
      write(*,*)'3cry'
      call mpifft3DCR(vy)
      write(*,*)'3crz'
      call mpifft3DCR(vz)

      ux = vx(1:lx,1:ly,:)
      uy = vy(1:lx,1:ly,:)
      uz = vz(1:lx,1:ly,:)

      END SUBROUTINE initvel


!-----------------------------------------------------------------
! The sign function
      REAL FUNCTION SGN(X)
      !$acc routine
      real x
      if(x.eq.0.0)then
      SGN = 0.0
      else if(x.lt.0.0)then
      SGN = -1.0
      else
      SGN = 1.0
      end if
      RETURN
      END FUNCTION SGN
!-----------------------------------------------------------------
!=========================================================
!contrain velocity and density
      SUBROUTINE vel_den_constrain
      use var_inc
      implicit none

      real T_factor,mu
      real qx,qy,qz,cux,cuy,cuz,tau
      integer  ix,iy,iz,ip,k
      real,dimension(0:npop-1) :: g_Sh,h_Sh
      real rh,uh,vh,wh,tt,Qgh

! constraining
      do iz=1,lz
      do iy=1,ly
      do ix=1,lx
      rh = rho(ix,iy,iz)
      uh = ux(ix,iy,iz)
      vh = uy(ix,iy,iz)
      wh = uz(ix,iy,iz)
      tt = te(ix,iy,iz)

      T_factor=exp(omega*log(tt/T_ref))
      mu=mu_ref*T_factor
      tau=mu/(rh*R*tt)
! Compute heat flux
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

      call Shakhov(rh,tt,uh,vh,wh,qx,qy,qz,g_Sh,h_Sh)
      g(:,ix,iy,iz) = g(:,ix,iy,iz) - g_sh(:)
      h(:,ix,iy,iz) = h(:,ix,iy,iz) - h_sh(:)

      rh = rhop(ix,iy,iz)
      uh = uxp(ix,iy,iz)
      vh = uyp(ix,iy,iz)
      wh = uzp(ix,iy,iz)

      call Shakhov(rh,tt,uh,vh,wh,qx,qy,qz,g_Sh,h_Sh)

      g(:,ix,iy,iz) = g(:,ix,iy,iz) + g_sh(:)
      h(:,ix,iy,iz) = h(:,ix,iy,iz) + h_sh(:)

      rho(ix,iy,iz) = rh
      ux(ix,iy,iz) = uh
      uy(ix,iy,iz) = vh
      uz(ix,iy,iz) = wh
 
      end do
      end do
      end do

      END SUBROUTINE vel_den_constrain
!=========================================================
!========================only to constrain velocity
      SUBROUTINE vel_constrain
      use var_inc
      implicit none

      real T_factor,mu
      real qx,qy,qz,cux,cuy,cuz,tau
      integer  ix,iy,iz,ip,k
      real,dimension(0:npop-1) :: g_Sh,h_Sh
      real rh,uh,vh,wh,tt,Qgh

! constraining
      do iz=1,lz
      do iy=1,ly
      do ix=1,lx
      uh = ux(ix,iy,iz)
      vh = uy(ix,iy,iz)
      wh = uz(ix,iy,iz)
      rh = rho(ix,iy,iz)
      tt = te(ix,iy,iz)

      T_factor=exp(omega*log(tt/T_ref))
      mu=mu_ref*T_factor
      tau=mu/(rh*R*tt)
! Compute heat flux
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

      call Shakhov(rh,tt,uh,vh,wh,qx,qy,qz,g_Sh,h_Sh)
      g(:,ix,iy,iz) = g(:,ix,iy,iz) - g_sh(:)
      h(:,ix,iy,iz) = h(:,ix,iy,iz) - h_sh(:)

      uh = uxp(ix,iy,iz)
      vh = uyp(ix,iy,iz)
      wh = uzp(ix,iy,iz)

      call Shakhov(rh,tt,uh,vh,wh,qx,qy,qz,g_Sh,h_Sh)

      g(:,ix,iy,iz) = g(:,ix,iy,iz) + g_sh(:)
      h(:,ix,iy,iz) = h(:,ix,iy,iz) + h_sh(:)

      ux(ix,iy,iz) = uh
      uy(ix,iy,iz) = vh
      uz(ix,iy,iz) = wh
 
      end do
      end do
      end do

      END SUBROUTINE vel_constrain

