      subroutine para
      use mpi
      use var_inc
      implicit none

      real,dimension(ny)   :: yc9,dy9
      real,dimension(nz)   :: zc9,dz9
      integer i,j,k,ip,jg,kg

      pi = 4.0*atan(1.0) 
      pi2 = 2.0*pi

      R = 0.5!gas contant
      !CFL number
      CFL = 0.5
!      CFL=0.4
      !D   \xi
      D = 3.0
      !L=3-D  \eta
      cL=0.0
      !K    \cta
      cK = 2.0
      !3-D+K  \eta \cta
      cLK=cL+cK
      !gam=cp/cv
      gam = (D+cL+cK+2.0)/(D+cL+cK)
      !Pr number
      Pr = 0.70
      !in model mu 
      omega = 0.76 
      T_ref=15.89
      mu_ref=0.00262   
      !reference density and kinematic viscosity
      rho0 = 1.0
      visc = mu_ref/rho0
      
       ttt = 0.0
!      tmassin=0.0 
!      tkein=0.0
!      tengin=0.0
!==========================================================
!set discreted velocities in Shan 2006
      call setVelocity

!==========================================================
!MPI
      nprocY =8  !!!!THIS IS MEANT TO BE CHANGED WITH NPROC
      nprocZ = nproc/nprocY
      ly = ny/nprocY         !local division of dist for procs in y dir
      lz = nz/nprocZ         !local division of dist for procs in z dir
      indy = mod(myid,nprocY)
      indz = int(myid/nprocY)
!========================================================
      allocate(yf(0:ly+2))
      allocate(dy(0:ly+1))
      allocate(yc(-1:ly+2))
      allocate(zf(0:lz+2))
      allocate(dz(0:lz+1))
      allocate(zc(-1:lz+2))
      allocate(ycg(1:ny))
      allocate(dyg(1:ny))
      allocate(zcg(1:nz))
      allocate(dzg(1:nz))

      tfluxx1=0.0
      tfluxx2=0.0
      tfluxx3=0.0

!========================================================
!grid  uniform
      dxyz = 2.0*pi/real(nx)
      do i=1,lx+1
      xf(i) = real(i-1)*dxyz
      end do

      do j=0,ly+2
      yf(j) = real(j+indy*ly-1)*dxyz
      end do

      do k=0,lz+2
      zf(k) = real(k+indz*lz-1)*dxyz
      end do
! Scale the maximum grid spacing to 1.0
      dxmax = xf(lx/2+1)-xf(lx/2)
      cwidth = real(nx)
      dymax = dxmax
      dzmax = dxmax
      ywidth = real(ny)
      zwidth = real(nz)
      
      if(myid.eq.0)write(*,*)'xf=',xf
      if(myid.eq.0)write(*,*)'yf=',yf
      if(myid.eq.0)write(*,*)'zf=',zf
!===============================================
      do i=1,lx
      dx(i) = xf(i+1)-xf(i)
      xc(i) = xf(i) + dx(i)/2.
      end do
      dx(0) = dx(1)
      dx(lx+1) = dx(lx)
      xc(0) = xf(1) - dx(0)/2.
      xc(lx+1) = xf(lx+1) + dx(lx+1)/2.
! the next two lines must be modified for nonuniform mesh
      xc(-1) = xc(0) - dx(0)
      xc(lx+2) = xc(lx+1)+dx(lx+1)!correction - to +

      if(myid.eq.0)write(*,*)'dx=',dx
      if(myid.eq.0)write(*,*)'xc=',xc

!===============================================
      do j=0,ly+1
      dy(j) = yf(j+1)-yf(j)
      yc(j) = yf(j)+dy(j)/2
      end do

      if(indy.eq.0)then
      dy(0) = dy(1)
      yc(0) = yf(1) - dy(0)/2.0
      end if

      if(indy.eq. (nprocY-1) )then
      dy(ly+1) = dy(ly)
      yc(ly+1) = yf(ly+1) + dy(ly+1)/2.0
      end if

!   need to reconsider for nonuniform mesh
      yc(-1) = yc(0) - dy(0)
      yc(ly+2) = yc(ly+1) + dy(ly+1)

!===============================================
      do k=0,lz+1
      dz(k) = zf(k+1)-zf(k)
      zc(k) = zf(k)+dz(k)/2.0
      end do
      if(indz.eq.0)then
      dz(0) = dz(1)
      zc(0) = zf(1) - dz(0)/2.0
      end if

      if(indz.eq. (nprocZ-1) )then
      dz(lz+1) = dz(lz)
      zc(lz+1) = zf(lz+1) + dz(lz+1)/2.0
      end if

!   need to reconsider for nonuniform mesh
      zc(-1) = zc(0) - dz(0)
      zc(lz+2) = zc(lz+1) + dz(lz+1)

!===============================================
       do j=1,ly
       jg = j+indy*ly
       dy9(jg)=dy(j)
       yc9(jg)=yc(j)
       end do

       do k=1,lz
       kg = k+indz*lz
       zc9(kg)=zc(k)
       dz9(kg)=dz(k)
       end do

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE (yc9,ycg,ny,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (dy9,dyg,ny,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (zc9,zcg,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (dz9,dzg,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       ycg = ycg/real(nprocZ)
       dyg = dyg/real(nprocZ)
       zcg = zcg/real(nprocY)
       dzg = dzg/real(nprocY)

      iseedf = 232300
!scales////////////////////////////////////////////////
      vscale = 1.0
      escale = 1.0
      dscale = 1.0
      tscale = 1.0
      istep0 = 0
      istep00 = 1
      istat = 0
      imovie = 0
!////////////////////////////////////////////////////////
      istpload =0    !!!!THIS NEEDS TO BE CHANGED WHEN STARTING NEW RUNS
      nsteps = 80000
      nforcing = 0

      nek = int(nx/2 - 1.5)

      newrun = .true.
!       newrun = .false.
      newinitflow = .true.
!      newinitflow = .false.
 
! find the neighboring nodes to each processor

      if((mod(myid,nprocY)+1).eq.nprocY) then
        lyext=2
      else
        lyext=0
      endif
      lly = ly+lyext

!******create index for proccessors*******

      mzp = mod(indz+1,nprocZ) * nprocY + indy  !top
      mzm = mod(indz + nprocZ - 1,nprocZ) * nprocY + indy !bottom

      myp = indz*nprocY + mod(indy+1,nprocY) !right
      mym = indz*nprocY + mod(indy+nprocY-1,nprocY) !left


      mypzp = mod(indz+1,nprocZ)*nprocY + mod(indy+1,nprocY) 
      mypzm = mod(indz+nprocZ-1,nprocZ)*nprocY + mod(indy+1,nprocY)
      mymzp = mod(indz+1,nprocZ)*nprocY + mod(indy+nprocY-1,nprocY)
      mymzm = mod(indz+nprocZ-1,nprocZ)*nprocY + mod(indy+nprocY-1,nprocY)


      teepsl = 1.e-07

! saving and loading directories relevant to flow
      dirgenr = '/home/zhangyq/dugks-pgi/'
      dirdiag = trim(dirgenr)//''
      dirstat = trim(dirgenr)//''
      dirinitflow = trim(dirgenr)//''
      dircntdflow =trim(dirgenr)//''
      dirflowout = trim(dirgenr)//''
      dirmoviedata = trim(dirgenr)//''

! Computing wave numbers and other values for the forcing
      do i = 1, nx+2, 2
       kxr(i)   = int(i/2)
       kxr(i+1) = kxr(i)
      end do

      do j = 1, ny
       if ( j.lt.ny/2+2 ) then
        kyr(j) = j - 1
       else
        kyr(j) = -(ny+1-j)
       endif
      end do

      do k = 1, nz
       if ( k.lt.nz/2+2 ) then
        kzr(k) = k - 1
       else
        kzr(k) = -(nz+1-k)
       endif
      end do

      iseedf = -iseedf
      ivf(:)  = 0
      iyf     = 0
      b1r = 0.0
      b2r = 0.0
      b3r = 0.0

      end subroutine para
!==================================================================

      subroutine allocarray
      use var_inc
      implicit none

      allocate (g(0:npop-1,lx,ly,lz))
!      allocate (geq0(0:npop-1,lx,ly,lz))
      allocate (gp(0:npop-1,-1:lx+2,-1:ly+2,-1:lz+2))

      allocate (h(0:npop-1,lx,ly,lz))
!      allocate (heq0(0:npop-1,lx,ly,lz))
      allocate (hp(0:npop-1,-1:lx+2,-1:ly+2,-1:lz+2))

      allocate (rho(lx,ly,lz))
      allocate (rhop(lx,ly,lz))

      allocate (Jx(lx,ly,lz))
      allocate (Jy(lx,ly,lz))
      allocate (Jz(lx,ly,lz))

      allocate (E(lx,ly,lz))

      allocate (udiv(lx,ly,lz))

      allocate (ux(lx,ly,lz))
      allocate (uy(lx,ly,lz))
      allocate (uz(lx,ly,lz))

      allocate (te(lx,ly,lz))

      allocate (uxp(lx,ly,lz))
      allocate (uyp(lx,ly,lz))
      allocate (uzp(lx,ly,lz))

      allocate (tep(lx,ly,lz))

      allocate (ox(lx,ly,lz))
      allocate (oy(lx,ly,lz))
      allocate (oz(lx,ly,lz))

      allocate (kx(lx+2,ly+lyext,lz))
      allocate (ky(lx+2,ly+lyext,lz))
      allocate (kz(lx+2,ly+lyext,lz))
      allocate (k2(lx+2,ly+lyext,lz))
      allocate (ik2(lx+2,ly+lyext,lz))

      allocate (ddiv(lx+2,ly+lyext,lz))
      allocate (prs(lx+2,ly+lyext,lz))

      allocate (vx(lx+2,ly+lyext,lz))
      allocate (vy(lx+2,ly+lyext,lz))
      allocate (vz(lx+2,ly+lyext,lz))
      allocate (wx(lx+2,ly+lyext,lz))
      allocate (wy(lx+2,ly+lyext,lz))
      allocate (wz(lx+2,ly+lyext,lz))

      allocate(force_realx(lx,ly,lz))
      allocate(force_realy(lx,ly,lz))
      allocate(force_realz(lx,ly,lz))


      end subroutine allocarray
!==================================================================



