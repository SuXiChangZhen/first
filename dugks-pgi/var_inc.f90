      module var_inc
      implicit none

      integer,parameter:: FFTW_FORWARD = -1, FFTW_BACKWARD = 1
      integer,parameter:: FFTW_REAL_TO_COMPLEX = -1,                   &
                          FFTW_COMPLEX_TO_REAL = 1 
      integer,parameter:: FFTW_ESTIMATE = 0, FFTW_MEASURE = 1
      integer,parameter:: FFTW_OUT_OF_PLACE = 0, FFTW_IN_PLACE = 8,    &
                          FFTW_USE_WISDOM = 16
      integer,parameter:: FFTW_THREADSAFE = 128 
      integer(kind = 8):: plan_RC, plan_CR, plan_F, plan_B  
 
      integer,parameter:: nx = 256, ny = 256, nz = 256
      integer,parameter:: lx = nx,nxz=nx*nz,nyz=ny*nz

      integer,parameter:: lxh = lx/2, lyh = ny/2
      integer,parameter:: nxh = nx/2, nyh = ny/2, nzh = nz/2
      integer,parameter:: npop = 27, npart = 1000
      integer,parameter:: iflowseed = 54321
      integer,parameter:: NTAB = 32 
      integer,parameter:: ndiag=10, nstat=10, nspec=500
      integer,parameter:: nflowout = 500
      integer,parameter:: nnuout = 8000
      integer,parameter:: nmovieout = 20000000, nsij = 4000     
      integer,parameter:: irelease = 1
      integer,parameter:: iprocrate = 2  

      real,parameter:: Tpd = 2000.0,Tpdp = 1500.0,beta9 = 3.0,gamma9=2.0
      real xl0,xls,phase9

      real rho0
      real tfluxx,tfluxx1,tfluxx2,tfluxx3

      real dt,hdt,CFL,grid_p,cwidth,dxmax,dymax,dzmax,ywidth,zwidth,dxyz
      integer ierr, myid, nproc
!**********temperature parameters**************
      real Th,Tc,T0,grav,Ra,DTE
!*********changes
      integer nprocY, nprocZ
     
      integer istep, istep0, istep00, nsteps, nforcing, istpload, imovie   
      integer lz, ly, lyext, lly, nek, MRTtype, mzp, mzm, istat
!*********changes
      integer indy, indz, myp, mym, mypzp, mypzm, mymzp, mymzm

      integer iseedp, msize, nps, iyp
      logical newrun, ipart, newinitflow, I_Limiter

      real force_in_y,force_in_x,ustar,u0,Rstar,ystar
      real tmassin, tkein, tengin
      real R,D,cL,cK,cLK,Pr,gam,omega,mu_ref,T_ref,Vmax,ttt, TEND, c0
      real dlv,cx,cy,cz   !cx,cy,cz represent wx,wy,wz in Wang Peng's code
      real pi, pi2, anu, visc, u0in 
      real vscale, escale, dscale, tscale 
      real omegepsl, omegepslj, omegxx 
      real s1, s2, s4, s9, s10, s13, s16
      real coef1, coef2, coef3, coef4, coef5, coef3i, coef4i
      real val1, val2, val3, val4, val5, val6, val7, val8, val9,       &
           val1i, val2i, val3i, val4i, val5i, val6i, val7i, val8i, val9i
      real ww0, ww1, ww2  
      real teerr,teerrm, rhoerr,rhomax,rhomin, rhoerrm, teepsl
      real volp, amp, aip, gscale, rhog
      real stf0, stf1
      real time_start, time_end, time_diff, time_max,time_sum,      &
                       time_lmt, time_buff, time_bond,     &
                       time1in,time2in
      real time1,time2,time_coll,time_strea,time_macro,   &
           time_collmax,time_streamax,time_macromax,      &
           time_collmin,time_streamin,time_macromin,      &
           time_stXp,time_stXm,time_stYp,time_stYm, &
           time_stZp,time_stZm,    &
           time_stXpmax,time_stXmmax,time_stYpmax,time_stYmmax, &
           time_stZpmax,time_stZmmax,    &
           time_stXpmin,time_stXmmin,time_stYpmin,time_stYmmin, &
           time_stZpmin,time_stZmmin
       real time_stream_start,time_stream_end,time_stream,  &
            time_stream_start2,time_stream_end2, time_stream2
      real time_stream_max, time_stream_avg, time_stream_sum
 
      integer,dimension(0:npop-1):: ipopp
      real,dimension(0:npop-1,0:nx+1):: sgx,sgy,sgz,shx,shy,shz
      real,dimension(0:npop-1,0:nx+1):: xg_face,yg_face,zg_face,xh_face,yh_face,zh_face
      !I think 1:nx+1 is enough 
      real,dimension(0:npop-1):: cix, ciy, ciz
      real,dimension(0:npop-1):: tp
      integer,dimension(NTAB):: ivp 
      real,dimension(npop-1):: wwp 

      INTEGER                :: iseedf, iyf
      integer,dimension(NTAB):: ivf
      real,dimension(6,5,5)  :: a1r,a2r,a3r,b1r,b2r,b3r
      real,dimension(nx+2)   :: kxr
      real,dimension(1:nx+1) :: xf
      real,dimension(0:nx+1) :: dx
      real,dimension(-1:nx+2) :: xc
      real,dimension(ny)     :: kyr,kzr
      real,allocatable,dimension(:,:,:):: force_realx,force_realy, &
                                            force_realz

      integer,allocatable,dimension(:,:,:):: ik2 
      integer,allocatable,dimension(:):: icouples 
      integer,allocatable,dimension(:):: ipglb

! The main memory array
      real,allocatable,dimension(:):: yf,dy,yc,zf,dz,zc,ycg,zcg,dzg,dyg
      real,allocatable,dimension(:,:,:,:):: g,h
      real,allocatable,dimension(:,:,:,:):: geq0,heq0
      real,allocatable,dimension(:,:,:,:):: gp,hp

      real,allocatable,dimension(:,:,:):: rho, rhop
      real,allocatable,dimension(:,:,:):: ux, uy, uz, te, udiv
      real,allocatable,dimension(:,:,:):: uxp, uyp, uzp, tep
      real,allocatable,dimension(:,:,:):: Jx,Jy,Jz,E
      real,allocatable,dimension(:,:,:):: ox, oy, oz
      real,allocatable,dimension(:,:,:):: kx, ky, kz, k2

! note that to make use of FFTW library on bluefire, the complex numbers 
! must have a size of at least complex*16, or even higher. For current 
! real*4, double complex is equivalent to complex*16. 
      real,allocatable,dimension(:,:,:):: vx, vy, vz, prs, ddiv
      real,allocatable,dimension(:,:,:):: wx, wy, wz


      character(len=150):: dirgenr, dirinitflow
      character(len=150):: dirdiag, dirstat   
      character(len=150):: dircntdflow
      character(len=150):: dirflowout
      character(len=150):: dirmoviedata


      end module var_inc
