      subroutine saveprerelax
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      fnm = trim(dirinitflow)//'prerelax_01/finit.'                  &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(10, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(10) istep
      write(10) g, rho
      write(10) ux, uy, uz

      close(10)

      end subroutine saveprerelax
!===========================================================================

      subroutine loadprerelax
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      fnm = trim(dirinitflow)//'prerelax_01/finit.'                  &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(10, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(10) istep
      read(10) g, rho
      read(10) ux, uy, uz

      close(10)

      end subroutine loadprerelax
!===========================================================================

      subroutine saveinitflow
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10) 

      fnm = trim(dirinitflow)//'finit.'                                &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(10, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(10) g
      write(10) h

      close(10)

      end subroutine saveinitflow
!===========================================================================

      subroutine saveinitflow2

      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

      fnm = trim(dirinitflow)//'finituvw.'                                &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(10, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      ux = ux / vscale
      uy = uy / vscale
      uz = uz / vscale
      write(10) ux,uy,uz

      close(10)

      end subroutine saveinitflow2

!===========================================================================

      subroutine saveinitflowascii
      use mpi
      use var_inc
      implicit none

      integer ip, ilen,ipp, indy9, indz9
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2

      real, dimension(lx,ny,nz):: rho9
      real, dimension(lx,ly,lz):: rho1,rho2

      character (len = 120):: fnm

      ilen = lx*ly*lz

      DO ipp = 0,npop-1 
      rho2 = g(ipp,:,:,:)


      if(myid == 0)then

        rho9(:,:,1:lz) = g(ipp,:,:,:)

        do ip = 1,nproc-1
          call MPI_RECV(rho1,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
          rho9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = rho1            
        end do

      else
        call MPI_SEND(rho2,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if


      if(myid == 0) write(*,*)'I am in step 2'

      if(myid == 0)then

        istp1 = ipp / 10
        istp2 = mod(ipp,10)

        fnm = trim(dirflowout)//'finit.'                                 &
            //char(istp1 + 48)//char(istp2 + 48)//'.dat'

       open(20, file = trim(fnm), status = 'unknown',                 &
                form = 'formatted')

        write(20,200) rho9
      end if
      
       if(myid == 0) write(*,*)'I am in step 3 ipp=', ipp

      END DO

200   format(2x,8(1pe22.13))

      end subroutine saveinitflowascii


!===========================================================================
      subroutine loadinitflow
      use var_inc
      implicit none
      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm
      real, allocatable, dimension(:,:,:,:):: g9
      real, allocatable, dimension(:,:,:,:):: h9
      allocate (g9(0:npop-1,lx,ly,lz))
      allocate (h9(0:npop-1,lx,ly,lz))

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10)

! binary format
      fnm = trim(dirinitflow)//'finit.'                                &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(10, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(10) g9
      read(10) h9
      close(10)
      
      g=g9
      h=h9 

! output double precision f9 to ascii format
!      fnm = trim(dirinitflow)//'real4_02/finit.'                       &
!          //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

!      open(10, file = trim(fnm), status = 'unknown',                   &
!               form = 'formatted')
!      write(10,105) f9 

!      close(10)

      deallocate (g9)
      deallocate (h9)


!      istat = 1  
! read in single precision f in ascii format
!      fnm = trim(dirinitflow)//'real4_02/finit.'                       &
!          //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

!       open(10, file = trim(fnm), status = 'unknown',                   &
!                form = 'formatted')
!       read(10,105) f

!       close(10)

! output single precision f in ascii format for checking
!       fnm = trim(dirinitflow)//'real4_02chk/finit.'                    &
!          //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)
 
!       open(10, file = trim(fnm), status = 'unknown',                   &
!               form = 'formatted')
!       write(10,105) f

!       close(10) 
     
      end subroutine loadinitflow
!===========================================================================

      subroutine savecntdflow
      use var_inc
      implicit none
      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5,istp6,istp7
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10) 

      istep = istep0+nsteps
      istp1 = istep / 1000000
      istp2 = mod(istep,1000000) / 100000
      istp3 = mod(istep,100000) / 10000
      istp4 = mod(istep,10000) / 1000
      istp5 = mod(istep,1000) / 100
      istp6 = mod(istep,100) / 10
      istp7 = mod(istep,10)

      fnm = trim(dircntdflow)//'endrunflow2D16x8.'                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6+48)//char(istp7 + 48)&
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(12, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

!      write(12) istep, istat, imovie 
      write(12) g
      write(12) h

      close(12)
      end subroutine savecntdflow      
!!===========================================================================

      SUBROUTINE input_outputf(idirec)

      use mpi
      use var_inc
      IMPLICIT NONE

      integer idump,idirec
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      IF(myid.eq.0) THEN
        if(idirec.eq.1) then
          idump = istpload
        elseif(idirec.eq.2) then
          idump = istep
        endif

        istp1 = idump / 100000
        istp2 = mod(idump,100000) / 10000
        istp3 = mod(idump,10000) / 1000
        istp4 = mod(idump,1000) / 100
        istp5 = mod(idump,100) / 10
        istp6 = mod(idump,10)

        fnm = trim(dircntdflow)//'force.'                              &
              //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)   &
              //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)

        if(idirec.eq.1) then
          open(10,file=trim(fnm),status='unknown',form='unformatted')
          read(10)iseedf
          read(10)ivf
          read(10)iyf
          read(10)b1r,b2r,b3r
          close(10)
        elseif(idirec.eq.2) then
          open(10,file=trim(fnm),status='unknown',form='unformatted')
          write(10)iseedf
          write(10)ivf
          write(10)iyf
          write(10)b1r,b2r,b3r
          close(10)
        endif
      ENDIF

      if(idirec.eq.1) then
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST (b1r,150,MPI_REAL8,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (b2r,150,MPI_REAL8,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (b3r,150,MPI_REAL8,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (iseedf,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (ivf,NTAB,MPI_INTEGER,0,mpi_comm_world,ierr)
        CALL MPI_BCAST (iyf,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      endif


      RETURN
      END SUBROUTINE input_outputf

!===========================================================================
      subroutine savecntdflowHR
      use var_inc
      implicit none
      integer iprc1, iprc2, iprc3, lyh2,jh,j,myida,myidb
      integer istp1, istp2, istp3, istp4, istp5,istp6,istp7
      character (len = 100):: fnm
      real, dimension(0:npop-1,lx,ly,lz):: f11,f22
! interpolation
      
      f11(0:npop-1,lx,1,lz) = 2.5/2.*g(0:npop-1,lx,1,lz) - 0.5/2.0*g(0:npop-1,lx,2,lz)
      do j=2,ly-1,2
      jh = j/2
      f11(0:npop-1,lx,j,lz) = 1.5/2.*g(0:npop-1,lx,jh,lz) + 0.5/2.0*g(0:npop-1,lx,jh+1,lz)
      f11(0:npop-1,lx,j+1,lz) = 0.5/2.*g(0:npop-1,lx,jh,lz) + 1.5/2.0*g(0:npop-1,lx,jh+1,lz)
      end do
      lyh2 = ly/2
      f11(0:npop-1,lx,ly,lz) = 1.5/2.*g(0:npop-1,lx,lyh2,lz) + 0.5/2.0*g(0:npop-1,lx,lyh2+1,lz)

      f22(0:npop-1,lx,1,lz) = 0.5/2.*g(0:npop-1,lx,lyh2,lz) + 1.5/2.0*g(0:npop-1,lx,lyh2+1,lz)
      do j=2,ly-1,2
      jh = j/2 + lyh2
      f22(0:npop-1,lx,j,lz) = 1.5/2.*g(0:npop-1,lx,jh,lz) + 0.5/2.0*g(0:npop-1,lx,jh+1,lz)
      f22(0:npop-1,lx,j+1,lz) = 0.5/2.*g(0:npop-1,lx,jh,lz) + 1.5/2.0*g(0:npop-1,lx,jh+1,lz)
      end do
      f22(0:npop-1,lx,ly,lz) = -0.5/2.*g(0:npop-1,lx,ly-1,lz) + 2.5/2.0*g(0:npop-1,lx,ly,lz)

      istep = istpload
      myida = 2*myid
      myidb = 2*myid + 1
      iprc1 = myida / 100
      iprc2 = mod(myida,100) / 10
      iprc3 = mod(myida,10)

      istp1 = istep / 1000000
      istp2 = mod(istep,1000000) / 100000
      istp3 = mod(istep,100000) / 10000
      istp4 = mod(istep,10000) / 1000
      istp5 = mod(istep,1000) / 100
      istp6 = mod(istep,100) / 10
      istp7 = mod(istep,10)

      fnm = trim(dircntdflow)//'endrunflow2D16x8.'                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6+48)//char(istp7 + 48)&
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(12, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(12) istep, istat, imovie
      write(12) f11

      close(12)

      iprc1 = myidb / 100
      iprc2 = mod(myidb,100) / 10
      iprc3 = mod(myidb,10)

      fnm = trim(dircntdflow)//'endrunflow2D16x8.'                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6+48)//char(istp7 + 48)&
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(12, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      write(12) istep, istat, imovie
      write(12) f22

      close(12)

      end subroutine savecntdflowHR

!===========================================================================
      
      subroutine loadcntdflow
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1,istp2,istp3,istp4,istp5,istp6,istp7
      character (len = 100):: fnm

      iprc1 = myid / 100
      iprc2 = mod(myid,100) / 10
      iprc3 = mod(myid,10) 

      istp1 = istpload / 1000000
      istp2 = mod(istpload,1000000) / 100000
      istp3 = mod(istpload,100000) / 10000
      istp4 = mod(istpload,10000) / 1000
      istp5 = mod(istpload,1000) / 100
      istp6 = mod(istpload,100) / 10
      istp7 = mod(istpload,10)

      fnm = trim(dircntdflow)//'endrunflow2D16x8.'                 &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)  &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)  &
            //char(istp7 + 48)//'.'                                 &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(12, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

!      read(12) istep0, istat, imovie 
      read(12) g
      read(12) h

      close(12)

      end subroutine loadcntdflow      
!===========================================================================
! load data from more number of processes
! e.g., (#_proc_current) = (1/2, 1/4, ..., etc.)*(#_proc_previous)
      subroutine loadcntdflow_frmmore
      use mpi
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm
      
      integer ii, myiid, lz9   
      real, allocatable, dimension(:,:,:,:):: f9

      lz9 = lz / iprocrate

      allocate (f9(0:npop-1,lx,ly,lz9))

      istp1 = istpload / 100000
      istp2 = mod(istpload,100000) / 10000
      istp3 = mod(istpload,10000) / 1000
      istp4 = mod(istpload,1000) / 100
      istp5 = mod(istpload,100) / 10
      istp6 = mod(istpload,10)

      do ii = 0, iprocrate-1 
        myiid = myid*iprocrate + ii 

        iprc1 = myiid / 100
        iprc2 = mod(myiid,100) / 10
        iprc3 = mod(myiid,10)

        fnm = trim(dircntdflow)//'endrunflow.'                         &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

        open(12, file = trim(fnm), status = 'unknown',                 &
               form = 'unformatted')

        read(12) istep0, istat, imovie
        read(12) f9
        close(12)

        g(:,:,:,1+ii*lz9:(ii+1)*lz9) = f9
      end do

      deallocate (f9)

      end subroutine loadcntdflow_frmmore
!===========================================================================
! load data from less number of processes
! e.g., (#_proc_current) = (2, 4, ..., etc.)*(#_proc_previous)
      subroutine loadcntdflow_frmless
      use mpi
      use var_inc
      implicit none

      integer iprc1, iprc2, iprc3
      integer istp1, istp2, istp3, istp4, istp5, istp6
      character (len = 100):: fnm

      integer ii, myiid, lz9  
      real, allocatable, dimension(:,:,:,:):: f9

      myiid = myid / iprocrate
      ii = mod(myid,iprocrate)   
      lz9 = lz * iprocrate 

      allocate (f9(0:npop-1,lx,ly,lz9))

      istp1 = istpload / 100000
      istp2 = mod(istpload,100000) / 10000
      istp3 = mod(istpload,10000) / 1000
      istp4 = mod(istpload,1000) / 100
      istp5 = mod(istpload,100) / 10
      istp6 = mod(istpload,10)

      iprc1 = myiid / 100
      iprc2 = mod(myiid,100) / 10
      iprc3 = mod(myiid,10)

      fnm = trim(dircntdflow)//'endrunflow.'                           &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //'.'                                                      &
            //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

      open(12, file = trim(fnm), status = 'unknown',                   &
               form = 'unformatted')

      read(12) istep0, istat, imovie
      read(12) f9
      close(12)

      g = f9(:,:,:,1+ii*lz:(ii+1)*lz)

      deallocate (f9)

      end subroutine loadcntdflow_frmless   
!===========================================================================
      subroutine outputflow   
      use var_inc
      implicit none

      call outputMa

      call outputux

      call outputuy

      call outputuz

      call outputte   

      call outputrho

      call outputheatflux

!      call outputvort

      end subroutine outputflow   
!===========================================================================
!theta=du_i/d_i
!in statistc
      subroutine outputdiv
      use mpi
      use var_inc
      implicit none

      integer ip, ilen, indy9, indz9
      integer i,j,k
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6, istp7

      real, dimension(lx,ny,nz):: ux9
      real, dimension(lx,ly,lz):: ux0

      character (len = 120):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

       ux9(:,1:ly,1:lz) = udiv

       do ip = 1,nproc-1
         call MPI_RECV(ux0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
      ux9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = ux0          
      end do

      else
       call MPI_SEND(udiv,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

! zero out the velocity inside particle for plot purpose
!        where(ibnodes9 > 0) ux9 = 0.0

        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'div'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(16, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(16,160) ux9

        close(16)

      end if

160   format(2x,8(1pe16.6))

      end subroutine outputdiv
!===========================================================================
! This has been modified for 2D DD
      subroutine outputux
      use mpi
      use var_inc
      implicit none

      integer ip, ilen, indy9, indz9
      integer i,j,k
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6, istp7

      real, dimension(lx,ny,nz):: ux9
      real, dimension(lx,ly,lz):: ux0

      character (len = 120):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

       ux9(:,1:ly,1:lz) = ux

       do ip = 1,nproc-1
         call MPI_RECV(ux0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
      ux9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = ux0          
      end do

      else
       call MPI_SEND(ux,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

! zero out the velocity inside particle for plot purpose
!        where(ibnodes9 > 0) ux9 = 0.0

        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'ux'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(16, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(16,160) ux9

        close(16)

      end if

160   format(2x,8(1pe16.6))

      end subroutine outputux
!===========================================================================

      subroutine outputuy
      use mpi
      use var_inc
      implicit none

      integer ip, ilen, indy9, indz9
      integer i,j,k
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6, istp7

      real, dimension(lx,ny,nz):: uy9
      real, dimension(lx,ly,lz):: uy0

      character (len = 120):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

        uy9(:,1:ly,1:lz) = uy

        do ip = 1,nproc-1
          call MPI_RECV(uy0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
          uy9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = uy0          
        end do

      else
        call MPI_SEND(uy,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if
      if(myid == 0)then

! zero out the velocity inside particle for plot purpose
!        where(ibnodes9 > 0) uy9 = 0.0

        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'uy'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(17, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(17,170) uy9

        close(17)

      end if

170   format(2x,8(1pe16.6))

      end subroutine outputuy

!===========================================================================

      subroutine outputuz
      use mpi
      use var_inc
      implicit none


      integer ip, ilen, indy9, indz9
      integer i,j,k
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6, istp7

      real, dimension(lx,ny,nz):: uz9
      real, dimension(lx,ly,lz):: uz0

      character (len = 120):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

        uz9(:,1:ly,1:lz) = uz

        do ip = 1,nproc-1
          call MPI_RECV(uz0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
          uz9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = uz0          
        end do

      else
        call MPI_SEND(uz,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then
! zero out the velocity inside particle for plot purpose
!        where(ibnodes9 > 0) uz9 = 0.0

        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'uz'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(18, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(18,180) uz9

        close(18)

      end if

180   format(2x,8(1pe16.6))

      end subroutine outputuz
!===========================================================================

      subroutine outputheatflux
      use mpi
      use var_inc
      implicit none


      integer ip, ilen, indy9, indz9
      integer ix,iy,iz
      integer i,j,k
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6, istp7
      real T_factor,mu,tau,cux,cuy,cuz,Qgh

      real, dimension(lx,ly,lz):: qx,qy,qz
      real, dimension(lx,ny,nz):: uz9
      real, dimension(lx,ly,lz):: uz0

      character (len = 120):: fnm

! heat flux
      do iz=1,lz
      do iy=1,ly
      do ix=1,lx
      T_factor=exp(omega*log(te(ix,iy,iz)/T_ref))
      mu=mu_ref*T_factor
      tau=mu/(rho(ix,iy,iz)*R*te(ix,iy,iz))

      qx(ix,iy,iz)=0.0
      qy(ix,iy,iz)=0.0
      qz(ix,iy,iz)=0.0
           do k=0,npop-1
           cux=cix(k)-ux(ix,iy,iz) 
           cuy=ciy(k)-uy(ix,iy,iz)
           cuz=ciz(k)-uz(ix,iy,iz)
           Qgh=(cux*cux+cuy*cuy+cuz*cuz)*g(k,ix,iy,iz)+h(k,ix,iy,iz)
           qx(ix,iy,iz)=qx(ix,iy,iz)+cux*Qgh
           qy(ix,iy,iz)=qy(ix,iy,iz)+cuy*Qgh
           qz(ix,iy,iz)=qz(ix,iy,iz)+cuz*Qgh
           end do
      qx=(tau*0.5*qx)/(tau+0.5*dt*Pr)
      qy=(tau*0.5*qy)/(tau+0.5*dt*Pr)
      qz=(tau*0.5*qz)/(tau+0.5*dt*Pr)
      end do
      end do
      end do

      ilen = lx*ly*lz

! output  qx
      IF(myid == 0)then
        uz9(:,1:ly,1:lz) = qx
        do ip = 1,nproc-1
          call MPI_RECV(uz0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
          uz9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = uz0          
        end do
      ELSE
        call MPI_SEND(qx,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      END IF

      IF(myid == 0)then
        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'qx'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(18, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(18,180) uz9

        close(18)

      END IF

! output  qy
      IF(myid == 0)then
        uz9(:,1:ly,1:lz) = qy
        do ip = 1,nproc-1
          call MPI_RECV(uz0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
          uz9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = uz0
        end do
      ELSE
        call MPI_SEND(qy,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      END IF

      IF(myid == 0)then
        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'qy'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(18, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(18,180) uz9

        close(18)

      END IF

! output qz 
      IF(myid == 0)then
        uz9(:,1:ly,1:lz) = qz
        do ip = 1,nproc-1
          call MPI_RECV(uz0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
          uz9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = uz0
        end do
      ELSE
        call MPI_SEND(qz,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      END IF

      IF(myid == 0)then
        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'qz'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(18, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(18,180) uz9

        close(18)

      END IF

180   format(2x,8(1pe16.6))

      end subroutine outputheatflux

!===========================================================================
      subroutine outputte
      use mpi
      use var_inc
      implicit none

      integer ip, ilen, indy9, indz9
      integer i,j,k
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6, istp7

      real, dimension(lx,ny,nz):: te9
      real, dimension(lx,ly,lz):: te0

      character (len = 120):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

       te9(:,1:ly,1:lz) = te

       do ip = 1,nproc-1
         call MPI_RECV(te0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
      te9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = te0          
      end do

      else
       call MPI_SEND(te,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'te'                                   &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(166, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(166,160) te9

        close(166)

      end if

160   format(2x,8(1pe16.6))



      end subroutine outputte
!===========================================================================
      subroutine outputrho
      use mpi
      use var_inc
      implicit none

!     integer, dimension(lx,ly,nz):: ibnodes9

      integer ip, ilen, indy9, indz9
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6, istp7

      real, dimension(lx,ny,nz):: rho9
      real, dimension(lx,ly,lz):: rho1

      character (len = 120):: fnm

      ilen = lx*ly*lz

      if(myid == 0)then

        rho9(:,1:ly,1:lz) = rho

        do ip = 1,nproc-1
          call MPI_RECV(rho1,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
          rho9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = rho1            
        end do

      else
        call MPI_SEND(rho,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

! zero out the density inside particle for plot purpose
!        where(ibnodes9 > 0) rho9 = 0.0

        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'rho'                                &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(19, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(19,190) rho9

        close(19)

      end if

190   format(2x,8(1pe16.6))

      end subroutine outputrho


!===========================================================================
      subroutine outputMa !local Mach number
      use mpi 
      use var_inc
      implicit none

      integer ip, ilen, indy9, indz9
      integer status(MPI_STATUS_SIZE)
      integer istp1, istp2, istp3, istp4, istp5, istp6, istp7
      
      real, dimension(lx,ny,nz):: vort9    
      real, dimension(lx,ly,lz):: vort, vort0    

      character (len = 100):: fnm

      vort = sqrt( (ux*ux+uy*uy+uz*uz)/(gam*R*Te ))
   
      ilen = lx*ly*lz

      if(myid == 0)then

        vort9(:,1:ly,1:lz) = vort         

       do ip = 1,nproc-1
         call MPI_RECV(vort0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

      indy9 = mod(ip,nprocY)
      indz9 = int(ip/nprocY)
      vort9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = vort0          
      end do

      else
       call MPI_SEND(vort,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      if(myid == 0)then

! zero out the vorticity inside particle for plot purpose
!        where(ibnodes9 > 0) vort9 = 0.0

        istp1 = istep / 1000000
        istp2 = mod(istep,1000000) / 100000
        istp3 = mod(istep,100000) / 10000
        istp4 = mod(istep,10000) / 1000
        istp5 = mod(istep,1000) / 100
        istp6 = mod(istep,100) / 10
        istp7 = mod(istep,10)

        fnm = trim(dirflowout)//'Ma'                                 &
            //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)     &
            //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)     &
            //char(istp7 + 48)//'.dat'

        open(20, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted')

        write(20,200) vort9 

        close(20)

      end if

200   format(2x,8(1pe16.6))

      end subroutine outputMa      
!===========================================================================
      subroutine outputvort
      use mpi
      use var_inc
      implicit none

! The vorticity calculation is done by each processor, to output, we
! need to sum all the information in processor0

! Here I just take the derivative to calculate the vorticity, when
! better method is found, I will improve this part

! As Dr. Wang suggested, here use allocatable arrays

      real, allocatable, dimension(:,:,:):: ox9, oy9, oz9, vort9
      real, allocatable, dimension(:,:,:):: ox0, oy0, oz0, vort0, vort
      integer, dimension(lx,ny,nz):: ibnodes9
      integer ip, ilen, indy9, indz9
      integer istp1, istp2, istp3, istp4, istp5, istp6, istp7
      integer status(MPI_STATUS_SIZE)
      character (len = 120):: fnm

      ilen = lx*ly*lz

       call vortcalc

      allocate (ox9(lx,ny,nz))
      allocate (oy9(lx,ny,nz))
      allocate (oz9(lx,ny,nz))
      allocate (vort9(lx,ny,nz))
      allocate (ox0(lx,ly,lz))
      allocate (oy0(lx,ly,lz))
      allocate (oz0(lx,ly,lz))
      allocate (vort0(lx,ly,lz))
      allocate (vort(lx,ly,lz))


      vort = sqrt(ox*ox + oy*oy + oz*oz)

      if(myid == 0 ) then

      ox9(:,1:ly,1:lz) = ox
      oy9(:,1:ly,1:lz) = oy
      oz9(:,1:ly,1:lz) = oz
      vort(:,1:ly,1:lz) = vort

      do ip = 1,nproc - 1

       call MPI_RECV(ox0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)
       call MPI_RECV(oy0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)
       call MPI_RECV(oz0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)
       call MPI_RECV(vort0,ilen,MPI_REAL8,ip,1,MPI_COMM_WORLD,status,ierr)

       indy9 = mod(ip,nprocY)
       indz9 = int(ip/nprocY)

       ox9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = ox0
       oy9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = oy0
       oz9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = oz0
       vort9(:,indy9*ly+1:indy9*ly+ly,indz9*lz+1:indz9*lz+lz) = vort0
      end do

      else
      call MPI_SEND(ox,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      call MPI_SEND(oy,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      call MPI_SEND(oz,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      call MPI_SEND(vort,ilen,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
      end if

      deallocate (ox0)
      deallocate (oy0)
      deallocate (oz0)
      deallocate (vort0)
      deallocate (vort)
! we can output the vorticity in each direction or the magnitude, here I
! just output the magnitude

! Again the writing is processed by processor0

      if (myid == 0 ) then

       istp1 = istep / 1000000
       istp2 = mod(istep,1000000) / 100000
       istp3 = mod(istep,100000) / 10000
       istp4 = mod(istep,10000) / 1000
       istp5 = mod(istep,1000) / 100
       istp6 = mod(istep,100) / 10
       istp7 = mod(istep,10)

       fnm = trim(dirflowout)//'oz'                               &
           //char(istp1 + 48)//char(istp2 + 48)//char(istp3 + 48)   &
           //char(istp4 + 48)//char(istp5 + 48)//char(istp6 + 48)   &
           //char(istp7 + 48)//'.dat'

       open(20, file = trim(fnm), status = 'unknown', form = 'formatted') 
       write(20,270) oz9

       close(20)
       end if

       deallocate (ox9)
       deallocate (oy9)
       deallocate (oz9)
       deallocate (vort9)

270    format(2x,8(1pe16.6))

      end subroutine outputvort
!===========================================================================
! to compute kinetic energy spectrum, dissipation rate spectrum
! also the skewness and flatness

      subroutine statistc
      use mpi
      use var_inc
      implicit none

      integer ik, nxyz,k,kglb,idZ
      character (len = 100):: fnm1, fnm2, fnm3, fnm4
      real, dimension(lx+2,lly,lz) :: tmph,tmp1h, tmpp, tmpd
      real ek,e_t,ekp,ep_t,eks, ekd, ed_t, e_ts, dissp, qq, qqp
      real qqqp, qqq, dissppp, eta, vk, tk, uprm
      real tmse, Re, xl, et, kmxeta, xintls, vskew, vflat
      real cc2, cc2t, cc3, cc3t, cc4, cc4t, cc5, cc5t

      REAL,DIMENSION     (lx,ly)  ::  tmp2D
      REAL,DIMENSION  (nz) :: vxave,vyave,vzave,vxsq,vysq,vzsq,stress_xz,stress_xy,stress_yz
      REAL,DIMENSION  (nz) :: vxavet,vyavet,vzavet,vxsqt,vysqt,vzsqt,stress_xzt,stress_xyt,stress_yzt

      vx = 0.0
      vy = 0.0
      vz = 0.0

      vx(1:lx,1:ly,:) = ux
      vy(1:lx,1:ly,:) = uy
      vz(1:lx,1:ly,:) = uz
      prs(1:lx,1:ly,:) = rho*R*te

      write(*,*)'3rcx'
      call mpifft3DRC(vx)
      call mpifft3DRC(vy)
      call mpifft3DRC(vz)
      call mpifft3DRC(prs)
      write(*,*)'3rcx done'


      if(myid == 0)then
        qqp = 0.0
        qq = 0.0
        dissp = 0.0
        xintls = 0.0

        fnm1 = trim(dirstat)//'spectrum.dat'
        fnm2 = trim(dirstat)//'monitora.dat'
        fnm3 = trim(dirstat)//'monitorb.dat'
        fnm4 = trim(dirstat)//'profiles.dat'

        open(24, file = trim(fnm1), status = 'unknown',                &
                 form = 'formatted', position = 'append')
        open(25, file = trim(fnm2), status = 'unknown',                &
                 form = 'formatted', position = 'append')
        open(26, file = trim(fnm3), status = 'unknown',                &
                 form = 'formatted', position = 'append')
        open(27, file = trim(fnm4), status = 'unknown',                &
                 form = 'formatted', position = 'append')

      end if

! skewness and flatness
      wx = kx * vx
      wy = ky * vy
      wz = kz * vz
      tmph = wx
      wx(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wx(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
      tmph = wy
      wy(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wy(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)
      tmph = wz
      wz(1:lx , 1:lly-1:2 , :) = -tmph(1:lx , 2:lly:2   , :)
      wz(1:lx , 2:lly:2   , :) =  tmph(1:lx , 1:lly-1:2 , :)

      call mpifft3DCR(wx)
      call mpifft3DCR(wy)
      call mpifft3DCR(wz)

! Method 1
      udiv(1:lx,1:ly,1:lz) = wx(1:lx,1:ly,1:lz) + wy(1:lx,1:ly,1:lz) &
                          + wz(1:lx,1:ly,1:lz)
      ddiv(1:lx,1:ly,:) = udiv

      call mpifft3DRC(ddiv)

!      if( mod(istep,nflowout) .eq. 0 ) call outputdiv

      tmph = wx**2 + wy**2 + wz**2
      tmph = tmph / 3.0

      cc2 = sum(tmph(1:lx,1:ly,:))
      call MPI_REDUCE(cc2,cc2t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      tmph = wx**3 + wy**3 + wz**3
      tmph = tmph / 3.0

      cc3 = sum(tmph(1:lx,1:ly,:))
      call MPI_REDUCE(cc3,cc3t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      tmph = wx**4 + wy**4 + wz**4
      tmph = tmph / 3.0

      cc4 = sum(tmph(1:lx,1:ly,:))
      call MPI_REDUCE(cc4,cc4t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      cc5 = sum(udiv(1:lx,1:ly,:)*udiv(1:lx,1:ly,:))
      call MPI_REDUCE(cc5,cc5t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid == 0)then
        nxyz = nx*ny*nz
        vskew = cc3t / real(nxyz) / (cc2t/real(nxyz))**1.5
        vflat = cc4t / real(nxyz) / (cc2t/real(nxyz))**2
        cc5t = sqrt(cc5t / real(nxyz))
      end if

! Profiles
       IF ( mod(istep,nspec).eq.0 ) THEN
        vxave = 0.0
        vyave = 0.0
        vzave = 0.0
        vxsq = 0.0
        vysq = 0.0
        vzsq = 0.0
        stress_xz = 0.0
        stress_xy = 0.0
        stress_yz = 0.0

       idz = int(myid/nprocY)

       do k=1,lz
       kglb = lz*idz+k
       tmp2D = ux(:,:,k)
       vxave(kglb) = sum (tmp2D(:,:) )
       tmp2D = uy(:,:,k)
       vyave(kglb) = sum (tmp2D(:,:) )
       tmp2D = uz(:,:,k)
       vzave(kglb) = sum (tmp2D(:,:) )

       tmp2D = ux(:,:,k)*uz(:,:,k)
       stress_xz(kglb) = sum ( tmp2D(:,:) )
       tmp2D = ux(:,:,k)*uy(:,:,k)
       stress_xy(kglb) = sum ( tmp2D(:,:) )
       tmp2D = uy(:,:,k)*uz(:,:,k)
       stress_yz(kglb) = sum ( tmp2D(:,:) )

       tmp2D = (ux(:,:,k))**2
       vxsq(kglb) = sum ( tmp2D(:,:) )
       tmp2D = (uy(:,:,k))**2
       vysq(kglb) = sum ( tmp2D(:,:) )
       tmp2D = (uz(:,:,k))**2
       vzsq(kglb) = sum ( tmp2D(:,:) )

       end do

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       CALL MPI_ALLREDUCE(vxave,vxavet,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE(vyave,vyavet,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE(vzave,vzavet,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE(vxsq,vxsqt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE(vysq,vysqt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE(vzsq,vzsqt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE(stress_xz,stress_xzt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE(stress_xy,stress_xyt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE(stress_yz,stress_yzt,nz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

       if (myid.eq.0) then
       vxavet = vxavet / (nx*ny) / vscale
       vyavet = vyavet / (nx*ny) /vscale
       vzavet = vzavet / (nx*ny) /vscale
       vxsqt = vxsqt / (nx*ny) /vscale**2
       vysqt = vysqt / (nx*ny) /vscale**2
       vzsqt = vzsqt / (nx*ny) /vscale**2
       stress_xzt = stress_xzt / (nx*ny) /vscale**2
       stress_xyt = stress_xyt / (nx*ny) /vscale**2
       stress_yzt = stress_yzt / (nx*ny) /vscale**2

       do k=1,nz
       write(27,460)k-1,vxavet(k),vyavet(k),vzavet(k), &
             vxsqt(k),vysqt(k),vzsqt(k), stress_xzt(k),stress_xyt(k),stress_yzt(k)
       end do
460    format(2x,i5,9(1pe15.6))
        close(27)
       end if

       ENDIF
! kinetic energy and dissipation rate
! add pressure
! add divergence
      tmph = vx*vx + vy*vy + vz*vz
      tmpp = 2.*prs*prs
      tmpd = 2.*ddiv*ddiv
      if(indy.eq.0) then
        tmph(:,1,:) = 0.5 * tmph(:,1,:) ! note: this is already1/2*(u')^2 = tke
        tmph(:,2,:) = 0.5 * tmph(:,2,:)
        tmpp(:,1,:) = 0.5 * tmpp(:,1,:) ! note: this is already1/2*(u')^2 = tke
        tmpp(:,2,:) = 0.5 * tmpp(:,2,:)
        tmpd(:,1,:) = 0.5 * tmpd(:,1,:) ! note: this is already1/2*(u')^2 = tke
        tmpd(:,2,:) = 0.5 * tmpd(:,2,:)
      endif
      tmp1h = 2. * visc * tmph * k2

      do ik = 1,nek
        ekp = 0.0
        ekd = 0.0
        ek = 0.0
        e_t = 0.0
        ed_t=0.0
        ep_t = 0.0
        ek = sum(tmph(1:lx,:,:), mask=(ik2(1:lx,:,:) == ik))
        ekp = sum(tmpp(1:lx,:,:), mask=(ik2(1:lx,:,:) == ik))
        ekd = sum(tmpd(1:lx,:,:), mask=(ik2(1:lx,:,:) == ik))

        call MPI_REDUCE(ek,e_t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ekp,ep_t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ekd,ed_t,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        eks = 0.0
        e_ts = 0.0
        eks = sum(tmp1h(1:lx,:,:), mask=(ik2(1:lx,:,:) == ik))

        call MPI_REDUCE(eks,e_ts,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(myid == 0 )then
          qq = qq + e_t
          qqp = qqp + ep_t
          dissp = dissp + e_ts
          xintls = xintls + e_t / real(ik)
          if( mod(istep,nspec)==0) write(24,240) ik, e_t*escale, e_ts*dscale, &
               ep_t*escale, ed_t
        end if
      end do

      if(myid == 0)then
        qqqp = qqp*escale
        qqq = qq*escale
        xintls = xintls*escale
        dissppp = dissp*dscale

        eta    = (visc**3 / dissppp)**0.25 ! Kolmogorov lengh scale
        vk     = (visc * dissppp)**0.25 ! Kolmogorov velocity scale
        tk     = (visc / dissppp)**0.5 ! Kolmogorov time scale
        uprm   = sqrt(2.0 * qqq / 3.0) ! u_prime = u_rms
        tmse   = sqrt(15.0 * visc * uprm**2 / dissppp) ! Taylor microscale
        Re     = uprm * tmse / visc ! Taylor microscale Reynolds number
        xl     = uprm**3 / dissppp ! large eddy lengthscale
        et     = xl / uprm ! large eddy turn-over time
        kmxeta = (real(lxh) - 1.5) * eta ! kmax * Kolmogorov lengthscale
        xintls = xintls * pi / 2.0 / uprm**2 !longitudinal integral lengthscale

        write(25,250) ttt, eta, tmse, xintls, vk, tk, et, xl
        write(26,250) ttt, uprm, qqq, dissppp,      &
                      Re, kmxeta, vskew, vflat, qqp, cc5t
        close(24)
        close(25)
        close(26)
      end if

      istat = istat + 1

240   format(2x,i5,4(E16.6e4))
250   format(2x,15(E16.6e4))

      end subroutine statistc
!===========================================================================

      subroutine gcheck
      use mpi
      use var_inc
      implicit none
      character (len = 120):: fnm
      real a1,a2,a3,Masum,Masumt
      real ama_max,ama_local,vel,velt,csl,cst,turbmt
      real tmassint,tkeint,tengint,tmasst,tket,tengt, amat
      integer ix,iy,iz,k

        ama_max=0.0
        vel = 0.0
        csl = 0.0
        Masum=0.0
        velt=0.0
        cst=0.0
        amat=0.0
        Masumt=0.0

        do iz=1,lz
        do iy=1,ly
        do ix=1,lx
        vel = vel + ux(ix,iy,iz)*ux(ix,iy,iz) + uy(ix,iy,iz)*uy(ix,iy,iz) &
                         + uz(ix,iy,iz)*uz(ix,iy,iz)
        csl = csl + sqrt(gam*R*te(ix,iy,iz))
        ama_local = ux(ix,iy,iz)*ux(ix,iy,iz) + uy(ix,iy,iz)*uy(ix,iy,iz) &
                         + uz(ix,iy,iz)*uz(ix,iy,iz)
        ama_local = sqrt(ama_local/(gam*R*te(ix,iy,iz)))
        Masum=Masum+ama_local**2
        if(ama_local.gt.ama_max)ama_max = ama_local
        end do
        end do
        end do

      CALL MPI_REDUCE(vel,velt,1,MPI_REAL8,MPI_SUM,0, &
                                   mpi_comm_world,ierr)
      CALL MPI_REDUCE(Masum,Masumt,1,MPI_REAL8,MPI_SUM,0, &
                                   mpi_comm_world,ierr)
      CALL MPI_REDUCE(csl,cst,1,MPI_REAL8,MPI_SUM,0, &
                                   mpi_comm_world,ierr)
      CALL MPI_REDUCE(ama_max,amat,1,MPI_REAL8,MPI_MAX,0, &
                                   mpi_comm_world,ierr)


!time , Ma_localmax, Ma_turb, <u^2>^(1/2),<(gamma*R*T)^(1/2)>,<Ma_local^2>^0.5
!====================================================================
      if(myid == 0)then
        turbmt = sqrt(velt/real(nx*ny*nz))/(cst/real(nx*ny*nz))!Mt
        a1=sqrt(velt/real(nx*ny*nz))!<u^2>^(1/2)
        a2=cst/real(nx*ny*nz)!<(gamma*R*T)^(1/2)>
        a3=sqrt(Masumt/real(nx*ny*nz))!<Ma_local^2>^0.5
        fnm = trim(dirstat)//'Mach_turb.dat'
        open(26, file = trim(fnm), status = 'unknown',                &
                 form = 'formatted', position = 'append')
       write(26,460)ttt,amat,turbmt,a1,a2,a3
       close(26)
      end if
460    format(2x,6(1pe15.5))


      end subroutine gcheck

!===========================================================================
! this is to monitor the mean and maximum flow velocity and particle
! velocity
      subroutine diag
      use mpi
      use var_inc
      implicit none

!********THIS IS CHANGED*******************************
!      integer, dimension(2):: idwp, idomgp  

      integer, dimension(1):: idwp, idomgp
      integer kg,jg

      real wpmax, omgpmax
      real, dimension(lx,ly,lz):: vel
      real, dimension(npart):: wpmag, omgpmag
      real, dimension(nproc):: vmax0,vmax0t
      integer, dimension(nproc):: im,jm,km,im0,jm0,km0
      real umean,vmean,wmean,pmean,prms,tmean,trms
      real flux,fluxt
      real pmeant,prmst,tmeant,trmst
      real umeant,vmeant,wmeant,umeantn,vmeantn,wmeantn
      real tnuba,tnuta
      REAL,DIMENSION     (ny,nz)  ::  nub,nut
      REAL,DIMENSION     (ny,nz)  ::  tnub,tnut
      integer i,j,k,imout,jmout,kmout

      character (len = 100):: fnm
!!!!!!!!!!!
      umean = 0.0
      vmean = 0.0
      wmean = 0.0
      pmean = 0.0
      tmean = 0.0
      prms = 0.0
      trms = 0.0
      flux = 0.0

      do k=1,lz
      do j=1,ly
      do i=1,lx
      flux = flux + ux(i,j,k)*te(i,j,k)*dx(i)*dy(j)*dz(k)
      umean = umean + ux(i,j,k)**2*dx(i)*dy(j)*dz(k)
      vmean = vmean + uy(i,j,k)**2*dx(i)*dy(j)*dz(k)
      wmean = wmean + uz(i,j,k)**2*dx(i)*dy(j)*dz(k)
      pmean = pmean + rho(i,j,k)*R*Te(i,j,k)*dx(i)*dy(j)*dz(k)
      prms = prms + (rho(i,j,k)*R*Te(i,j,k))**2*dx(i)*dy(j)*dz(k)
      tmean = tmean + te(i,j,k)*dx(i)*dy(j)*dz(k)
      trms = trms + (te(i,j,k))**2*dx(i)*dy(j)*dz(k)
      
      end do
      end do
      end do


       nub = 0.0
       nut = 0.0

       do k=1,lz
       kg = k+indz*lz
       do j=1,ly
       jg = j+indy*ly
        nub(jg,kg) = nx*(Tc - Te(1,j,k) )/(0.5*dx(1))
        nut(jg,kg) = nx*(Te(lx,j,k)-Th)/(0.5*dx(lx))
       end do
       end do


      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      CALL MPI_ALLREDUCE(flux,fluxt,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(umean,umeant,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(vmean,vmeant,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(wmean,wmeant,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(pmean,pmeant,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(tmean,tmeant,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(prms,prmst,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      CALL MPI_ALLREDUCE(trms,trmst,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
!
       CALL MPI_ALLREDUCE (nub,tnub,nyz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_ALLREDUCE (nut,tnut,nyz,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

      vel = sqrt(ux*ux + uy*uy + uz*uz)

!     vmax0 = maxval(vel) 
      im = 0
      jm = 0
      km = 0
      vmax0 = 0.0
      do k=1,lz
      do j=1,ly
      do i=1,lx
      if(vel(i,j,k).gt.vmax0(myid) ) then
      vmax0(myid) = vel(i,j,k)
      im(myid)=i
      jm(myid)=j +  indy*ly
      km(myid)=k +  indz*lz
      end if
      end do
      end do
      end do
!
!     call
!     MPI_REDUCE(vmax0,vmax,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(im,im0,nproc,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(jm,jm0,nproc,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(km,km0,nproc,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(vmax0,vmax0t,nproc,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!     if(myid == 0)then
!       wpmax = 0.0
!       omgpmax = 0.0
!       idwp = 0
!       idomgp = 0 

!       if(ipart)then
!         wpmag = sqrt(wp(1,:)**2 + wp(2,:)**2 + wp(3,:)**2)
!         wpmax = maxval(wpmag)
!         idwp = maxloc(wpmag)

!         omgpmag = sqrt(omgp(1,:)**2 + omgp(2,:)**2 + omgp(3,:)**2)
!         omgpmax = maxval(omgpmag)
!         idomgp = maxloc(omgpmag)
!       end if

!     end if

          rhoerr = maxval(rho)
          call MPI_ALLREDUCE(rhoerr,rhomax,1,MPI_REAL8,MPI_MAX, &
                             MPI_COMM_WORLD,ierr)
          rhoerr = minval(rho)
          call MPI_ALLREDUCE(rhoerr,rhomin,1,MPI_REAL8,MPI_MIN, &
                             MPI_COMM_WORLD,ierr)
          if(myid == 0 ) write(*,*)istep, rhomax, rhomin

      if(myid == 0)then
       fluxt = fluxt/real(nx*ny*nz)
       fluxt = fluxt*real(nx)/(visc/Pr)
       fluxt = fluxt + 1.d0
       umeant = umeant/real(nx*ny*nz)
       vmeant = vmeant/real(nx*ny*nz)
       wmeant = wmeant/real(nx*ny*nz)
       pmeant = pmeant/real(nx*ny*nz)
       tmeant = tmeant/real(nx*ny*nz)
       prmst = prmst/real(nx*ny*nz)
       trmst = trmst/real(nx*ny*nz)

       
       prmst = (prmst - pmeant*pmeant)
       if(prmst.gt.0.0) then
         prmst = sqrt(prmst)
       else
         prmst = 0.0
       end if

       trmst = sqrt(trmst - tmeant*tmeant)

       umeantn = umeant/u0**2
       vmeantn = vmeant/u0**2
       wmeantn = wmeant/u0**2
       pmeant = pmeant/u0**2
       prmst = prmst/u0**2

        fnm = trim(dirdiag)//'diag.dat'

        open(26, file = trim(fnm), status = 'unknown',                 &
                 form = 'formatted', position = 'append')

        vmax = 0.0
        imout = 0
        jmout = 0
        kmout = 0
        do i = 1, nproc
        if(vmax0t(i).gt.vmax)then
        vmax = vmax0t(i)
        imout = im0(i)
        jmout = jm0(i)
        kmout = km0(i)
        end if

        end do

       tnuba = 0.0
       tnuta = 0.0
       do k=1,nz
       do j=1,ny
       tnuba = tnuba + tnub(j,k)*dzg(k)*dyg(j)
       tnuta = tnuta + tnut(j,k)*dzg(k)*dyg(j)
       end do
       end do
       tnuba = tnuba/real(nz)/real(ny)
       tnuta = tnuta/real(nz)/real(ny)

       write(26,260) ttt,vmax,imout,jmout,kmout,umeantn,vmeantn,wmeantn,&
        pmeant,prmst,tmeant,trmst,rhomax,rhomin,tnuba,tnuta,fluxt
       write(*,*)'ttt,u^2,v^2,w^2,pmeant,prmst,tmeant,trmst,rhomax,rhomin=', &
       ttt,umeantn,vmeantn,wmeantn,pmeant,prmst,tmeant,trmst,rhomax,rhomin
        close(26)

      end if

260   format(2x, 2(1pe16.6),3I6,12(1pe16.6) ) 
 
      end subroutine diag    

!==========================================================================
      subroutine vortcalc
      use mpi
      use var_inc
      implicit none
! Since the vorticity is only decided by the derivative in x, y and z
! directions, there is no need to consider the conner points
      real, dimension(lx,ly) :: tmpuxF, tmpuyF, tmpuzF
      real, dimension(lx,ly) :: tmpuxB, tmpuyB, tmpuzB
      real, dimension(lx,lz) :: tmpuxL, tmpuyL, tmpuzL
      real, dimension(lx,lz) :: tmpuxR, tmpuyR, tmpuzR
      real, dimension(lx,ly,lz):: pwy,pvz,puz,pwx,pvx,puy
      integer i, j, k, id
      real omg1, omg2, omg3

      integer, dimension(lx,ly):: tmpiF, tmpiB
      integer, dimension(lx,lz):: tmpiL, tmpiR

      call exchng8(ux(:,:,1),tmpuxB,ux(:,:,lz),tmpuxF,ux(:,1,:),tmpuxR &
           ,ux(:,ly,:),tmpuxL)
      call exchng8(uy(:,:,1),tmpuyB,uy(:,:,lz),tmpuyF,uy(:,1,:),tmpuyR &
           ,uy(:,ly,:),tmpuyL)
      call exchng8(uz(:,:,1),tmpuzB,uz(:,:,lz),tmpuzF,uz(:,1,:),tmpuzR &
           ,uz(:,ly,:),tmpuzL)
!      call exchng8i(ibnodes(:,:,1),tmpiB,ibnodes(:,:,lz),tmpiF,
!      &
!           ibnodes(:,1,:),tmpiR,ibnodes(:,ly,:),tmpiL)

      ox = 0.0
      oy = 0.0
      oz = 0.0

      do k = 1,lz
      do j = 1,ly
      do i = 1,lx
! Here note that inside the particle, we also have calculated the
! velocity, so we dont need to particularly deal with the situation that
! some point is inside the particle. The only situation that needs to be
! particularly treated is near the wall 

      if(i == 1) then
       pwx(i,j,k) = 2.*dx(2)/dx(1)/(dx(1)+dx(2))*uz(i,j,k)   &
                 +  dx(1)/(dx(2)+2.*dx(1))/(xc(2)-xc(1))*uz(i+1,j,k)
       pvx(i,j,k) = 2.*dx(2)/dx(1)/(dx(1)+dx(2))*uy(i,j,k)   &
                 +  dx(1)/(dx(2)+2.*dx(1))/(xc(2)-xc(1))*uy(i+1,j,k)
      else if(i == lx) then
       pwx(i,j,k) = - 2.*dx(lx-1)/dx(lx)/(dx(lx)+dx(lx-1))*uz(i,j,k)   &
                -  dx(lx)/(dx(lx-1)+2.*dx(lx))/(xc(lx-1)-xc(lx))*uz(i-1,j,k)
       pvx(i,j,k) = - 2.*dx(lx-1)/dx(lx)/(dx(lx)+dx(lx-1))*uy(i,j,k)   &
                -  dx(lx)/(dx(lx-1)+2.*dx(lx))/(xc(lx-1)-xc(lx))*uy(i-1,j,k)
      else
       pwx(i,j,k) = (uz(i+1,j,k)-uz(i-1,j,k))/(xc(i+1)-xc(i-1))
       pvx(i,j,k) = (uy(i+1,j,k)-uy(i-1,j,k))/(xc(i+1)-xc(i-1))
      end if

      if(j == 1) then
       pwy(i,j,k) = (uz(i,j+1,k)-tmpuzL(i,k))/2.d0/(yc(j+1)-yc(j))
       puy(i,j,k) = (ux(i,j+1,k)-tmpuxL(i,k))/2.d0/(yc(j+1)-yc(j))
      else if(j == ly) then
       pwy(i,j,k) = (tmpuzR(i,k)-uz(i,j-1,k))/2.d0/(yc(j)-yc(j-1))
       puy(i,j,k) = (tmpuxR(i,k)-ux(i,j-1,k))/2.d0/(yc(j)-yc(j-1))
      else
       pwy(i,j,k) = (uz(i,j+1,k)-uz(i,j-1,k))/2.d0/(yc(j+1)-yc(j-1))
       puy(i,j,k) = (ux(i,j+1,k)-ux(i,j-1,k))/2.d0/(yc(j+1)-yc(j-1))
      end if

      if(k == 1) then
       puz(i,j,k) = (ux(i,j,k+1)-tmpuxF(i,j))/2.d0/(zc(k+1)-zc(k))
       pvz(i,j,k) = (uy(i,j,k+1)-tmpuyF(i,j))/2.d0/(zc(k+1)-zc(k))
      else if(k == lz) then
       puz(i,j,k) = (tmpuxB(i,j)-ux(i,j,k-1))/2.d0/(zc(k)-zc(k-1))
       pvz(i,j,k) = (tmpuyB(i,j)-uy(i,j,k-1))/2.d0/(zc(k)-zc(k-1))
      else
       puz(i,j,k) = (ux(i,j,k+1)-ux(i,j,k-1))/2.d0/(zc(k+1)-zc(k-1))
       pvz(i,j,k) = (uy(i,j,k+1)-uy(i,j,k-1))/2.d0/(zc(k+1)-zc(k-1))
      end if


      ox(i,j,k) = pwy(i,j,k) - pvz(i,j,k)
      oy(i,j,k) = puz(i,j,k) - pwx(i,j,k)
      oz(i,j,k) = pvx(i,j,k) - puy(i,j,k)

      end do
      end do
      end do


      end subroutine vortcalc

!===========================================================================

!==========================================================================
! This subroutine is for writing out a slice of the instantaneous flow field
      subroutine writeflowfieldstart
      use mpi
      use var_inc
      implicit none

      real,dimension(ly,lz)::uxplane
      integer i,k,j
      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm1

        iprc1 = myid / 100
        iprc2 = mod(myid,100) / 10
        iprc3 = mod(myid,10)

        fnm1 = trim(dirstat)//'flowfield2Dstart.dat.'         &
              //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

        open(44, file = trim(fnm1), status = 'unknown',                 &
                 form = 'unformatted')

          uxplane = ux(lx/2,:,:)

          write(44) uxplane

       close(44)

      end subroutine writeflowfieldstart
!===========================================================================

!==========================================================================
! This subroutine is for writing out a slice of the instantaneous flow field
      subroutine writeflowfield
      use mpi
      use var_inc
      implicit none

      real,dimension(ly,lz)::uxplane
      integer i,k,j
      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm1

        iprc1 = myid / 100
        iprc2 = mod(myid,100) / 10
        iprc3 = mod(myid,10)

        fnm1 = trim(dirstat)//'flowfield2D.dat.'         &
              //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

        open(66, file = trim(fnm1), status = 'unknown',                 &
                 form = 'unformatted')

          uxplane = ux(lx/2,:,:)

          write(66) uxplane

          close(66)

      end subroutine writeflowfield
!===========================================================================

!==========================================================================
! This subroutine is for writing out a slice of the instantaneous flow field
      subroutine writeflowfieldend
      use mpi
      use var_inc
      implicit none

      real,dimension(ly,lz)::uxplane
      integer i,k,j
      integer iprc1, iprc2, iprc3
      character (len = 100):: fnm1

        iprc1 = myid / 100
        iprc2 = mod(myid,100) / 10
        iprc3 = mod(myid,10)

        fnm1 = trim(dirstat)//'flowfield2Dend.dat.'         &
              //char(iprc1 + 48)//char(iprc2 + 48)//char(iprc3 + 48)

        open(55, file = trim(fnm1), status = 'unknown',                 &
                 form = 'unformatted')

          uxplane = ux(lx/2,:,:)

          write(55) uxplane

          close(55)

      end subroutine writeflowfieldend
!===========================================================================
      subroutine exchng8(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)
      use mpi
      use var_inc
      implicit none
      real, dimension(lx,ly):: tmp1,tmp2,tmp3,tmp4
      real, dimension(lx,lz):: tmp5,tmp6,tmp7,tmp8

      integer ileny, ilenz
      integer error, status_array(MPI_STATUS_SIZE,4),req(4)

      ilenz = lx*ly
      ileny = lx*lz

      call MPI_IRECV(tmp2,ilenz,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp4,ilenz,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp1,ilenz,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmp3,ilenz,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

      call MPI_IRECV(tmp6,ileny,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmp8,ileny,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(2),ierr)
      call MPI_ISEND(tmp5,ileny,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(3),ierr) 
      call MPI_ISEND(tmp7,ileny,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

      end subroutine exchng8
