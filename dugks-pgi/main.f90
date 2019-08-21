! --------------------------------------------------------------------------------
!       MPI F90 Code for Simulating 3D Decaying Compressible Homogeneous
!       Isotropic
!       Turbulence (DCHIT)
!
!       Using E_3,7^27 Gauss Hermite Guadrature from Shan et al. (2006
!       JFM)
!
!       The code is originally written by Dr. Lian-Ping Wang, June 2017.
!
! --------------------------------------------------------------------------------
      PROGRAM main
     use mpi
      use var_inc
      implicit none
      real pmean,pmeant,prms,prmst,rkin,rkinm
      real time1a,time2a,time3a,time4a,time5a,time6a
      character(len=200)::fiter
      character(len=200)::timerecord

      time1a= 0.0
      time2a= 0.0
      time3a= 0.0
      time4a= 0.0
      time5a= 0.0
      time6a= 0.0
  

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)      


! Create FFTW plans, store in plan1DF, plan1DB, plan2DRC, plan2DCR


      call para
     ! write(*,*)'para done'
      call allocarray
     ! write(*,*)'alloc done'
      call wave
     ! write(*,*)'wave done'

      IF(newrun)THEN   

      if(newinitflow)then  

      fiter=trim(dirdiag)//'iteration.dat'
      write(*,*)'pre_initrand'
    !    call initrand(iflowseed)
        write(*,*)'pre_initvel'
        call initvel
        write(*,*)'pre_initpop'
        call initpop
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        write(*,*)'pre_interation'
        GO TO 111
! interation stage
! keep velocity ux uy uz and density rho unchanged
! evolve f~ g~ T

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%store the initial velocity and density
        uxp = ux
        uyp = uy
        uzp = uz
        rhop = rho
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%istep=0 step number
        istep = 0
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% interation stage do not use
!limiter will not be used in interation
        I_Limiter = .false.
!        I_Limiter = .true.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start interation
        write(*,*)'iterstart'
      do
      !save T in last step
        tep = te

      !update rho u T f~ g~
        call interface 
        !$acc kernels
        call evol
!$acc end kernels
      !contrain:rho,(ux,uy,uz); evolve: T g~ h~  
      !another choice is only to constrain velocity
        call vel_den_constrain

      !call vel_constrain
        istep = istep + 1

      !monitor  kinetic energy
      rkin = sum(ux*ux+uy*uy+uz*uz)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(rkin,rkinm,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
      if(myid == 0) rkinm = rkinm/2.0/real(nx*ny*nz)

      !monitor  temperature 
      teerr = sum(abs(te - tep))
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(teerr,teerrm,1,MPI_REAL8,MPI_SUM,     &
                             MPI_COMM_WORLD,ierr)
      if(myid == 0)  teerrm = teerrm/real(nx*ny*nz)/T_ref
 
      !Output temperature norm
        open(24, file = trim(fiter), status = 'unknown',&
                 form = 'formatted', position = 'append')

        if(myid == 0 .and. mod(istep,1) == 0)                        &
            write(24,240)istep,dt,rkinm,teerrm
        close(24)
240   format(2x,i5,3(1pe16.6))

      !criterion for exit 
write(*,*)'here!'
      if(teerrm <= teepsl .or. istep > 15000)then
            if(myid == 0) write(*,*)'final relaxation => ',            &
                                    istep,teerrm,teepsl
            exit
      end if

        end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
111     continue
        write(*,*)'here'
        I_Limiter = .true.
        istep=0
        ttt=0.0!the total time
! save population for input of next run
!save g~ and h~ 
!        call saveinitflow
!        call outputflow
        write(*,*)'pre stat'
        call statistc!the skewness at 0
        write(*,*)'stat done'
        call gcheck
        write(*,*)'gchecked'

      else     
! readin population from previous saving
        call loadinitflow     
  
        istep = 0
        call initrand(iflowseed)
        call initvel
        call outputflow
!
      end if     

      ELSE
! load data from same number of processes  
      call loadcntdflow      
      END IF
!==========================================================
       timerecord=trim(dirdiag)//'timerecord.dat'
!==========================================================
      time_start = MPI_WTIME()
      ttt = 0.0!must 0

!==============start main loop
        do istep=istep0+1,istep0+nsteps
        !==================
        time_start = MPI_WTIME()
        call interface
        time_end = MPI_WTIME()
        time1a=time1a+(time_end-time_start)

        !==================
        time_start = MPI_WTIME()
        call evol
        time_end = MPI_WTIME()
        time3a=time3a+(time_end-time_start)

        !==================
        ttt = ttt + dt

        !%%%%%%%%%%time step and time used at the moment
        open(29, file = trim(timerecord), status = 'unknown',&
                 form = 'formatted', position = 'append')
        if(myid == 0 .and. mod(istep,10) == 0)                        &
            write(29,299)istep,dt,ttt
        close(29)
299   format(2x,i5,2(1pe16.6))

 
        !output the result
        !500 
!       if(mod(istep,nflowout) == 0) call outputflow 
        !10
       if(mod(istep,nstat) == 0) call  statistc
        !10
       if(mod(istep,ndiag) == 0)  call gcheck

      end do
!=================================================================
!=================================================================
! main loop ends here

       time_end = MPI_WTIME()
       time_diff = time_end - time_start

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(time_diff,time_max,1,MPI_REAL8,  &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

!test interface subroutie time
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(time1a,time_sum,1,MPI_REAL8,  &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
       if(myid == 0) write(*,*) 'interface',time_sum/istep/nproc
!test flux subroutie time
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(time2a,time_sum,1,MPI_REAL8,  &
                           MPI_SUM,MPI_COMM_WORLD,ierr)

       if(myid == 0) write(*,*) 'post-processing',time_sum/istep/nproc
!test evol subroutie time
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_ALLREDUCE(time3a,time_sum,1,MPI_REAL8,  &
                           MPI_SUM,MPI_COMM_WORLD,ierr)
       if(myid == 0) write(*,*) 'evol',time_sum/istep/nproc

! save data for continue run
      istep = istep0+nsteps
      call savecntdflow      

! Destroy plans
      call destroy_plan1D_RC
      call destroy_plan1D_CC

      call MPI_FINALIZE(ierr)

      END PROGRAM main
