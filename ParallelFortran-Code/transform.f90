MODULE transform
  USE runparms
  USE gridparms
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_fft,             finalize_fft,                  &
            xz_fft_real_to_complex,     xz_transform,                  &
            xz_fft_complex_to_real,     xz_fft_complex_to_real_padded, &
            dxz_transform,              dxz_transform_padded,          &
            small_fft_complex_to_real,  small_fft_complex_to_real_pad, &
            small_fft_real_to_complex,  compute_divergence_of_fluxes,  &
            output_fft_real_to_complex, output_fft_complex_to_real,    &
            FFTW_REAL_TO_COMPLEX,       FFTW_COMPLEX_TO_REAL


  INCLUDE  'mpif.h'

  COMPLEX(KIND=prec), DIMENSION(:,:), ALLOCATABLE, SAVE :: cf, coutput
  COMPLEX(KIND=prec), DIMENSION(:), ALLOCATABLE, SAVE   :: ctmp1, ctmp2, ctmp3
  REAL(KIND=prec), DIMENSION(:,:), ALLOCATABLE, SAVE    :: rtmp, rpad, rsmall,&
       rsmallpad

  REAL(KIND=prec) :: scale
  REAL(KIND=prec), PARAMETER :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

! fftw plans
  INTEGER*8,SAVE :: &
       x_p2f_plan, z_p2f_plan, x_f2p_plan, z_f2p_plan, x_padded_f2p_plan, &
       xz_p2f_plan, xz_f2p_plan, output_p2f_plan, output_f2p_plan, &
       x_small_p2f_plan, x_small_f2p_plan, z_small_p2f_plan,  &
       z_small_f2p_plan, x_small_pad_plan, x_p2f_plan3, z_p2f_plan3

! fftw and mpi_fftw constants
  INTEGER, PARAMETER :: &
       FFTW_FORWARD=-1,               FFTW_BACKWARD=+1, &
       FFTW_REAL_TO_COMPLEX=-1,       FFTW_COMPLEX_TO_REAL=1, &
       FFTW_DESTROY_INPUT=1, &
       FFTW_UNALIGNED=2,              FFTW_CONSERVE_MEMORY=4, &
       FFTW_EXHAUSTIVE=8,             FFTW_PRESERVE_INPUT=16, &
       FFTW_PATIENT=32,               FFTW_ESTIMATE=64, &
       FFTW_ESTIMATE_PATIENT=128,     FFTW_BELIEVE_PCOST=256, &
       FFTW_DFT_R2HC_ICKY=512,        FFTW_NONTHREADED_ICKY=1024, &
       FFTW_NO_BUFFERING=2048,        FFTW_NO_INDIRECT_OP=4096,   &
       FFTW_ALLOW_LARGE_GENERIC=8192, FFTW_NO_RANK_SPLITS=16384,   &
       FFTW_NO_VRANK_SPLITS=32768,    FFTW_NO_VRECURSE=65536, &
       FFTW_R2HC=0,                   FFTW_HC2R=1,       &
       FFTW_DHT=2,                    FFTW_REDFT00=3,    &
       FFTW_REDFT01=4,                FFTW_REDFT10=5,    &
       FFTW_REDFT11=6,                FFTW_RODFT00=7,    &
       FFTW_RODFT01=8,                FFTW_RODFT10=9,    &
       FFTW_RODFT11=10

  INTEGER n, x, z, jj, index, index2, offset, stride, mpierr, ierr
  INTEGER, DIMENSION(2) :: inembed(2), onembed(2)

CONTAINS
!=================================================================
!==========================INITIALIZE FFT=========================
!=================================================================
  SUBROUTINE initialize_fft
    IMPLICIT NONE

    ! ALLOCATE ARRAYS TO BE USED IN TRANSFORMS.
    ALLOCATE(cf(nx2/2+1,nz2), coutput(nx/2+1,nz), &
         ctmp1(padded_local_size), ctmp2(padded_local_size), &
         ctmp3(padded_local_size), rtmp(nx2,local_nz), &
         rpad(1:nx2,0:local_nz+1), rsmall(nx,local_nz_small), &
         rsmallpad(1:nx,0:local_nz_small+1), STAT=ierr)

    ! SET UP PHYSICAL TO FOURIER TRANSFORM FOR X DIRECTION.
    inembed(1) = nx2+2
    inembed(2) = local_nz
    onembed(1) = nx2/2+1
    onembed(2) = local_nz
    CALL dfftw_plan_many_dft_r2c(x_p2f_plan, 1, nx2, local_nz, &
         ctmp1,inembed,1,nx2+2, ctmp1,onembed,1,nx2/2+1, FFTW_PATIENT)
    CALL dfftw_plan_many_dft_r2c(x_p2f_plan3, 1, nx2, local_nz, &
         ctmp3,inembed,1,nx2+2, ctmp3,onembed,1,nx2/2+1, FFTW_PATIENT)

    ! SET UP FOURIER TO PHYSICAL TRANSFORM FOR X DIRECTION.
    inembed(1) = nx2/2+1
    inembed(2) = local_nz
    onembed(1) = nx2
    onembed(2) = local_nz
    CALL dfftw_plan_many_dft_c2r(x_f2p_plan, 1, nx2, local_nz, &
         ctmp1,inembed,1,nx2/2+1, rtmp,onembed,1,nx2, FFTW_PATIENT)

    ! SET UP PADDED FOURIER TO PHYSICAL TRANSFORM FOR X DIRECTION.
    inembed(1) = nx2/2+1
    inembed(2) = local_nz+2
    onembed(1) = nx2
    onembed(2) = local_nz+2
    CALL dfftw_plan_many_dft_c2r(x_padded_f2p_plan, 1, nx2, local_nz+2, &
         ctmp1,inembed,1,nx2/2+1, rpad,onembed,1,nx2, FFTW_PATIENT)

    ! SET UP PHYSICAL TO FOURIER TRANSFORM FOR Z DIRECTION.
    inembed(1) = nz2
    inembed(2) = local_nx
    onembed(1) = nz2
    onembed(2) = local_nx
    CALL dfftw_plan_many_dft(z_p2f_plan, 1, nz2, local_nx, &
         ctmp2,inembed,1,nz2, ctmp2,onembed,1,nz2, FFTW_FORWARD, FFTW_PATIENT)
    CALL dfftw_plan_many_dft(z_p2f_plan3, 1, nz2, local_nx, &
         ctmp3,inembed,1,nz2, ctmp3,onembed,1,nz2, FFTW_FORWARD, FFTW_PATIENT)
    ! SET UP FOURIER TO PHYSICAL TRANSFORM FOR Z DIRECTION.
    CALL dfftw_plan_many_dft(z_f2p_plan, 1, nz2, local_nx,           &
         ctmp2,inembed,1,nz2, ctmp2,onembed,1,nz2, FFTW_BACKWARD, FFTW_PATIENT)

    !! INITIALIZE 2D DE-ALIASED (i.e. NX2 x NZ2) TRANSFORMS
    ! INITIALIZE PLAN FOR 2D PHYSICAL TO FOURIER TRANSFORM
    CALL dfftw_plan_dft_r2c_2d(xz_p2f_plan,nx2,nz2,cf,cf, FFTW_PATIENT)
    ! INITIALIZE PLAN FOR 2D FOURIER TO PHYSICAL TRANSFORM
    CALL dfftw_plan_dft_c2r_2d(xz_f2p_plan,nx2,nz2,cf,cf, FFTW_PATIENT)

    !! INITIALIZE 2D SMALL (i.e. NX x NZ) TRANSFORMS
    ! INITIALIZE PLAN FOR 2D PHYSICAL TO FOURIER TRANSFORM
    CALL dfftw_plan_dft_r2c_2d(output_p2f_plan,nx,nz,coutput,coutput, &
         FFTW_PATIENT)
    ! INITIALIZE PLAN FOR 2D FOURIER TO PHYSICAL TRANSFORM
    CALL dfftw_plan_dft_c2r_2d(output_f2p_plan,nx,nz,coutput,coutput, &
         FFTW_PATIENT)

    !! INITIALIZE 1D SMALL (i.e. NX x NZ) TRANSFORMS
    ! SET UP PHYSICAL TO FOURIER TRANSFORM FOR X DIRECTION.
    inembed(1) = nx+2
    inembed(2) = local_nz_small
    onembed(1) = nx/2+1
    onembed(2) = local_nz_small
    CALL dfftw_plan_many_dft_r2c(x_small_p2f_plan, 1, nx, local_nz_small, &
         ctmp1,inembed,1,nx+2,ctmp1,onembed,1,nx/2+1, FFTW_PATIENT)

    ! SET UP FOURIER TO PHYSICAL TRANSFORM FOR X DIRECTION.
    inembed(1) = nx/2+1
    inembed(2) = local_nz_small
    onembed(1) = nx
    onembed(2) = local_nz_small
    CALL dfftw_plan_many_dft_c2r(x_small_f2p_plan, 1, nx, local_nz_small, &
         ctmp1,inembed,1,nx/2+1, rsmall,onembed,1,nx, FFTW_PATIENT)

    ! SET UP PADDED FOURIER TO PHYSICAL TRANSFORM FOR X DIRECTION.
    inembed(1) = nx/2+1
    inembed(2) = local_nz_small+2
    onembed(1) = nx
    onembed(2) = local_nz_small+2
    CALL dfftw_plan_many_dft_c2r(x_small_pad_plan, 1, nx, local_nz_small+2, &
         ctmp1,inembed,1,nx/2+1, rsmallpad,onembed,1,nx, FFTW_PATIENT)

    ! SET UP PHYSICAL TO FOURIER TRANSFORM FOR Z DIRECTION.
    inembed(1) = nz
    inembed(2) = local_nx
    onembed(1) = nz
    onembed(2) = local_nx
    CALL dfftw_plan_many_dft(z_small_p2f_plan, 1, nz, local_nx, &
         ctmp2,inembed,1,nz, ctmp2,onembed,1,nz, FFTW_FORWARD, FFTW_PATIENT)

    ! SET UP FOURIER TO PHYSICAL TRANSFORM FOR Z DIRECTION.
    CALL dfftw_plan_many_dft(z_small_f2p_plan, 1, nz, local_nx,           &
         ctmp2,inembed,1,nz, ctmp2,onembed,1,nz, FFTW_BACKWARD, FFTW_PATIENT)

  END SUBROUTINE initialize_fft
!=================================================================
!===========================FINALIZE FFT==========================
!=================================================================
  SUBROUTINE finalize_fft
    IMPLICIT NONE

    CALL dfftw_cleanup()

    CALL dfftw_destroy_plan(x_p2f_plan)
    CALL dfftw_destroy_plan(xz_p2f_plan)
    CALL dfftw_destroy_plan(z_p2f_plan)
    CALL dfftw_destroy_plan(output_p2f_plan)
    CALL dfftw_destroy_plan(x_f2p_plan)
    CALL dfftw_destroy_plan(x_p2f_plan3)
    CALL dfftw_destroy_plan(z_p2f_plan3)
    CALL dfftw_destroy_plan(xz_f2p_plan)
    CALL dfftw_destroy_plan(z_f2p_plan)
    CALL dfftw_destroy_plan(output_f2p_plan)
    CALL dfftw_destroy_plan(x_padded_f2p_plan)
    CALL dfftw_destroy_plan(x_small_p2f_plan)
    CALL dfftw_destroy_plan(x_small_f2p_plan)
    CALL dfftw_destroy_plan(z_small_p2f_plan)
    CALL dfftw_destroy_plan(z_small_f2p_plan)
    CALL dfftw_destroy_plan(x_small_pad_plan)

    IF (me == 0) write(*,*) 'fftw plans destroyed'

    DEALLOCATE(cf, coutput, ctmp1, ctmp2, ctmp3, &
         rtmp, rpad, rsmall, rsmallpad, STAT=ierr)

  END SUBROUTINE finalize_fft
!=====================================================================
!=================XZ FFT REAL TO COMPLEX==============================
!=====================================================================
  SUBROUTINE xz_fft_real_to_complex(data,cdata)
    IMPLICIT NONE

    REAL(KIND=prec), INTENT(IN)     :: data(1:nx2,1:local_nz)
    COMPLEX(KIND=prec), INTENT(OUT) :: cdata(1:nz,1:local_nx)
    COMPLEX(KIND=prec), DIMENSION(nx/2+1,nz) :: ctmp
  
    IF (nproc == 0) THEN
       CALL xz_transform(ctmp,data,FFTW_REAL_TO_COMPLEX)
       cdata = TRANSPOSE(ctmp)
    ELSE
       CALL xz_fft_real_to_complex_parallel(data,cdata)
    END IF

  END SUBROUTINE xz_fft_real_to_complex
!=====================================================================
!=================XZ FFT COMPLEX TO REAL==============================
!=====================================================================
  SUBROUTINE xz_fft_complex_to_real(cdata,data)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(IN), DIMENSION(1:nz,1:local_nx) :: cdata
    REAL(KIND=prec), INTENT(OUT), DIMENSION(1:nx2,1:local_nz)  :: data
    COMPLEX(KIND=prec), DIMENSION(nx/2+1,nz) :: ctmp
  
    IF (nproc == 0) THEN
       ctmp = TRANSPOSE(cdata)
       CALL xz_transform(ctmp,data,FFTW_COMPLEX_TO_REAL)
    ELSE
       CALL xz_fft_complex_to_real_parallel(cdata,data)
    END IF

  END SUBROUTINE xz_fft_complex_to_real
!=====================================================================
!=================COMPUTE DIVERGENCE OF FLUXES========================
!=====================================================================
  SUBROUTINE compute_divergence_of_fluxes(datax,datay,dataz,cdata)
    IMPLICIT NONE

    REAL(KIND=prec),INTENT(IN),DIMENSION(1:nx2,1:local_nz) :: datax,datay,dataz
    COMPLEX(KIND=prec), INTENT(OUT), DIMENSION(1:nz,1:local_nx) :: cdata

    IF (nproc == 0) THEN
!!$       CALL xz_transform(ctmp,datay,FFTW_REAL_TO_COMPLEX)
!!$       DO x = 1,nx/2
!!$          cdata(1:nz,x) = ctmp(x,1:nz)
!!$       END DO
!!$       CALL xz_transform(ctmp,datax,FFTW_REAL_TO_COMPLEX)
!!$       DO x = 1,nx/2
!!$          cdata(1:nz,x) = cdata(1:nz,x) + ikx(x)*ctmp(x,1:nz)
!!$       END DO
!!$       CALL xz_transform(ctmp,dataz,FFTW_REAL_TO_COMPLEX)
!!$       DO x = 1,nx/2
!!$          cdata(1:nz,x) = cdata(1:nz,x) + ikz(1:nz)*ctmp(x,1:nz)
!!$       END DO
    ELSE
       CALL parallel_divergence_of_fluxes(datax,datay,dataz,cdata)
    END IF

  END SUBROUTINE compute_divergence_of_fluxes
!=====================================================================
!=================XZ FFT REAL TO COMPLEX PARALLEL=====================
!=====================================================================
  SUBROUTINE xz_fft_real_to_complex_parallel(data,cdata)
    IMPLICIT NONE

    REAL(KIND=prec), INTENT(IN)     :: data(1:nx2,1:local_nz)
    COMPLEX(KIND=prec), INTENT(OUT) :: cdata(1:nz,1:local_nx)

    ! COPY data INTO ctmp1
    index = 1
    DO z = 1,local_nz
       DO x = 1,nx2/2
          ctmp1(index) = CMPLX(data(2*x-1,z),data(2*x,z),KIND=prec)
          index = index + 1
       END DO
       ctmp1(index) = zero
       index = index + 1
    END DO

    ! FFT IN X-DIRECTION, PERFORMING IN-PLACE TRANSFORM IN ctmp1.
    CALL dfftw_execute(x_p2f_plan)
!!$    CALL rfftwnd_f77_real_to_complex(x_p2f_plan, local_nz,            &
!!$                                     data, 1, nx2, ctmp1, 1, nx2/2+1)

    ! SORT DATA INTO ARRAY ctmp2 IN PREPARATION FOR THE DISTRIBUTION 
    ! OF THE DATA AMONGST THE PROCESSORS.  THE DATA IS NOW DIVIDED 
    ! ALONG THE Z-DIRECTION INTO CHUNKS OF SIZE NX2*LOCAL_NZ.  
    ! AFTER PASSING, IT WILL BE DIVIDED ALONG THE X-DIRECTION 
    ! INTO CHUNKS OF SIZE local_nx.
    stride = nx2/2+1
    DO n = 0,nproc-1
       index = 1 + sdispl(n)
       DO x = x_start(n)+1,x_start(n)+nx_proc(n)
          DO z = 1,local_nz
             ctmp2(index) = ctmp1(x+(z-1)*stride)
             index = index + 1
          END DO
       END DO
    END DO

    ! PASS THE DATA AMONG THE PROCESSORS.  ALTHOUGH THE CHUNKS MAY NOT BE
    ! OF UNIFORM SIZE (LOCAL_NX AND LOCAL_NZ MAY VARY AMONG THE PROCESSORS),
    ! WE HAVE CHOSEN TO PASS EQUALLY-SIZED CHUNKS TO TAKE ADVANTAGE OF THE
    ! MPI_ALLTOALL COMMAND FOR SIMPLER CODING AND IN THE HOPE THAT THIS
    ! ROUTINE WILL BE BETTER OPTIMIZED THAN MPI_ALLTOALLV.
    CALL MPI_Alltoallv( &
         ctmp2, scount, sdispl, MPI_DOUBLE_COMPLEX, &
         ctmp1, rcount, rdispl, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, mpierr)

    ! SORT THE DATA ON THE NEW PROCESSOR TO COMPLETE THE Z-DIRECTION
    ! ALONG THE FASTEST-VARYING DIRECTION OF THE ARRAY.
    index = 0
    DO x = 1,local_nx
       DO n = 0,nproc-1
          index2 = rdispl(n) + (x-1)*nz_proc(n)
          ctmp2(index+1:index+nz_proc(n))=ctmp1(index2+1:index2+nz_proc(n))
          index = index + nz_proc(n)
       END DO
    END DO

    ! FFT IN THE Z-DIRECTION.
    CALL dfftw_execute(z_p2f_plan)
!!$    CALL fftwnd_f77(z_p2f_plan,local_nx,ctmp2,1,nz2,ctmp1,1,nz2)

    ! SCALE THE TRANSFORMED VARIABLE, REMOVE PADDING, ADD COEFFICIENTS 
    ! FOR TWO DIFFERENT PIECES OF SAWTOOTH MODE, AND PUT INTO cdata
    scale = one/DBLE(nx2*nz2)
    offset = nz2 - nz
    DO x = 1,local_nx
       index = nz2*(x-1)
       cdata(1:nz/2,x)    = scale*ctmp2(index+1:index+nz/2)
       cdata(nz/2+1,x)    =  &
            scale*(ctmp2(index+nz/2+1) + CONJG(ctmp2(index+offset+nz/2+1)))
       cdata(nz/2+2:nz,x) = scale*ctmp2(index+offset+nz/2+2:index+nz2)
    END DO
    IF (local_x_start + local_nx == nx/2+1) THEN
       cdata(1:nz,local_nx) = two*cdata(1:nz,local_nx)
    END IF

  END SUBROUTINE xz_fft_real_to_complex_parallel
!=====================================================================
!=================XZ FFT COMPLEX TO REAL PARALLEL=====================
!=====================================================================
  SUBROUTINE xz_fft_complex_to_real_parallel(cdata,data)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(IN), DIMENSION(1:nz,1:local_nx) :: cdata
    REAL(KIND=prec), INTENT(OUT), DIMENSION(1:nx2,1:local_nz)  :: data

    IF (nz2 > nz) THEN
       ! PAD VARIABLE IN Z-DIRECTION (TO LENGTH NZ2) AND PUT INTO ctmp2
       DO x = 1,local_nx
          index = nz2*(x-1)
          ctmp2(index+1:index+nz/2)          = cdata(1:nz/2,x)
          ctmp2(index+nz/2+1)                = half*cdata(nz/2+1,x)
          ctmp2(index+nz/2+2:index+nz2-nz/2) = zero
          ctmp2(index+nz2-nz/2+1)            = half*CONJG(cdata(nz/2+1,x))
          ctmp2(index+nz2-nz/2+2:index+nz2)  = cdata(nz/2+2:nz,x)
       END DO
       IF ((nx2 > nx) .and. (local_x_start + local_nx == nx/2+1)) THEN
          ctmp2(index+1:index+nz2) = half*ctmp2(index+1:index+nz2)
       END IF
    ELSE
       DO x = 1,local_nx
          index = nz2*(x-1)
          ctmp2(index+1:index+nz/2+1)          = cdata(1:nz/2+1,x)
          ctmp2(index+nz/2+2:index+nz2-nz/2+1) = zero
          ctmp2(index+nz2-nz/2+2:index+nz2)    = cdata(nz/2+2:nz,x)
       END DO
    END IF

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN Z-DIRECTION.
    CALL dfftw_execute(z_f2p_plan)
!!$    CALL fftwnd_f77(z_f2p_plan,local_nx,ctmp2,1,nz2,ctmp2,1,nz2) ! FFT IN Z

    ! SORT DATA INTO ARRAY tmp1, DISTRIBUTING ARRAY INTO CHUNKS ACCORDING
    ! TO THEIR Z COORDINATE, WHICH WILL BE BROKEN UP ACROSS THE PROCESSORS.
    DO n = 0,nproc-1
       index = 1 + rdispl(n)
       DO z = z_start(n)+1,z_start(n)+nz_proc(n)
          DO x = 1,local_nx
             ctmp1(index) = ctmp2(z+(x-1)*nz2)
             index = index + 1                               
          END DO
       END DO
    END DO

    ! EXCHANGE DATA AMONG PROCESSORS TO COMPLETE x-DIRECTION
    CALL MPI_Alltoallv( &
         ctmp1,rcount,rdispl,MPI_DOUBLE_COMPLEX,            &
         ctmp2,scount,sdispl,MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, mpierr)

    ! SORT DATA TO MAKE X DIRECTION CONTIGUOUS IN ARRAY
    index = 0
    DO z = 1,local_nz
       DO n = 0,nproc-1
          index2 = sdispl(n) + (z-1)*nx_proc(n)
          ctmp1(index+1:index+nx_proc(n))=ctmp2(index2+1:index2+nx_proc(n))
          index = index + nx_proc(n)
       END DO
       !! PAD THE NX2/2+1 ELEMENT IN THE x DIRECTION.
       IF (nx2 > nx) THEN
          ctmp1(index+1:index+(nx2-nx)/2) = zero  
          index = index + (nx2-nx)/2
       END IF
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN THE X-DIRECTION
    CALL dfftw_execute(x_f2p_plan)
    data = rtmp
!!$    CALL rfftwnd_f77_complex_to_real(x_f2p_plan,local_nz,          &
!!$         ctmp1, 1, nx2/2+1, data, 1, nx2)

  END SUBROUTINE xz_fft_complex_to_real_parallel
!=====================================================================
!=================XZ FFT COMPLEX TO REAL PADDED=======================
!=====================================================================
  SUBROUTINE xz_fft_complex_to_real_padded(cdata,data)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(IN), DIMENSION(1:nz,1:local_nx)  :: cdata
    REAL(KIND=prec),INTENT(OUT),DIMENSION(0:nx2+1,0:local_nz+1) :: data

    IF (nz2 > nz) THEN
       ! PAD VARIABLE IN Z-DIRECTION (TO LENGTH NZ2) AND PUT INTO ctmp2
       DO x = 1,local_nx
          index = nz2*(x-1)
          ctmp2(index+1:index+nz/2)          = cdata(1:nz/2,x)
          ctmp2(index+nz/2+1)                = half*cdata(nz/2+1,x)
          ctmp2(index+nz/2+2:index+nz2-nz/2) = zero
          ctmp2(index+nz2-nz/2+1)            = half*CONJG(cdata(nz/2+1,x))
          ctmp2(index+nz2-nz/2+2:index+nz2)  = cdata(nz/2+2:nz,x)
       END DO
       IF ((nx2 > nx) .and. (local_x_start + local_nx == nx/2+1)) THEN
          ctmp2(index+1:index+nz2) = half*ctmp2(index+1:index+nz2)
       END IF
    ELSE
       DO x = 1,local_nx
          index = nz2*(x-1)
          ctmp2(index+1:index+nz/2+1)          = cdata(1:nz/2+1,x)
          ctmp2(index+nz/2+2:index+nz2-nz/2+1) = zero
          ctmp2(index+nz2-nz/2+2:index+nz2)    = cdata(nz/2+2:nz,x)
       END DO
    END IF

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN Z-DIRECTION.
    CALL dfftw_execute(z_f2p_plan)
!!$    CALL fftwnd_f77(z_f2p_plan,local_nx,ctmp2,1,nz2,ctmp1,1,nz2)

    ! SORT DATA INTO ARRAY tmp1, DISTRIBUTING ARRAY INTO CHUNKS ACCORDING
    ! TO THEIR Z COORDINATE, WHICH WILL BE BROKEN UP ACROSS THE PROCESSORS.
    ! SEND EACH PROCESSOR TWO EXTRA VECTORS, WHICH SIT ALONGSIDE THAT
    ! PROCESSOR'S POINTS IN THE Z-DIRECTION.  (THIS IS THE PADDING.)
    ! THIS WILL ALLOW THE PROCESSORS TO COMPUTE LOCAL AVERAGES (NEEDED WHEN
    ! COMPUTING DYNAMIC MODEL COEFFICIENTS) WITHOUT EXTRA MESSAGE PASSING.
    !    scount_padded = scount + 2*maxval(nx_proc)
    DO n = 0,nproc-1
       index = 1 + rdispl_pad(n)
       DO z = z_start(n),z_start(n)+nz_proc(n)+1
          IF (z == 0) THEN
             jj = nz2
          ELSEIF (z == nz2+1) THEN
             jj = 1
          ELSE
             jj = z
          END IF
          DO x = 1,local_nx
             ctmp1(index) = ctmp2(jj+(x-1)*nz2)
             index = index + 1                               
          END DO
       END DO
    END DO

    ! EXCHANGE DATA AMONG PROCESSORS TO COMPLETE x-DIRECTION
    CALL Mpi_alltoallv( &
         ctmp1, rcount_pad, rdispl_pad, MPI_DOUBLE_COMPLEX, & 
         ctmp2, scount_pad, sdispl_pad, MPI_DOUBLE_COMPLEX, &
         MPI_COMM_WORLD, mpierr)

    ! SORT DATA TO MAKE X DIRECTION CONTIGUOUS IN ARRAY
    index = 0
    DO z = 1,local_nz+2
       DO n = 0,nproc-1
          index2 = sdispl_pad(n) + (z-1)*nx_proc(n)
          ctmp1(index+1:index+nx_proc(n))=ctmp2(index2+1:index2+nx_proc(n))
          index = index + nx_proc(n)
       END DO
       !! PAD THE NX2/2+1 ELEMENT IN THE x DIRECTION.
       IF (nx2 > nx) THEN
          ctmp1(index+1:index+(nx2-nx)/2) = zero  
          index = index + (nx2-nx)/2
       END IF
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN THE X-DIRECTION
    CALL dfftw_execute(x_padded_f2p_plan)
!!$    CALL rfftwnd_f77_complex_to_real(x_f2p_plan,local_nz+2, &
!!$         ctmp1, 1, nx2/2+1, rtmp, 1, nx2)
    data(0,0:local_nz+1)     = rpad(nx2,0:local_nz+1)
    data(1:nx2,0:local_nz+1) = rpad
    data(nx2+1,0:local_nz+1) = rpad(1,0:local_nz+1)
    
  END SUBROUTINE xz_fft_complex_to_real_padded
!=====================================================================
!========================XZ TRANSFORM=================================
!=====================================================================
  SUBROUTINE XZ_TRANSFORM(CU,U,ISIGN)
    IMPLICIT NONE

    COMPLEX(KIND=prec) :: CU(NX/2+1,NZ)
    REAL(KIND=prec)    :: U(NX2,NZ2)
    INTEGER     ISIGN, offset

    IF (ISIGN .EQ. FFTW_COMPLEX_TO_REAL) THEN
       !---------------------------------------------------
       !  TRANSFORM FROM FOURIER TO PHYSICAL SPACE
       !---------------------------------------------------
       cf = zero
       IF (NZ2 > NZ) THEN
          offset = nz2 - nz
          ! PAD THE DATA BEFORE THE TRANSFORM TO ALLOW FOR DE-ALIASING.
          DO z = 1,nz/2
             cf(1:nx/2+1,z) = cu(1:nx/2+1,z)
          END DO

          cf(1:nx/2+1,nz/2+1)     = half*cu(1:nx/2+1,nz/2+1)
          cf(1:nx/2+1,nz2-offset+1) = half*CONJG(cu(1:nx/2+1,nz/2+1))

          DO z = nz/2+2,nz
             cf(1:nx/2+1,z+offset) = cu(1:nx/2+1,z)
          END DO
          IF (nx2 > nx) cf(nx/2+1,1:nz2) = half*cf(nx/2+1,1:nz2)
       ELSE
          cf(1:nx/2+1,1:nz) = cu
       END IF
       ! CALL THE FFTW SUBROUTINE WHICH PERFORMS IN-PLACE 2D TRANSFORM IN cf.
       CALL dfftw_execute(xz_f2p_plan)
!!$          CALL rfftwnd_f77_one_complex_to_real(xz_f2p_plan, cf, u)
       
       DO z = 1,nz2
          DO x = 1,nx2/2
             u(2*x-1,z) = DBLE(cf(x,z))
             u(2*x,z)   = AIMAG(cf(x,z))
          END DO
       END DO

    ELSEIF (ISIGN .EQ. FFTW_REAL_TO_COMPLEX) THEN
       !------------------------------------------
       !  TRANSFORM FROM PHYSICAL TO FOURIER SPACE.
       !------------------------------------------
       
       ! COPY INPUT DATA INTO cf.
       DO z = 1,nz2
          DO x = 1,nx2/2
             cf(x,z) = CMPLX(u(2*x-1,z),u(2*x,z),KIND=prec)
          END DO
          cf(nx2/2+1,z) = zero
       END DO

       ! CALL THE FFTW SUBROUTINE WHICH PERFORMS IN-PLACE 2D TRANSFORM IN cf.
       CALL dfftw_execute(xz_p2f_plan)
!!$          CALL RFFTWND_F77_ONE_REAL_TO_COMPLEX(XZ_P2F_PLAN, U, CF)

       ! SCALE THE DATA.
       scale = one/DBLE(nx2*nz2)
       offset = nz2 - nz
       cu(1:nx/2+1,1:nz/2)    = scale*cf(1:nx/2+1,1:nz/2)
       cu(1:nx/2,nz/2+1)      = scale*(cf(1:nx/2,nz/2+1) &
                                     + CONJG(cf(1:nx/2,nz2-nz/2+1)))
       cu(1:nx/2+1,nz/2+2:nz) = scale*cf(1:nx/2+1,nz/2+2+offset:nz2)
       cu(nx/2+1,1:nz)        = two*cu(nx/2+1,1:nz)
    END IF

  END SUBROUTINE xz_transform
!=====================================================================
!=================COMPUTE DIVERGENCE OF FLUXES========================
!=====================================================================
  SUBROUTINE parallel_divergence_of_fluxes(datax,datay,dataz,cdata)
    IMPLICIT NONE

    REAL(KIND=prec),INTENT(IN),DIMENSION(1:nx2,1:local_nz) :: datax,datay,dataz
    COMPLEX(KIND=prec), INTENT(OUT), DIMENSION(1:nz,1:local_nx) :: cdata

    ! COPY datay INTO ctmp1 AND datax INTO ctmp3.
    index = 1
    DO z = 1,local_nz
       DO x = 1,nx2/2
          ctmp1(index) = CMPLX(datay(2*x-1,z),datay(2*x,z),KIND=prec)
          ctmp3(index) = CMPLX(datax(2*x-1,z),datax(2*x,z),KIND=prec)
          index = index + 1
       END DO
       ctmp1(index) = zero
       ctmp3(index) = zero
       index = index + 1
    END DO

    ! NOTE THAT THE X AND Z FLUXES ARE INPUT DIRECTLY,
    ! AND THEIR DERIVATIVES ARE COMPUTED IN THIS SUBROUTINE. 
    ! THE Y-FLUX IS ALREADY DIFFERENTIATED BEFORE BEING INPUT TO THIS ROUTINE.
    ! FFT X FLUX AND ALREADY COMPUTED Y-DERIVATIVE OF Y-FLUX IN X-DIRECTION.
    CALL dfftw_execute(x_p2f_plan)
    CALL dfftw_execute(x_p2f_plan3)
!!$       CALL rfftwnd_f77_real_to_complex(                               &
!!$            x_p2f_plan, local_nz, datay, 1, nx2, ctmp1, 1, nx2/2+1)
!!$       CALL rfftwnd_f77_real_to_complex(                               &
!!$            x_p2f_plan, local_nz, datax, 1, nx2, ctmp3, 1, nx2/2+1)

    ! COMPUTE D/DX(X-FLUX) + D/DY(Y-FLUX) AND SORT INTO ARRAY ctmp2
    ! IN PREPARATION FOR THE DISTRIBUTION OF THE DATA AMONGST THE
    ! PROCESSORS.  THE DATA IS NOW DIVIDED ALONG THE Z-DIRECTION INTO
    ! CHUNKS OF SIZE NX2*LOCAL_NZ.  AFTER PASSING, IT WILL BE DIVIDED
    ! ALONG THE X-DIRECTION INTO CHUNKS OF SIZE local_nx.
    stride = nx2/2+1                                  
    DO n = 0,nproc-1                                  
       index = 1 + sdispl(n)
       DO x = x_start(n)+1,x_start(n)+nx_proc(n)       
          DO z = 1,local_nz                             
             ctmp2(index)=ctmp1(x+(z-1)*stride)+ikx(x)*ctmp3(x+(z-1)*stride)
             index = index + 1                           
          END DO
       END DO
    END DO

    ! PASS THE DATA AMONG THE PROCESSORS.  ALTHOUGH THE CHUNKS MAY NOT BE
    ! OF UNIFORM SIZE (LOCAL_NX AND LOCAL_NZ MAY VARY AMONG THE PROCESSORS),
    ! WE HAVE CHOSEN TO PASS EQUALLY-SIZED CHUNKS TO TAKE ADVANTAGE OF THE
    ! MPI_ALLTOALL COMMAND FOR SIMPLER CODING AND IN THE HOPE THAT THIS
    ! ROUTINE WILL BE BETTER OPTIMIZED THAN MPI_ALLTOALLV.
    CALL Mpi_alltoallv( &
         ctmp2, scount, sdispl, MPI_DOUBLE_COMPLEX, &
         ctmp1, rcount, rdispl, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, mpierr)

    ! SORT THE DATA ON THE NEW PROCESSOR TO COMPLETE THE Z-DIRECTION
    ! ALONG THE FASTEST-VARYING DIRECTION OF THE ARRAY.
    index = 0
    DO x = 1,local_nx                              ! SORT DATA TO
       DO n = 0,nproc-1                            ! MAKE Z DIRECTION
          index2 = rdispl(n) + (x-1)*nz_proc(n)    ! CONTIGUOUS IN ARRAY
          ctmp3(index+1:index+nz_proc(n))=ctmp1(index2+1:index2+nz_proc(n))
          index = index + nz_proc(n)
       END DO
    END DO

    ! FFT IN THE Z-DIRECTION.  THE FOURIER TRANSFORMED QUANTITY
    !  D/DX(X-FLUX) + D/DY(Y-FLUX) WILL NOW SIT IN THE ARRAY tmp3
    CALL dfftw_execute(z_p2f_plan3) ! FFT IN Z
!!$       CALL fftwnd_f77(z_p2f_plan,local_nx,ctmp3,1,nz2,ctmp3,1,nz2) 

    ! NOW DEAL WITH THE Z-FLUX.  CALL THE FFT ROUTINE TO TRANSFORM IT FROM
    ! PHYSICAL TO FOURIER SPACE.
    CALL xz_fft_real_to_complex(dataz,cdata)

    ! NOW, TAKE THE Z-DERIVATIVE OF THE Z-FLUX (SITTING IN cdata) AND ADD
    ! IT TO THE FOURIER TRANSFORMED QUANTITY D/DX(X-FLUX) + D/DY(Y-FLUX) 
    ! IN tmp3 SCALE THE QUANTITY IN tmp3 AND REMOVE ITS PADDING IN THE 
    ! PROCESS.
    scale = one/DBLE(nx2*nz2)
    offset = nz2 - nz
    IF (local_x_start + local_nx == nx/2+1) THEN
       DO x = 1,local_nx-1
          index = nz2*(x-1)
          cdata(1:nz/2,x)    = ikz(1:nz/2)*cdata(1:nz/2,x) &
               + scale*ctmp3(index+1:index+nz/2)
          ! ADD TOGETHER TWO DIFFERENT COMPONENTS OF SAWTOOTH MODE THAT
          ! APPEAR ON DE-ALIASED GRID.
          cdata(nz/2+1,x)    = ikz(nz/2+1)*cdata(nz/2+1,x) &
              + scale*(ctmp3(index+nz/2+1) + CONJG(ctmp3(index+offset+nz/2+1)))
          cdata(nz/2+2:nz,x) = ikz(nz/2+2:nz)*cdata(nz/2+2:nz,x) &
               + scale*ctmp3(index+offset+nz/2+2:index+nz2)
       END DO
       x = local_nx
       index = nz2*(x-1)
       cdata(1:nz/2,x)    = ikz(1:nz/2)*cdata(1:nz/2,x) &
            + two*scale*ctmp3(index+1:index+nz/2)
       ! ADD TOGETHER TWO DIFFERENT COMPONENTS OF SAWTOOTH MODE THAT
       ! APPEAR ON DE-ALIASED GRID.
       cdata(nz/2+1,x)    = ikz(nz/2+1)*cdata(nz/2+1,x) &
          + two*scale*(ctmp3(index+nz/2+1) + CONJG(ctmp3(index+offset+nz/2+1)))
       cdata(nz/2+2:nz,x) = ikz(nz/2+2:nz)*cdata(nz/2+2:nz,x) &
            + two*scale*ctmp3(index+offset+nz/2+2:index+nz2)
    ELSE
       DO x = 1,local_nx
          index = nz2*(x-1)
          cdata(1:nz/2,x)    = ikz(1:nz/2)*cdata(1:nz/2,x) &
               + scale*ctmp3(index+1:index+nz/2)
          ! ADD TOGETHER TWO DIFFERENT COMPONENTS OF SAWTOOTH MODE THAT
          ! APPEAR ON DE-ALIASED GRID.
          cdata(nz/2+1,x)    = ikz(nz/2+1)*cdata(nz/2+1,x) &
               + scale*(ctmp3(index+nz/2+1) + CONJG(ctmp3(index+offset+nz/2+1)))
          cdata(nz/2+2:nz,x) = ikz(nz/2+2:nz)*cdata(nz/2+2:nz,x) &
               + scale*ctmp3(index+offset+nz/2+2:index+nz2)
       END DO
    END IF

  END SUBROUTINE parallel_divergence_of_fluxes
!=====================================================================
!=================DXZ TRANSFORM COMPLEX TO REAL=======================
!=====================================================================
  SUBROUTINE dxz_transform(cdata,data,datax,dataz)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(IN), DIMENSION(1:nz,1:local_nx) :: cdata
    REAL(KIND=prec),INTENT(OUT),DIMENSION(1:nx2,1:local_nz) :: data,datax,dataz

! TEMPORARY ARRAYS AND VARIABLES.
    COMPLEX(KIND=prec), DIMENSION(1:nz,1:local_nx)            :: cdataz

    ! THIS SUBROUTINE TAKES A FIELD IN FOURIER SPACE AND RETURNS THAT FIELD
    ! ALONG WITH ITS X- AND Z-DERIVATES, ALSO IN PHYSICAL SPACE.  THIS IS
    ! USEFUL WHEN THE STRESS TENSOR IS NEEDED IN PHYSICAL SPACE, FOR EXAMPLE
    ! WHEN COMPUTING THE COEFFICIENTS ASSOCIATED WITH THE DYNAMIC MODEL
    ! FOR LARGE EDDY SIMULATION.

    ! PAD VARIABLE IN Z-DIRECTION (TO LENGTH NZ2) AND PUT INTO ctmp2
    IF (nz2 > nz) THEN
       DO x = 1,local_nx
          index = nz2*(x-1)
          ctmp2(index+1:index+nz/2)          = cdata(1:nz/2,x)
          ctmp2(index+nz/2+1)                = half*cdata(nz/2+1,x)
          ctmp2(index+nz/2+2:index+nz2-nz/2) = zero
          ctmp2(index+nz2-nz/2+1)            = half*CONJG(cdata(nz/2+1,x))
          ctmp2(index+nz2-nz/2+2:index+nz2)  = cdata(nz/2+2:nz,x)
       END DO
       IF ((nx2 > nx) .and. (local_x_start + local_nx == nx/2+1)) THEN
          ctmp2(index+1:index+nz2) = half*ctmp2(index+1:index+nz2)
       END IF
    ELSE
       DO x = 1,local_nx
          index = nz2*(x-1)
          ctmp2(index+1:index+nz/2+1)          = cdata(1:nz/2+1,x)
          ctmp2(index+nz/2+2:index+nz2-nz/2+1) = zero
          ctmp2(index+nz2-nz/2+2:index+nz2)    = cdata(nz/2+2:nz,x)
       END DO
    END IF

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN Z-DIRECTION.
    CALL dfftw_execute(z_f2p_plan)
!!$       CALL fftwnd_f77(z_f2p_plan,local_nx,ctmp2,1,nz2,ctmp2,1,nz2)

    ! SORT DATA INTO ARRAY tmp1, DISTRIBUTING ARRAY INTO CHUNKS ACCORDING
    ! TO THEIR Z COORDINATE, WHICH WILL BE BROKEN UP ACROSS THE PROCESSORS.
    DO n = 0,nproc-1
       index = 1 + rdispl(n)
       DO z = z_start(n)+1,z_start(n)+nz_proc(n)
          DO x = 1,local_nx
             ctmp1(index) = ctmp2(z+(x-1)*nz2)
             index = index + 1                               
          END DO
       END DO
    END DO

    ! EXCHANGE DATA AMONG PROCESSORS TO COMPLETE x-DIRECTION
    CALL Mpi_alltoallv( &
         ctmp1,rcount,rdispl,MPI_DOUBLE_COMPLEX,          &
         ctmp2,scount,sdispl,MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, mpierr)

    ! SORT DATA TO MAKE X DIRECTION CONTIGUOUS IN ARRAY
    index = 0
    DO z = 1,local_nz
       DO n = 0,nproc-1
          index2 = sdispl(n) + (z-1)*nx_proc(n)
          ctmp1(index+1:index+nx_proc(n))=ctmp2(index2+1:index2+nx_proc(n))
          index = index + nx_proc(n)
       END DO
       !! PAD THE NX2/2+1 ELEMENT IN THE x DIRECTION.
       IF (nx2 > nx) THEN
          ctmp1(index+1:index+(nx2-nx)/2) = zero  
          index = index + (nx2-nx)/2
       END IF
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN THE X-DIRECTION
    CALL dfftw_execute(x_f2p_plan)
    data = rtmp
!!$       CALL rfftwnd_f77_complex_to_real(x_f2p_plan,local_nz,          &
!!$            ctmp1, 1, nx2/2+1, data, 1, nx2)
    ! REPEAT THIS PROCESS FOR THE X-DERIVATIVE
    ! SORT DATA TO MAKE X DIRECTION CONTIGUOUS IN ARRAY
    index = 0
    DO z = 1,local_nz
       DO n = 0,nproc-1
          index2 = sdispl(n) + (z-1)*nx_proc(n)
          ctmp1(index+1:index+nx_proc(n)) =                        &
               ikx(x_start(n)+1:x_start(n)+nx_proc(n))             &
               *ctmp2(index2+1:index2+nx_proc(n))
          index = index + nx_proc(n)
       END DO
       ! PAD ARRAY IN X-DIRECTION TO LENGTH NX2/2+1 AND ZERO OUT 
       ! THE NX2/2+1 ELEMENT IN THE x DIRECTION (THE SAWTOOTH MODE).
       IF (nx2 > nx) THEN
          ctmp1(index+1:index+(nx2-nx)/2) = zero  
          index = index + (nx2-nx)/2
       END IF
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN THE X-DIRECTION
    CALL dfftw_execute(x_f2p_plan)
    datax = rtmp
!!$       CALL rfftwnd_f77_complex_to_real(x_f2p_plan,local_nz,          &
!!$            ctmp1, 1, nx2/2+1, datax, 1, nx2)

    ! NEXT, COMPUTE Z-DERIVATE OF ORIGINAL FIELD (STILL IN cdata)
    ! AND THEN TRANSFORM IT, PUTTING THE RESULT IN dataz.
    DO x = 1,local_nx
       cdataz(1:nz,x) = ikz(1:nz)*cdata(1:nz,x)
    END DO
    CALL xz_fft_complex_to_real(cdataz,dataz)

  END SUBROUTINE dxz_transform
!=====================================================================
!=================DXZ TRANSFORM COMPLEX TO REAL PADDED================
!=====================================================================
  SUBROUTINE dxz_transform_padded(cdata,data,datax,dataz)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(IN), DIMENSION(1:nz,1:local_nx) :: cdata
    REAL(KIND=prec),INTENT(OUT),DIMENSION(0:nx2+1,0:local_nz+1) :: &
         data,datax,dataz

! TEMPORARY ARRAYS AND VARIABLES.
    COMPLEX(KIND=prec), DIMENSION(1:nz,1:local_nx)     :: cdataz

    ! THIS SUBROUTINE TAKES A FIELD IN FOURIER SPACE AND RETURNS THAT FIELD
    ! ALONG WITH ITS X- AND Z-DERIVATES, ALSO IN PHYSICAL SPACE.  THIS IS
    ! USEFUL WHEN THE STRESS TENSOR IS NEEDED IN PHYSICAL SPACE, FOR EXAMPLE
    ! WHEN COMPUTING THE COEFFICIENTS ASSOCIATED WITH THE DYNAMIC MODEL
    ! FOR LARGE EDDY SIMULATION.

    ! PAD VARIABLE IN Z-DIRECTION (TO LENGTH NZ2) AND PUT INTO ctmp2
    IF (nz2 > nz) THEN
       DO x = 1,local_nx
          index = nz2*(x-1)
          ctmp2(index+1:index+nz/2)          = cdata(1:nz/2,x)
          ctmp2(index+nz/2+1)                = half*cdata(nz/2+1,x)
          ctmp2(index+nz/2+2:index+nz2-nz/2) = zero
          ctmp2(index+nz2-nz/2+1)            = half*CONJG(cdata(nz/2+1,x))
          ctmp2(index+nz2-nz/2+2:index+nz2)  = cdata(nz/2+2:nz,x)
       END DO
       IF ((nx2 > nx) .and. (local_x_start + local_nx == nx/2+1)) THEN
          ctmp2(index+1:index+nz2) = half*ctmp2(index+1:index+nz2)
       END IF
    ELSE
       DO x = 1,local_nx
          index = nz2*(x-1)
          ctmp2(index+1:index+nz/2+1)          = cdata(1:nz/2+1,x)
          ctmp2(index+nz/2+2:index+nz2-nz/2+1) = zero
          ctmp2(index+nz2-nz/2+2:index+nz2)    = cdata(nz/2+2:nz,x)
       END DO
    END IF

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN Z-DIRECTION.
    CALL dfftw_execute(z_f2p_plan)
!!$    CALL fftwnd_f77(z_f2p_plan,local_nx,ctmp2,1,nz2,ctmp2,1,nz2)

    ! SORT DATA INTO ARRAY tmp1, DISTRIBUTING ARRAY INTO CHUNKS ACCORDING
    ! TO THEIR Z COORDINATE, WHICH WILL BE BROKEN UP ACROSS THE PROCESSORS.
    ! SEND EACH PROCESSOR TWO EXTRA VECTORS, WHICH SIT ALONGSIDE THAT
    ! PROCESSOR'S POINTS IN THE Z-DIRECTION.  (THIS IS THE PADDING.)
    ! THIS WILL ALLOW THE PROCESSORS TO COMPUTE LOCAL AVERAGES (NEEDED WHEN
    ! COMPUTING DYNAMIC MODEL COEFFICIENTS) WITHOUT EXTRA MESSAGE PASSING.
    !    scount_padded = scount + 2*maxval(nx_proc)
    DO n = 0,nproc-1
       index = 1 + rdispl_pad(n)
       DO z = z_start(n),z_start(n)+nz_proc(n)+1
          IF (z == 0) THEN
             jj = nz2
          ELSEIF (z == nz2+1) THEN
             jj = 1
          ELSE
             jj = z
          END IF
          DO x = 1,local_nx
             ctmp1(index) = ctmp2(jj+(x-1)*nz2)
             index = index + 1                               
          END DO
       END DO
    END DO

    ! EXCHANGE DATA AMONG PROCESSORS TO COMPLETE x-DIRECTION
    CALL Mpi_alltoallv( &
         ctmp1, rcount_pad, rdispl_pad, MPI_DOUBLE_COMPLEX, & 
         ctmp2, scount_pad, sdispl_pad, MPI_DOUBLE_COMPLEX, &
         MPI_COMM_WORLD, mpierr)

    ! SORT DATA TO MAKE X DIRECTION CONTIGUOUS IN ARRAY
    index = 0
    DO z = 1,local_nz+2
       DO n = 0,nproc-1
          index2 = sdispl_pad(n) + (z-1)*nx_proc(n)
          ctmp1(index+1:index+nx_proc(n))=ctmp2(index2+1:index2+nx_proc(n))
          index = index + nx_proc(n)
       END DO
       !! PAD THE NX2/2+1 ELEMENT IN THE x DIRECTION.
       IF (nx2 > nx) THEN
          ctmp1(index+1:index+(nx2-nx)/2) = zero  
          index = index + (nx2-nx)/2
       END IF
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN THE X-DIRECTION
    !    WRITE(*,*) me, 'CALL FFT in X', x_f2p_plan, local_nz, nx2
    CALL dfftw_execute(x_padded_f2p_plan)
!!$       CALL rfftwnd_f77_complex_to_real(x_f2p_plan,local_nz+2, &
!!$            ctmp1, 1, nx2/2+1, rtmp, 1, nx2)

    ! FILL THE FIRST AND LAST ROWS OF THE MATRIX WITH VALUES FROM THE
    ! SECOND LAST AND SECOND ROWS.  (THESE ARE THE NEIGHBORING VALUES
    ! SINCE THE ARRAY IS PERIODIC IN X.)
    data(0,0:local_nz+1)     = rpad(nx2,0:local_nz+1)
    data(1:nx2,0:local_nz+1) = rpad
    data(nx2+1,0:local_nz+1) = rpad(1,0:local_nz+1)

    ! REPEAT THIS SORTING AND TRANSFORM FOR THE X-DERIVATIVE.
    ! SORT DATA TO MAKE X DIRECTION CONTIGUOUS IN ARRAY
    index = 0
    DO z = 1,local_nz+2
       DO n = 0,nproc-1
          index2 = sdispl_pad(n) + (z-1)*nx_proc(n)
          ctmp1(index+1:index+nx_proc(n)) =                        &
               ikx(x_start(n)+1:x_start(n)+nx_proc(n))             &
               *ctmp2(index2+1:index2+nx_proc(n))
          index = index + nx_proc(n)
       END DO
       !! PAD THE NX2/2+1 ELEMENT IN THE x DIRECTION.
       IF (nx2 > nx) THEN
          ctmp1(index+1:index+(nx2-nx)/2) = zero  
          index = index + (nx2-nx)/2
       END IF
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN THE X-DIRECTION
    !    WRITE(*,*) me, 'CALL FFT in X', x_f2p_plan, local_nz, nx2
    CALL dfftw_execute(x_padded_f2p_plan)
!!$       CALL rfftwnd_f77_complex_to_real(x_f2p_plan,local_nz+2, &
!!$            ctmp1, 1, nx2/2+1, rtmp, 1, nx2)

    ! FILL THE FIRST AND LAST ROWS OF THE MATRIX WITH VALUES FROM THE
    ! SECOND LAST AND SECOND ROWS.  (THESE ARE THE NEIGHBORING VALUES
    ! SINCE THE ARRAY IS PERIODIC IN X.)
    datax(0,0:local_nz+1)     = rpad(nx2,0:local_nz+1)
    datax(1:nx2,0:local_nz+1) = rpad
    datax(nx2+1,0:local_nz+1) = rpad(1,0:local_nz+1)

    ! NEXT, COMPUTE Z-DERIVATE OF ORIGINAL FIELD (STILL IN cdata)
    ! AND THEN TRANSFORM IT, PUTTING THE RESULT IN dataz.
    DO x = 1,local_nx
       cdataz(1:nz,x) = ikz(1:nz)*cdata(1:nz,x)
    END DO
    CALL xz_fft_complex_to_real_padded(cdataz,dataz)

  END SUBROUTINE dxz_transform_padded
!=====================================================================
!=================SMALL FFT COMPLEX TO REAL===========================
!=====================================================================
  SUBROUTINE small_fft_complex_to_real(cdata,data)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(IN), DIMENSION(1:nz,1:local_nx)     :: cdata
    REAL(KIND=prec), INTENT(OUT), DIMENSION(1:nx,1:local_nz_small) ::data

    ! PAD VARIABLE IN Z-DIRECTION (TO LENGTH NZ2) AND PUT INTO ctmp2
    DO x = 1,local_nx
       index = nz*(x-1)
       ctmp2(index+1:index+nz) = cdata(1:nz,x)
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN Z-DIRECTION.
    CALL dfftw_execute(z_small_f2p_plan)
!!$       CALL fftwnd_f77(z_small_f2p_plan,local_nx,ctmp2,1,nz,ctmp2,1,nz)

    ! SORT DATA INTO ARRAY tmp1, DISTRIBUTING ARRAY INTO CHUNKS ACCORDING
    ! TO THEIR Z COORDINATE, WHICH WILL BE BROKEN UP ACROSS THE PROCESSORS.
    DO n = 0,nproc-1
       index = 1 + rdispl_small(n)
       DO z = z_small_start(n)+1,z_small_start(n)+nz_small_proc(n)
          DO x = 1,local_nx
             ctmp1(index) = ctmp2(z+(x-1)*nz)
             index = index + 1                               
          END DO
       END DO
    END DO

    ! EXCHANGE DATA AMONG PROCESSORS TO COMPLETE x-DIRECTION
    CALL Mpi_alltoallv( &
         ctmp1,rcount_small,rdispl_small,MPI_DOUBLE_COMPLEX, &
         ctmp2,scount_small,sdispl_small,MPI_DOUBLE_COMPLEX, &
         MPI_COMM_WORLD, mpierr)

    ! SORT DATA TO MAKE X DIRECTION CONTIGUOUS IN ARRAY
    index = 0
    DO z = 1,local_nz_small
       DO n = 0,nproc-1
          index2 = sdispl_small(n) + (z-1)*nx_proc(n)
          ctmp1(index+1:index+nx_proc(n))=ctmp2(index2+1:index2+nx_proc(n))
          index = index + nx_proc(n)
       END DO
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN THE X-DIRECTION
    CALL dfftw_execute(x_small_f2p_plan)
    data = rsmall
!!$       CALL rfftwnd_f77_complex_to_real(x_small_f2p_plan,local_nz_small, &
!!$            ctmp1, 1, nx/2+1, data, 1, nx)

  END SUBROUTINE small_fft_complex_to_real
!=====================================================================
!=================SMALL FFT COMPLEX TO REAL PAD====================
!=====================================================================
  SUBROUTINE small_fft_complex_to_real_pad(cdata,data)
    IMPLICIT NONE

    COMPLEX(KIND=prec), INTENT(IN), DIMENSION(1:nz,1:local_nx)         :: cdata
    REAL(KIND=prec), INTENT(OUT), DIMENSION(0:nx+1,0:local_nz_small+1) :: data

    ! PAD VARIABLE IN Z-DIRECTION (TO LENGTH NZ2) AND PUT INTO ctmp2
    DO x = 1,local_nx
       index = nz*(x-1)
       ctmp2(index+1:index+nz) = cdata(1:nz,x)
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN Z-DIRECTION.
    CALL dfftw_execute(z_small_f2p_plan)
!!$       CALL fftwnd_f77(z_small_f2p_plan,local_nx,ctmp2,1,nz,ctmp1,1,nz)

    ! SORT DATA INTO ARRAY tmp1, DISTRIBUTING ARRAY INTO CHUNKS ACCORDING
    ! TO THEIR Z COORDINATE, WHICH WILL BE BROKEN UP ACROSS THE PROCESSORS.
    ! SEND EACH PROCESSOR TWO EXTRA VECTORS, WHICH SIT ALONGSIDE THAT
    ! PROCESSOR'S POINTS IN THE Z-DIRECTION.  (THIS IS THE PADDING.)
    ! THIS WILL ALLOW THE PROCESSORS TO COMPUTE LOCAL AVERAGES (NEEDED WHEN
    ! COMPUTING DYNAMIC MODEL COEFFICIENTS) WITHOUT ANY EXTRA MESSAGE PASSING.
    DO n = 0,nproc-1
       index = 1 + rdispl_small_pad(n)
       DO z = z_small_start(n),z_small_start(n)+nz_small_proc(n)+1
          IF (z == 0) THEN
             jj = nz
          ELSEIF (z == nz+1) THEN
             jj = 1
          ELSE
             jj = z
          END IF
          DO x = 1,local_nx
             ctmp1(index) = ctmp2(jj+(x-1)*nz)
             index = index + 1                               
          END DO
       END DO
    END DO

    ! EXCHANGE DATA AMONG PROCESSORS TO COMPLETE x-DIRECTION
    CALL Mpi_alltoallv( &
         ctmp1,rcount_small_pad,rdispl_small_pad,MPI_DOUBLE_COMPLEX, &
         ctmp2,scount_small_pad,sdispl_small_pad,MPI_DOUBLE_COMPLEX, &
         MPI_COMM_WORLD, mpierr)

    ! SORT DATA TO MAKE X DIRECTION CONTIGUOUS IN ARRAY
    index = 0
    DO z = 1,local_nz_small+2
       DO n = 0,nproc-1
          index2 = sdispl_small_pad(n) + (z-1)*nx_proc(n)
          ctmp1(index+1:index+nx_proc(n))=ctmp2(index2+1:index2+nx_proc(n))
          index = index + nx_proc(n)
       END DO
    END DO

    ! FOURIER TRANSFORM INTO PHYSICAL SPACE IN THE X-DIRECTION
    CALL dfftw_execute(x_small_pad_plan)
!!$       CALL rfftwnd_f77_complex_to_real(x_small_f2p_plan,local_nz_small+2, &
!!$            ctmp1, 1, nx/2+1, rsmallpad, 1, nx)
    data(0,   0:local_nz_small+1) = rsmallpad(nx,0:local_nz_small+1)
    data(1:nx,0:local_nz_small+1) = rsmallpad
    data(nx+1,0:local_nz_small+1) = rsmallpad(1, 0:local_nz_small+1)

  END SUBROUTINE small_fft_complex_to_real_pad
!=====================================================================
!=================SMALL FFT REAL TO COMPLEX===========================
!=====================================================================
  SUBROUTINE small_fft_real_to_complex(data,cdata)
    IMPLICIT NONE

    REAL(KIND=prec), INTENT(IN)     :: data(1:nx,1:local_nz_small)
    COMPLEX(KIND=prec), INTENT(OUT) :: cdata(1:nz,1:local_nx)

    ! COPY data INTO ctmp1
    index = 1
    DO z = 1,local_nz_small
       DO x = 1,nx/2
          ctmp1(index) = CMPLX(data(2*x-1,z),data(2*x,z),KIND=prec)
          index = index + 1
       END DO
       ctmp1(index) = zero
       index = index + 1
    END DO

    ! FFT IN X-DIRECTION, IN-PLACE TRANSFORM USING THE ARRAY ctmp1
    CALL dfftw_execute(x_small_p2f_plan)
!!$       CALL rfftwnd_f77_real_to_complex(x_small_p2f_plan, local_nz,   &
!!$            data, 1, nx, ctmp1, 1, nx/2+1)

    ! SORT DATA INTO ARRAY ctmp2 IN PREPARATION FOR THE DISTRIBUTION 
    ! OF THE DATA AMONGST THE PROCESSORS.  THE DATA IS NOW DIVIDED 
    ! ALONG THE Z-DIRECTION INTO CHUNKS OF SIZE NX2*LOCAL_NZ.  
    ! AFTER PASSING, IT WILL BE DIVIDED ALONG THE X-DIRECTION 
    ! INTO CHUNKS OF SIZE local_nx.
    stride = nx/2+1
    DO n = 0,nproc-1
       index = 1 + sdispl_small(n)
       DO x = x_start(n)+1,x_start(n)+nx_proc(n)
          DO z = 1,local_nz_small
             ctmp2(index) = ctmp1(x+(z-1)*stride)
             index = index + 1
          END DO
       END DO
    END DO

    ! PASS THE DATA AMONG THE PROCESSORS.  ALTHOUGH THE CHUNKS MAY NOT BE
    ! OF UNIFORM SIZE (LOCAL_NX AND LOCAL_NZ MAY VARY AMONG THE PROCESSORS),
    ! WE HAVE CHOSEN TO PASS EQUALLY-SIZED CHUNKS TO TAKE ADVANTAGE OF THE
    ! MPI_ALLTOALL COMMAND FOR SIMPLER CODING AND IN THE HOPE THAT THIS
    ! ROUTINE WILL BE BETTER OPTIMIZED THAN MPI_ALLTOALLV.
    CALL MPI_Alltoallv( &
         ctmp2,scount_small,sdispl_small,MPI_DOUBLE_COMPLEX, &
         ctmp1,rcount_small,rdispl_small,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, &
         mpierr)

    ! SORT THE DATA ON THE NEW PROCESSOR TO COMPLETE THE Z-DIRECTION
    ! ALONG THE FASTEST-VARYING DIRECTION OF THE ARRAY.
    index = 0
    DO x = 1,local_nx
       DO n = 0,nproc-1
          index2 = rdispl_small(n) + (x-1)*nz_small_proc(n)
          ctmp2(index+1:index+nz_small_proc(n)) = &
               ctmp1(index2+1:index2+nz_small_proc(n))
          index = index + nz_small_proc(n)
       END DO
    END DO

    ! FFT IN THE Z-DIRECTION.
    CALL dfftw_execute(z_small_p2f_plan)
!!$       CALL fftwnd_f77(z_small_p2f_plan,local_nx,ctmp2,1,nz,ctmp1,1,nz)

    ! SCALE THE TRANSFORMED VARIABLE AND PUT INTO cdata.
    scale = one/DBLE(nx*nz)
    DO x = 1,local_nx
       index = nz*(x-1)
       cdata(1:nz,x)  = scale*ctmp2(index+1:index+nz)
    END DO

  END SUBROUTINE small_fft_real_to_complex
!=================================================================
!================ OUTPUT FFT REAL TO COMPLEX =====================
!=================================================================
  SUBROUTINE output_fft_real_to_complex(in,cout)
    IMPLICIT NONE
    REAL(KIND=prec), DIMENSION(nx,nz), TARGET :: in
    COMPLEX(KIND=prec), DIMENSION(nx/2+1,nz), TARGET :: cout
    
    ! COPY INPUT ARRAY INTO coutput
    DO z = 1,nz
       DO x = 1,nx/2
          coutput(x,z) = CMPLX(in(2*x-1,z),in(2*x,z),KIND=prec)
       END DO
       coutput(nx/2+1,z) = zero
    END DO

    ! PERFORM 2D IN-PLACE TRANSFORM IN coutput.
    CALL dfftw_execute(output_p2f_plan)
!!$       CALL rfftwnd_f77_one_real_to_complex(output_p2f_plan, in, cout)

    ! SCALE AND COPY ARRAY INTO cout.
    cout = coutput/DBLE(nx*nz)

  END SUBROUTINE output_fft_real_to_complex
!=================================================================
!================ OUTPUT FFT REAL TO COMPLEX =====================
!=================================================================
  SUBROUTINE output_fft_complex_to_real(cin,out)
    IMPLICIT NONE
    REAL(KIND=prec), DIMENSION(nx,nz), TARGET :: out
    COMPLEX(KIND=prec), DIMENSION(nx/2+1,nz), TARGET :: cin

    ! COPY INPUT ARRAY INTO coutput
    coutput = cin
    
    ! PERFORM IN-PLACE 2D TRANSFORM IN coutput.
    CALL dfftw_execute(output_f2p_plan)
!!$       CALL rfftwnd_f77_one_complex_to_real(output_f2p_plan, cin, out)

    ! COPY INTO OUTPUT ARRAY out.
    DO z = 1,nz
       DO x = 1,nx/2
          out(2*x-1,z) = DBLE(coutput(x,z))
          out(2*x,z)   = AIMAG(coutput(x,z))
       END DO
    END DO

  END SUBROUTINE output_fft_complex_to_real
END MODULE transform



