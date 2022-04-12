MODULE pass_data
  USE runparms
  USE gridparms
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pass_neighbors, pass_slabs

  INCLUDE  'mpif.h'

CONTAINS
!=====================================================================
!==========================PASS NEIGHBORS=============================
!=====================================================================
  SUBROUTINE pass_neighbors(n1,n2,npad,f)
    IMPLICIT NONE

    INTEGER, INTENT(IN)                 :: n1, n2, npad
    REAL(KIND=prec), INTENT(INOUT)      :: f(1-npad:n1+npad,1-npad:n2+npad)
    REAL(KIND=prec), DIMENSION(n1)      :: passl, passr, recvl, recvr
    REAL(KIND=prec), DIMENSION(n1,npad) :: passl2, passr2, recvl2, recvr2

    INTEGER :: n, i, status(MPI_STATUS_SIZE,4), mpierr, req(4), nz_min
  
! ON CALLING pass_neighbors, f SHOULD HOLD DATA VALUES IN INTERIOR OF
! ARRAY, i.e. f(1:n1,1:n2) SHOULD HOLD NEW VALUES.  THIS ROUTINE WILL
! PAD THIS ARRAY WITH NEIGHBORING VALUES (WHICH MIGHT BE HELD ON OTHER
! PROCESSORS.

    nz_min = MINVAL(nz_proc)
    IF (npad <= nz_min) THEN

          CALL MPI_IRECV(recvr2,npad*n1,MPI_DOUBLE_PRECISION,&
               right,right,MPI_COMM_WORLD,req(1),mpierr)
          CALL MPI_IRECV(recvl2,npad*n1,MPI_DOUBLE_PRECISION,&
               left,nproc+left,MPI_COMM_WORLD,req(2),mpierr)
          passl2 = f(1:n1,1:npad)
          passr2 = f(1:n1,n2+1-npad:n2)
          CALL MPI_SEND(passl2,npad*n1,MPI_DOUBLE_PRECISION,&
               left,me,MPI_COMM_WORLD,mpierr)
          CALL MPI_SEND(passr2,npad*n1,MPI_DOUBLE_PRECISION,&
               right,nproc+me,MPI_COMM_WORLD,mpierr)
          DO i = 1,2
             CALL MPI_WAIT(req(i),status(1,i),mpierr)
          END DO
          f(1:n1,n2+1:n2+npad) = recvr2
          f(1:n1,1-npad:0)     = recvl2

    ELSE

       DO n = 1,npad
          CALL MPI_IRECV(recvr,  n1,MPI_DOUBLE_PRECISION,&
               right,right,MPI_COMM_WORLD,req(1),mpierr)
          CALL MPI_IRECV(recvl,   n1,MPI_DOUBLE_PRECISION,&
               left,nproc+left,MPI_COMM_WORLD,req(2),mpierr)
          passl = f(1:n1,n)
          passr = f(1:n1,n2+1-n)
          CALL MPI_SEND(passl,     n1,MPI_DOUBLE_PRECISION,&
               left,me,MPI_COMM_WORLD,mpierr)
          CALL MPI_SEND(passr,n1,MPI_DOUBLE_PRECISION,&
               right,nproc+me,MPI_COMM_WORLD,mpierr)
          DO i = 1,2
             CALL MPI_WAIT(req(i),status(1,i),mpierr)
          END DO
          f(1:n1,n2+n) = recvr
          f(1:n1,1-n)  = recvl
       END DO

    END IF

    DO n = 1,npad
       f(1-n,:)  = f(n1+1-n,:) ! PAD AT BEGINNING OF x-DIRECTION
       f(n1+n,:) = f(n,:)      ! PAD AT END OF x-DIRECTION
    END DO

  END SUBROUTINE pass_neighbors
!=====================================================================
!==========================PASS SLABS=============================
!=====================================================================
  SUBROUTINE pass_slabs(n1,n2,n3,npad,f)
    IMPLICIT NONE

    INTEGER, INTENT(IN)            :: n1, n2, n3, npad
    REAL(KIND=prec), DIMENSION(1-npad:n1+npad,1-npad:n2+npad,1:n3), &
         INTENT(INOUT) :: f
    REAL(KIND=prec), DIMENSION(n1,n3) :: passl, passr, recvl, recvr
    REAL(KIND=prec), DIMENSION(n1,n3,npad) :: passl2, passr2, recvl2, recvr2

    INTEGER :: m, n, i, status(MPI_STATUS_SIZE,4), mpierr, req(4), nz_min
  
! ON CALLING pass_slabs, f SHOULD HOLD DATA VALUES IN INTERIOR OF
! ARRAY, i.e. f(1:n1,1:n2) SHOULD HOLD NEW VALUES.  THIS ROUTINE WILL
! PAD THIS ARRAY WITH NEIGHBORING VALUES (WHICH MIGHT BE HELD ON OTHER
! PROCESSORS.

    nz_min = MINVAL(nz_proc)
    IF (npad <= nz_min) THEN
          DO n = 1,npad
             DO m = 1,n3
                passl2(:,m,n) = f(1:n1,n,m)
                passr2(:,m,n) = f(1:n1,n2+1-n,m)
             END DO
          END DO

          CALL MPI_IRECV(recvr2,npad*n1*n3,MPI_DOUBLE_PRECISION,&
               right,right,MPI_COMM_WORLD,req(1),mpierr)
          CALL MPI_IRECV(recvl2,npad*n1*n3,MPI_DOUBLE_PRECISION,&
               left,nproc+left,MPI_COMM_WORLD,req(2),mpierr)

          CALL MPI_SEND(passl2,npad*n1*n3,MPI_DOUBLE_PRECISION,&
               left,me,MPI_COMM_WORLD,mpierr)
          CALL MPI_SEND(passr2,npad*n1*n3,MPI_DOUBLE_PRECISION,&
               right,nproc+me,MPI_COMM_WORLD,mpierr)
          DO i = 1,2
             CALL MPI_WAIT(req(i),status(1,i),mpierr)
          END DO

          DO n = 1,npad
             DO m = 1,n3
                f(1:n1,n2+n,m) = recvr2(1:n1,m,n)
                f(1:n1,1-n, m) = recvl2(1:n1,m,n)
             END DO
          END DO

    ELSE

       DO n = 1,npad
          CALL MPI_IRECV(recvr,n1*n3,MPI_DOUBLE_PRECISION,&
               right,right,MPI_COMM_WORLD,req(1),mpierr)
          CALL MPI_IRECV(recvl,n1*n3,MPI_DOUBLE_PRECISION,&
               left,nproc+left,MPI_COMM_WORLD,req(2),mpierr)
          passl = f(1:n1,n,1:n3)
          passr = f(1:n1,n2+1-n,1:n3)
          CALL MPI_SEND(passl,n1*n3,MPI_DOUBLE_PRECISION,&
               left,me,MPI_COMM_WORLD,mpierr)
          CALL MPI_SEND(passr,n1*n3,MPI_DOUBLE_PRECISION,&
               right,nproc+me,MPI_COMM_WORLD,mpierr)
          DO i = 1,2
             CALL MPI_WAIT(req(i),status(1,i),mpierr)
          END DO
          f(1:n1,n2+n,1:n3) = recvr
          f(1:n1,1-n,1:n3)  = recvl
       END DO

    END IF

    DO n = 1,npad
       f(1-n,:,:)  = f(n1+1-n,:,:) ! PAD AT BEGINNING OF x-DIRECTION
       f(n1+n,:,:) = f(n,:,:)      ! PAD AT END OF x-DIRECTION
    END DO

  END SUBROUTINE pass_slabs
END MODULE pass_data





