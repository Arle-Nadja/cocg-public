PROGRAM MAIN 
    USE COCG_Kernel
    IMPLICIT NONE 
    COMPLEX(KIND(0d0)), ALLOCATABLE :: H(:,:), b(:), x(:)
    COMPLEX(KIND(0d0)), PARAMETER :: zI = (0d0, 1d0)
    INTEGER :: N, Nconv
    DOUBLE PRECISION :: delta_conv, dum
    INTEGER :: ix, jx

    N = 100
    delta_conv = 1d-14

    ALLOCATE(H(N,N), b(N), x(N))

    x = 0d0

    OPEN(999, file = "Hseed", status = "old")

    DO ix = 1, N 
        DO jx = 1, N 
            ! H(ix, jx) = DBLE(ix) + DBLE(jx) + (DBLE(ix) - DBLE(jx)) * zI
            ! CALL random_number(dum)
            ! H(ix,jx) = dum
            READ(999,*) dum
            H(ix,jx) = dum
        END DO
        ! H(ix,ix) = H(ix,ix) + N
        b(ix) = DBLE(ix)
    END DO

    CLOSE(999)

    ! H(:,:) = MATMUL(TRANSPOSE(H), H)


    CALL CG_Solve_Normal(N, H, b, delta_conv, x, Nconv)

    OPEN(999, file = "CGtest_output", status = "replace")

    DO ix = 1, N 
        WRITE(999,*) x(ix)
    END DO 

    CLOSE(999)

    PRINT *, "Calculation End."
    PRINT *, "Iteration Number", Nconv

    DEALLOCATE(H, b, x)

    STOP 
END PROGRAM