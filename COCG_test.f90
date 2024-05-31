PROGRAM MAIN 
    USE COCG_Kernel
    IMPLICIT NONE 
    COMPLEX(KIND(0d0)), ALLOCATABLE :: H(:,:), b(:), x(:,:), sigs(:)
    COMPLEX(KIND(0d0)), PARAMETER :: zI = (0d0, 1d0)
    INTEGER :: N, Nconv, Nsig
    DOUBLE PRECISION :: delta_conv, dum
    INTEGER :: ix, jx

    N = 100
    Nsig = 10
    delta_conv = 1d-12

    ALLOCATE(H(N,N), b(N), x(N,0:Nsig), sigs(0:Nsig))

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
        H(ix,ix) = H(ix,ix) 
        b(ix) = DBLE(ix)
    END DO

    CLOSE(999)

    !! Set sigma
    DO ix = 0, Nsig
        sigs(ix) = DBLE(ix) * 0.1d0
    END DO 

    ! H(:,:) = MATMUL(TRANSPOSE(H), H)

    CALL Shifted_COCG_Solve_Normal(N, H, Nsig, sigs, b, delta_conv, x, Nconv)
    ! CALL CG_Solve_Normal(N, H, b, delta_conv, x, Nconv)

    OPEN(999, file = "CGtest_output", status = "replace")
    DO jx = 0, Nsig
    DO ix = 1, N 
        WRITE(999,*) x(ix, jx)
    END DO
    WRITE(999,*) ; WRITE(999,*) ; WRITE(999,*) 
    END DO 

    CLOSE(999)

    PRINT *, "Calculation End."
    PRINT *, "Iteration Number", Nconv

    DEALLOCATE(H, b, x, sigs)

    STOP 
END PROGRAM