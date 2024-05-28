MODULE COCG_Kernel 
    USE COCG_Procedure
    IMPLICIT NONE
    CONTAINS 
    SUBROUTINE Shifted_COCG_Solve_Normal(N, H, Nsig, sigs, b, delta_conv, xs, Nconv)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N, Nsig
        COMPLEX(KIND(0d0)), INTENT(in) :: H(N,N)
        COMPLEX(KIND(0d0)), INTENT(in) :: sigs(Nsig)
        COMPLEX(KIND(0d0)), INTENT(in) :: b(N)
        DOUBLE PRECISION, INTENT(in) :: delta_conv 
        COMPLEX(KIND(0d0)), INTENT(out) :: xs(N, Nsig)
        INTEGER, INTENT(out) :: Nconv

        RETURN
    END SUBROUTINE

END MODULE