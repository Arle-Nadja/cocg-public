MODULE COCG_Cores
    USE COCG_Components
    IMPLICIT NONE 
    CONTAINS 
    SUBROUTINE COCG_Step(N, rn_1, pn_1, xn_1, A, alphan_1, betan_1, rn, pn, xn)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N 
        DOUBLE PRECISION, INTENT(in) :: rn_1(N), pn_1(N), xn_1(N), A(N,N)
        DOUBLE PRECISION, INTENT(out) :: alphan_1, betan_1
        DOUBLE PRECISION, INTENT(inout) :: rn(N), pn(N), xn(N)

        CALL Calc_Next_Alpha(N, rn_1, pn_1, A, alphan_1)
        CALL Calc_Next_x(N, pn_1, xn_1, alphan_1, xn)
        CALL Calc_Next_r(N, rn_1, pn_1, A, alphan_1, rn)
        CALL Calc_Next_beta(N, rn, rn_1, betan_1)
        CALL Calc_Next_p(N, betan_1, pn_1, rn, pn)

        RETURN
    END SUBROUTINE

    ! SUBROUTINE Shifted_COCG_Step(N, rn_1, pn_1, xn_1, A, alphan_1, betan_1, rn, pn, xn,&
    !                             )

    ! END SUBROUTINE
    
END MODULE