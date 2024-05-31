MODULE COCG_Procedure
    USE COCG_Components
    IMPLICIT NONE 
    CONTAINS 
    SUBROUTINE COCG_Step(N, rn_1, pn_1, xn_1, A, Ap, rr, alphan_1, betan_1, rn, pn, xn)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N 
        COMPLEX(KIND(0d0)), INTENT(in) :: rn_1(N), pn_1(N), xn_1(N), A(N,N)
        COMPLEX(KIND(0d0)), INTENT(out) :: alphan_1, betan_1
        COMPLEX(KIND(0d0)), INTENT(inout) :: rn(N), pn(N), xn(N), Ap(N), rr

        rr = dot_product(rn_1, rn_1)
        Ap(:) = MATMUL(A, pn_1)

        CALL Calc_Next_Alpha(N, rr, pn_1, Ap, alphan_1)
        CALL Calc_Next_x(N, pn_1, xn_1, alphan_1, xn)
        CALL Calc_Next_r(N, rn_1, Ap, alphan_1, rn)
        CALL Calc_Next_beta(N, rn, rr, betan_1)
        CALL Calc_Next_p(N, betan_1, pn_1, rn, pn)

        RETURN
    END SUBROUTINE

    SUBROUTINE Shifted_COCG_Step(N, rn, pn_1s, xn_1s, alphan_1, alphan_2, betan_1, betan_2, sigma, pin_1s, pin_2s, &
                                    rns, xns, pns, pins, alphan_1s, betan_1s)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: N 
        COMPLEX(KIND(0d0)), INTENT(in) :: rn(N), alphan_1, alphan_2, betan_1, betan_2, sigma, pin_1s, pin_2s
        COMPLEX(KIND(0d0)), INTENT(in) :: pn_1s(N), xn_1s(N)
        COMPLEX(KIND(0d0)), INTENT(inout) :: rns(N), xns(N), pns(N), pins, alphan_1s, betan_1s
        
        CALL Calc_Next_pis(pins, sigma, alphan_1, alphan_2, betan_2, pin_1s, pin_2s)
        CALL Calc_Shifted_Alpha(pin_1s, pins, alphan_1, alphan_1s)
        CALL Calc_Shifted_Beta(pin_1s, pins, betan_1, betan_1s)

        CALL Calc_Shifted_rn(N, rn, pins, rns)
        CALL Calc_Next_x(N, pn_1s, xn_1s, alphan_1s, xns)
        CALL Calc_Next_p(N, betan_1s, pn_1s, rns, pns)

        RETURN
    END SUBROUTINE
    
END MODULE