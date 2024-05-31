MODULE COCG_Components
    IMPLICIT NONE

    CONTAINS 
    SUBROUTINE Calc_Next_Alpha(N, rr, pn_1, Ap, alphan_1)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N 
        COMPLEX(KIND(0d0)), INTENT(in) :: rr, pn_1(N), AP(N)
        COMPLEX(KIND(0d0)), INTENT(out) :: alphan_1

        alphan_1 = rr / dot_product(pn_1, Ap)

        RETURN
    END SUBROUTINE

    SUBROUTINE Calc_Next_x(N, pn_1, xn_1, alphan_1, xn)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N 
        COMPLEX(KIND(0d0)), INTENT(in) :: alphan_1, pn_1(N), xn_1(N) 
        COMPLEX(KIND(0d0)), INTENT(out) :: xn(N)

        xn = xn_1 + alphan_1 * pn_1

        RETURN
    END SUBROUTINE

    SUBROUTINE Calc_Next_r(N, rn_1, Ap, alphan_1, rn)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N 
        COMPLEX(KIND(0d0)), INTENT(in) :: alphan_1, Ap(N), rn_1(N)
        COMPLEX(KIND(0d0)), INTENT(out) :: rn(N)

        rn = rn_1 - alphan_1 * Ap

        RETURN
    END SUBROUTINE

    SUBROUTINE Calc_Next_beta(N, rn, rr, betan_1)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N 
        COMPLEX(KIND(0d0)), INTENT(in) :: rr, rn(N)
        COMPLEX(KIND(0d0)), INTENT(out) :: betan_1

        betan_1 = dot_product(rn, rn) / rr

        RETURN
    END SUBROUTINE

    SUBROUTINE Calc_Next_p(N, betan_1, pn_1, rn, pn)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N 
        COMPLEX(KIND(0d0)), INTENT(in) :: betan_1, pn_1(N), rn(N) 
        COMPLEX(KIND(0d0)), INTENT(out) :: pn(N)

        pn = rn + betan_1 * pn_1

        RETURN 
    END SUBROUTINE 

    SUBROUTINE Calc_Shifted_rn(N, rn, pins, rns)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N
        COMPLEX(KIND(0d0)), INTENT(in) :: rn(N), pins 
        COMPLEX(KIND(0d0)), INTENT(out) :: rns(N)

        rns = rn / pins

        RETURN 
    END SUBROUTINE

    SUBROUTINE Calc_Next_pis(pins, sigma, alphan_1, alphan_2, betan_2, pin_1s, pin_2s)
        IMPLICIT NONE 
        COMPLEX(KIND(0d0)), INTENT(in) :: sigma, alphan_1, alphan_2, betan_2, pin_1s, pin_2s
        COMPLEX(KIND(0d0)), INTENT(out) :: pins 

        pins = (1d0 + alphan_1 * sigma + alphan_1 * betan_2 / alphan_2) * pin_1s &
                    - alphan_1 * betan_2 / alphan_2 * pin_2s
        
        RETURN
    END SUBROUTINE

    SUBROUTINE Calc_Shifted_Alpha(pin_1s, pins, alphan_1, alphan_1s)
        IMPLICIT NONE 
        COMPLEX(KIND(0d0)), INTENT(in) :: pin_1s, pins, alphan_1
        COMPLEX(KIND(0d0)), INTENT(out) :: alphan_1s

        alphan_1s = pin_1s * alphan_1 / pins

        RETURN 
    END SUBROUTINE 

    SUBROUTINE Calc_Shifted_Beta(pin_1s, pins, betan_1, betan_1s)
        IMPLICIT NONE 
        COMPLEX(KIND(0d0)), INTENT(in) :: pin_1s, pins, betan_1
        COMPLEX(KIND(0d0)), INTENT(out) :: betan_1s

        betan_1s = pin_1s / pins * pin_1s / pins * betan_1

        RETURN 
    END SUBROUTINE

END MODULE