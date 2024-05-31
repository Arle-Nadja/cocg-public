MODULE COCG_Kernel 
    USE COCG_Procedure
    IMPLICIT NONE
    CONTAINS 
    SUBROUTINE CG_Solve_Normal(N, H, b, delta_conv, x, Nconv)
        IMPLICIT NONE 
        INTEGER :: N 
        COMPLEX(KIND(0d0)), INTENT(in) :: H(N,N)
        COMPLEX(KIND(0d0)), INTENT(in) :: b(N)
        DOUBLE PRECISION, INTENT(in) :: delta_conv 
        COMPLEX(KIND(0d0)), INTENT(out) :: x(N)
        INTEGER, INTENT(out) :: Nconv
        COMPLEX(KIND(0d0)), ALLOCATABLE :: xn(:), xn_1(:), pn(:), pn_1(:), rn(:), rn_1(:), Ap(:)
        COMPLEX(KIND(0d0)) :: alphan_1, betan_1, rr
        INTEGER :: iter 
        DOUBLE PRECISION :: conv 

        ALLOCATE(xn(N), xn_1(N))
        ALLOCATE(pn(N), pn_1(N))
        ALLOCATE(rn(N), rn_1(N))
        ALLOCATE(Ap(N))

        !! INITIAL SETTINGS
        xn_1(:) = 0d0 
        rn_1(:) = b(:) ; pn_1(:) = b(:)
        alphan_1 = 1d0 ; betan_1 = 0d0
        iter = 0
        conv = 1d10

        DO WHILE(ABS(conv) > delta_conv)
            iter = iter + 1
            CALL COCG_Step(N, rn_1, pn_1, xn_1, H, Ap, rr, alphan_1, betan_1, rn, pn, xn)

            conv = SUM(ABS(rn(:))**2)

            IF(MOD(iter,100) == 0) THEN 
                PRINT *, "Step:", iter, "conv", conv
            END IF

            rn_1(:) = rn(:)
            xn_1(:) = xn(:)
            pn_1(:) = pn(:)
        END DO

        x(:) = xn(:)
        Nconv = iter

        DEALLOCATE(xn, xn_1, pn, pn_1, rn, rn_1, Ap)

        RETURN
    END SUBROUTINE


    SUBROUTINE Shifted_COCG_Solve_Normal(N, H, Nsig, sigs, b, delta_conv, xs, Nconv)
        IMPLICIT NONE 
        INTEGER, INTENT(in) :: N, Nsig
        COMPLEX(KIND(0d0)), INTENT(in) :: H(N,N)
        COMPLEX(KIND(0d0)), INTENT(in) :: sigs(0:Nsig)
        COMPLEX(KIND(0d0)), INTENT(in) :: b(N)
        DOUBLE PRECISION, INTENT(in) :: delta_conv 
        COMPLEX(KIND(0d0)), INTENT(out) :: xs(N, 0:Nsig)
        INTEGER, INTENT(out) :: Nconv
        COMPLEX(KIND(0d0)), ALLOCATABLE :: xn_1s(:,:), pn_1s(:,:), rns(:,:), rn_1(:)
        COMPLEX(KIND(0d0)), ALLOCATABLE :: pns(:,:), xns(:,:)
        COMPLEX(KIND(0d0)), ALLOCATABLE :: alphan_1s(:), alphan_2s(:), betan_1s(:), betan_2s(:)
        COMPLEX(KIND(0d0)), ALLOCATABLE :: pins(:), pin_1s(:), pin_2s(:)
        COMPLEX(KIND(0d0)), ALLOCATABLE :: Ap(:)
        COMPLEX(KIND(0d0)) :: rr
        INTEGER :: isig, iter
        DOUBLE PRECISION, ALLOCATABLE :: convs(:)
        DOUBLE PRECISION :: conv

        ALLOCATE(xn_1s(N,0:Nsig), xns(N,0:Nsig))
        ALLOCATE(pn_1s(N,0:Nsig), pns(N,0:Nsig))
        ALLOCATE(rns(N,0:Nsig), rn_1(N))
        ALLOCATE(Ap(N))
        ALLOCATE(alphan_1s(0:Nsig), alphan_2s(0:Nsig), betan_1s(0:Nsig), betan_2s(0:Nsig))
        ALLOCATE(pins(0:Nsig), pin_1s(0:Nsig), pin_2s(0:Nsig))
        ALLOCATE(convs(0:Nsig))
        

        !! Initial Condition
        xn_1s(:,:) = 0d0 
        rn_1(:) = b(:)
        DO isig = 0, Nsig
            pn_1s(:,isig) = b(:)
        END DO
        alphan_2s(:) = 1d0 ; betan_2s(:) = 0d0
        pin_1s(:) = 1d0 ; pin_2s(:) = 1d0
        conv = 1d10

        iter = 0

        DO WHILE(ABS(conv) > ABS(delta_conv))
            iter = iter + 1
            CALL COCG_Step(N, rn_1(:), pn_1s(:,0), xn_1s(:,0), H, Ap(:), rr, alphan_1s(0), betan_1s(0), &
                            rns(:,0), pns(:,0), xns(:,0))
            
            DO isig = 1, Nsig 
                CALL Shifted_COCG_Step(N, rns(:,0), pn_1s(:,isig), xn_1s(:,isig), alphan_1s(0), alphan_2s(0),&
                        betan_1s(0), betan_2s(0), sigs(isig) - sigs(0), pin_1s(isig), pin_2s(isig),&
                        rns(:,isig), xns(:,isig), pns(:,isig), pins(isig), alphan_1s(isig), betan_1s(isig))
            END DO

            !! Convergence Check
            DO isig = 0, Nsig 
                convs(isig) = SUM(ABS(rns(:,isig)))
            END DO

            conv = SUM(convs)

            IF(MOD(iter,100) == 0) THEN 
                PRINT *, "Step:", iter, "conv", conv
            END IF

            xn_1s(:,:) = xns(:,:)
            pn_1s(:,:) = pns(:,:)
            rn_1(:) = rns(:,0)

            alphan_2s(:) = alphan_1s(:)
            betan_2s(:) = betan_1s(:)

            pin_2s(:) = pin_1s(:)
            pin_1s(:) = pins(:)

        END DO

        xs(:,:) = xns(:,:)
        Nconv = iter

        DEALLOCATE(xn_1s, pn_1s, rns, rn_1)
        DEALLOCATE(alphan_1s, alphan_2s, betan_1s, betan_2s)
        DEALLOCATE(pins, pin_1s, pin_2s)
        DEALLOCATE(Ap)

        RETURN
    END SUBROUTINE

END MODULE