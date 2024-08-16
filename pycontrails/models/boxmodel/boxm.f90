MODULE BOXM
    USE NETCDF
    IMPLICIT NONE

    INTEGER :: IOSTAT
    REAL :: PI, TIME1, TSTORE, DTS

    INTEGER :: NTS, NCELL, NS
    INTEGER :: NTC, NPC, NPP, NFL, NEMI
    INTEGER :: S, CELL, TS

    INTEGER :: NCID, DIMID_TIME, DIMID_CELL, DIMID_NS, DIMID_NEMI
    INTEGER :: DIMID_NPP, DIMID_NPC, DIMID_NTC, DIMID_NFL
    INTEGER, PRIVATE :: VARID_TIME, VARID_LEVEL, VARID_LON, VARID_LAT, VARID_PRESSURE, VARID_ALT
    INTEGER, PRIVATE :: VARID_SPECIES, VARID_EMI_SPECIES, VARID_BG_CHEM
    INTEGER, PRIVATE :: VARID_TEMP, VARID_M, VARID_H2O, VARID_O2, VARID_N2, VARID_SZA, VARID_EMI
    
    INTEGER :: VARID_Y, VARID_J, VARID_DJ, VARID_RC, VARID_FL
    INTEGER, PRIVATE :: IERR

    ! DEFINE MET INPUTS
    DOUBLE PRECISION, ALLOCATABLE :: TIME(:), LEVEL(:), LON(:), LAT(:), TEMP(:), PRESSURE(:), ALT(:)
    CHARACTER(LEN=80), ALLOCATABLE :: SPECIES(:), EMI_SPECIES(:)
    DOUBLE PRECISION, ALLOCATABLE :: EMI(:,:), EMIP(:,:), BG_CHEM(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: M(:), H2O(:), O2(:), N2(:), SZA(:)

    ! DEFINE CHEM VARIABLES
    DOUBLE PRECISION, ALLOCATABLE :: Y(:,:), YP(:,:), RC(:,:), J(:,:), DJ(:,:), FL(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: SOA(:), MOM(:), BR01(:), RO2(:), P(:), L(:), Y_PPB(:,:), EMI_PPB(:,:)

    ! CHEMCO
    ! SIMPLE RATE COEFFICIENT SCALARS
    DOUBLE PRECISION, PRIVATE :: KRO2NO3,KDEC

    ! COMPLEX RATE COEFFICIENT SCALARS
    DOUBLE PRECISION, PRIVATE :: FCC,FCD,FC1,K2I,FC2,FC7,FC8,K9I,FC9,FC10,KI,K13I
    DOUBLE PRECISION, PRIVATE :: FC13,FC14,FC15,FC16,FCX

    ! SIMPLE RATE COEFFICIENT CELLS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: KRO2NO,KAPNO,KRO2HO2,KAPHO2,KNO3AL,KALKOXY
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: KALKPXY,KIN,KOUT2604,KOUT4608,KOUT2631
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: KOUT2635,KOUT4610,KOUT2605,KOUT2630,KOUT2629
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: KOUT2632,KOUT2637,KOUT3612,KOUT3613,KOUT3442

    ! COMPLEX RATE COEFFICIENT CELLS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: KC0,KCI,KRC,FC,KFPAN,KD0,KDI,KRD,FD,KBPAN,K10
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: K1I,KR1,F1,KMT01,K20,KR2,Fa2,KMT02,K30,K3I
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: KR3,FC3,F3,KMT03,K40,K4I,KR4,FC4,Fa4,KMT04
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: KMT05,KMT06,K70,K7I,KR7,F7,KMT07,K80,K8I,KR8
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: F8,KMT08,K90,KR9,F9,KMT09,K100,K10I,KR10,F10
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: KMT10,K1,K3,K4,K2,KMT11,K0,F,KMT12,K130
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: KR13,F13,KMT13,K140,K14I,KR14,F14,KMT14,K150
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: K15I,KR15,F15,KMT15,K160,K16I,KR16,F16,KMT16
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), PRIVATE :: K170,K17I,KR17,FC17,F17,KMT17

CONTAINS 
    ! PRE INTEGRATION

    SUBROUTINE CHECK(IERR, MSG)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: IERR
        CHARACTER(LEN=*), INTENT(IN) :: MSG
        IF (IERR /= NF90_NOERR) THEN
            PRINT *, MSG
            STOP
        END IF
    END SUBROUTINE CHECK

    SUBROUTINE OPEN_NC
        IMPLICIT NONE
        ! OPEN BOXM INPUT NC
        IERR = NF90_OPEN('/home/ktait98/pycontrails_kt/pycontrails/models/files/boxm_ds.nc', NF90_WRITE, NCID)
        IF (IERR /= NF90_NOERR) THEN
            PRINT *, NF90_STRERROR(IERR)
        END IF
        CALL CHECK(IERR, 'ERROR: CANNOT OPEN BOXM.NC')
    END SUBROUTINE OPEN_NC

    SUBROUTINE GET_DIMS
        IMPLICIT NONE
        ! GET DIMIDS AND LENGTHS
        IERR = NF90_INQ_DIMID(NCID, 'time', DIMID_TIME)
        CALL CHECK(IERR, 'ERROR: CANNOT GET TIME DIMID')
        IERR = NF90_INQUIRE_DIMENSION(NCID, DIMID_TIME, LEN=NTS)
        CALL CHECK(IERR, 'ERROR: CANNOT GET TIME LENGTH')
        IERR = NF90_INQ_DIMID(NCID, 'cell', DIMID_CELL)
        CALL CHECK(IERR, 'ERROR: CANNOT GET CELL DIMID')
        IERR = NF90_INQUIRE_DIMENSION(NCID, DIMID_CELL, LEN=NCELL)
        CALL CHECK(IERR, 'ERROR: CANNOT GET CELL LENGTH')
        IERR = NF90_INQ_DIMID(NCID, 'species', DIMID_NS)
        CALL CHECK(IERR, 'ERROR: CANNOT GET NS DIMID')
        IERR = NF90_INQUIRE_DIMENSION(NCID, DIMID_NS, LEN=NS)
        CALL CHECK(IERR, 'ERROR: CANNOT GET NS LENGTH')
        IERR = NF90_INQ_DIMID(NCID, 'emi_species', DIMID_NEMI)
        CALL CHECK(IERR, 'ERROR: CANNOT GET NEMI DIMID')
        IERR = NF90_INQUIRE_DIMENSION(NCID, DIMID_NEMI, LEN=NEMI)
        CALL CHECK(IERR, 'ERROR: CANNOT GET NEMI LENGTH')
        ! IERR = NF90_INQ_DIMID(NCID, 'photol_params', DIMID_NPP) 
        ! CALL CHECK(IERR, 'ERROR: CANNOT DEFINE NPP DIM')
        ! IERR = NF90_INQUIRE_DIMENSION(NCID, DIMID_NPP, LEN=NPP)
        ! CALL CHECK(IERR, 'ERROR: CANNOT DEFINE NPP DIM')
        ! IERR = NF90_INQ_DIMID(NCID, 'photol_coeffs', DIMID_NPC)
        ! CALL CHECK(IERR, 'ERROR: CANNOT DEFINE NPC DIM')
        ! IERR = NF90_INQUIRE_DIMENSION(NCID, DIMID_NPC, LEN=NPC)
        ! CALL CHECK(IERR, 'ERROR: CANNOT DEFINE NPC DIM')
        ! IERR = NF90_INQ_DIMID(NCID, 'therm_coeffs', DIMID_NTC)
        ! CALL CHECK(IERR, 'ERROR: CANNOT DEFINE NTC DIM')
        ! IERR = NF90_INQUIRE_DIMENSION(NCID, DIMID_NTC, LEN=NTC)
        ! CALL CHECK(IERR, 'ERROR: CANNOT DEFINE NTC DIM')
        ! IERR = NF90_INQ_DIMID(NCID, 'flux_rates', DIMID_NFL)
        ! CALL CHECK(IERR, 'ERROR: CANNOT DEFINE NFL DIM')
        ! IERR = NF90_INQUIRE_DIMENSION(NCID, DIMID_NFL, LEN=NFL)
        ! CALL CHECK(IERR, 'ERROR: CANNOT DEFINE NFL DIM')

        NPP = 57
        NPC = 96
        NTC = 512
        NFL = 130
        
        TIME1 = 0.0
        TSTORE = 0.0
    END SUBROUTINE GET_DIMS

    SUBROUTINE INIT_VARS
        IMPLICIT NONE
        ! ALLOCATE CHEM DATA STRUCTURE
        IF (.NOT. ALLOCATED(Y)) THEN
            ALLOCATE(Y(NCELL,NS))
            Y(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(YP)) THEN
            ALLOCATE(YP(NCELL,NS))
            YP(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(Y_PPB)) THEN
            ALLOCATE(Y_PPB(NCELL,NS))
            Y_PPB(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(RC)) THEN
            ALLOCATE(RC(NCELL,NTC))
            RC(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(J)) THEN
            ALLOCATE(J(NCELL,NPP))
            J(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(DJ)) THEN
            ALLOCATE(DJ(NCELL,NPC))
            DJ(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(FL)) THEN
            ALLOCATE(FL(NCELL,NFL))
            FL(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(EMI)) THEN
            ALLOCATE(EMI(NCELL,NEMI))
            EMI(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(EMIP)) THEN
            ALLOCATE(EMIP(NCELL,NEMI))
            EMIP(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(EMI_PPB)) THEN
            ALLOCATE(EMI_PPB(NCELL,NEMI))
            EMI_PPB(:,:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(SOA)) THEN
            ALLOCATE(SOA(NCELL))
            SOA(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(MOM)) THEN
            ALLOCATE(MOM(NCELL))
            MOM(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(BR01)) THEN
            ALLOCATE(BR01(NCELL))
            BR01(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(RO2)) THEN
            ALLOCATE(RO2(NCELL))
            RO2(:) = 0.0
        END IF

        IF (.NOT. ALLOCATED(P)) THEN
            ALLOCATE(P(NCELL))
            P(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(L)) THEN
            ALLOCATE(L(NCELL))
            L(:) = 0.0
        END IF
                
        IF (.NOT. ALLOCATED(KRO2NO)) THEN
            ALLOCATE(KRO2NO(NCELL))
            KRO2NO(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KAPNO)) THEN
            ALLOCATE(KAPNO(NCELL))
            KAPNO(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KRO2HO2)) THEN
            ALLOCATE(KRO2HO2(NCELL))
            KRO2HO2(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KAPHO2)) THEN
            ALLOCATE(KAPHO2(NCELL))
            KAPHO2(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KNO3AL)) THEN
            ALLOCATE(KNO3AL(NCELL))
            KNO3AL(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KALKOXY)) THEN
            ALLOCATE(KALKOXY(NCELL))
            KALKOXY(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KALKPXY)) THEN
            ALLOCATE(KALKPXY(NCELL))
            KALKPXY(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KIN)) THEN
            ALLOCATE(KIN(NCELL))
            KIN(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT2604)) THEN
            ALLOCATE(KOUT2604(NCELL))
            KOUT2604(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT4608)) THEN
            ALLOCATE(KOUT4608(NCELL))
            KOUT4608(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT2631)) THEN
            ALLOCATE(KOUT2631(NCELL))
            KOUT2631(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT2635)) THEN
            ALLOCATE(KOUT2635(NCELL))
            KOUT2635(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT4610)) THEN
            ALLOCATE(KOUT4610(NCELL))
            KOUT4610(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT2605)) THEN
            ALLOCATE(KOUT2605(NCELL))
            KOUT2605(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT2630)) THEN
            ALLOCATE(KOUT2630(NCELL))
            KOUT2630(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT2629)) THEN
            ALLOCATE(KOUT2629(NCELL))
            KOUT2629(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT2632)) THEN
            ALLOCATE(KOUT2632(NCELL))
            KOUT2632(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT2637)) THEN
            ALLOCATE(KOUT2637(NCELL))
            KOUT2637(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT3612)) THEN
            ALLOCATE(KOUT3612(NCELL))
            KOUT3612(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT3613)) THEN
            ALLOCATE(KOUT3613(NCELL))
            KOUT3613(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KOUT3442)) THEN
            ALLOCATE(KOUT3442(NCELL))
            KOUT3442(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KC0)) THEN
            ALLOCATE(KC0(NCELL))
            KC0(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KCI)) THEN
            ALLOCATE(KCI(NCELL))
            KCI(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KRC)) THEN
            ALLOCATE(KRC(NCELL))
            KRC(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(FC)) THEN
            ALLOCATE(FC(NCELL))
            FC(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KFPAN)) THEN
            ALLOCATE(KFPAN(NCELL))
            KFPAN(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KD0)) THEN
            ALLOCATE(KD0(NCELL))
            KD0(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KDI)) THEN
            ALLOCATE(KDI(NCELL))
            KDI(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KRD)) THEN
            ALLOCATE(KRD(NCELL))
            KRD(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(FD)) THEN
            ALLOCATE(FD(NCELL))
            FD(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KBPAN)) THEN
            ALLOCATE(KBPAN(NCELL))
            KBPAN(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K10)) THEN
            ALLOCATE(K10(NCELL))
            K10(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K1I)) THEN
            ALLOCATE(K1I(NCELL))
            K1I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR1)) THEN
            ALLOCATE(KR1(NCELL))
            KR1(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F1)) THEN
            ALLOCATE(F1(NCELL))
            F1(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT01)) THEN
            ALLOCATE(KMT01(NCELL))
            KMT01(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K20)) THEN
            ALLOCATE(K20(NCELL))
            K20(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR2)) THEN
            ALLOCATE(KR2(NCELL))
            KR2(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(Fa2)) THEN
            ALLOCATE(Fa2(NCELL))
            Fa2(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT02)) THEN
            ALLOCATE(KMT02(NCELL))
            KMT02(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K30)) THEN
            ALLOCATE(K30(NCELL))
            K30(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K3I)) THEN
            ALLOCATE(K3I(NCELL))
            K3I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR3)) THEN
            ALLOCATE(KR3(NCELL))
            KR3(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(FC3)) THEN
            ALLOCATE(FC3(NCELL))
            FC3(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F3)) THEN
            ALLOCATE(F3(NCELL))
            F3(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT03)) THEN
            ALLOCATE(KMT03(NCELL))
            KMT03(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K40)) THEN
            ALLOCATE(K40(NCELL))
            K40(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K4I)) THEN
            ALLOCATE(K4I(NCELL))
            K4I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR4)) THEN
            ALLOCATE(KR4(NCELL))
            KR4(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(FC4)) THEN
            ALLOCATE(FC4(NCELL))
            FC4(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(Fa4)) THEN
            ALLOCATE(Fa4(NCELL))
            Fa4(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT04)) THEN
            ALLOCATE(KMT04(NCELL))
            KMT04(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT05)) THEN
            ALLOCATE(KMT05(NCELL))
            KMT05(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT06)) THEN
            ALLOCATE(KMT06(NCELL))
            KMT06(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K70)) THEN
            ALLOCATE(K70(NCELL))
            K70(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K7I)) THEN
            ALLOCATE(K7I(NCELL))
            K7I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR7)) THEN
            ALLOCATE(KR7(NCELL))
            KR7(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F7)) THEN
            ALLOCATE(F7(NCELL))
            F7(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT07)) THEN
            ALLOCATE(KMT07(NCELL))
            KMT07(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K80)) THEN
            ALLOCATE(K80(NCELL))
            K80(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K8I)) THEN
            ALLOCATE(K8I(NCELL))
            K8I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR8)) THEN
            ALLOCATE(KR8(NCELL))
            KR8(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F8)) THEN
            ALLOCATE(F8(NCELL))
            F8(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT08)) THEN
            ALLOCATE(KMT08(NCELL))
            KMT08(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K90)) THEN
            ALLOCATE(K90(NCELL))
            K90(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR9)) THEN
            ALLOCATE(KR9(NCELL))
            KR9(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F9)) THEN
            ALLOCATE(F9(NCELL))
            F9(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT09)) THEN
            ALLOCATE(KMT09(NCELL))
            KMT09(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K100)) THEN
            ALLOCATE(K100(NCELL))
            K100(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K10I)) THEN
            ALLOCATE(K10I(NCELL))
            K10I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR10)) THEN
            ALLOCATE(KR10(NCELL))
            KR10(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F10)) THEN
            ALLOCATE(F10(NCELL))
            F10(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT10)) THEN
            ALLOCATE(KMT10(NCELL))
            KMT10(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K1)) THEN
            ALLOCATE(K1(NCELL))
            K1(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K3)) THEN
            ALLOCATE(K3(NCELL))
            K3(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K4)) THEN
            ALLOCATE(K4(NCELL))
            K4(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K2)) THEN
            ALLOCATE(K2(NCELL))
            K2(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT11)) THEN
            ALLOCATE(KMT11(NCELL))
            KMT11(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K0)) THEN
            ALLOCATE(K0(NCELL))
            K0(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F)) THEN
            ALLOCATE(F(NCELL))
            F(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT12)) THEN
            ALLOCATE(KMT12(NCELL))
            KMT12(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K130)) THEN
            ALLOCATE(K130(NCELL))
            K130(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR13)) THEN
            ALLOCATE(KR13(NCELL))
            KR13(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F13)) THEN
            ALLOCATE(F13(NCELL))
            F13(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT13)) THEN
            ALLOCATE(KMT13(NCELL))
            KMT13(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K140)) THEN
            ALLOCATE(K140(NCELL))
            K140(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K14I)) THEN
            ALLOCATE(K14I(NCELL))
            K14I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR14)) THEN
            ALLOCATE(KR14(NCELL))
            KR14(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F14)) THEN
            ALLOCATE(F14(NCELL))
            F14(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT14)) THEN
            ALLOCATE(KMT14(NCELL))
            KMT14(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K150)) THEN
            ALLOCATE(K150(NCELL))
            K150(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K15I)) THEN
            ALLOCATE(K15I(NCELL))
            K15I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR15)) THEN
            ALLOCATE(KR15(NCELL))
            KR15(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F15)) THEN
            ALLOCATE(F15(NCELL))
            F15(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT15)) THEN
            ALLOCATE(KMT15(NCELL))
            KMT15(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K160)) THEN
            ALLOCATE(K160(NCELL))
            K160(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K16I)) THEN
            ALLOCATE(K16I(NCELL))
            K16I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR16)) THEN
            ALLOCATE(KR16(NCELL))
            KR16(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F16)) THEN
            ALLOCATE(F16(NCELL))
            F16(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT16)) THEN
            ALLOCATE(KMT16(NCELL))
            KMT16(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K170)) THEN
            ALLOCATE(K170(NCELL))
            K170(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(K17I)) THEN
            ALLOCATE(K17I(NCELL))
            K17I(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KR17)) THEN
            ALLOCATE(KR17(NCELL))
            KR17(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(FC17)) THEN
            ALLOCATE(FC17(NCELL))
            FC17(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(F17)) THEN
            ALLOCATE(F17(NCELL))
            F17(:) = 0.0
        END IF
        IF (.NOT. ALLOCATED(KMT17)) THEN
            ALLOCATE(KMT17(NCELL))
            KMT17(:) = 0.0
        END IF
    
    END SUBROUTINE INIT_VARS

    SUBROUTINE GET_PRE_INT_VARIDS
        IMPLICIT NONE
        ! TIME
        IERR = NF90_INQ_VARID(NCID, 'time', VARID_TIME)
        CALL CHECK(IERR, 'ERROR: CANNOT GET TIME VARID')
        IF (.NOT. ALLOCATED(TIME)) THEN
            ALLOCATE(TIME(NTS))
        END IF
        ! LEVEL
        IERR = NF90_INQ_VARID(NCID, 'level', VARID_LEVEL)
        CALL CHECK(IERR, 'ERROR: CANNOT GET LEVEL VARID')
        IF (.NOT. ALLOCATED(LEVEL)) THEN
            ALLOCATE(LEVEL(NCELL))
        END IF
        ! LON
        IERR = NF90_INQ_VARID(NCID, 'longitude', VARID_LON)
        CALL CHECK(IERR, 'ERROR: CANNOT GET LON VARID')
        IF (.NOT. ALLOCATED(LON)) THEN
            ALLOCATE(LON(NCELL))
        END IF
        ! LAT
        IERR = NF90_INQ_VARID(NCID, 'latitude', VARID_LAT)
        CALL CHECK(IERR, 'ERROR: CANNOT GET LAT VARID')
        IF (.NOT. ALLOCATED(LAT)) THEN
            ALLOCATE(LAT(NCELL))
        END IF
        ! PRESSURE
        IERR = NF90_INQ_VARID(NCID, 'air_pressure', VARID_PRESSURE)
        CALL CHECK(IERR, 'ERROR: CANNOT GET PRESSURE VARID')
        IF (.NOT. ALLOCATED(PRESSURE)) THEN
            ALLOCATE(PRESSURE(NCELL))
        END IF
        ! ALTITUDE
        IERR = NF90_INQ_VARID(NCID, 'altitude', VARID_ALT)
        CALL CHECK(IERR, 'ERROR: CANNOT GET ALT VARID')
        IF (.NOT. ALLOCATED(ALT)) THEN
            ALLOCATE(ALT(NCELL))
        END IF
        ! ! SPECIES
        ! IERR = NF90_INQ_VARID(NCID, 'species', VARID_SPECIES)
        ! CALL CHECK(IERR, 'ERROR: CANNOT GET SPECIES VARID')
        ! IF (.NOT. ALLOCATED(SPECIES)) THEN
        !     ALLOCATE(SPECIES(NS))
        ! END IF
        ! ! EMI_SPECIES
        ! IERR = NF90_INQ_VARID(NCID, 'emi_species', VARID_EMI_SPECIES)
        ! CALL CHECK(IERR, 'ERROR: CANNOT GET EMI_SPECIES VARID')
        ! IF (.NOT. ALLOCATED(EMI_SPECIES)) THEN
        !     ALLOCATE(EMI_SPECIES(NEMI))
        ! END IF

        ! BG_CHEM
        IERR = NF90_INQ_VARID(NCID, 'bg_chem', VARID_BG_CHEM)
        CALL CHECK(IERR, 'ERROR: CANNOT GET BG_CHEM VARID')
        IF (.NOT. ALLOCATED(BG_CHEM)) THEN
            ALLOCATE(BG_CHEM(NCELL,NS))
        END IF
    END SUBROUTINE GET_PRE_INT_VARIDS

    SUBROUTINE GET_INT_VARIDS
        IMPLICIT NONE
        ! TEMP
        IERR = NF90_INQ_VARID(NCID, 'air_temperature', VARID_TEMP)
        CALL CHECK(IERR, 'ERROR: CANNOT GET TEMP VARID')
        IF (.NOT. ALLOCATED(TEMP)) THEN
            ALLOCATE(TEMP(NCELL))
        END IF
        ! M
        IERR = NF90_INQ_VARID(NCID, 'M', VARID_M)
        CALL CHECK(IERR, 'ERROR: CANNOT GET M VARID')
        IF (.NOT. ALLOCATED(M)) THEN
            ALLOCATE(M(NCELL))
        END IF
        ! H2O
        IERR = NF90_INQ_VARID(NCID, 'H2O', VARID_H2O)
        CALL CHECK(IERR, 'ERROR: CANNOT GET H2O VARID')
        IF (.NOT. ALLOCATED(H2O)) THEN
            ALLOCATE(H2O(NCELL))
        END IF
        ! O2
        IERR = NF90_INQ_VARID(NCID, 'O2', VARID_O2)
        CALL CHECK(IERR, 'ERROR: CANNOT GET O2 VARID')
        IF (.NOT. ALLOCATED(O2)) THEN
            ALLOCATE(O2(NCELL))
        END IF
        ! N2
        IERR = NF90_INQ_VARID(NCID, 'N2', VARID_N2)
        CALL CHECK(IERR, 'ERROR: CANNOT GET N2 VARID')
        IF (.NOT. ALLOCATED(N2)) THEN
            ALLOCATE(N2(NCELL))
        END IF
        ! SZA
        IERR = NF90_INQ_VARID(NCID, 'sza', VARID_SZA)
        CALL CHECK(IERR, 'ERROR: CANNOT GET SZA VARID')
        IF (.NOT. ALLOCATED(SZA)) THEN
            ALLOCATE(SZA(NCELL))
        END IF
        ! EMI
        IERR = NF90_INQ_VARID(NCID, 'emi', VARID_EMI)
        CALL CHECK(IERR, 'ERROR: CANNOT GET EMI VARID')
        IF (.NOT. ALLOCATED(EMI)) THEN
            ALLOCATE(EMI(NCELL,NEMI))
        END IF
        ! Y
        IERR = NF90_INQ_VARID(NCID, 'Y', VARID_Y)
        CALL CHECK(IERR, 'ERROR: CANNOT GET Y VARID')
        IF (.NOT. ALLOCATED(Y)) THEN
            ALLOCATE(Y(NCELL,NS))
        END IF
        ! J
        IERR = NF90_INQ_VARID(NCID, 'J', VARID_J)
        CALL CHECK(IERR, 'ERROR: CANNOT GET J VARID')
        IF (.NOT. ALLOCATED(J)) THEN
            ALLOCATE(J(NCELL,NPP))
        END IF
        ! DJ
        IERR = NF90_INQ_VARID(NCID, 'DJ', VARID_DJ)
        CALL CHECK(IERR, 'ERROR: CANNOT GET DJ VARID')
        IF (.NOT. ALLOCATED(DJ)) THEN
            ALLOCATE(DJ(NCELL,NPC))
        END IF
        ! RC
        IERR = NF90_INQ_VARID(NCID, 'RC', VARID_RC)
        CALL CHECK(IERR, 'ERROR: CANNOT GET RC VARID')
        IF (.NOT. ALLOCATED(RC)) THEN
            ALLOCATE(RC(NCELL,NTC))
        END IF
        ! FL
        ! IERR = NF90_INQ_VARID(NCID, 'FL', VARID_FL)
        ! CALL CHECK(IERR, 'ERROR: CANNOT GET FL VARID')
        ! IF (.NOT. ALLOCATED(FL)) THEN
        !     ALLOCATE(FL(NCELL,NFL))
        ! END IF

    END SUBROUTINE GET_INT_VARIDS

    SUBROUTINE GET_PRE_INT_VARS
        IMPLICIT NONE
        ! ALLOCATE AND POPULATE CONSTANT INPUT VARIABLES
        
        IERR = NF90_GET_VAR(NCID, VARID_TIME, TIME)
        CALL CHECK(IERR, 'ERROR: CANNOT GET TIME ARRAY')

        IERR = NF90_GET_VAR(NCID, VARID_LEVEL, LEVEL)
        CALL CHECK(IERR, 'ERROR: CANNOT GET LEVEL ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_LON, LON)
        CALL CHECK(IERR, 'ERROR: CANNOT GET LON ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_LAT, LAT)
        CALL CHECK(IERR, 'ERROR: CANNOT GET LAT ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_PRESSURE, PRESSURE)
        CALL CHECK(IERR, 'ERROR: CANNOT GET PRESSURE ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_ALT, ALT)
        CALL CHECK(IERR, 'ERROR: CANNOT GET ALT ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_BG_CHEM, BG_CHEM, (/1, 1/), (/NCELL, NS/))
        CALL CHECK(IERR, 'ERROR: CANNOT GET BG_CHEM ARRAY')

    END SUBROUTINE GET_PRE_INT_VARS

    ! INTEGRATION
    SUBROUTINE GET_INT_VARS(TS)
        IMPLICIT NONE
        INTEGER :: TS
        INTEGER, DIMENSION(5) :: DIMIDS
        ! ALLOCATE AND POPULATE TIME-DEPENDENT INPUT VARIABLES
        IERR = NF90_GET_VAR(NCID, VARID_TEMP, TEMP, (/1, TS/), (/NCELL, 1/))
        CALL CHECK(IERR, 'ERROR: CANNOT GET TEMP ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_M, M, (/1, TS/), (/NCELL, 1/))
        CALL CHECK(IERR, 'ERROR: CANNOT GET M ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_H2O, H2O, (/1, TS/), (/NCELL, 1/))
        CALL CHECK(IERR, 'ERROR: CANNOT GET H2O ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_O2, O2, (/1, TS/), (/NCELL, 1/)) 
        CALL CHECK(IERR, 'ERROR: CANNOT GET O2 ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_N2, N2, (/1, TS/), (/NCELL, 1/))
        CALL CHECK(IERR, 'ERROR: CANNOT GET N2 ARRAY')
        
        IERR = NF90_GET_VAR(NCID, VARID_SZA, SZA, (/1, TS/), (/NCELL, 1/))
        CALL CHECK(IERR, 'ERROR: CANNOT GET SZA ARRAY')
        
        ! ASSIGN BG_CHEM TO FIRST TIMESTEP OF Y
        IF (TS == 1) THEN
            Y(:,:) = BG_CHEM(:,:) 
            DO S = 1, NS
                Y(:,S) = Y(:,S) * M(:)/1E+09 ! CONVERT FROM PPB TO MOL/CM3 FOR KINETICS
            END DO
        END IF

        ! SET PREVIOUS EMI TO EMIP
        EMIP(:,:) = EMI(:,:)
        IERR = NF90_GET_VAR(NCID, VARID_EMI, EMI, (/1, 1, TS/), (/NCELL, NEMI, 1/))
        CALL CHECK(IERR, 'ERROR: CANNOT GET EMI ARRAY')

        

    END SUBROUTINE GET_INT_VARS

    SUBROUTINE CALC_AEROSOL
        IMPLICIT NONE
        DOUBLE PRECISION :: BGOAM

        BGOAM = 0.7

        SOA(:) = Y(:,204)*3.574E-10 + Y(:,205)*3.574E-10 + &
        Y(:,206)*3.059E-10 + Y(:,207)*3.126E-10 + Y(:,208)*3.093E-10 + &
        Y(:,209)*3.093E-10 + Y(:,210)*3.325E-10 + Y(:,211)*4.072E-10 + &
        Y(:,212)*2.860E-10 + Y(:,213)*3.391E-10 + Y(:,214)*2.310E-10 + &
        Y(:,215)*2.543E-10 + Y(:,216)*1.628E-10 + Y(:,219)*2.493E-10
        
        MOM(:) = Y(:,218) + BGOAM + SOA(:)

    END SUBROUTINE CALC_AEROSOL

    SUBROUTINE CHEMCO
        IMPLICIT NONE
        DOUBLE PRECISION :: R

        R = 8.314
      
        !     SIMPLE RATE COEFFICIENTS                     
                                                                  
        KRO2NO(:)  = 2.70D-12*EXP(360/TEMP(:)) 
        KAPNO(:)   = 7.50D-12*EXP(290/TEMP(:)) 
        KRO2NO3    = 2.30D-12 
        KRO2HO2(:) = 2.91D-13*EXP(1300/TEMP(:)) 
        KAPHO2(:)  = 5.20D-13*EXP(980/TEMP(:)) 
        KNO3AL(:)  = 1.44D-12*EXP(-1862/TEMP(:)) 
        KDEC       = 1.0D+06
        KALKOXY(:) = 3.70D-14*EXP(-460/TEMP(:))*O2(:) 
        KALKPXY(:) = 1.80D-14*EXP(-260/TEMP(:))*O2(:) 
        BR01(:)    = (0.156 + 9.77D+08*EXP(-6415/TEMP(:))) 

        KIN(:)      = 6.2E-03*MOM(:)
        KOUT2604(:) = 4.34*EXP(-7776/(R*TEMP(:)))
        KOUT4608(:) = 4.34*EXP(-9765/(R*TEMP(:)))
        KOUT2631(:) = 4.34*EXP(-14500/(R*TEMP(:)))
        KOUT2635(:) = 4.34*EXP(-12541/(R*TEMP(:)))
        KOUT4610(:) = 4.34*EXP(-10513/(R*TEMP(:)))
        KOUT2605(:) = 4.34*EXP(-8879/(R*TEMP(:)))
        KOUT2630(:) = 4.34*EXP(-12639/(R*TEMP(:)))
        KOUT2629(:) = 4.34*EXP(-4954/(R*TEMP(:)))
        KOUT2632(:) = 4.34*EXP(-3801/(R*TEMP(:)))
        KOUT2637(:) = 4.34*EXP(-16752/(R*TEMP(:)))
        KOUT3612(:) = 4.34*EXP(-8362/(R*TEMP(:)))
        KOUT3613(:) = 4.34*EXP(-11003/(R*TEMP(:)))
        KOUT3442(:) = 4.34*EXP(-12763/(R*TEMP(:)))

        !    COMPLEX RATE COEFFICIENTS                    
                                                                     
        !    KFPAN                                                          
        KC0(:)     = 3.28D-28*M(:)*(TEMP(:)/300)**-7.1 
        KCI(:)     = 1.125D-11*(TEMP(:)/300)**-1.105  
        KRC(:)     = KC0(:)/KCI(:)    
        FCC        = 0.30       
        FC(:)      = 10**(LOG10(FCC)/(1+(LOG10(KRC(:)))**2)) 
        KFPAN(:)   = (KC0(:)*KCI(:))*FC(:)/(KC0(:)+KCI(:)) 

        !    KBPAN                                                   
        KD0(:)     = 1.1D-05*M(:)*EXP(-10100/TEMP(:)) 
        KDI(:)     = 1.9D+17*EXP(-14100/TEMP(:))  
        KRD(:)     = KD0(:)/KDI(:)    
        FCD        = 0.30       
        FD(:)      = 10**(LOG10(FCD)/(1+(LOG10(KRD(:)))**2)) 
        KBPAN(:)   = (KD0(:)*KDI(:))*FD(:)/(KD0(:)+KDI(:)) 
                                                             
        !     KMT01                                                   
        K10(:)     = 9.00D-32*M(:)*(TEMP(:)/300)**-1.5 
        K1I(:)     = 3.00D-11*(TEMP(:)/300)**0.3    
        KR1(:)     = K10(:)/K1I(:)    
        FC1        = 0.6 
        F1(:)      = 10**(LOG10(FC1)/(1+(LOG10(KR1(:)))**2)) 
        KMT01(:)   = (K10(:)*K1I(:))*F1(:)/(K10+K1I(:)) 
                                                                     
        !     KMT02                                                   
        K20(:)     = 9.00D-32*((TEMP(:)/300)**-2.0)*M(:) 
        K2I        = 2.20D-11
        KR2(:)     = K20(:)/K2I    
        FC2        = 0.6 
        Fa2(:)     = 10**(LOG10(FC2)/(1+(LOG10(KR2(:)))**2)) 
        KMT02(:)   = (K20(:)*K2I)*Fa2(:)/(K20(:)+K2I) 
                                                                     
        !      KMT03  : NO2      + NO3     = N2O5                               
        !    IUPAC 2001                                                       
        K30(:)     = 2.70D-30*M(:)*(TEMP(:)/300)**-3.4 
        K3I(:)     = 2.00D-12*(TEMP(:)/300)**0.2    
        KR3(:)     = K30(:)/K3I(:)   
        FC3(:)     = (EXP(-TEMP(:)/250) + EXP(-1050/TEMP(:))) 
        F3(:)      = 10**(LOG10(FC3(:))/(1+(LOG10(KR3(:)))**2)) 
        KMT03(:)   = (K30(:)*K3I(:))*F3(:)/(K30(:)+K3I(:)) 
                                                             
        !     KMT04  : N2O5               = NO2     + NO3                     
        ! IUPAC 1997/2001                                                 
        K40(:)     = (2.20D-03*M(:)*(TEMP(:)/300)**-4.34)*(EXP(-11080/TEMP(:)))
        K4I(:)     = (9.70D+14*(TEMP(:)/300)**0.1)*EXP(-11080/TEMP(:))    
        KR4(:)     = K40(:)/K4I(:)    
        FC4(:)     = (EXP(-TEMP(:)/250) + EXP(-1050/TEMP(:)))
        Fa4(:)     = 10**(LOG10(FC4(:))/(1+(LOG10(KR4(:)))**2)) 
        KMT04(:)   = (K40(:)*K4I(:))*Fa4(:)/(K40(:)+K4I(:))       
                                                           
        !       KMT05                                                   
        KMT05(:)   = 1 + ((0.6*M(:))/(2.687D+19*(273/TEMP(:)))) 
                                                                     
        !    KMT06                                                   
        KMT06(:)   = 1 + (1.40D-21*EXP(2200/TEMP(:))*H2O(:)) 
                                                                     
        !    KMT07  : OH       + NO      = HONO                              
        !    IUPAC 2001                                                      
        K70(:)     = 7.00D-31*M(:)*(TEMP(:)/300)**-2.6 
        K7I(:)     = 3.60D-11*(TEMP(:)/300)**0.1    
        KR7(:)     = K70(:)/K7I(:)  
        FC7        = 0.6  
        F7(:)      = 10**(LOG10(FC7)/(1+(LOG10(KR7(:)))**2)) 
        KMT07(:)   = (K70(:)*K7I(:))*F7(:)/(K70(:)+K7I(:)) 
                                                                     
        ! NASA 2000                                                           
        !    KMT08                                                    
        K80(:)     = 2.50D-30*((TEMP(:)/300)**-4.4)*M(:) 
        K8I(:)     = 1.60D-11 
        KR8(:)     = K80(:)/K8I(:)
        FC8        = 0.6 
        F8(:)      = 10**(LOG10(FC8)/(1+(LOG10(KR8(:)))**2)) 
        KMT08(:)   = (K80(:)*K8I(:))*F8(:)/(K80(:)+K8I(:)) 
                                                                     
        !    KMT09  : HO2      + NO2     = HO2NO2                            
        !    IUPAC 1997/2001                                                 
        K90(:)     = 1.80D-31*M(:)*(TEMP(:)/300)**-3.2 
        K9I        = 4.70D-12    
        KR9(:)     = K90(:)/K9I    
        FC9        = 0.6 
        F9(:)      = 10**(LOG10(FC9)/(1+(LOG10(KR9(:)))**2)) 
        KMT09(:)   = (K90(:)*K9I)*F9(:)/(K90(:)+K9I) 
                                                                     
        ! KMT10  : HO2NO2             = HO2     + NO2                     
        ! IUPAC 2001                                                      

        K100(:)    = 4.10D-05*M(:)*EXP(-10650/TEMP(:)) 
        K10I(:)    = 5.70D+15*EXP(-11170/TEMP(:))   
        KR10(:)    = K100(:)/K10I(:)    
        FC10       = 0.5 
        F10(:)     = 10**(LOG10(FC10)/(1+(LOG10(KR10(:)))**2)) 
        KMT10(:)   = (K100(:)*K10I(:))*F10(:)/(K100(:)+K10I(:)) 
                                                                     
        !   KMT11  : OH       + HNO3    = H2O + NO3                     
        !   IUPAC 2001                                                      
        K1(:)      = 7.20D-15*EXP(785/TEMP(:)) 
        K3(:)      = 1.90D-33*EXP(725/TEMP(:)) 
        K4(:)      = 4.10D-16*EXP(1440/TEMP(:)) 
        K2(:)      = (K3(:)*M(:))/(1+(K3(:)*M(:)/K4(:))) 
        KMT11(:)   = K1(:) + K2(:) 
                                                              
        ! KMT12 : OH    +   SO2  =  HSO3                                  
        ! IUPAC 2003                                                      
        K0(:)      = 3.0D-31*((TEMP(:)/300)**-3.3)*M(:) 
        KI         = 1.5D-12 
        KR1(:)     = K0(:)/KI 
        FCX        = 0.6 
        F(:)       = 10**(LOG10(FCX)/(1+(LOG10(KR1(:)))**2)) 
        KMT12(:)   = (K0(:)*KI*F(:))/(K0(:)+KI) 
                                                              
        ! KMT13  : CH3O2    + NO2     = CH3O2NO2                           
        ! IUPAC 2003                                                       
        K130(:)     = 1.20D-30*((TEMP(:)/300)**-6.9)*M(:) 
        K13I        = 1.80D-11 
        KR13(:)     = K130(:)/K13I 
        FC13        = 0.36 
        F13(:)      = 10**(LOG10(FC13)/(1+(LOG10(KR13(:)))**2)) 
        KMT13(:)    = (K130(:)*K13I)*F13(:)/(K130(:)+K13I) 
                                                                     
        !  KMT14  : CH3O2NO2           = CH3O2   + NO2                      
        !  IUPAC 2001                                                       
        K140(:)     = 9.00D-05*EXP(-9690/TEMP(:))*M(:) 
        K14I(:)     = 1.10D+16*EXP(-10560/TEMP(:)) 
        KR14(:)     = K140(:)/K14I(:)
        FC14        = 0.36 
        F14(:)      = 10**(LOG10(FC14)/(1+(LOG10(KR14(:)))**2)) 
        KMT14(:)    = (K140(:)*K14I(:))*F14(:)/(K140(:)+K14I(:)) 
                                                          
        ! KMT15  :    OH  +  C2H4  =                                       
        ! IUPAC 2001                                                      
        K150(:)     = 6.00D-29*((TEMP(:)/298)**-4.0)*M(:) 
        K15I(:)     = 9.00D-12*((TEMP(:)/298)**-1.1) 
        KR15(:)     = K150(:)/K15I(:)
        FC15        = 0.7
        F15(:)      = 10**(LOG10(FC15)/(1+(LOG10(KR15(:)))**2)) 
        KMT15(:)    = (K150(:)*K15I(:))*F15(:)/(K150(:)+K15I(:)) 
                                                             
        ! KMT16  :  OH  +  C3H6         
        ! IUPAC 2003                                                     
        K160(:)     = 3.00D-27*((TEMP(:)/298)**-3.0)*M(:) 
        K16I(:)     = 2.80D-11*((TEMP(:)/298)**-1.3) 
        KR16(:)     = K160(:)/K16I(:) 
        FC16        = 0.5 
        F16(:)      = 10**(LOG10(FC16)/(1+(LOG10(KR16(:)))**2)) 
        KMT16(:)    = (K160(:)*K16I(:))*F16(:)/(K160(:)+K16I(:)) 

        ! KMT17                                                   
        K170(:)     = 5.00D-30*((TEMP(:)/298)**-1.5)*M(:) 
        K17I(:)     = 9.40D-12*EXP(-700/TEMP(:)) 
        KR17(:)     = K170(:)/K17I(:)
        FC17(:)     = (EXP(-TEMP(:)/580) + EXP(-2320/TEMP(:))) 
        F17(:)      = 10**(LOG10(FC17(:))/(1+(LOG10(KR17(:)))**2)) 
        KMT17(:)    = (K170(:) * K17I(:)) * F17(:)/(K170(:) + K17I(:)) 

        ! LIST OF ALL REACTIONS
        ! Reaction (1) O = O3                                                             
        RC(:,1) = 5.60D-34*O2(:)*N2(:)*((TEMP(:)/300)**-2.6)

        ! Reaction (2) O = O3                                                             
        RC(:,2) = 6.00D-34*O2(:)*O2(:)*((TEMP(:)/300)**-2.6)

        ! Reaction (3) O + O3 =                                                           
        RC(:,3) = 8.00D-12*EXP(-2060/TEMP(:))         

        ! Reaction (4) O + NO = NO2                                                       
        RC(:,4) = KMT01(:)                          

        ! Reaction (5) O + NO2 = NO                                                       
        RC(:,5) = 5.50D-12*EXP(188/TEMP(:))           

        ! Reaction (6) O + NO2 = NO3                                                      
        RC(:,6) = KMT02(:)                            

        ! Reaction (7) O1D = O                                                            
        RC(:,7) = 3.20D-11*O2(:)*EXP(67/TEMP(:))         

        ! Reaction (8) O1D = O                                                            
        RC(:,8) = 1.80D-11*N2(:)*EXP(107/TEMP(:))        

        ! Reaction (9) NO + O3 = NO2                                                      
        RC(:,9) = 1.40D-12*EXP(-1310/TEMP(:))         

        ! Reaction (10) NO2 + O3 = NO3                                                     
        RC(:,10) = 1.40D-13*EXP(-2470/TEMP(:))         

        ! Reaction (11) NO + NO = NO2 + NO2                                                
        RC(:,11) = 3.30D-39*EXP(530/TEMP(:))*O2(:)        

        ! Reaction (12) NO + NO3 = NO2 + NO2                                               
        RC(:,12) = 1.80D-11*EXP(110/TEMP(:))           

        ! Reaction (13) NO2 + NO3 = NO + NO2                                               
        RC(:,13) = 4.50D-14*EXP(-1260/TEMP(:))         

        ! Reaction (14) NO2 + NO3 = N2O5                                                   
        RC(:,14) = KMT03(:)                            

        ! Reaction (15) N2O5 = NO2 + NO3                                                   
        RC(:,15) = KMT04(:)                            

        ! Reaction (16) O1D = OH + OH                                                      
        RC(:,16) = 2.20D-10                     

        ! Reaction (17) OH + O3 = HO2                                                      
        RC(:,17) = 1.70D-12*EXP(-940/TEMP(:))          

        ! Reaction (18) OH + H2 = HO2                                                      
        RC(:,18) = 7.70D-12*EXP(-2100/TEMP(:))         

        ! Reaction (19) OH + CO = HO2                                                      
        RC(:,19) = 1.30D-13*KMT05(:)                   

        ! Reaction (20) OH + H2O2 = HO2                                                    
        RC(:,20) = 2.90D-12*EXP(-160/TEMP(:))          

        ! Reaction (21) HO2 + O3 = OH                                                      
        RC(:,21) = 2.03D-16*((TEMP(:)/300)**4.57)*EXP(693/TEMP(:))  

        ! Reaction (22) OH + HO2 =                                                         
        RC(:,22) = 4.80D-11*EXP(250/TEMP(:))           

        ! Reaction (23) HO2 + HO2 = H2O2                                                   
        RC(:,23) = 2.20D-13*KMT06(:)*EXP(600/TEMP(:))     

        ! Reaction (24) HO2 + HO2 = H2O2                                                   
        RC(:,24) = 1.90D-33*M(:)*KMT06(:)*EXP(980/TEMP(:))   

        ! Reaction (25) OH + NO = HONO                                                     
        RC(:,25) = KMT07(:)                            

        ! Reaction (26) NO2 = HONO                                                         
        RC(:,26) = 5.0D-07                          

        ! Reaction (27) OH + NO2 = HNO3                                                    
        RC(:,27) = KMT08(:)                            

        ! Reaction (28) OH + NO3 = HO2 + NO2                                               
        RC(:,28) = 2.00D-11                         

        ! Reaction (29) HO2 + NO = OH + NO2                                                
        RC(:,29) = 3.60D-12*EXP(270/TEMP(:))           

        ! Reaction (30) HO2 + NO2 = HO2NO2                                                 
        RC(:,30) = KMT09(:)                            

        ! Reaction (31) HO2NO2 = HO2 + NO2                                                 
        RC(:,31) = KMT10(:)                            

        ! Reaction (32) OH + HO2NO2 = NO2                                                  
        RC(:,32) = 1.90D-12*EXP(270/TEMP(:))           

        ! Reaction (33) HO2 + NO3 = OH + NO2                                               
        RC(:,33) = 4.00D-12                         

        ! Reaction (34) OH + HONO = NO2                                                    
        RC(:,34) = 2.50D-12*EXP(-260/TEMP(:))          

        ! Reaction (35) OH + HNO3 = NO3                                                    
        RC(:,35) = KMT11(:)                            

        ! Reaction (36) O + SO2 = SO3                                                      
        RC(:,36) = 4.00D-32*EXP(-1000/TEMP(:))*M(:)       

        ! Reaction (37) OH + SO2 = HSO3                                                    
        RC(:,37) = KMT12(:)                            

        ! Reaction (38) HSO3 = HO2 + SO3                                                   
        RC(:,38) = 1.30D-12*EXP(-330/TEMP(:))*O2(:)       

        ! Reaction (39) HNO3 = NA                                                          
        RC(:,39) = 6.00D-06                         

        ! Reaction (40) N2O5 = NA + NA                                                     
        RC(:,40) = 4.00D-05                       

        ! Reaction (41) SO3 = SA                                                           
        RC(:,41) = 1.20D-15*H2O(:)                     

        ! Reaction (42) OH + CH4 = CH3O2                                                   
        RC(:,42) = 9.65D-20*TEMP(:)**2.58*EXP(-1082/TEMP(:)) 

        ! Reaction (43) OH + C2H6 = C2H5O2                                                 
        RC(:,43) = 1.52D-17*TEMP(:)**2*EXP(-498/TEMP(:)) 

        ! Reaction (44) OH + C3H8 = IC3H7O2                                                
        RC(:,44) = 1.55D-17*TEMP(:)**2*EXP(-61/TEMP(:))*0.736  

        ! Reaction (45) OH + C3H8 = RN10O2                                                 
        RC(:,45) = 1.55D-17*TEMP(:)**2*EXP(-61/TEMP(:))*0.264  

        ! Reaction (46) OH + NC4H10 = RN13O2                                               
        RC(:,46) = 1.69D-17*TEMP(:)**2*EXP(145/TEMP(:))  

        ! Reaction (47) OH + C2H4 = HOCH2CH2O2                                             
        RC(:,47) = KMT15(:)                        

        ! Reaction (48) OH + C3H6 = RN9O2                                                  
        RC(:,48) = KMT16(:)                        

        ! Reaction (49) OH + TBUT2ENE = RN12O2                                             
        RC(:,49) = 1.01D-11*EXP(550/TEMP(:))       

        ! Reaction (50) NO3 + C2H4 = NRN6O2                                                
        RC(:,50) = 2.10D-16                     

        ! Reaction (51) NO3 + C3H6 = NRN9O2                                                
        RC(:,51) = 9.40D-15                     

        ! Reaction (52) NO3 + TBUT2ENE = NRN12O2                                           
        RC(:,52) = 3.90D-13                     

        ! Reaction (53) O3 + C2H4 = HCHO + CO + HO2 + OH                                   
        RC(:,53) = 9.14D-15*EXP(-2580/TEMP(:))*0.13  

        ! Reaction (54) O3 + C2H4 = HCHO + HCOOH                                           
        RC(:,54) = 9.14D-15*EXP(-2580/TEMP(:))*0.87  

        ! Reaction (55) O3 + C3H6 = HCHO + CO + CH3O2 + OH                                 
        RC(:,55) = 5.51D-15*EXP(-1878/TEMP(:))*0.36  

        ! Reaction (56) O3 + C3H6 = HCHO + CH3CO2H                                         
        RC(:,56) = 5.51D-15*EXP(-1878/TEMP(:))*0.64  

        ! Reaction (57) O3 + TBUT2ENE = CH3CHO + CO + CH3O2 + OH                           
        RC(:,57) = 6.64D-15*EXP(-1059/TEMP(:))*0.69 

        ! Reaction (58) O3 + TBUT2ENE = CH3CHO + CH3CO2H                                   
        RC(:,58) = 6.64D-15*EXP(-1059/TEMP(:))*0.31 

        ! Reaction (59) OH + C5H8 = RU14O2                                                 
        RC(:,59) = 2.70D-11*EXP(390/TEMP(:))       

        ! Reaction (60) NO3 + C5H8 = NRU14O2                                               
        RC(:,60) = 3.15D-12*EXP(-450/TEMP(:))      

        ! Reaction (61) O3 + C5H8 = UCARB10 + CO + HO2 + OH                                
        RC(:,61) = 1.03D-14*EXP(-1995/TEMP(:))*0.27 

        ! Reaction (62) O3 + C5H8 = UCARB10 + HCOOH                                        
        RC(:,62) = 1.03D-14*EXP(-1995/TEMP(:))*0.73 

        ! Reaction (63) APINENE + OH = RTN28O2                                             
        RC(:,63) = 1.20D-11*EXP(444/TEMP(:))           

        ! Reaction (64) APINENE + NO3 = NRTN28O2                                           
        RC(:,64) = 1.19D-12*EXP(490/TEMP(:))           

        ! Reaction (65) APINENE + O3 = OH + RTN26O2                                        
        RC(:,65) = 1.01D-15*EXP(-732/TEMP(:))*0.80  

        ! Reaction (66) APINENE + O3 = TNCARB26 + H2O2                                     
        RC(:,66) = 1.01D-15*EXP(-732/TEMP(:))*0.075  

        ! Reaction (67) APINENE + O3 = RCOOH25                                             
        RC(:,67) = 1.01D-15*EXP(-732/TEMP(:))*0.125  

        ! Reaction (68) BPINENE + OH = RTX28O2                                             
        RC(:,68) = 2.38D-11*EXP(357/TEMP(:)) 

        ! Reaction (69) BPINENE + NO3 = NRTX28O2                                           
        RC(:,69) = 2.51D-12 

        ! Reaction (70) BPINENE + O3 =  RTX24O2 + OH                                       
        RC(:,70) = 1.50D-17*0.35 

        ! Reaction (71) BPINENE + O3 =  HCHO + TXCARB24 + H2O2                             
        RC(:,71) = 1.50D-17*0.20 

        ! Reaction (72) BPINENE + O3 =  HCHO + TXCARB22                                    
        RC(:,72) = 1.50D-17*0.25 

        ! Reaction (73) BPINENE + O3 =  TXCARB24 + CO                                      
        RC(:,73) = 1.50D-17*0.20 

        ! Reaction (74) C2H2 + OH = HCOOH + CO + HO2                                       
        RC(:,74) = KMT17(:)*0.364 

        ! Reaction (75) C2H2 + OH = CARB3 + OH                                             
        RC(:,75) = KMT17(:)*0.636 

        ! Reaction (76) BENZENE + OH = RA13O2                                              
        RC(:,76) = 2.33D-12*EXP(-193/TEMP(:))*0.47 

        ! Reaction (77) BENZENE + OH = AROH14 + HO2                                        
        RC(:,77) = 2.33D-12*EXP(-193/TEMP(:))*0.53 

        ! Reaction (78) TOLUENE + OH = RA16O2                                              
        RC(:,78) = 1.81D-12*EXP(338/TEMP(:))*0.82 

        ! Reaction (79) TOLUENE + OH = AROH17 + HO2                                        
        RC(:,79) = 1.81D-12*EXP(338/TEMP(:))*0.18 

        ! Reaction (80) OXYL + OH = RA19AO2                                                
        RC(:,80) = 1.36D-11*0.70 

        ! Reaction (81) OXYL + OH = RA19CO2                                                
        RC(:,81) = 1.36D-11*0.30 

        ! Reaction (82) OH + HCHO = HO2 + CO                                               
        RC(:,82) = 1.20D-14*TEMP(:)*EXP(287/TEMP(:))  

        ! Reaction (83) OH + CH3CHO = CH3CO3                                               
        RC(:,83) = 5.55D-12*EXP(311/TEMP(:))             

        ! Reaction (84) OH + C2H5CHO = C2H5CO3                                             
        RC(:,84) = 1.96D-11                                

        ! Reaction (85) NO3 + HCHO = HO2 + CO + HNO3                                       
        RC(:,85) = 5.80D-16                  

        ! Reaction (86) NO3 + CH3CHO = CH3CO3 + HNO3                                       
        RC(:,86) = KNO3AL(:)                   

        ! Reaction (87) NO3 + C2H5CHO = C2H5CO3 + HNO3                                     
        RC(:,87) = KNO3AL(:)*2.4             

        ! Reaction (88) OH + CH3COCH3 = RN8O2                                              
        RC(:,88) = 5.34D-18*TEMP(:)**2*EXP(-230/TEMP(:)) 

        ! Reaction (89) MEK + OH = RN11O2                                                  
        RC(:,89) = 3.24D-18*TEMP(:)**2*EXP(414/TEMP(:))

        ! Reaction (90) OH + CH3OH = HO2 + HCHO                                            
        RC(:,90) = 6.01D-18*TEMP(:)**2*EXP(170/TEMP(:))  

        ! Reaction (91) OH + C2H5OH = CH3CHO + HO2                                         
        RC(:,91) = 6.18D-18*TEMP(:)**2*EXP(532/TEMP(:))*0.887 

        ! Reaction (92) OH + C2H5OH = HOCH2CH2O2                                           
        RC(:,92) = 6.18D-18*TEMP(:)**2*EXP(532/TEMP(:))*0.113 

        ! Reaction (93) NPROPOL + OH = C2H5CHO + HO2                                       
        RC(:,93) = 5.53D-12*0.49 

        ! Reaction (94) NPROPOL + OH = RN9O2                                               
        RC(:,94) = 5.53D-12*0.51 

        ! Reaction (95) OH + IPROPOL = CH3COCH3 + HO2                                      
        RC(:,95) = 4.06D-18*TEMP(:)**2*EXP(788/TEMP(:))*0.86 

        ! Reaction (96) OH + IPROPOL = RN9O2                                               
        RC(:,96) = 4.06D-18*TEMP(:)**2*EXP(788/TEMP(:))*0.14 

        ! Reaction (97) HCOOH + OH = HO2                                                   
        RC(:,97) = 4.50D-13 

        ! Reaction (98) CH3CO2H + OH = CH3O2                                               
        RC(:,98) = 8.00D-13 

        ! Reaction (99) OH + CH3CL = CH3O2                                                 
        RC(:,99) = 7.33D-18*TEMP(:)**2*EXP(-809/TEMP(:))   

        ! Reaction (100) OH + CH2CL2 = CH3O2                                                
        RC(:,100) = 6.14D-18*TEMP(:)**2*EXP(-389/TEMP(:))   

        ! Reaction (101) OH + CHCL3 = CH3O2                                                 
        RC(:,101) = 1.80D-18*TEMP(:)**2*EXP(-129/TEMP(:))   

        ! Reaction (102) OH + CH3CCL3 = C2H5O2                                              
        RC(:,102) = 2.25D-18*TEMP(:)**2*EXP(-910/TEMP(:))   

        ! Reaction (103) OH + TCE = HOCH2CH2O2                                              
        RC(:,103) = 9.64D-12*EXP(-1209/TEMP(:))         

        ! Reaction (104) OH + TRICLETH = HOCH2CH2O2                                         
        RC(:,104) = 5.63D-13*EXP(427/TEMP(:))            

        ! Reaction (105) OH + CDICLETH = HOCH2CH2O2                                         
        RC(:,105) = 1.94D-12*EXP(90/TEMP(:))            

        ! Reaction (106) OH + TDICLETH = HOCH2CH2O2                                         
        RC(:,106) = 1.01D-12*EXP(250/TEMP(:))           

        ! Reaction (107) CH3O2 + NO = HCHO + HO2 + NO2                                      
        RC(:,107) = 3.00D-12*EXP(280/TEMP(:))*0.999 

        ! Reaction (108) C2H5O2 + NO = CH3CHO + HO2 + NO2                                   
        RC(:,108) = 2.60D-12*EXP(365/TEMP(:))*0.991 

        ! Reaction (109) RN10O2 + NO = C2H5CHO + HO2 + NO2                                  
        RC(:,109) = 2.80D-12*EXP(360/TEMP(:))*0.980 

        ! Reaction (110) IC3H7O2 + NO = CH3COCH3 + HO2 + NO2                                
        RC(:,110) = 2.70D-12*EXP(360/TEMP(:))*0.958 

        ! Reaction (111) RN13O2 + NO = CH3CHO + C2H5O2 + NO2                                
        RC(:,111) = KRO2NO(:)*0.917*BR01(:)       

        ! Reaction (112) RN13O2 + NO = CARB11A + HO2 + NO2                                  
        RC(:,112) = KRO2NO(:)*0.917*(1-BR01(:))   

        ! Reaction (113) RN16O2 + NO = RN15AO2 + NO2                                        
        RC(:,113) = KRO2NO(:)*0.877                 

        ! Reaction (114) RN19O2 + NO = RN18AO2 + NO2                                        
        RC(:,114) = KRO2NO(:)*0.788                 

        ! Reaction (115) RN13AO2 + NO = RN12O2 + NO2                                        
        RC(:,115) = KRO2NO(:)                       

        ! Reaction (116) RN16AO2 + NO = RN15O2 + NO2                                        
        RC(:,116) = KRO2NO(:)                       

        ! Reaction (117) RA13O2 + NO = CARB3 + UDCARB8 + HO2 + NO2                          
        RC(:,117) = KRO2NO(:)*0.918       

        ! Reaction (118) RA16O2 + NO = CARB3 + UDCARB11 + HO2 + NO2                         
        RC(:,118) = KRO2NO(:)*0.889*0.7 

        ! Reaction (119) RA16O2 + NO = CARB6 + UDCARB8 + HO2 + NO2                          
        RC(:,119) = KRO2NO(:)*0.889*0.3 

        ! Reaction (120) RA19AO2 + NO = CARB3 + UDCARB14 + HO2 + NO2                        
        RC(:,120) = KRO2NO(:)*0.862       

        ! Reaction (121) RA19CO2 + NO = CARB9 + UDCARB8 + HO2 + NO2                         
        RC(:,121) = KRO2NO(:)*0.862       

        ! Reaction (122) HOCH2CH2O2 + NO = HCHO + HCHO + HO2 + NO2                          
        RC(:,122) = KRO2NO(:)*0.995*0.776  

        ! Reaction (123) HOCH2CH2O2 + NO = HOCH2CHO + HO2 + NO2                             
        RC(:,123) = KRO2NO(:)*0.995*0.224  

        ! Reaction (124) RN9O2 + NO = CH3CHO + HCHO + HO2 + NO2                             
        RC(:,124) = KRO2NO(:)*0.979     

        ! Reaction (125) RN12O2 + NO = CH3CHO + CH3CHO + HO2 + NO2                          
        RC(:,125) = KRO2NO(:)*0.959     

        ! Reaction (126) RN15O2 + NO = C2H5CHO + CH3CHO + HO2 + NO2                         
        RC(:,126) = KRO2NO(:)*0.936     

        ! Reaction (127) RN18O2 + NO = C2H5CHO + C2H5CHO + HO2 + NO2                        
        RC(:,127) = KRO2NO(:)*0.903     

        ! Reaction (128) RN15AO2 + NO = CARB13 + HO2 + NO2                                  
        RC(:,128) = KRO2NO(:)*0.975     

        ! Reaction (129) RN18AO2 + NO = CARB16 + HO2 + NO2                                  
        RC(:,129) = KRO2NO(:)*0.946     

        ! Reaction (130) CH3CO3 + NO = CH3O2 + NO2                                          
        RC(:,130) = KAPNO(:)                      

        ! Reaction (131) C2H5CO3 + NO = C2H5O2 + NO2                                        
        RC(:,131) = KAPNO(:)                      

        ! Reaction (132) HOCH2CO3 + NO = HO2 + HCHO + NO2                                   
        RC(:,132) = KAPNO(:)                      

        ! Reaction (133) RN8O2 + NO = CH3CO3 + HCHO + NO2                                   
        RC(:,133) = KRO2NO(:)                     

        ! Reaction (134) RN11O2 + NO = CH3CO3 + CH3CHO + NO2                                
        RC(:,134) = KRO2NO(:)                     

        ! Reaction (135) RN14O2 + NO = C2H5CO3 + CH3CHO + NO2                               
        RC(:,135) = KRO2NO(:)                     

        ! Reaction (136) RN17O2 + NO = RN16AO2 + NO2                                        
        RC(:,136) = KRO2NO(:)                     

        ! Reaction (137) RU14O2 + NO = UCARB12 + HO2 +  NO2                                 
        RC(:,137) = KRO2NO(:)*0.900*0.252  

        ! Reaction (138) RU14O2 + NO = UCARB10 + HCHO + HO2 + NO2                           
        RC(:,138) = KRO2NO(:)*0.900*0.748 

        ! Reaction (139) RU12O2 + NO = CH3CO3 + HOCH2CHO + NO2                              
        RC(:,139) = KRO2NO(:)*0.7         

        ! Reaction (140) RU12O2 + NO = CARB7 + CO + HO2 + NO2                               
        RC(:,140) = KRO2NO(:)*0.3         

        ! Reaction (141) RU10O2 + NO = CH3CO3 + HOCH2CHO + NO2                              
        RC(:,141) = KRO2NO(:)*0.670      

        ! Reaction (142) RU10O2 + NO = CARB6 + HCHO + HO2 + NO2                             
        RC(:,142) = KRO2NO(:)*0.295         

        ! Reaction (143) RU10O2 + NO = CARB7 + HCHO + HO2 + NO2                             
        RC(:,143) = KRO2NO(:)*0.035          

        ! Reaction (144) NRN6O2 + NO = HCHO + HCHO + NO2 + NO2                              
        RC(:,144) = KRO2NO(:)                 

        ! Reaction (145) NRN9O2 + NO = CH3CHO + HCHO + NO2 + NO2                            
        RC(:,145) = KRO2NO(:)                 

        ! Reaction (146) NRN12O2 + NO = CH3CHO + CH3CHO + NO2 + NO2                         
        RC(:,146) = KRO2NO(:)                 

        ! Reaction (147) NRU14O2 + NO = NUCARB12 + HO2 + NO2                                
        RC(:,147) = KRO2NO(:)                 

        ! Reaction (148) NRU12O2 + NO = NOA + CO + HO2 + NO2                                
        RC(:,148) = KRO2NO(:)                 

        ! Reaction (149) RTN28O2 + NO = TNCARB26 + HO2 + NO2                                
        RC(:,149) = KRO2NO(:)*0.767*0.915  

        ! Reaction (150) RTN28O2 + NO = CH3COCH3 + RN19O2 + NO2                             
        RC(:,150) = KRO2NO(:)*0.767*0.085  

        ! Reaction (151) NRTN28O2 + NO = TNCARB26 + NO2 + NO2                               
        RC(:,151) = KRO2NO(:)                  

        ! Reaction (152) RTN26O2 + NO = RTN25O2 + NO2                                       
        RC(:,152) = KAPNO(:)                   

        ! Reaction (153) RTN25O2 + NO = RTN24O2 + NO2                                       
        RC(:,153) = KRO2NO(:)*0.840        

        ! Reaction (154) RTN24O2 + NO = RTN23O2 + NO2                                       
        RC(:,154) = KRO2NO(:)                   

        ! Reaction (155) RTN23O2 + NO = CH3COCH3 + RTN14O2 + NO2                            
        RC(:,155) = KRO2NO(:)                  

        ! Reaction (156) RTN14O2 + NO = HCHO + TNCARB10 + HO2 + NO2                         
        RC(:,156) = KRO2NO(:)               

        ! Reaction (157) RTN10O2 + NO = RN8O2 + CO + NO2                                    
        RC(:,157) = KRO2NO(:)               

        ! Reaction (158) RTX28O2 + NO = TXCARB24 + HCHO + HO2 + NO2                         
        RC(:,158) = KRO2NO(:)*0.767*0.915  

        ! Reaction (159) RTX28O2 + NO = CH3COCH3 + RN19O2 + NO2                             
        RC(:,159) = KRO2NO(:)*0.767*0.085  

        ! Reaction (160) NRTX28O2 + NO = TXCARB24 + HCHO + NO2 + NO2                        
        RC(:,160) = KRO2NO(:)            

        ! Reaction (161) RTX24O2 + NO = TXCARB22 + HO2 + NO2                                
        RC(:,161) = KRO2NO(:)*0.843*0.6  

        ! Reaction (162) RTX24O2 + NO = CH3COCH3 + RN13AO2 + HCHO + NO2                     
        RC(:,162) = KRO2NO(:)*0.843*0.4  

        ! Reaction (163) RTX22O2 + NO = CH3COCH3 + RN13O2 + NO2                             
        RC(:,163) = KRO2NO(:)*0.700         

        ! Reaction (164) CH3O2    + NO2     = CH3O2NO2                                      
        RC(:,164) = KMT13(:)         

        ! Reaction (165) CH3O2NO2           = CH3O2   + NO2                                 
        RC(:,165) = KMT14(:)         

        ! Reaction (166) CH3O2 + NO = CH3NO3                                                
        RC(:,166) = 3.00D-12*EXP(280/TEMP(:))*0.001 

        ! Reaction (167) C2H5O2 + NO = C2H5NO3                                              
        RC(:,167) = 2.60D-12*EXP(365/TEMP(:))*0.009 

        ! Reaction (168) RN10O2 + NO = RN10NO3                                              
        RC(:,168) = 2.80D-12*EXP(360/TEMP(:))*0.020 

        ! Reaction (169) IC3H7O2 + NO = IC3H7NO3                                            
        RC(:,169) = 2.70D-12*EXP(360/TEMP(:))*0.042 

        ! Reaction (170) RN13O2 + NO = RN13NO3                                              
        RC(:,170) = KRO2NO(:)*0.083                 

        ! Reaction (171) RN16O2 + NO = RN16NO3                                              
        RC(:,171) = KRO2NO(:)*0.123                 

        ! Reaction (172) RN19O2 + NO = RN19NO3                                              
        RC(:,172) = KRO2NO(:)*0.212                 

        ! Reaction (173) HOCH2CH2O2 + NO = HOC2H4NO3                                        
        RC(:,173) = KRO2NO(:)*0.005                 

        ! Reaction (174) RN9O2 + NO = RN9NO3                                                
        RC(:,174) = KRO2NO(:)*0.021                 

        ! Reaction (175) RN12O2 + NO = RN12NO3                                              
        RC(:,175) = KRO2NO(:)*0.041                 

        ! Reaction (176) RN15O2 + NO = RN15NO3                                              
        RC(:,176) = KRO2NO(:)*0.064                 

        ! Reaction (177) RN18O2 + NO = RN18NO3                                              
        RC(:,177) = KRO2NO(:)*0.097                 

        ! Reaction (178) RN15AO2 + NO = RN15NO3                                             
        RC(:,178) = KRO2NO(:)*0.025                 

        ! Reaction (179) RN18AO2 + NO = RN18NO3                                             
        RC(:,179) = KRO2NO(:)*0.054                 

        ! Reaction (180) RU14O2 + NO = RU14NO3                                              
        RC(:,180) = KRO2NO(:)*0.100                 

        ! Reaction (181) RA13O2 + NO = RA13NO3                                              
        RC(:,181) = KRO2NO(:)*0.082                 

        ! Reaction (182) RA16O2 + NO = RA16NO3                                              
        RC(:,182) = KRO2NO(:)*0.111                 

        ! Reaction (183) RA19AO2 + NO = RA19NO3                                             
        RC(:,183) = KRO2NO(:)*0.138                 

        ! Reaction (184) RA19CO2 + NO = RA19NO3                                             
        RC(:,184) = KRO2NO(:)*0.138                 

        ! Reaction (185) RTN28O2 + NO = RTN28NO3                                            
        RC(:,185) = KRO2NO(:)*0.233        

        ! Reaction (186) RTN25O2 + NO = RTN25NO3                                            
        RC(:,186) = KRO2NO(:)*0.160        

        ! Reaction (187) RTX28O2 + NO = RTX28NO3                                            
        RC(:,187) = KRO2NO(:)*0.233        

        ! Reaction (188) RTX24O2 + NO = RTX24NO3                                            
        RC(:,188) = KRO2NO(:)*0.157        

        ! Reaction (189) RTX22O2 + NO = RTX22NO3                                            
        RC(:,189) = KRO2NO(:)*0.300        

        ! Reaction (190) CH3O2 + NO3 = HCHO + HO2 + NO2                                     
        RC(:,190) = KRO2NO3*0.40          

        ! Reaction (191) C2H5O2 + NO3 = CH3CHO + HO2 + NO2                                  
        RC(:,191) = KRO2NO3               

        ! Reaction (192) RN10O2 + NO3 = C2H5CHO + HO2 + NO2                                 
        RC(:,192) = KRO2NO3               

        ! Reaction (193) IC3H7O2 + NO3 = CH3COCH3 + HO2 + NO2                               
        RC(:,193) = KRO2NO3               

        ! Reaction (194) RN13O2 + NO3 = CH3CHO + C2H5O2 + NO2                               
        RC(:,194) = KRO2NO3*BR01(:)     

        ! Reaction (195) RN13O2 + NO3 = CARB11A + HO2 + NO2                                 
        RC(:,195) = KRO2NO3*(1-BR01(:)) 

        ! Reaction (196) RN16O2 + NO3 = RN15AO2 + NO2                                       
        RC(:,196) = KRO2NO3               

        ! Reaction (197) RN19O2 + NO3 = RN18AO2 + NO2                                       
        RC(:,197) = KRO2NO3               

        ! Reaction (198) RN13AO2 + NO3 = RN12O2 + NO2                                       
        RC(:,198) = KRO2NO3                      

        ! Reaction (199) RN16AO2 + NO3 = RN15O2 + NO2                                       
        RC(:,199) = KRO2NO3                      

        ! Reaction (200) RA13O2 + NO3 = CARB3 + UDCARB8 + HO2 + NO2                         
        RC(:,200) = KRO2NO3            

        ! Reaction (201) RA16O2 + NO3 = CARB3 + UDCARB11 + HO2 + NO2                        
        RC(:,201) = KRO2NO3*0.7     

        ! Reaction (202) RA16O2 + NO3 = CARB6 + UDCARB8 + HO2 + NO2                         
        RC(:,202) = KRO2NO3*0.3     

        ! Reaction (203) RA19AO2 + NO3 = CARB3 + UDCARB14 + HO2 + NO2                       
        RC(:,203) = KRO2NO3           

        ! Reaction (204) RA19CO2 + NO3 = CARB9 + UDCARB8 + HO2 + NO2                        
        RC(:,204) = KRO2NO3           

        ! Reaction (205) HOCH2CH2O2 + NO3 = HCHO + HCHO + HO2 + NO2                         
        RC(:,205) = KRO2NO3*0.776  

        ! Reaction (206) HOCH2CH2O2 + NO3 = HOCH2CHO + HO2 + NO2                            
        RC(:,206) = KRO2NO3*0.224  

        ! Reaction (207) RN9O2 + NO3 = CH3CHO + HCHO + HO2 + NO2                            
        RC(:,207) = KRO2NO3         

        ! Reaction (208) RN12O2 + NO3 = CH3CHO + CH3CHO + HO2 + NO2                         
        RC(:,208) = KRO2NO3         

        ! Reaction (209) RN15O2 + NO3 = C2H5CHO + CH3CHO + HO2 + NO2                        
        RC(:,209) = KRO2NO3         

        ! Reaction (210) RN18O2 + NO3 = C2H5CHO + C2H5CHO + HO2 + NO2                       
        RC(:,210) = KRO2NO3         

        ! Reaction (211) RN15AO2 + NO3 = CARB13 + HO2 + NO2                                 
        RC(:,211) = KRO2NO3         

        ! Reaction (212) RN18AO2 + NO3 = CARB16 + HO2 + NO2                                 
        RC(:,212) = KRO2NO3         

        ! Reaction (213) CH3CO3 + NO3 = CH3O2 + NO2                                         
        RC(:,213) = KRO2NO3*1.60          

        ! Reaction (214) C2H5CO3 + NO3 = C2H5O2 + NO2                                       
        RC(:,214) = KRO2NO3*1.60          

        ! Reaction (215) HOCH2CO3 + NO3 = HO2 + HCHO + NO2                                  
        RC(:,215) = KRO2NO3*1.60         

        ! Reaction (216) RN8O2 + NO3 = CH3CO3 + HCHO + NO2                                  
        RC(:,216) = KRO2NO3               

        ! Reaction (217) RN11O2 + NO3 = CH3CO3 + CH3CHO + NO2                               
        RC(:,217) = KRO2NO3               

        ! Reaction (218) RN14O2 + NO3 = C2H5CO3 + CH3CHO + NO2                              
        RC(:,218) = KRO2NO3               

        ! Reaction (219) RN17O2 + NO3 = RN16AO2 + NO2                                       
        RC(:,219) = KRO2NO3               

        ! Reaction (220) RU14O2 + NO3 = UCARB12 + HO2 + NO2                                 
        RC(:,220) = KRO2NO3*0.032     

        ! Reaction (221) RU14O2 + NO3 = UCARB10 + HCHO + HO2 + NO2                          
        RC(:,221) = KRO2NO3*0.968     

        ! Reaction (222) RU12O2 + NO3 = CH3CO3 + HOCH2CHO + NO2                             
        RC(:,222) = KRO2NO3*0.7         

        ! Reaction (223) RU12O2 + NO3 = CARB7 + CO + HO2 + NO2                              
        RC(:,223) = KRO2NO3*0.3         

        ! Reaction (224) RU10O2 + NO3 = CH3CO3 + HOCH2CHO + NO2                             
        RC(:,224) = KRO2NO3*0.7         

        ! Reaction (225) RU10O2 + NO3 = CARB6 + HCHO + HO2 + NO2                            
        RC(:,225) = KRO2NO3*0.3         

        ! Reaction (226) RU10O2 + NO3 = CARB7 + HCHO + HO2 + NO2                            
        RC(:,226) = KRO2NO3*0.0        

        ! Reaction (227) NRN6O2 + NO3 = HCHO + HCHO + NO2 + NO2                             
        RC(:,227) = KRO2NO3               

        ! Reaction (228) NRN9O2 + NO3 = CH3CHO + HCHO + NO2 + NO2                           
        RC(:,228) = KRO2NO3               

        ! Reaction (229) NRN12O2 + NO3 = CH3CHO + CH3CHO + NO2 + NO2                        
        RC(:,229) = KRO2NO3               

        ! Reaction (230) NRU14O2 + NO3 = NUCARB12 + HO2 + NO2                               
        RC(:,230) = KRO2NO3               

        ! Reaction (231) NRU12O2 + NO3 = NOA + CO + HO2 + NO2                               
        RC(:,231) = KRO2NO3               

        ! Reaction (232) RTN28O2 + NO3 = TNCARB26 + HO2 + NO2                               
        RC(:,232) = KRO2NO3                

        ! Reaction (233) NRTN28O2 + NO3 = TNCARB26 + NO2 + NO2                              
        RC(:,233) = KRO2NO3                

        ! Reaction (234) RTN26O2 + NO3 = RTN25O2 + NO2                                      
        RC(:,234) = KRO2NO3*1.60                   

        ! Reaction (235) RTN25O2 + NO3 = RTN24O2 + NO2                                      
        RC(:,235) = KRO2NO3                 

        ! Reaction (236) RTN24O2 + NO3 = RTN23O2 + NO2                                      
        RC(:,236) = KRO2NO3                   

        ! Reaction (237) RTN23O2 + NO3 = CH3COCH3 + RTN14O2 + NO2                           
        RC(:,237) = KRO2NO3                 

        ! Reaction (238) RTN14O2 + NO3 = HCHO + TNCARB10 + HO2 + NO2                        
        RC(:,238) = KRO2NO3             

        ! Reaction (239) RTN10O2 + NO3 = RN8O2 + CO + NO2                                   
        RC(:,239) = KRO2NO3               

        ! Reaction (240) RTX28O2 + NO3 = TXCARB24 + HCHO + HO2 + NO2                        
        RC(:,240) = KRO2NO3             

        ! Reaction (241) RTX24O2 + NO3 = TXCARB22 + HO2 + NO2                               
        RC(:,241) = KRO2NO3             

        ! Reaction (242) RTX22O2 + NO3 = CH3COCH3 + RN13O2 + NO2                            
        RC(:,242) = KRO2NO3             

        ! Reaction (243) NRTX28O2 + NO3 = TXCARB24 + HCHO + NO2 + NO2                       
        RC(:,243) = KRO2NO3            

        ! Reaction (244) CH3O2 + HO2 = CH3OOH                                               
        RC(:,244) = 4.10D-13*EXP(790/TEMP(:))  

        ! Reaction (245) C2H5O2 + HO2 = C2H5OOH                                             
        RC(:,245) = 7.50D-13*EXP(700/TEMP(:))  

        ! Reaction (246) RN10O2 + HO2 = RN10OOH                                             
        RC(:,246) = KRO2HO2(:)*0.520           

        ! Reaction (247) IC3H7O2 + HO2 = IC3H7OOH                                           
        RC(:,247) = KRO2HO2(:)*0.520           

        ! Reaction (248) RN13O2 + HO2 = RN13OOH                                             
        RC(:,248) = KRO2HO2(:)*0.625           

        ! Reaction (249) RN16O2 + HO2 = RN16OOH                                             
        RC(:,249) = KRO2HO2(:)*0.706           

        ! Reaction (250) RN19O2 + HO2 = RN19OOH                                             
        RC(:,250) = KRO2HO2(:)*0.770           

        ! Reaction (251) RN13AO2 + HO2 = RN13OOH                                            
        RC(:,251) = KRO2HO2(:)*0.625           

        ! Reaction (252) RN16AO2 + HO2 = RN16OOH                                            
        RC(:,252) = KRO2HO2(:)*0.706           

        ! Reaction (253) RA13O2 + HO2 = RA13OOH                                             
        RC(:,253) = KRO2HO2(:)*0.770           

        ! Reaction (254) RA16O2 + HO2 = RA16OOH                                             
        RC(:,254) = KRO2HO2(:)*0.820           

        ! Reaction (255) RA19AO2 + HO2 = RA19OOH                                            
        RC(:,255) = KRO2HO2(:)*0.859           

        ! Reaction (256) RA19CO2 + HO2 = RA19OOH                                            
        RC(:,256) = KRO2HO2(:)*0.859           

        ! Reaction (257) HOCH2CH2O2 + HO2 = HOC2H4OOH                                       
        RC(:,257) = 2.03D-13*EXP(1250/TEMP(:)) 

        ! Reaction (258) RN9O2 + HO2 = RN9OOH                                               
        RC(:,258) = KRO2HO2(:)*0.520           

        ! Reaction (259) RN12O2 + HO2 = RN12OOH                                             
        RC(:,259) = KRO2HO2(:)*0.625           

        ! Reaction (260) RN15O2 + HO2 = RN15OOH                                             
        RC(:,260) = KRO2HO2(:)*0.706           

        ! Reaction (261) RN18O2 + HO2 = RN18OOH                                             
        RC(:,261) = KRO2HO2(:)*0.770           

        ! Reaction (262) RN15AO2 + HO2 = RN15OOH                                            
        RC(:,262) = KRO2HO2(:)*0.706           

        ! Reaction (263) RN18AO2 + HO2 = RN18OOH                                            
        RC(:,263) = KRO2HO2(:)*0.770           

        ! Reaction (264) CH3CO3 + HO2 = CH3CO3H                                             
        RC(:,264) = KAPHO2(:)*0.560                  

        ! Reaction (265) C2H5CO3 + HO2 = C2H5CO3H                                           
        RC(:,265) = KAPHO2(:)*0.560                  

        ! Reaction (266) HOCH2CO3 + HO2 = HOCH2CO3H                                         
        RC(:,266) = KAPHO2(:)*0.560                 

        ! Reaction (267) RN8O2 + HO2 = RN8OOH                                               
        RC(:,267) = KRO2HO2(:)*0.520           

        ! Reaction (268) RN11O2 + HO2 = RN11OOH                                             
        RC(:,268) = KRO2HO2(:)*0.625           

        ! Reaction (269) RN14O2 + HO2 = RN14OOH                                             
        RC(:,269) = KRO2HO2(:)*0.706           

        ! Reaction (270) RN17O2 + HO2 = RN17OOH                                             
        RC(:,270) = KRO2HO2(:)*0.770           

        ! Reaction (271) RU14O2 + HO2 = RU14OOH                                             
        RC(:,271) = KRO2HO2(:)*0.770           

        ! Reaction (272) RU12O2 + HO2 = RU12OOH                                             
        RC(:,272) = KRO2HO2(:)*0.706           

        ! Reaction (273) RU10O2 + HO2 = RU10OOH                                             
        RC(:,273) = KRO2HO2(:)*0.625           

        ! Reaction (274) NRN6O2 + HO2 = NRN6OOH                                             
        RC(:,274) = KRO2HO2(:)*0.387         

        ! Reaction (275) NRN9O2 + HO2 = NRN9OOH                                             
        RC(:,275) = KRO2HO2(:)*0.520         

        ! Reaction (276) NRN12O2 + HO2 = NRN12OOH                                           
        RC(:,276) = KRO2HO2(:)*0.625         

        ! Reaction (277) NRU14O2 + HO2 = NRU14OOH                                           
        RC(:,277) = KRO2HO2(:)*0.706         

        ! Reaction (278) NRU12O2 + HO2 = NRU12OOH                                           
        RC(:,278) = KRO2HO2(:)*0.706        

        ! Reaction (279) RTN28O2 + HO2 = RTN28OOH                                           
        RC(:,279) = KRO2HO2(:)*0.914         

        ! Reaction (280) NRTN28O2 + HO2 = NRTN28OOH                                         
        RC(:,280) = KRO2HO2(:)*0.914         

        ! Reaction (281) RTN26O2 + HO2 = RTN26OOH                                           
        RC(:,281) = KAPHO2(:)*0.56                     

        ! Reaction (282) RTN25O2 + HO2 = RTN25OOH                                           
        RC(:,282) = KRO2HO2(:)*0.890       

        ! Reaction (283) RTN24O2 + HO2 = RTN24OOH                                           
        RC(:,283) = KRO2HO2(:)*0.890       

        ! Reaction (284) RTN23O2 + HO2 = RTN23OOH                                           
        RC(:,284) = KRO2HO2(:)*0.890       

        ! Reaction (285) RTN14O2 + HO2 = RTN14OOH                                           
        RC(:,285) = KRO2HO2(:)*0.770       

        ! Reaction (286) RTN10O2 + HO2 = RTN10OOH                                           
        RC(:,286) = KRO2HO2(:)*0.706       

        ! Reaction (287) RTX28O2 + HO2 = RTX28OOH                                           
        RC(:,287) = KRO2HO2(:)*0.914       

        ! Reaction (288) RTX24O2 + HO2 = RTX24OOH                                           
        RC(:,288) = KRO2HO2(:)*0.890       

        ! Reaction (289) RTX22O2 + HO2 = RTX22OOH                                           
        RC(:,289) = KRO2HO2(:)*0.890       

        ! Reaction (290) NRTX28O2 + HO2 = NRTX28OOH                                         
        RC(:,290) = KRO2HO2(:)*0.914       

        ! Reaction (291) CH3O2 = HCHO + HO2                                                 
        RC(:,291) = 1.82D-13*EXP(416/TEMP(:))*0.33*RO2(:)  

        ! Reaction (292) CH3O2 = HCHO                                                       
        RC(:,292) = 1.82D-13*EXP(416/TEMP(:))*0.335*RO2(:) 

        ! Reaction (293) CH3O2 = CH3OH                                                      
        RC(:,293) = 1.82D-13*EXP(416/TEMP(:))*0.335*RO2(:) 

        ! Reaction (294) C2H5O2 = CH3CHO + HO2                                              
        RC(:,294) = 3.10D-13*0.6*RO2(:)             

        ! Reaction (295) C2H5O2 = CH3CHO                                                    
        RC(:,295) = 3.10D-13*0.2*RO2(:)             

        ! Reaction (296) C2H5O2 = C2H5OH                                                    
        RC(:,296) = 3.10D-13*0.2*RO2(:)             

        ! Reaction (297) RN10O2 = C2H5CHO + HO2                                             
        RC(:,297) = 6.00D-13*0.6*RO2(:)             

        ! Reaction (298) RN10O2 = C2H5CHO                                                   
        RC(:,298) = 6.00D-13*0.2*RO2(:)             

        ! Reaction (299) RN10O2 = NPROPOL                                                   
        RC(:,299) = 6.00D-13*0.2*RO2(:)             

        ! Reaction (300) IC3H7O2 = CH3COCH3 + HO2                                           
        RC(:,300) = 4.00D-14*0.6*RO2(:)             

        ! Reaction (301) IC3H7O2 = CH3COCH3                                                 
        RC(:,301) = 4.00D-14*0.2*RO2(:)             

        ! Reaction (302) IC3H7O2 = IPROPOL                                                  
        RC(:,302) = 4.00D-14*0.2*RO2(:)             

        ! Reaction (303) RN13O2 = CH3CHO + C2H5O2                                           
        RC(:,303) = 2.50D-13*RO2(:)*BR01(:)       

        ! Reaction (304) RN13O2 = CARB11A + HO2                                             
        RC(:,304) = 2.50D-13*RO2(:)*(1-BR01(:))   

        ! Reaction (305) RN13AO2 = RN12O2                                                   
        RC(:,305) = 8.80D-13*RO2(:)                 

        ! Reaction (306) RN16AO2 = RN15O2                                                   
        RC(:,306) = 8.80D-13*RO2(:)                 

        ! Reaction (307) RA13O2 = CARB3 + UDCARB8 + HO2                                     
        RC(:,307) = 8.80D-13*RO2(:)                 

        ! Reaction (308) RA16O2 = CARB3 + UDCARB11 + HO2                                    
        RC(:,308) = 8.80D-13*RO2(:)*0.7          

        ! Reaction (309) RA16O2 = CARB6 + UDCARB8 + HO2                                     
        RC(:,309) = 8.80D-13*RO2(:)*0.3          

        ! Reaction (310) RA19AO2 = CARB3 + UDCARB14 + HO2                                   
        RC(:,310) = 8.80D-13*RO2(:)                 

        ! Reaction (311) RA19CO2 = CARB3 + UDCARB14 + HO2                                   
        RC(:,311) = 8.80D-13*RO2(:)                 

        ! Reaction (312) RN16O2 = RN15AO2                                                   
        RC(:,312) = 2.50D-13*RO2(:)                 

        ! Reaction (313) RN19O2 = RN18AO2                                                   
        RC(:,313) = 2.50D-13*RO2(:)                 

        ! Reaction (314) HOCH2CH2O2 = HCHO + HCHO + HO2                                     
        RC(:,314) = 2.00D-12*RO2(:)*0.776       

        ! Reaction (315) HOCH2CH2O2 = HOCH2CHO + HO2                                        
        RC(:,315) = 2.00D-12*RO2(:)*0.224       

        ! Reaction (316) RN9O2 = CH3CHO + HCHO + HO2                                        
        RC(:,316) = 8.80D-13*RO2(:)                 

        ! Reaction (317) RN12O2 = CH3CHO + CH3CHO + HO2                                     
        RC(:,317) = 8.80D-13*RO2(:)                 

        ! Reaction (318) RN15O2 = C2H5CHO + CH3CHO + HO2                                    
        RC(:,318) = 8.80D-13*RO2(:)                 

        ! Reaction (319) RN18O2 = C2H5CHO + C2H5CHO + HO2                                   
        RC(:,319) = 8.80D-13*RO2(:)                 

        ! Reaction (320) RN15AO2 = CARB13 + HO2                                             
        RC(:,320) = 8.80D-13*RO2(:)                 

        ! Reaction (321) RN18AO2 = CARB16 + HO2                                             
        RC(:,321) = 8.80D-13*RO2(:)                 

        ! Reaction (322) CH3CO3 = CH3O2                                                     
        RC(:,322) = 1.00D-11*RO2(:)                 

        ! Reaction (323) C2H5CO3 = C2H5O2                                                   
        RC(:,323) = 1.00D-11*RO2(:)                 

        ! Reaction (324) HOCH2CO3 = HCHO + HO2                                              
        RC(:,324) = 1.00D-11*RO2(:)                 

        ! Reaction (325) RN8O2 = CH3CO3 + HCHO                                              
        RC(:,325) = 1.40D-12*RO2(:)                 

        ! Reaction (326) RN11O2 = CH3CO3 + CH3CHO                                           
        RC(:,326) = 1.40D-12*RO2(:)                 

        ! Reaction (327) RN14O2 = C2H5CO3 + CH3CHO                                          
        RC(:,327) = 1.40D-12*RO2(:)                 

        ! Reaction (328) RN17O2 = RN16AO2                                                   
        RC(:,328) = 1.40D-12*RO2(:)                 

        ! Reaction (329) RU14O2 = UCARB12 + HO2                                             
        RC(:,329) = 1.26D-12*RO2(:)*0.1        

        ! Reaction (330) RU14O2 = UCARB10 + HCHO + HO2                                      
        RC(:,330) = 1.26D-12*RO2(:)*0.9       

        ! Reaction (331) RU12O2 = CH3CO3 + HOCH2CHO                                         
        RC(:,331) = 4.20D-13*RO2(:)*0.7           

        ! Reaction (332) RU12O2 = CARB7 + HOCH2CHO + HO2                                    
        RC(:,332) = 4.20D-13*RO2(:)*0.3            

        ! Reaction (333) RU10O2 = CH3CO3 + HOCH2CHO                                         
        RC(:,333) = 1.83D-12*RO2(:)*0.7            

        ! Reaction (334) RU10O2 = CARB6 + HCHO + HO2                                        
        RC(:,334) = 1.83D-12*RO2(:)*0.3            

        ! Reaction (335) RU10O2 = CARB7 + HCHO + HO2                                        
        RC(:,335) = 1.83D-12*RO2(:)*0.0             

        ! Reaction (336) NRN6O2 = HCHO + HCHO + NO2                                         
        RC(:,336) = 6.00D-13*RO2(:)                 

        ! Reaction (337) NRN9O2 = CH3CHO + HCHO + NO2                                       
        RC(:,337) = 2.30D-13*RO2(:)                 

        ! Reaction (338) NRN12O2 = CH3CHO + CH3CHO + NO2                                    
        RC(:,338) = 2.50D-13*RO2(:)                 

        ! Reaction (339) NRU14O2 = NUCARB12 + HO2                                           
        RC(:,339) = 1.30D-12*RO2(:)                 

        ! Reaction (340) NRU12O2 = NOA + CO + HO2                                           
        RC(:,340) = 9.60D-13*RO2(:)                 

        ! Reaction (341) RTN28O2 = TNCARB26 + HO2                                           
        RC(:,341) = 2.85D-13*RO2(:)                 

        ! Reaction (342) NRTN28O2 = TNCARB26 + NO2                                          
        RC(:,342) = 1.00D-13*RO2(:)                 

        ! Reaction (343) RTN26O2 = RTN25O2                                                  
        RC(:,343) = 1.00D-11*RO2(:)                   

        ! Reaction (344) RTN25O2 = RTN24O2                                                  
        RC(:,344) = 1.30D-12*RO2(:)           

        ! Reaction (345) RTN24O2 = RTN23O2                                                  
        RC(:,345) = 6.70D-15*RO2(:)             

        ! Reaction (346) RTN23O2 = CH3COCH3 + RTN14O2                                       
        RC(:,346) = 6.70D-15*RO2(:)            

        ! Reaction (347) RTN14O2 = HCHO + TNCARB10 + HO2                                    
        RC(:,347) = 8.80D-13*RO2(:)        

        ! Reaction (348) RTN10O2 = RN8O2 + CO                                               
        RC(:,348) = 2.00D-12*RO2(:)        

        ! Reaction (349) RTX28O2 = TXCARB24 + HCHO + HO2                                    
        RC(:,349) = 2.00D-12*RO2(:)       

        ! Reaction (350) RTX24O2 = TXCARB22 + HO2                                           
        RC(:,350) = 2.50D-13*RO2(:)       

        ! Reaction (351) RTX22O2 = CH3COCH3 + RN13O2                                        
        RC(:,351) = 2.50D-13*RO2(:)       

        ! Reaction (352) NRTX28O2 = TXCARB24 + HCHO + NO2                                   
        RC(:,352) = 9.20D-14*RO2(:)       

        ! Reaction (353) OH + CARB14 = RN14O2                                               
        RC(:,353) = 1.87D-11       

        ! Reaction (354) OH + CARB17 = RN17O2                                               
        RC(:,354) = 4.36D-12       

        ! Reaction (355) OH + CARB11A = RN11O2                                              
        RC(:,355) = 3.24D-18*TEMP(:)**2*EXP(414/TEMP(:))

        ! Reaction (356) OH + CARB7 = CARB6 + HO2                                           
        RC(:,356) = 1.60D-12*EXP(305/TEMP(:))      

        ! Reaction (357) OH + CARB10 = CARB9 + HO2                                          
        RC(:,357) = 5.86D-12       

        ! Reaction (358) OH + CARB13 = RN13O2                                               
        RC(:,358) = 1.65D-11       

        ! Reaction (359) OH + CARB16 = RN16O2                                               
        RC(:,359) = 1.25D-11       

        ! Reaction (360) OH + UCARB10 = RU10O2                                              
        RC(:,360) = 3.84D-12*EXP(533/TEMP(:))*0.693     

        ! Reaction (361) NO3 + UCARB10 = RU10O2 + HNO3                                      
        RC(:,361) = KNO3AL(:)*0.415    

        ! Reaction (362) O3 + UCARB10 = HCHO + CH3CO3 + CO + OH                             
        RC(:,362) = 1.20D-15*EXP(-1710/TEMP)*0.32       

        ! Reaction (363) O3 + UCARB10 = HCHO + CARB6 + H2O2                                 
        RC(:,363) = 1.20D-15*EXP(-1710/TEMP)*0.68       

        ! Reaction (364) OH + HOCH2CHO = HOCH2CO3                                           
        RC(:,364) = 1.00D-11       

        ! Reaction (365) NO3 + HOCH2CHO = HOCH2CO3 + HNO3                                   
        RC(:,365) = KNO3AL(:)        

        ! Reaction (366) OH + CARB3 = CO + CO + HO2                                         
        RC(:,366) = 3.10D-12*EXP(340/TEMP(:))*0.8       

        ! Reaction (367) OH + CARB6 = CH3CO3 + CO                                           
        RC(:,367) = 1.90D-12*EXP(575/TEMP(:))      

        ! Reaction (368) OH + CARB9 = RN9O2                                                 
        RC(:,368) = 2.40D-13       

        ! Reaction (369) OH + CARB12 = RN12O2                                               
        RC(:,369) = 1.38D-12       

        ! Reaction (370) OH + CARB15 = RN15O2                                               
        RC(:,370) = 4.81D-12       

        ! Reaction (371) OH + CCARB12 = RN12O2                                              
        RC(:,371) = 4.79D-12       

        ! Reaction (372) OH + UCARB12 = RU12O2                                              
        RC(:,372) = 4.52D-11            

        ! Reaction (373) NO3 + UCARB12 = RU12O2 + HNO3                                      
        RC(:,373) = KNO3AL(:)*4.25    

        ! Reaction (374) O3 + UCARB12 = HOCH2CHO + CH3CO3 + CO + OH                         
        RC(:,374) = 2.40D-17*0.89   

        ! Reaction (375) O3 + UCARB12 = HOCH2CHO + CARB6 + H2O2                             
        RC(:,375) = 2.40D-17*0.11   

        ! Reaction (376) OH + NUCARB12 = NRU12O2                                            
        RC(:,376) = 4.16D-11            

        ! Reaction (377) OH + NOA = CARB6 + NO2                                             
        RC(:,377) = 6.70D-13             

        ! Reaction (378) OH + UDCARB8 = C2H5O2                                              
        RC(:,378) = 5.20D-11*0.50        

        ! Reaction (379) OH + UDCARB8 = ANHY + HO2                                          
        RC(:,379) = 5.20D-11*0.50        

        ! Reaction (380) OH + UDCARB11 = RN10O2                                             
        RC(:,380) = 5.58D-11*0.55     

        ! Reaction (381) OH + UDCARB11 = ANHY + CH3O2                                       
        RC(:,381) = 5.58D-11*0.45     

        ! Reaction (382) OH + UDCARB14 = RN13O2                                             
        RC(:,382) = 7.00D-11*0.55     

        ! Reaction (383) OH + UDCARB14 = ANHY + C2H5O2                                      
        RC(:,383) = 7.00D-11*0.45     

        ! Reaction (384) OH + TNCARB26 = RTN26O2                                            
        RC(:,384) = 4.20D-11           

        ! Reaction (385) OH + TNCARB15 = RN15AO2                                            
        RC(:,385) = 1.00D-12           

        ! Reaction (386) OH + TNCARB10 = RTN10O2                                            
        RC(:,386) = 1.00D-10           

        ! Reaction (387) NO3 + TNCARB26 = RTN26O2 + HNO3                                    
        RC(:,387) = 3.80D-14            

        ! Reaction (388) NO3 + TNCARB10 = RTN10O2 + HNO3                                    
        RC(:,388) = KNO3AL(:)*5.5      

        ! Reaction (389) OH + RCOOH25 = RTN25O2                                             
        RC(:,389) = 6.65D-12            

        ! Reaction (390) OH + TXCARB24 = RTX24O2                                            
        RC(:,390) = 1.55D-11           

        ! Reaction (391) OH + TXCARB22 = RTX22O2                                            
        RC(:,391) = 4.55D-12           

        ! Reaction (392) OH + CH3NO3 = HCHO + NO2                                           
        RC(:,392) = 1.00D-14*EXP(1060/TEMP(:))      

        ! Reaction (393) OH + C2H5NO3 = CH3CHO + NO2                                        
        RC(:,393) = 4.40D-14*EXP(720/TEMP(:))       

        ! Reaction (394) OH + RN10NO3 = C2H5CHO + NO2                                       
        RC(:,394) = 7.30D-13                     

        ! Reaction (395) OH + IC3H7NO3 = CH3COCH3 + NO2                                     
        RC(:,395) = 4.90D-13                     

        ! Reaction (396) OH + RN13NO3 = CARB11A + NO2                                       
        RC(:,396) = 9.20D-13                     

        ! Reaction (397) OH + RN16NO3 = CARB14 + NO2                                        
        RC(:,397) = 1.85D-12                     

        ! Reaction (398) OH + RN19NO3 = CARB17 + NO2                                        
        RC(:,398) = 3.02D-12                     

        ! Reaction (399) OH + HOC2H4NO3 = HOCH2CHO + NO2                                    
        RC(:,399) = 1.09D-12               

        ! Reaction (400) OH + RN9NO3 = CARB7 + NO2                                          
        RC(:,400) = 1.31D-12               

        ! Reaction (401) OH + RN12NO3 = CARB10 + NO2                                        
        RC(:,401) = 1.79D-12               

        ! Reaction (402) OH + RN15NO3 = CARB13 + NO2                                        
        RC(:,402) = 1.03D-11               

        ! Reaction (403) OH + RN18NO3 = CARB16 + NO2                                        
        RC(:,403) = 1.34D-11               

        ! Reaction (404) OH + RU14NO3 = UCARB12 + NO2                                       
        RC(:,404) = 3.00D-11*0.34               

        ! Reaction (405) OH + RA13NO3 = CARB3 + UDCARB8 + NO2                               
        RC(:,405) = 7.30D-11               

        ! Reaction (406) OH + RA16NO3 = CARB3 + UDCARB11 + NO2                              
        RC(:,406) = 7.16D-11               

        ! Reaction (407) OH + RA19NO3 = CARB6 + UDCARB11 + NO2                              
        RC(:,407) = 8.31D-11               

        ! Reaction (408) OH + RTN28NO3 = TNCARB26 + NO2                                     
        RC(:,408) = 4.35D-12               

        ! Reaction (409) OH + RTN25NO3 = CH3COCH3 + TNCARB15 + NO2                          
        RC(:,409) = 2.88D-12               

        ! Reaction (410) OH + RTX28NO3 = TXCARB24 + HCHO + NO2                              
        RC(:,410) = 3.53D-12                  

        ! Reaction (411) OH + RTX24NO3 = TXCARB22 + NO2                                     
        RC(:,411) = 6.48D-12                  

        ! Reaction (412) OH + RTX22NO3 = CH3COCH3 + CCARB12 + NO2                           
        RC(:,412) = 4.74D-12                  

        ! Reaction (413) OH + AROH14 = RAROH14                                              
        RC(:,413) = 2.63D-11             

        ! Reaction (414) NO3 + AROH14 = RAROH14 + HNO3                                      
        RC(:,414) = 3.78D-12               

        ! Reaction (415) RAROH14 + NO2 = ARNOH14                                            
        RC(:,415) = 2.08D-12               

        ! Reaction (416) OH + ARNOH14 = CARB13 + NO2                                        
        RC(:,416) = 9.00D-13               

        ! Reaction (417) NO3 + ARNOH14 = CARB13 + NO2 + HNO3                                
        RC(:,417) = 9.00D-14               

        ! Reaction (418) OH + AROH17 = RAROH17                                              
        RC(:,418) = 4.65D-11               

        ! Reaction (419) NO3 + AROH17 = RAROH17 + HNO3                                      
        RC(:,419) = 1.25D-11               

        ! Reaction (420) RAROH17 + NO2 = ARNOH17                                            
        RC(:,420) = 2.08D-12               

        ! Reaction (421) OH + ARNOH17 = CARB16 + NO2                                        
        RC(:,421) = 1.53D-12               

        ! Reaction (422) NO3 + ARNOH17 = CARB16 + NO2 + HNO3                                
        RC(:,422) = 3.13D-13               

        ! Reaction (423) OH + CH3OOH = CH3O2                                                
        RC(:,423) = 1.90D-11*EXP(190/TEMP(:))       

        ! Reaction (424) OH + CH3OOH = HCHO + OH                                            
        RC(:,424) = 1.00D-11*EXP(190/TEMP(:))       

        ! Reaction (425) OH + C2H5OOH = CH3CHO + OH                                         
        RC(:,425) = 1.36D-11               

        ! Reaction (426) OH + RN10OOH = C2H5CHO + OH                                        
        RC(:,426) = 1.89D-11               

        ! Reaction (427) OH + IC3H7OOH = CH3COCH3 + OH                                      
        RC(:,427) = 2.78D-11               

        ! Reaction (428) OH + RN13OOH = CARB11A + OH                                        
        RC(:,428) = 3.57D-11               

        ! Reaction (429) OH + RN16OOH = CARB14 + OH                                         
        RC(:,429) = 4.21D-11               

        ! Reaction (430) OH + RN19OOH = CARB17 + OH                                         
        RC(:,430) = 4.71D-11               

        ! Reaction (431) OH + CH3CO3H = CH3CO3                                              
        RC(:,431) = 3.70D-12                     

        ! Reaction (432) OH + C2H5CO3H = C2H5CO3                                            
        RC(:,432) = 4.42D-12                     

        ! Reaction (433) OH + HOCH2CO3H = HOCH2CO3                                          
        RC(:,433) = 6.19D-12                     

        ! Reaction (434) OH + RN8OOH = CARB6 + OH                                           
        RC(:,434) = 1.2D-11                     

        ! Reaction (435) OH + RN11OOH = CARB9 + OH                                          
        RC(:,435) = 2.50D-11                     

        ! Reaction (436) OH + RN14OOH = CARB12 + OH                                         
        RC(:,436) = 3.20D-11                     

        ! Reaction (437) OH + RN17OOH = CARB15 + OH                                         
        RC(:,437) = 3.35D-11                     

        ! Reaction (438) OH + RU14OOH = UCARB12 + OH                                        
        RC(:,438) = 7.51D-11                     

        ! Reaction (439) OH + RU12OOH = RU12O2                                              
        RC(:,439) = 3.50D-11                     

        ! Reaction (440) OH + RU10OOH = RU10O2                                              
        RC(:,440) = 3.50D-11                     

        ! Reaction (441) OH + NRU14OOH = NUCARB12 + OH                                      
        RC(:,441) = 1.03D-10                     

        ! Reaction (442) OH + NRU12OOH = NOA + CO + OH                                      
        RC(:,442) = 2.65D-11                     

        ! Reaction (443) OH + HOC2H4OOH = HOCH2CHO + OH                                     
        RC(:,443) = 2.13D-11               

        ! Reaction (444) OH + RN9OOH = CARB7 + OH                                           
        RC(:,444) = 2.50D-11               

        ! Reaction (445) OH + RN12OOH = CARB10 + OH                                         
        RC(:,445) = 3.25D-11               

        ! Reaction (446) OH + RN15OOH = CARB13 + OH                                         
        RC(:,446) = 3.74D-11               

        ! Reaction (447) OH + RN18OOH = CARB16 + OH                                         
        RC(:,447) = 3.83D-11               

        ! Reaction (448) OH + NRN6OOH = HCHO + HCHO + NO2 + OH                              
        RC(:,448) = 5.22D-12               

        ! Reaction (449) OH + NRN9OOH = CH3CHO + HCHO + NO2 + OH                            
        RC(:,449) = 6.50D-12               

        ! Reaction (450) OH + NRN12OOH = CH3CHO + CH3CHO + NO2 + OH                         
        RC(:,450) = 7.15D-12               

        ! Reaction (451) OH + RA13OOH = CARB3 + UDCARB8 + OH                                
        RC(:,451) = 9.77D-11               

        ! Reaction (452) OH + RA16OOH = CARB3 + UDCARB11 + OH                               
        RC(:,452) = 9.64D-11               

        ! Reaction (453) OH + RA19OOH = CARB6 + UDCARB11 + OH                               
        RC(:,453) = 1.12D-10               

        ! Reaction (454) OH + RTN28OOH = TNCARB26 + OH                                      
        RC(:,454) = 2.38D-11               

        ! Reaction (455) OH + RTN26OOH = RTN26O2                                            
        RC(:,455) = 1.20D-11               

        ! Reaction (456) OH + NRTN28OOH = TNCARB26 + NO2 + OH                               
        RC(:,456) = 9.50D-12               

        ! Reaction (457) OH + RTN25OOH = RTN25O2                                            
        RC(:,457) = 1.66D-11               

        ! Reaction (458) OH + RTN24OOH = RTN24O2                                            
        RC(:,458) = 1.05D-11               

        ! Reaction (459) OH + RTN23OOH = RTN23O2                                            
        RC(:,459) = 2.05D-11               

        ! Reaction (460) OH + RTN14OOH = RTN14O2                                            
        RC(:,460) = 8.69D-11               

        ! Reaction (461) OH + RTN10OOH = RTN10O2                                            
        RC(:,461) = 4.23D-12               

        ! Reaction (462) OH + RTX28OOH = RTX28O2                                            
        RC(:,462) = 2.00D-11               

        ! Reaction (463) OH + RTX24OOH = TXCARB22 + OH                                      
        RC(:,463) = 8.59D-11               

        ! Reaction (464) OH + RTX22OOH = CH3COCH3 + CCARB12 + OH                            
        RC(:,464) = 7.50D-11               

        ! Reaction (465) OH + NRTX28OOH = NRTX28O2                                          
        RC(:,465) = 9.58D-12               

        ! Reaction (466) OH + ANHY = HOCH2CH2O2                                             
        RC(:,466) = 1.50D-12        

        ! Reaction (467) CH3CO3 + NO2 = PAN                                                 
        RC(:,467) = KFPAN(:)                        

        ! Reaction (468) PAN = CH3CO3 + NO2                                                 
        RC(:,468) = KBPAN(:)                        

        ! Reaction (469) C2H5CO3 + NO2 = PPN                                                
        RC(:,469) = KFPAN(:)                        

        ! Reaction (470) PPN = C2H5CO3 + NO2                                                
        RC(:,470) = KBPAN(:)                        

        ! Reaction (471) HOCH2CO3 + NO2 = PHAN                                              
        RC(:,471) = 1.125D-11*(TEMP(:)/300)**(-1.105)                            

        ! Reaction (472) PHAN = HOCH2CO3 + NO2                                              
        RC(:,472) = 5.2D16*EXP(-13850/TEMP(:))                       

        ! Reaction (473) OH + PAN = HCHO + CO + NO2                                         
        RC(:,473) = 3.00D-14    

        ! Reaction (474) OH + PPN = CH3CHO + CO + NO2                                       
        RC(:,474) = 1.27D-12                       

        ! Reaction (475) OH + PHAN = HCHO + CO + NO2                                        
        RC(:,475) = 1.12D-12                       

        ! Reaction (476) RU12O2 + NO2 = RU12PAN                                             
        RC(:,476) = KFPAN(:)*0.061             

        ! Reaction (477) RU12PAN = RU12O2 + NO2                                             
        RC(:,477) = KBPAN(:)                   

        ! Reaction (478) RU10O2 + NO2 = MPAN                                                
        RC(:,478) = KFPAN(:)*0.041             

        ! Reaction (479) MPAN = RU10O2 + NO2                                                
        RC(:,479) = KBPAN(:)                  

        ! Reaction (480) OH + MPAN = CARB7 + CO + NO2                                       
        RC(:,480) = 2.90D-11*0.22 

        ! Reaction (481) OH + RU12PAN = UCARB10 + NO2                                       
        RC(:,481) = 2.52D-11 

        ! Reaction (482) RTN26O2 + NO2 = RTN26PAN                                           
        RC(:,482) = 1.125D-11*(TEMP(:)/300)**(-1.105)*0.722      

        ! Reaction (483) RTN26PAN = RTN26O2 + NO2                                           
        RC(:,483) = 5.2D16*EXP(-13850/TEMP(:))                   

        ! Reaction (484) OH + RTN26PAN = CH3COCH3 + CARB16 + NO2                            
        RC(:,484) = 3.66D-12  

        ! Reaction (485) RTN28NO3 = P2604                                                   
        RC(:,485) = KIN(:)  		

        ! Reaction (486) P2604 = RTN28NO3                                                   
        RC(:,486) = KOUT2604(:)

        ! Reaction (487) RTX28NO3 = P4608                                                   
        RC(:,487) = KIN(:)	

        ! Reaction (488) P4608 = RTX28NO3                                                   
        RC(:,488) = KOUT4608(:) 	

        ! Reaction (489) RCOOH25 = P2631                                                    
        RC(:,489) = KIN(:)  		

        ! Reaction (490) P2631 = RCOOH25                                                    
        RC(:,490) = KOUT2631(:) 	

        ! Reaction (491) RTN24OOH = P2635                                                   
        RC(:,491) = KIN(:)  		

        ! Reaction (492) P2635 = RTN24OOH                                                   
        RC(:,492) = KOUT2635(:) 	

        ! Reaction (493) RTX28OOH = P4610                                                   
        RC(:,493) = KIN(:)  		

        ! Reaction (494) P4610 = RTX28OOH                                                   
        RC(:,494) = KOUT4610(:) 	

        ! Reaction (495) RTN28OOH = P2605                                                   
        RC(:,495) = KIN(:)  		

        ! Reaction (496) P2605 = RTN28OOH                                                   
        RC(:,496) = KOUT2605(:) 	

        ! Reaction (497) RTN26OOH = P2630                                                   
        RC(:,497) = KIN(:)		

        ! Reaction (498) P2630 = RTN26OOH                                                   
        RC(:,498) = KOUT2630(:) 	

        ! Reaction (499) RTN26PAN = P2629                                                   
        RC(:,499) = KIN(:) 		

        ! Reaction (500) P2629 = RTN26PAN                                                   
        RC(:,500) = KOUT2629(:) 	

        ! Reaction (501) RTN25OOH = P2632                                                   
        RC(:,501) = KIN(:) 		

        ! Reaction (502) P2632 = RTN25OOH                                                   
        RC(:,502) = KOUT2632(:) 	

        ! Reaction (503) RTN23OOH = P2637                                                   
        RC(:,503) = KIN(:)  		

        ! Reaction (504) P2637 = RTN23OOH                                                   
        RC(:,504) = KOUT2637(:) 	

        ! Reaction (505) ARNOH14 = P3612                                                    
        RC(:,505) = KIN(:)  		

        ! Reaction (506) P3612 = ARNOH14                                                    
        RC(:,506) = KOUT3612(:) 	

        ! Reaction (507) ARNOH17 = P3613                                                    
        RC(:,507) = KIN(:) 		

        ! Reaction (508) P3613 = ARNOH17                                                    
        RC(:,508) = KOUT3613(:) 	

        ! Reaction (509) ANHY = P3442                                                       
        RC(:,509) = KIN(:)  		

        ! Reaction (510) P3442 = ANHY                                                       
        RC(:,510) = KOUT3442(:)

    END SUBROUTINE CHEMCO

    SUBROUTINE CALC_J
        IMPLICIT NONE
        REAL :: PI, COSX, SECX

        PI = 4.00E+00*ATAN(1.00E+00)

        DO CELL = 1,NCELL
            IF((SZA(CELL) - 0.5*PI).LT.0) THEN
                COSX = COS(SZA(CELL))
                SECX = 1.0/COSX
                J(CELL,1)=6.073D-05*(COSX**(1.743))*EXP(-0.474*SECX) 
                J(CELL,2)=4.775D-04*(COSX**(0.298))*EXP(-0.080*SECX) 
                J(CELL,3)=1.041D-05*(COSX**(0.723))*EXP(-0.279*SECX) 
                J(CELL,4)=1.165D-02*(COSX**(0.244))*EXP(-0.267*SECX) 
                J(CELL,5)=2.485D-02*(COSX**(0.168))*EXP(-0.108*SECX) 
                J(CELL,6)=1.747D-01*(COSX**(0.155))*EXP(-0.125*SECX) 
                J(CELL,7)=2.644D-03*(COSX**(0.261))*EXP(-0.288*SECX) 
                J(CELL,8)=9.312D-07*(COSX**(1.230))*EXP(-0.307*SECX) 
                J(CELL,11)=4.642D-05*(COSX**(0.762))*EXP(-0.353*SECX)
                J(CELL,12)=6.853D-05*(COSX**(0.477))*EXP(-0.323*SECX) 
                J(CELL,13)=7.344D-06*(COSX**(1.202))*EXP(-0.417*SECX) 
                J(CELL,14)=2.879D-05*(COSX**(1.067))*EXP(-0.358*SECX) 
                J(CELL,15)=2.792D-05*(COSX**(0.805))*EXP(-0.338*SECX) 
                J(CELL,16)=1.675D-05*(COSX**(0.805))*EXP(-0.338*SECX) 
                J(CELL,17)=7.914D-05*(COSX**(0.764))*EXP(-0.364*SECX) 
                J(CELL,18)=1.140D-05*(COSX**(0.396))*EXP(-0.298*SECX) 
                J(CELL,19)=1.140D-05*(COSX**(0.396))*EXP(-0.298*SECX) 
                J(CELL,21)=7.992D-07*(COSX**(1.578))*EXP(-0.271*SECX) 
                J(CELL,22)=5.804D-06*(COSX**(1.092))*EXP(-0.377*SECX) 
                J(CELL,23)=1.836D-05*(COSX**(0.395))*EXP(-0.296*SECX) 
                J(CELL,24)=1.836D-05*(COSX**(0.395))*EXP(-0.296*SECX) 
                J(CELL,31)=6.845D-05*(COSX**(0.130))*EXP(-0.201*SECX) 
                J(CELL,32)=1.032D-05*(COSX**(0.130))*EXP(-0.201*SECX) 
                J(CELL,33)=3.802D-05*(COSX**(0.644))*EXP(-0.312*SECX) 
                J(CELL,34)=1.537D-04*(COSX**(0.170))*EXP(-0.208*SECX) 
                J(CELL,35)=3.326D-04*(COSX**(0.148))*EXP(-0.215*SECX) 
                J(CELL,41)=7.649D-06*(COSX**(0.682))*EXP(-0.279*SECX) 
                J(CELL,51)=1.588D-06*(COSX**(1.154))*EXP(-0.318*SECX) 
                J(CELL,52)=1.907D-06*(COSX**(1.244))*EXP(-0.335*SECX) 
                J(CELL,53)=2.485D-06*(COSX**(1.196))*EXP(-0.328*SECX) 
                J(CELL,54)=4.095D-06*(COSX**(1.111))*EXP(-0.316*SECX) 
                J(CELL,55)=1.135D-05*(COSX**(0.974))*EXP(-0.309*SECX) 
                J(CELL,56)=7.549D-06*(COSX**(1.015))*EXP(-0.324*SECX) 
                J(CELL,57)=3.363D-06*(COSX**(1.296))*EXP(-0.322*SECX)
            ELSE
                J(CELL,:) = 1.0E-30
            END IF
            
        END DO

    END SUBROUTINE CALC_J

    SUBROUTINE PHOTOL
        IMPLICIT NONE
        ! Photol Reaction (1) O3 = O1D                                                           
        DJ(:,1) = J(:,1)                             

        ! Photol Reaction (2) O3 = O                                                             
        DJ(:,2) = J(:,2)                             

        ! Photol Reaction (3) H2O2 = OH + OH                                                     
        DJ(:,3) = J(:,3)                             

        ! Photol Reaction (4) NO2 = NO + O                                                       
        DJ(:,4) = J(:,4)                             

        ! Photol Reaction (5) NO3 = NO                                                           
        DJ(:,5) = J(:,5)                             

        ! Photol Reaction (6) NO3 = NO2 + O                                                      
        DJ(:,6) = J(:,6)                             

        ! Photol Reaction (7) HONO = OH + NO                                                     
        DJ(:,7) = J(:,7)                             

        ! Photol Reaction (8) HNO3 = OH + NO2                                                    
        DJ(:,8) = J(:,8)                             

        ! Photol Reaction (9) HCHO = CO + HO2 + HO2                                              
        DJ(:,9) = J(:,11)                        

        ! Photol Reaction (10) HCHO = H2 + CO                                                     
        DJ(:,10) = J(:,12)                        

        ! Photol Reaction (11) CH3CHO = CH3O2 + HO2 + CO                                          
        DJ(:,11) = J(:,13)                        

        ! Photol Reaction (12) C2H5CHO = C2H5O2 + CO + HO2                                        
        DJ(:,12) = J(:,14)                        

        ! Photol Reaction (13) CH3COCH3 = CH3CO3 + CH3O2                                          
        DJ(:,13) = J(:,21)                        

        ! Photol Reaction (14) MEK = CH3CO3 + C2H5O2                                              
        DJ(:,14) = J(:,22)                        

        ! Photol Reaction (15) CARB14 = CH3CO3 + RN10O2                                           
        DJ(:,15) = J(:,22)*4.74               

        ! Photol Reaction (16) CARB17 = RN8O2 + RN10O2                                            
        DJ(:,16) = J(:,22)*1.33               

        ! Photol Reaction (17) CARB11A = CH3CO3 + C2H5O2                                          
        DJ(:,17) = J(:,22)                        

        ! Photol Reaction (18) CARB7 = CH3CO3 + HCHO + HO2                                        
        DJ(:,18) = J(:,22)                        

        ! Photol Reaction (19) CARB10 = CH3CO3 + CH3CHO + HO2                                     
        DJ(:,19) = J(:,22)                        

        ! Photol Reaction (20) CARB13 = RN8O2 + CH3CHO + HO2                                      
        DJ(:,20) = J(:,22)*3.00               

        ! Photol Reaction (21) CARB16 = RN8O2 + C2H5CHO + HO2                                     
        DJ(:,21) = J(:,22)*3.35               

        ! Photol Reaction (22) HOCH2CHO = HCHO + CO + HO2 + HO2                                   
        DJ(:,22) = J(:,15)                        

        ! Photol Reaction (23) UCARB10 = CH3CO3 + HCHO + HO2                                      
        DJ(:,23) = J(:,18)*2                       

        ! Photol Reaction (24) CARB3 = CO + CO + HO2 + HO2                                        
        DJ(:,24) = J(:,33)                        

        ! Photol Reaction (25) CARB6 = CH3CO3 + CO + HO2                                          
        DJ(:,25) = J(:,34)                        

        ! Photol Reaction (26) CARB9 = CH3CO3 + CH3CO3                                            
        DJ(:,26) = J(:,35)                        

        ! Photol Reaction (27) CARB12 = CH3CO3 + RN8O2                                            
        DJ(:,27) = J(:,35)                        

        ! Photol Reaction (28) CARB15 = RN8O2 + RN8O2                                             
        DJ(:,28) = J(:,35)                        

        ! Photol Reaction (29) UCARB12 = CH3CO3 + HOCH2CHO + CO + HO2                             
        DJ(:,29) = J(:,18)*2           

        ! Photol Reaction (30) NUCARB12 = NOA + CO + CO + HO2 + HO2                               
        DJ(:,30) = J(:,18)             

        ! Photol Reaction (31) NOA = CH3CO3 + HCHO + NO2                                          
        DJ(:,31) = J(:,56)             

        ! Photol Reaction (32) NOA = CH3CO3 + HCHO + NO2                                          
        DJ(:,32) = J(:,57)             

        ! Photol Reaction (33) UDCARB8 = C2H5O2 + HO2                                             
        DJ(:,33) = J(:,4)*0.02*0.64   

        ! Photol Reaction (34) UDCARB8 = ANHY + HO2 + HO2                                         
        DJ(:,34) = J(:,4)*0.02*0.36   

        ! Photol Reaction (35) UDCARB11 = RN10O2 + HO2                                            
        DJ(:,35) = J(:,4)*0.02*0.55   

        ! Photol Reaction (36) UDCARB11 = ANHY + HO2 + CH3O2                                      
        DJ(:,36) = J(:,4)*0.02*0.45   

        ! Photol Reaction (37) UDCARB14 = RN13O2 + HO2                                            
        DJ(:,37) = J(:,4)*0.02*0.55   

        ! Photol Reaction (38) UDCARB14 = ANHY + HO2 + C2H5O2                                     
        DJ(:,38) = J(:,4)*0.02*0.45   

        ! Photol Reaction (39) TNCARB26 = RTN26O2 + HO2                                           
        DJ(:,39) = J(:,15)             

        ! Photol Reaction (40) TNCARB10 = CH3CO3 + CH3CO3 + CO                                    
        DJ(:,40) = J(:,35)*0.5        

        ! Photol Reaction (41) CH3NO3 = HCHO + HO2 + NO2                                          
        DJ(:,41) = J(:,51)                        

        ! Photol Reaction (42) C2H5NO3 = CH3CHO + HO2 + NO2                                       
        DJ(:,42) = J(:,52)                        

        ! Photol Reaction (43) RN10NO3 = C2H5CHO + HO2 + NO2                                      
        DJ(:,43) = J(:,53)                        

        ! Photol Reaction (44) IC3H7NO3 = CH3COCH3 + HO2 + NO2                                    
        DJ(:,44) = J(:,54)                        

        ! Photol Reaction (45) RN13NO3 =  CH3CHO + C2H5O2 + NO2                                   
        DJ(:,45) = J(:,53)*BR01(:)               

        ! Photol Reaction (46) RN13NO3 =  CARB11A + HO2 + NO2                                     
        DJ(:,46) = J(:,53)*(1-BR01(:))           

        ! Photol Reaction (47) RN16NO3 = RN15O2 + NO2                                             
        DJ(:,47) = J(:,53)                        

        ! Photol Reaction (48) RN19NO3 = RN18O2 + NO2                                             
        DJ(:,48) = J(:,53)                        

        ! Photol Reaction (49) RA13NO3 = CARB3 + UDCARB8 + HO2 + NO2                              
        DJ(:,49) = J(:,54)                    

        ! Photol Reaction (50) RA16NO3 = CARB3 + UDCARB11 + HO2 + NO2                             
        DJ(:,50) = J(:,54)                    

        ! Photol Reaction (51) RA19NO3 = CARB6 + UDCARB11 + HO2 + NO2                             
        DJ(:,51) = J(:,54)                    

        ! Photol Reaction (52) RTX24NO3 = TXCARB22 + HO2 + NO2                                    
        DJ(:,52) = J(:,54)                    

        ! Photol Reaction (53) CH3OOH = HCHO + HO2 + OH                                           
        DJ(:,53) = J(:,41)                        

        ! Photol Reaction (54) C2H5OOH = CH3CHO + HO2 + OH                                        
        DJ(:,54) = J(:,41)                        

        ! Photol Reaction (55) RN10OOH = C2H5CHO + HO2 + OH                                       
        DJ(:,55) = J(:,41)                        

        ! Photol Reaction (56) IC3H7OOH = CH3COCH3 + HO2 + OH                                     
        DJ(:,56) = J(:,41)                        

        ! Photol Reaction (57) RN13OOH =  CH3CHO + C2H5O2 + OH                                    
        DJ(:,57) = J(:,41)*BR01(:)         

        ! Photol Reaction (58) RN13OOH =  CARB11A + HO2 + OH                                      
        DJ(:,58) = J(:,41)*(1-BR01(:))     

        ! Photol Reaction (59) RN16OOH = RN15AO2 + OH                                             
        DJ(:,59) = J(:,41)                        

        ! Photol Reaction (60) RN19OOH = RN18AO2 + OH                                             
        DJ(:,60) = J(:,41)                        

        ! Photol Reaction (61) CH3CO3H = CH3O2 + OH                                               
        DJ(:,61) = J(:,41)                        

        ! Photol Reaction (62) C2H5CO3H = C2H5O2 + OH                                             
        DJ(:,62) = J(:,41)                        

        ! Photol Reaction (63) HOCH2CO3H = HCHO + HO2 + OH                                        
        DJ(:,63) = J(:,41)                        

        ! Photol Reaction (64) RN8OOH = C2H5O2 + OH                                               
        DJ(:,64) = J(:,41)                        

        ! Photol Reaction (65) RN11OOH = RN10O2 + OH                                              
        DJ(:,65) = J(:,41)                        

        ! Photol Reaction (66) RN14OOH = RN13O2 + OH                                              
        DJ(:,66) = J(:,41)                        

        ! Photol Reaction (67) RN17OOH = RN16O2 + OH                                              
        DJ(:,67) = J(:,41)                        

        ! Photol Reaction (68) RU14OOH = UCARB12 + HO2 + OH                                       
        DJ(:,68) = J(:,41)*0.252              

        ! Photol Reaction (69) RU14OOH = UCARB10 + HCHO + HO2 + OH                                
        DJ(:,69) = J(:,41)*0.748              

        ! Photol Reaction (70) RU12OOH = CARB6 + HOCH2CHO + HO2 + OH                              
        DJ(:,70) = J(:,41)                   

        ! Photol Reaction (71) RU10OOH = CH3CO3 + HOCH2CHO + OH                                   
        DJ(:,71) = J(:,41)                   

        ! Photol Reaction (72) NRU14OOH = NUCARB12 + HO2 + OH                                     
        DJ(:,72) = J(:,41)                   

        ! Photol Reaction (73) NRU12OOH = NOA + CO + HO2 + OH                                     
        DJ(:,73) = J(:,41)                   

        ! Photol Reaction (74) HOC2H4OOH = HCHO + HCHO + HO2 + OH                                 
        DJ(:,74) = J(:,41)                  

        ! Photol Reaction (75) RN9OOH = CH3CHO + HCHO + HO2 + OH                                  
        DJ(:,75) = J(:,41)                  

        ! Photol Reaction (76) RN12OOH = CH3CHO + CH3CHO + HO2 + OH                               
        DJ(:,76) = J(:,41)                  

        ! Photol Reaction (77) RN15OOH = C2H5CHO + CH3CHO + HO2 + OH                              
        DJ(:,77) = J(:,41)                  

        ! Photol Reaction (78) RN18OOH = C2H5CHO + C2H5CHO + HO2 + OH                             
        DJ(:,78) = J(:,41)                 

        ! Photol Reaction (79) NRN6OOH = HCHO + HCHO + NO2 + OH                                   
        DJ(:,79) = J(:,41)                  

        ! Photol Reaction (80) NRN9OOH = CH3CHO + HCHO + NO2 + OH                                 
        DJ(:,80) = J(:,41)                  

        ! Photol Reaction (81) NRN12OOH = CH3CHO + CH3CHO + NO2 + OH                              
        DJ(:,81) = J(:,41)                  

        ! Photol Reaction (82) RA13OOH = CARB3 + UDCARB8 + HO2 + OH                               
        DJ(:,82) = J(:,41)                  

        ! Photol Reaction (83) RA16OOH = CARB3 + UDCARB11 + HO2 + OH                              
        DJ(:,83) = J(:,41)                  

        ! Photol Reaction (84) RA19OOH = CARB6 + UDCARB11 + HO2 + OH                              
        DJ(:,84) = J(:,41)                  

        ! Photol Reaction (85) RTN28OOH = TNCARB26 + HO2 + OH                                     
        DJ(:,85) = J(:,41)                  

        ! Photol Reaction (86) NRTN28OOH = TNCARB26 + NO2 + OH                                    
        DJ(:,86) = J(:,41)                  

        ! Photol Reaction (87) RTN26OOH = RTN25O2 + OH                                            
        DJ(:,87) = J(:,41)             

        ! Photol Reaction (88) RTN25OOH = RTN24O2 + OH                                            
        DJ(:,88) = J(:,41)             

        ! Photol Reaction (89) RTN24OOH = RTN23O2 + OH                                            
        DJ(:,89) = J(:,41)             

        ! Photol Reaction (90) RTN23OOH = CH3COCH3 + RTN14O2 + OH                                 
        DJ(:,90) = J(:,41)             

        ! Photol Reaction (91) RTN14OOH = TNCARB10 + HCHO + HO2 + OH                              
        DJ(:,91) = J(:,41)             

        ! Photol Reaction (92) RTN10OOH = RN8O2 + CO + OH                                         
        DJ(:,92) = J(:,41)             

        ! Photol Reaction (93) RTX28OOH = TXCARB24 + HCHO + HO2 + OH                              
        DJ(:,93) = J(:,41)                  

        ! Photol Reaction (94) RTX24OOH = TXCARB22 + HO2 + OH                                     
        DJ(:,94) = J(:,41)                  

        ! Photol Reaction (95) RTX22OOH = CH3COCH3 + RN13O2 + OH                                  
        DJ(:,95) = J(:,41)                  

        ! Photol Reaction (96) NRTX28OOH = TXCARB24 + HCHO + NO2 + OH                             
        DJ(:,96) = J(:,41)

    END SUBROUTINE PHOTOL

    SUBROUTINE DERIV(DTS)
        IMPLICIT NONE
        REAL, INTENT(IN) :: DTS
        ! YP IS PREVIOUS Y LESS PREVIOUS EMI PLUS CURRENT EMI
        YP(:,:) = Y(:,:)

        DO CELL = 1,NCELL
            IF (EMI(CELL,1) == 0) THEN
                EMIP(CELL,1) = 0
            ENDIF

            YP(CELL,8) = Y(CELL,8) - EMIP(CELL,1) + EMI(CELL,1)

            IF (EMI(CELL,2) == 0) THEN
                EMIP(CELL,2) = 0
            ENDIF

            YP(CELL,4) = Y(CELL,4) - EMIP(CELL,2) + EMI(CELL,2)

            IF (EMI(CELL,3) == 0) THEN
                EMIP(CELL,3) = 0
            ENDIF

            YP(CELL,11) = Y(CELL,11) - EMIP(CELL,3) + EMI(CELL,3)
        END DO

        ! PRINT *, Y(182,8), YP(182,8), EMIP(182,1), EMI(182,1)
        
        ! YP(:,39) = Y(:,39) - EMIP(:,3) + EMI(:,3)
        ! YP(:,42) = Y(:,42) - EMIP(:,4) + EMI(:,4)
        ! YP(:,73) = Y(:,73) - EMIP(:,5) + EMI(:,5)
        ! YP(:,23) = Y(:,23) - EMIP(:,6) + EMI(:,6)
        ! YP(:,30) = Y(:,30) - EMIP(:,7) + EMI(:,7)
        ! YP(:,25) = Y(:,25) - EMIP(:,8) + EMI(:,8)
        ! YP(:,32) = Y(:,32) - EMIP(:,9) + EMI(:,9)
        ! YP(:,59) = Y(:,59) - EMIP(:,10) + EMI(:,10)
        ! YP(:,61) = Y(:,61) - EMIP(:,11) + EMI(:,11)
        ! YP(:,64) = Y(:,64) - EMIP(:,12) + EMI(:,12)
        ! YP(:,71) = Y(:,71) - EMIP(:,13) + EMI(:,13) 

        !          O1D              Y(:,1) 
        P(:) = DJ(:,1) * Y(:,6)
        L(:) = 0.0 + (RC(:,7)) + (RC(:,8)) + (RC(:,16) *  H2O(:)) 
        Y(:,1) = P(:)/L(:)
        
        !          O                Y(:,2)
        P(:) = (DJ(:,6) * Y(:,5)) &
        + (DJ(:,2) * Y(:,6)) + (DJ(:,4) * Y(:,4)) &
        + (RC(:,7) * Y(:,1)) + (RC(:,8) * Y(:,1)) 
        L(:) = 0.0 &
        + (RC(:,36) * Y(:,16)) &
        + (RC(:,4) * Y(:,8)) + (RC(:,5) * Y(:,4)) + (RC(:,6) * Y(:,4)) &
        + (RC(:,1)) + (RC(:,2)) + (RC(:,3) * Y(:,6)) 
        Y(:,2) = P(:)/L(:)
        
        !          OH               Y(:,3) 
        P(:) = (DJ(:,95) * Y(:,184)) + (DJ(:,96) * Y(:,185)) &
        + (DJ(:,93) * Y(:,182)) + (DJ(:,94) * Y(:,183)) &
        + (DJ(:,91) * Y(:,180)) + (DJ(:,92) * Y(:,181)) &
        + (DJ(:,89) * Y(:,178)) + (DJ(:,90) * Y(:,179)) &
        + (DJ(:,87) * Y(:,176)) + (DJ(:,88) * Y(:,177)) &
        + (DJ(:,85) * Y(:,174)) + (DJ(:,86) * Y(:,175)) &
        + (DJ(:,83) * Y(:,152)) + (DJ(:,84) * Y(:,153)) &
        + (DJ(:,81) * Y(:,171)) + (DJ(:,82) * Y(:,151)) &
        + (DJ(:,79) * Y(:,169)) + (DJ(:,80) * Y(:,170)) &
        + (DJ(:,77) * Y(:,157)) + (DJ(:,78) * Y(:,158)) &
        + (DJ(:,75) * Y(:,155)) + (DJ(:,76) * Y(:,156)) &
        + (DJ(:,73) * Y(:,173)) + (DJ(:,74) * Y(:,154)) &
        + (DJ(:,71) * Y(:,168)) + (DJ(:,72) * Y(:,172)) &
        + (DJ(:,69) * Y(:,166)) + (DJ(:,70) * Y(:,167)) &
        + (DJ(:,67) * Y(:,165)) + (DJ(:,68) * Y(:,166)) &
        + (DJ(:,65) * Y(:,163)) + (DJ(:,66) * Y(:,164)) &
        + (DJ(:,63) * Y(:,161)) + (DJ(:,64) * Y(:,162)) &
        + (DJ(:,61) * Y(:,159)) + (DJ(:,62) * Y(:,160)) &
        + (DJ(:,59) * Y(:,149)) + (DJ(:,60) * Y(:,150)) &
        + (DJ(:,57) * Y(:,148)) + (DJ(:,58) * Y(:,148)) &
        + (DJ(:,55) * Y(:,146)) + (DJ(:,56) * Y(:,147)) &
        + (DJ(:,53) * Y(:,144)) + (DJ(:,54) * Y(:,145)) &
        + (DJ(:,7) * Y(:,13)) + (DJ(:,8) * Y(:,14)) &
        + (RC(:,464) * Y(:,3) * Y(:,184)) + (DJ(:,3) * Y(:,12) * 2.00) &
        + (RC(:,456) * Y(:,3) * Y(:,175)) + (RC(:,463) * Y(:,3) * Y(:,183)) &
        + (RC(:,453) * Y(:,3) * Y(:,153)) + (RC(:,454) * Y(:,3) * Y(:,174)) &
        + (RC(:,451) * Y(:,3) * Y(:,151)) + (RC(:,452) * Y(:,3) * Y(:,152)) &
        + (RC(:,449) * Y(:,3) * Y(:,170)) + (RC(:,450) * Y(:,3) * Y(:,171)) &
        + (RC(:,447) * Y(:,3) * Y(:,158)) + (RC(:,448) * Y(:,3) * Y(:,169)) &
        + (RC(:,445) * Y(:,3) * Y(:,156)) + (RC(:,446) * Y(:,3) * Y(:,157)) &
        + (RC(:,443) * Y(:,3) * Y(:,154)) + (RC(:,444) * Y(:,3) * Y(:,155)) &
        + (RC(:,441) * Y(:,3) * Y(:,172)) + (RC(:,442) * Y(:,3) * Y(:,173)) &
        + (RC(:,437) * Y(:,3) * Y(:,165)) + (RC(:,438) * Y(:,3) * Y(:,166)) &
        + (RC(:,435) * Y(:,3) * Y(:,163)) + (RC(:,436) * Y(:,3) * Y(:,164)) &
        + (RC(:,430) * Y(:,3) * Y(:,150)) + (RC(:,434) * Y(:,3) * Y(:,162)) &
        + (RC(:,428) * Y(:,3) * Y(:,148)) + (RC(:,429) * Y(:,3) * Y(:,149)) &
        + (RC(:,426) * Y(:,3) * Y(:,146)) + (RC(:,427) * Y(:,3) * Y(:,147)) &
        + (RC(:,424) * Y(:,3) * Y(:,144)) + (RC(:,425) * Y(:,3) * Y(:,145)) &
        + (RC(:,362) * Y(:,6) * Y(:,46)) + (RC(:,374) * Y(:,6) * Y(:,109)) &
        + (RC(:,70) * Y(:,53) * Y(:,6)) + (RC(:,75) * Y(:,59) * Y(:,3)) &
        + (RC(:,61) * Y(:,6) * Y(:,43)) + (RC(:,65) * Y(:,47) * Y(:,6)) &
        + (RC(:,55) * Y(:,6) * Y(:,32)) + (RC(:,57) * Y(:,6) * Y(:,34)) &
        + (RC(:,33) * Y(:,9) * Y(:,5)) + (RC(:,53) * Y(:,6) * Y(:,30)) &
        + (RC(:,21) * Y(:,9) * Y(:,6)) + (RC(:,29) * Y(:,9) * Y(:,8)) &
        + (RC(:,16) * Y(:,1) *  H2O(:)*2.00) 
        L(:) = 0.0 &
        + (RC(:,481) * Y(:,201)) + (RC(:,484) * Y(:,203)) &
        + (RC(:,474) * Y(:,199)) + (RC(:,475) * Y(:,200)) + (RC(:,480) * Y(:,202)) &
        + (RC(:,465) * Y(:,185)) + (RC(:,466) * Y(:,192)) + (RC(:,473) * Y(:,198)) &
        + (RC(:,462) * Y(:,182)) + (RC(:,463) * Y(:,183)) + (RC(:,464) * Y(:,184)) &
        + (RC(:,459) * Y(:,179)) + (RC(:,460) * Y(:,180)) + (RC(:,461) * Y(:,181)) &
        + (RC(:,456) * Y(:,175)) + (RC(:,457) * Y(:,177)) + (RC(:,458) * Y(:,178)) &
        + (RC(:,453) * Y(:,153)) + (RC(:,454) * Y(:,174)) + (RC(:,455) * Y(:,176)) &
        + (RC(:,450) * Y(:,171)) + (RC(:,451) * Y(:,151)) + (RC(:,452) * Y(:,152)) &
        + (RC(:,447) * Y(:,158)) + (RC(:,448) * Y(:,169)) + (RC(:,449) * Y(:,170)) &
        + (RC(:,444) * Y(:,155)) + (RC(:,445) * Y(:,156)) + (RC(:,446) * Y(:,157)) &
        + (RC(:,441) * Y(:,172)) + (RC(:,442) * Y(:,173)) + (RC(:,443) * Y(:,154)) &
        + (RC(:,438) * Y(:,166)) + (RC(:,439) * Y(:,167)) + (RC(:,440) * Y(:,168)) &
        + (RC(:,435) * Y(:,163)) + (RC(:,436) * Y(:,164)) + (RC(:,437) * Y(:,165)) &
        + (RC(:,432) * Y(:,160)) + (RC(:,433) * Y(:,161)) + (RC(:,434) * Y(:,162)) &
        + (RC(:,429) * Y(:,149)) + (RC(:,430) * Y(:,150)) + (RC(:,431) * Y(:,159)) &
        + (RC(:,426) * Y(:,146)) + (RC(:,427) * Y(:,147)) + (RC(:,428) * Y(:,148)) &
        + (RC(:,423) * Y(:,144)) + (RC(:,424) * Y(:,144)) + (RC(:,425) * Y(:,145)) &
        + (RC(:,416) * Y(:,195)) + (RC(:,418) * Y(:,66)) + (RC(:,421) * Y(:,197)) &
        + (RC(:,411) * Y(:,142)) + (RC(:,412) * Y(:,143)) + (RC(:,413) * Y(:,63)) &
        + (RC(:,408) * Y(:,139)) + (RC(:,409) * Y(:,140)) + (RC(:,410) * Y(:,141)) &
        + (RC(:,405) * Y(:,136)) + (RC(:,406) * Y(:,137)) + (RC(:,407) * Y(:,138)) &
        + (RC(:,402) * Y(:,133)) + (RC(:,403) * Y(:,134)) + (RC(:,404) * Y(:,135)) &
        + (RC(:,399) * Y(:,130)) + (RC(:,400) * Y(:,131)) + (RC(:,401) * Y(:,132)) &
        + (RC(:,396) * Y(:,127)) + (RC(:,397) * Y(:,128)) + (RC(:,398) * Y(:,129)) &
        + (RC(:,393) * Y(:,124)) + (RC(:,394) * Y(:,125)) + (RC(:,395) * Y(:,126)) &
        + (RC(:,390) * Y(:,57)) + (RC(:,391) * Y(:,58)) + (RC(:,392) * Y(:,123)) &
        + (RC(:,385) * Y(:,193)) + (RC(:,386) * Y(:,120)) + (RC(:,389) * Y(:,52)) &
        + (RC(:,382) * Y(:,99)) + (RC(:,383) * Y(:,99)) + (RC(:,384) * Y(:,51)) &
        + (RC(:,379) * Y(:,96)) + (RC(:,380) * Y(:,97)) + (RC(:,381) * Y(:,97)) &
        + (RC(:,376) * Y(:,113)) + (RC(:,377) * Y(:,115)) + (RC(:,378) * Y(:,96)) &
        + (RC(:,370) * Y(:,190)) + (RC(:,371) * Y(:,191)) + (RC(:,372) * Y(:,109)) &
        + (RC(:,367) * Y(:,98)) + (RC(:,368) * Y(:,100)) + (RC(:,369) * Y(:,189)) &
        + (RC(:,360) * Y(:,46)) + (RC(:,364) * Y(:,102)) + (RC(:,366) * Y(:,60)) &
        + (RC(:,357) * Y(:,188)) + (RC(:,358) * Y(:,104)) + (RC(:,359) * Y(:,105)) &
        + (RC(:,354) * Y(:,187)) + (RC(:,355) * Y(:,88)) + (RC(:,356) * Y(:,111)) &
        + (RC(:,105) * Y(:,86)) + (RC(:,106) * Y(:,87)) + (RC(:,353) * Y(:,186)) &
        + (RC(:,102) * Y(:,83)) + (RC(:,103) * Y(:,84)) + (RC(:,104) * Y(:,85)) &
        + (RC(:,99) * Y(:,80)) + (RC(:,100) * Y(:,81)) + (RC(:,101) * Y(:,82)) &
        + (RC(:,96) * Y(:,79)) + (RC(:,97) * Y(:,40)) + (RC(:,98) * Y(:,41)) &
        + (RC(:,93) * Y(:,78)) + (RC(:,94) * Y(:,78)) + (RC(:,95) * Y(:,79)) &
        + (RC(:,90) * Y(:,76)) + (RC(:,91) * Y(:,77)) + (RC(:,92) * Y(:,77)) &
        + (RC(:,84) * Y(:,71)) + (RC(:,88) * Y(:,73)) + (RC(:,89) * Y(:,101)) &
        + (RC(:,81) * Y(:,67)) + (RC(:,82) * Y(:,39)) + (RC(:,83) * Y(:,42)) &
        + (RC(:,78) * Y(:,64)) + (RC(:,79) * Y(:,64)) + (RC(:,80) * Y(:,67)) &
        + (RC(:,75) * Y(:,59)) + (RC(:,76) * Y(:,61)) + (RC(:,77) * Y(:,61)) &
        + (RC(:,63) * Y(:,47)) + (RC(:,68) * Y(:,53)) + (RC(:,74) * Y(:,59)) &
        + (RC(:,48) * Y(:,32)) + (RC(:,49) * Y(:,34)) + (RC(:,59) * Y(:,43)) &
        + (RC(:,45) * Y(:,25)) + (RC(:,46) * Y(:,28)) + (RC(:,47) * Y(:,30)) &
        + (RC(:,42) * Y(:,21)) + (RC(:,43) * Y(:,23)) + (RC(:,44) * Y(:,25)) &
        + (RC(:,34) * Y(:,13)) + (RC(:,35) * Y(:,14)) + (RC(:,37) * Y(:,16)) &
        + (RC(:,27) * Y(:,4)) + (RC(:,28) * Y(:,5)) + (RC(:,32) * Y(:,15)) &
        + (RC(:,20) * Y(:,12)) + (RC(:,22) * Y(:,9)) + (RC(:,25) * Y(:,8)) &
        + (RC(:,17) * Y(:,6)) + (RC(:,18) * Y(:,10)) + (RC(:,19) * Y(:,11)) 
        Y(:,3) = P(:)/L(:)

        !          NO2              Y(:,4) 
        P(:) = (DJ(:,86) * Y(:,175)) + (DJ(:,96) * Y(:,185)) &
        + (DJ(:,80) * Y(:,170)) + (DJ(:,81) * Y(:,171)) &
        + (DJ(:,52) * Y(:,142)) + (DJ(:,79) * Y(:,169)) &
        + (DJ(:,50) * Y(:,137)) + (DJ(:,51) * Y(:,138)) &
        + (DJ(:,48) * Y(:,129)) + (DJ(:,49) * Y(:,136)) &
        + (DJ(:,46) * Y(:,127)) + (DJ(:,47) * Y(:,128)) &
        + (DJ(:,44) * Y(:,126)) + (DJ(:,45) * Y(:,127)) &
        + (DJ(:,42) * Y(:,124)) + (DJ(:,43) * Y(:,125)) &
        + (DJ(:,32) * Y(:,115)) + (DJ(:,41) * Y(:,123)) &
        + (DJ(:,8) * Y(:,14)) + (DJ(:,31) * Y(:,115)) &
        + (RC(:,484) * Y(:,3) * Y(:,203)) + (DJ(:,6) * Y(:,5)) &
        + (RC(:,481) * Y(:,3) * Y(:,201)) + (RC(:,483) * Y(:,203)) &
        + (RC(:,479) * Y(:,202)) + (RC(:,480) * Y(:,3) * Y(:,202)) &
        + (RC(:,475) * Y(:,3) * Y(:,200)) + (RC(:,477) * Y(:,201)) &
        + (RC(:,473) * Y(:,3) * Y(:,198)) + (RC(:,474) * Y(:,3) * Y(:,199)) &
        + (RC(:,470) * Y(:,199)) + (RC(:,472) * Y(:,200)) &
        + (RC(:,456) * Y(:,3) * Y(:,175)) + (RC(:,468) * Y(:,198)) &
        + (RC(:,449) * Y(:,3) * Y(:,170)) + (RC(:,450) * Y(:,3) * Y(:,171)) &
        + (RC(:,422) * Y(:,5) * Y(:,197)) + (RC(:,448) * Y(:,3) * Y(:,169)) &
        + (RC(:,417) * Y(:,5) * Y(:,195)) + (RC(:,421) * Y(:,3) * Y(:,197)) &
        + (RC(:,412) * Y(:,3) * Y(:,143)) + (RC(:,416) * Y(:,3) * Y(:,195)) &
        + (RC(:,410) * Y(:,3) * Y(:,141)) + (RC(:,411) * Y(:,3) * Y(:,142)) &
        + (RC(:,408) * Y(:,3) * Y(:,139)) + (RC(:,409) * Y(:,3) * Y(:,140)) &
        + (RC(:,406) * Y(:,3) * Y(:,137)) + (RC(:,407) * Y(:,3) * Y(:,138)) &
        + (RC(:,404) * Y(:,3) * Y(:,135)) + (RC(:,405) * Y(:,3) * Y(:,136)) &
        + (RC(:,402) * Y(:,3) * Y(:,133)) + (RC(:,403) * Y(:,3) * Y(:,134)) &
        + (RC(:,400) * Y(:,3) * Y(:,131)) + (RC(:,401) * Y(:,3) * Y(:,132)) &
        + (RC(:,398) * Y(:,3) * Y(:,129)) + (RC(:,399) * Y(:,3) * Y(:,130)) &
        + (RC(:,396) * Y(:,3) * Y(:,127)) + (RC(:,397) * Y(:,3) * Y(:,128)) &
        + (RC(:,394) * Y(:,3) * Y(:,125)) + (RC(:,395) * Y(:,3) * Y(:,126)) &
        + (RC(:,392) * Y(:,3) * Y(:,123)) + (RC(:,393) * Y(:,3) * Y(:,124)) &
        + (RC(:,352) * Y(:,55)) + (RC(:,377) * Y(:,3) * Y(:,115)) &
        + (RC(:,338) * Y(:,38)) + (RC(:,342) * Y(:,49)) &
        + (RC(:,336) * Y(:,36)) + (RC(:,337) * Y(:,37)) &
        + (RC(:,242) * Y(:,122) * Y(:,5)) &
        + (RC(:,243) * Y(:,55) * Y(:,5) * 2.00) &
        + (RC(:,240) * Y(:,54) * Y(:,5)) + (RC(:,241) * Y(:,56) * Y(:,5)) &
        + (RC(:,238) * Y(:,119) * Y(:,5)) + (RC(:,239) * Y(:,121) * Y(:,5)) &
        + (RC(:,236) * Y(:,117) * Y(:,5)) + (RC(:,237) * Y(:,118) * Y(:,5)) &
        + (RC(:,234) * Y(:,50) * Y(:,5)) + (RC(:,235) * Y(:,116) * Y(:,5)) &
        + (RC(:,232) * Y(:,48) * Y(:,5)) + (RC(:,233) * Y(:,49) * Y(:,5) * 2.00) &
        + (RC(:,230) * Y(:,45) * Y(:,5)) + (RC(:,231) * Y(:,114) * Y(:,5)) &
        + (RC(:,229) * Y(:,38) * Y(:,5) * 2.00) &
        + (RC(:,228) * Y(:,37) * Y(:,5) * 2.00) &
        + (RC(:,226) * Y(:,112) * Y(:,5)) &
        + (RC(:,227) * Y(:,36) * Y(:,5) * 2.00) &
        + (RC(:,224) * Y(:,112) * Y(:,5)) + (RC(:,225) * Y(:,112) * Y(:,5)) &
        + (RC(:,222) * Y(:,110) * Y(:,5)) + (RC(:,223) * Y(:,110) * Y(:,5)) &
        + (RC(:,220) * Y(:,44) * Y(:,5)) + (RC(:,221) * Y(:,44) * Y(:,5)) &
        + (RC(:,218) * Y(:,107) * Y(:,5)) + (RC(:,219) * Y(:,108) * Y(:,5)) &
        + (RC(:,216) * Y(:,74) * Y(:,5)) + (RC(:,217) * Y(:,75) * Y(:,5)) &
        + (RC(:,214) * Y(:,72) * Y(:,5)) + (RC(:,215) * Y(:,106) * Y(:,5)) &
        + (RC(:,212) * Y(:,92) * Y(:,5)) + (RC(:,213) * Y(:,70) * Y(:,5)) &
        + (RC(:,210) * Y(:,103) * Y(:,5)) + (RC(:,211) * Y(:,90) * Y(:,5)) &
        + (RC(:,208) * Y(:,35) * Y(:,5)) + (RC(:,209) * Y(:,95) * Y(:,5)) &
        + (RC(:,206) * Y(:,31) * Y(:,5)) + (RC(:,207) * Y(:,33) * Y(:,5)) &
        + (RC(:,204) * Y(:,69) * Y(:,5)) + (RC(:,205) * Y(:,31) * Y(:,5)) &
        + (RC(:,202) * Y(:,65) * Y(:,5)) + (RC(:,203) * Y(:,68) * Y(:,5)) &
        + (RC(:,200) * Y(:,62) * Y(:,5)) + (RC(:,201) * Y(:,65) * Y(:,5)) &
        + (RC(:,198) * Y(:,93) * Y(:,5)) + (RC(:,199) * Y(:,94) * Y(:,5)) &
        + (RC(:,196) * Y(:,89) * Y(:,5)) + (RC(:,197) * Y(:,91) * Y(:,5)) &
        + (RC(:,194) * Y(:,29) * Y(:,5)) + (RC(:,195) * Y(:,29) * Y(:,5)) &
        + (RC(:,192) * Y(:,27) * Y(:,5)) + (RC(:,193) * Y(:,26) * Y(:,5)) &
        + (RC(:,190) * Y(:,22) * Y(:,5)) + (RC(:,191) * Y(:,24) * Y(:,5)) &
        + (RC(:,163) * Y(:,122) * Y(:,8)) + (RC(:,165) * Y(:,217)) &
        + (RC(:,161) * Y(:,56) * Y(:,8)) + (RC(:,162) * Y(:,56) * Y(:,8)) &
        + (RC(:,160) * Y(:,55) * Y(:,8) * 2.00) &
        + (RC(:,158) * Y(:,54) * Y(:,8)) + (RC(:,159) * Y(:,54) * Y(:,8)) &
        + (RC(:,156) * Y(:,119) * Y(:,8)) + (RC(:,157) * Y(:,121) * Y(:,8)) &
        + (RC(:,154) * Y(:,117) * Y(:,8)) + (RC(:,155) * Y(:,118) * Y(:,8)) &
        + (RC(:,152) * Y(:,50) * Y(:,8)) + (RC(:,153) * Y(:,116) * Y(:,8)) &
        + (RC(:,151) * Y(:,49) * Y(:,8) * 2.00) &
        + (RC(:,149) * Y(:,48) * Y(:,8)) + (RC(:,150) * Y(:,48) * Y(:,8)) &
        + (RC(:,147) * Y(:,45) * Y(:,8)) + (RC(:,148) * Y(:,114) * Y(:,8)) &
        + (RC(:,146) * Y(:,38) * Y(:,8) * 2.00) &
        + (RC(:,145) * Y(:,37) * Y(:,8) * 2.00) &
        + (RC(:,143) * Y(:,112) * Y(:,8)) &
        + (RC(:,144) * Y(:,36) * Y(:,8) * 2.00) &
        + (RC(:,141) * Y(:,112) * Y(:,8)) + (RC(:,142) * Y(:,112) * Y(:,8)) &
        + (RC(:,139) * Y(:,110) * Y(:,8)) + (RC(:,140) * Y(:,110) * Y(:,8)) &
        + (RC(:,137) * Y(:,44) * Y(:,8)) + (RC(:,138) * Y(:,44) * Y(:,8)) &
        + (RC(:,135) * Y(:,107) * Y(:,8)) + (RC(:,136) * Y(:,108) * Y(:,8)) &
        + (RC(:,133) * Y(:,74) * Y(:,8)) + (RC(:,134) * Y(:,75) * Y(:,8)) &
        + (RC(:,131) * Y(:,72) * Y(:,8)) + (RC(:,132) * Y(:,106) * Y(:,8)) &
        + (RC(:,129) * Y(:,92) * Y(:,8)) + (RC(:,130) * Y(:,70) * Y(:,8)) &
        + (RC(:,127) * Y(:,103) * Y(:,8)) + (RC(:,128) * Y(:,90) * Y(:,8)) &
        + (RC(:,125) * Y(:,35) * Y(:,8)) + (RC(:,126) * Y(:,95) * Y(:,8)) &
        + (RC(:,123) * Y(:,31) * Y(:,8)) + (RC(:,124) * Y(:,33) * Y(:,8)) &
        + (RC(:,121) * Y(:,69) * Y(:,8)) + (RC(:,122) * Y(:,31) * Y(:,8)) &
        + (RC(:,119) * Y(:,65) * Y(:,8)) + (RC(:,120) * Y(:,68) * Y(:,8)) &
        + (RC(:,117) * Y(:,62) * Y(:,8)) + (RC(:,118) * Y(:,65) * Y(:,8)) &
        + (RC(:,115) * Y(:,93) * Y(:,8)) + (RC(:,116) * Y(:,94) * Y(:,8)) &
        + (RC(:,113) * Y(:,89) * Y(:,8)) + (RC(:,114) * Y(:,91) * Y(:,8)) &
        + (RC(:,111) * Y(:,29) * Y(:,8)) + (RC(:,112) * Y(:,29) * Y(:,8)) &
        + (RC(:,109) * Y(:,27) * Y(:,8)) + (RC(:,110) * Y(:,26) * Y(:,8)) &
        + (RC(:,107) * Y(:,22) * Y(:,8)) + (RC(:,108) * Y(:,24) * Y(:,8)) &
        + (RC(:,33) * Y(:,9) * Y(:,5)) + (RC(:,34) * Y(:,3) * Y(:,13)) &
        + (RC(:,31) * Y(:,15)) + (RC(:,32) * Y(:,3) * Y(:,15)) &
        + (RC(:,28) * Y(:,3) * Y(:,5)) + (RC(:,29) * Y(:,9) * Y(:,8)) &
        + (RC(:,13) * Y(:,4) * Y(:,5)) + (RC(:,15) * Y(:,7)) &
        + (RC(:,12) * Y(:,8) * Y(:,5) * 2.00) &
        + (RC(:,11) * Y(:,8) * Y(:,8) * 2.00) &
        + (RC(:,4) * Y(:,2) * Y(:,8)) + (RC(:,9) * Y(:,8) * Y(:,6)) 
        !     
        !          
        L(:) = 0.0 &
        + (RC(:,478) * Y(:,112)) + (RC(:,482) * Y(:,50)) + (DJ(:,4)) &
        + (RC(:,469) * Y(:,72)) + (RC(:,471) * Y(:,106)) + (RC(:,476) * Y(:,110)) &
        + (RC(:,415) * Y(:,194)) + (RC(:,420) * Y(:,196)) + (RC(:,467) * Y(:,70)) &
        + (RC(:,27) * Y(:,3)) + (RC(:,30) * Y(:,9)) + (RC(:,164) * Y(:,22)) &
        + (RC(:,13) * Y(:,5)) + (RC(:,14) * Y(:,5)) + (RC(:,26)) &
        + (RC(:,5) * Y(:,2)) + (RC(:,6) * Y(:,2)) + (RC(:,10) * Y(:,6)) 
        Y(:,4) = (YP(:,4) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NO3              Y(:,5) 
        P(:) = (RC(:,15) * Y(:,7)) + (RC(:,35) * Y(:,3) * Y(:,14)) &
        + (RC(:,6) * Y(:,2) * Y(:,4)) + (RC(:,10) * Y(:,4) * Y(:,6)) 
        L(:) = 0.0 &
        + (DJ(:,6)) &
        + (RC(:,419) * Y(:,66)) + (RC(:,422) * Y(:,197)) + (DJ(:,5)) &
        + (RC(:,388) * Y(:,120)) + (RC(:,414) * Y(:,63)) + (RC(:,417) * Y(:,195)) &
        + (RC(:,365) * Y(:,102)) + (RC(:,373) * Y(:,109)) + (RC(:,387) * Y(:,51)) &
        + (RC(:,242) * Y(:,122)) + (RC(:,243) * Y(:,55)) + (RC(:,361) * Y(:,46)) &
        + (RC(:,239) * Y(:,121)) + (RC(:,240) * Y(:,54)) + (RC(:,241) * Y(:,56)) &
        + (RC(:,236) * Y(:,117)) + (RC(:,237) * Y(:,118)) + (RC(:,238) * Y(:,119)) &
        + (RC(:,233) * Y(:,49)) + (RC(:,234) * Y(:,50)) + (RC(:,235) * Y(:,116)) &
        + (RC(:,230) * Y(:,45)) + (RC(:,231) * Y(:,114)) + (RC(:,232) * Y(:,48)) &
        + (RC(:,227) * Y(:,36)) + (RC(:,228) * Y(:,37)) + (RC(:,229) * Y(:,38)) &
        + (RC(:,224) * Y(:,112)) + (RC(:,225) * Y(:,112)) + (RC(:,226) * Y(:,112)) &
        + (RC(:,221) * Y(:,44)) + (RC(:,222) * Y(:,110)) + (RC(:,223) * Y(:,110)) &
        + (RC(:,218) * Y(:,107)) + (RC(:,219) * Y(:,108)) + (RC(:,220) * Y(:,44)) &
        + (RC(:,215) * Y(:,106)) + (RC(:,216) * Y(:,74)) + (RC(:,217) * Y(:,75)) &
        + (RC(:,212) * Y(:,92)) + (RC(:,213) * Y(:,70)) + (RC(:,214) * Y(:,72)) &
        + (RC(:,209) * Y(:,95)) + (RC(:,210) * Y(:,103)) + (RC(:,211) * Y(:,90)) &
        + (RC(:,206) * Y(:,31)) + (RC(:,207) * Y(:,33)) + (RC(:,208) * Y(:,35)) &
        + (RC(:,203) * Y(:,68)) + (RC(:,204) * Y(:,69)) + (RC(:,205) * Y(:,31)) &
        + (RC(:,200) * Y(:,62)) + (RC(:,201) * Y(:,65)) + (RC(:,202) * Y(:,65)) &
        + (RC(:,197) * Y(:,91)) + (RC(:,198) * Y(:,93)) + (RC(:,199) * Y(:,94)) &
        + (RC(:,194) * Y(:,29)) + (RC(:,195) * Y(:,29)) + (RC(:,196) * Y(:,89)) &
        + (RC(:,191) * Y(:,24)) + (RC(:,192) * Y(:,27)) + (RC(:,193) * Y(:,26)) &
        + (RC(:,86) * Y(:,42)) + (RC(:,87) * Y(:,71)) + (RC(:,190) * Y(:,22)) &
        + (RC(:,64) * Y(:,47)) + (RC(:,69) * Y(:,53)) + (RC(:,85) * Y(:,39)) &
        + (RC(:,51) * Y(:,32)) + (RC(:,52) * Y(:,34)) + (RC(:,60) * Y(:,43)) &
        + (RC(:,28) * Y(:,3)) + (RC(:,33) * Y(:,9)) + (RC(:,50) * Y(:,30)) &
        + (RC(:,12) * Y(:,8)) + (RC(:,13) * Y(:,4)) + (RC(:,14) * Y(:,4)) 
        Y(:,5) = (YP(:,5) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          O3               Y(:,6) 
        P(:) = (RC(:,1) * Y(:,2)) + (RC(:,2) * Y(:,2)) 
        !     

        L(:) = 0.0 &
        + (DJ(:,1)) + (DJ(:,2)) &
        + (RC(:,363) * Y(:,46)) + (RC(:,374) * Y(:,109)) + (RC(:,375) * Y(:,109)) &
        + (RC(:,72) * Y(:,53)) + (RC(:,73) * Y(:,53)) + (RC(:,362) * Y(:,46)) &
        + (RC(:,67) * Y(:,47)) + (RC(:,70) * Y(:,53)) + (RC(:,71) * Y(:,53)) &
        + (RC(:,62) * Y(:,43)) + (RC(:,65) * Y(:,47)) + (RC(:,66) * Y(:,47)) &
        + (RC(:,57) * Y(:,34)) + (RC(:,58) * Y(:,34)) + (RC(:,61) * Y(:,43)) &
        + (RC(:,54) * Y(:,30)) + (RC(:,55) * Y(:,32)) + (RC(:,56) * Y(:,32)) &
        + (RC(:,17) * Y(:,3)) + (RC(:,21) * Y(:,9)) + (RC(:,53) * Y(:,30)) &
        + (RC(:,3) * Y(:,2)) + (RC(:,9) * Y(:,8)) + (RC(:,10) * Y(:,4)) 
        Y(:,6) = (YP(:,6) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          N2O5             Y(:,7) 
        P(:) = (RC(:,14) * Y(:,4) * Y(:,5)) 
        !   

        L(:) = 0.0 &
        + (RC(:,15)) + (RC(:,40)) 
        Y(:,7) = (YP(:,7) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NO               Y(:,8) 
        P(:) = (DJ(:,7) * Y(:,13)) &
        + (DJ(:,4) * Y(:,4)) + (DJ(:,5) * Y(:,5)) &
        + (RC(:,5) * Y(:,2) * Y(:,4)) + (RC(:,13) * Y(:,4) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,189) * Y(:,122)) &
        + (RC(:,186) * Y(:,116)) + (RC(:,187) * Y(:,54)) + (RC(:,188) * Y(:,56)) &
        + (RC(:,183) * Y(:,68)) + (RC(:,184) * Y(:,69)) + (RC(:,185) * Y(:,48)) &
        + (RC(:,180) * Y(:,44)) + (RC(:,181) * Y(:,62)) + (RC(:,182) * Y(:,65)) &
        + (RC(:,177) * Y(:,103)) + (RC(:,178) * Y(:,90)) + (RC(:,179) * Y(:,92)) &
        + (RC(:,174) * Y(:,33)) + (RC(:,175) * Y(:,35)) + (RC(:,176) * Y(:,95)) &
        + (RC(:,171) * Y(:,89)) + (RC(:,172) * Y(:,91)) + (RC(:,173) * Y(:,31)) &
        + (RC(:,168) * Y(:,27)) + (RC(:,169) * Y(:,26)) + (RC(:,170) * Y(:,29)) &
        + (RC(:,163) * Y(:,122)) + (RC(:,166) * Y(:,22)) + (RC(:,167) * Y(:,24)) &
        + (RC(:,160) * Y(:,55)) + (RC(:,161) * Y(:,56)) + (RC(:,162) * Y(:,56)) &
        + (RC(:,157) * Y(:,121)) + (RC(:,158) * Y(:,54)) + (RC(:,159) * Y(:,54)) &
        + (RC(:,154) * Y(:,117)) + (RC(:,155) * Y(:,118)) + (RC(:,156) * Y(:,119)) &
        + (RC(:,151) * Y(:,49)) + (RC(:,152) * Y(:,50)) + (RC(:,153) * Y(:,116)) &
        + (RC(:,148) * Y(:,114)) + (RC(:,149) * Y(:,48)) + (RC(:,150) * Y(:,48)) &
        + (RC(:,145) * Y(:,37)) + (RC(:,146) * Y(:,38)) + (RC(:,147) * Y(:,45)) &
        + (RC(:,142) * Y(:,112)) + (RC(:,143) * Y(:,112)) + (RC(:,144) * Y(:,36)) &
        + (RC(:,139) * Y(:,110)) + (RC(:,140) * Y(:,110)) + (RC(:,141) * Y(:,112)) &
        + (RC(:,136) * Y(:,108)) + (RC(:,137) * Y(:,44)) + (RC(:,138) * Y(:,44)) &
        + (RC(:,133) * Y(:,74)) + (RC(:,134) * Y(:,75)) + (RC(:,135) * Y(:,107)) &
        + (RC(:,130) * Y(:,70)) + (RC(:,131) * Y(:,72)) + (RC(:,132) * Y(:,106)) &
        + (RC(:,127) * Y(:,103)) + (RC(:,128) * Y(:,90)) + (RC(:,129) * Y(:,92)) &
        + (RC(:,124) * Y(:,33)) + (RC(:,125) * Y(:,35)) + (RC(:,126) * Y(:,95)) &
        + (RC(:,121) * Y(:,69)) + (RC(:,122) * Y(:,31)) + (RC(:,123) * Y(:,31)) &
        + (RC(:,118) * Y(:,65)) + (RC(:,119) * Y(:,65)) + (RC(:,120) * Y(:,68)) &
        + (RC(:,115) * Y(:,93)) + (RC(:,116) * Y(:,94)) + (RC(:,117) * Y(:,62)) &
        + (RC(:,112) * Y(:,29)) + (RC(:,113) * Y(:,89)) + (RC(:,114) * Y(:,91)) &
        + (RC(:,109) * Y(:,27)) + (RC(:,110) * Y(:,26)) + (RC(:,111) * Y(:,29)) &
        + (RC(:,29) * Y(:,9)) + (RC(:,107) * Y(:,22)) + (RC(:,108) * Y(:,24)) &
        + (RC(:,11) * Y(:,8)) + (RC(:,12) * Y(:,5)) + (RC(:,25) * Y(:,3)) &
        + (RC(:,4) * Y(:,2)) + (RC(:,9) * Y(:,6)) + (RC(:,11) * Y(:,8)) 
        Y(:,8) = (YP(:,8) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HO2              Y(:,9) 
        P(:) = (DJ(:,94) * Y(:,183)) &
        + (DJ(:,91) * Y(:,180)) + (DJ(:,93) * Y(:,182)) &
        + (DJ(:,84) * Y(:,153)) + (DJ(:,85) * Y(:,174)) &
        + (DJ(:,82) * Y(:,151)) + (DJ(:,83) * Y(:,152)) &
        + (DJ(:,77) * Y(:,157)) + (DJ(:,78) * Y(:,158)) &
        + (DJ(:,75) * Y(:,155)) + (DJ(:,76) * Y(:,156)) &
        + (DJ(:,73) * Y(:,173)) + (DJ(:,74) * Y(:,154)) &
        + (DJ(:,70) * Y(:,167)) + (DJ(:,72) * Y(:,172)) &
        + (DJ(:,68) * Y(:,166)) + (DJ(:,69) * Y(:,166)) &
        + (DJ(:,58) * Y(:,148)) + (DJ(:,63) * Y(:,161)) &
        + (DJ(:,55) * Y(:,146)) + (DJ(:,56) * Y(:,147)) &
        + (DJ(:,53) * Y(:,144)) + (DJ(:,54) * Y(:,145)) &
        + (DJ(:,51) * Y(:,138)) + (DJ(:,52) * Y(:,142)) &
        + (DJ(:,49) * Y(:,136)) + (DJ(:,50) * Y(:,137)) &
        + (DJ(:,44) * Y(:,126)) + (DJ(:,46) * Y(:,127)) &
        + (DJ(:,42) * Y(:,124)) + (DJ(:,43) * Y(:,125)) &
        + (DJ(:,39) * Y(:,51)) + (DJ(:,41) * Y(:,123)) &
        + (DJ(:,37) * Y(:,99)) + (DJ(:,38) * Y(:,99)) &
        + (DJ(:,35) * Y(:,97)) + (DJ(:,36) * Y(:,97)) &
        + (DJ(:,33) * Y(:,96)) + (DJ(:,34) * Y(:,96) * 2.00) &
        + (DJ(:,30) * Y(:,113) * 2.00) &
        + (DJ(:,25) * Y(:,98)) + (DJ(:,29) * Y(:,109)) &
        + (DJ(:,23) * Y(:,46)) + (DJ(:,24) * Y(:,60) * 2.00) &
        + (DJ(:,22) * Y(:,102) * 2.00) &
        + (DJ(:,20) * Y(:,104)) + (DJ(:,21) * Y(:,105)) &
        + (DJ(:,18) * Y(:,111)) + (DJ(:,19) * Y(:,188)) &
        + (DJ(:,11) * Y(:,42)) + (DJ(:,12) * Y(:,71)) &
        + (RC(:,379) * Y(:,3) * Y(:,96)) + (DJ(:,9) * Y(:,39) * 2.00) &
        + (RC(:,357) * Y(:,3) * Y(:,188)) + (RC(:,366) * Y(:,3) * Y(:,60)) &
        + (RC(:,350) * Y(:,56)) + (RC(:,356) * Y(:,3) * Y(:,111)) &
        + (RC(:,347) * Y(:,119)) + (RC(:,349) * Y(:,54)) &
        + (RC(:,340) * Y(:,114)) + (RC(:,341) * Y(:,48)) &
        + (RC(:,335) * Y(:,112)) + (RC(:,339) * Y(:,45)) &
        + (RC(:,332) * Y(:,110)) + (RC(:,334) * Y(:,112)) &
        + (RC(:,329) * Y(:,44)) + (RC(:,330) * Y(:,44)) &
        + (RC(:,321) * Y(:,92)) + (RC(:,324) * Y(:,106)) &
        + (RC(:,319) * Y(:,103)) + (RC(:,320) * Y(:,90)) &
        + (RC(:,317) * Y(:,35)) + (RC(:,318) * Y(:,95)) &
        + (RC(:,315) * Y(:,31)) + (RC(:,316) * Y(:,33)) &
        + (RC(:,311) * Y(:,69)) + (RC(:,314) * Y(:,31)) &
        + (RC(:,309) * Y(:,65)) + (RC(:,310) * Y(:,68)) &
        + (RC(:,307) * Y(:,62)) + (RC(:,308) * Y(:,65)) &
        + (RC(:,300) * Y(:,26)) + (RC(:,304) * Y(:,29)) &
        + (RC(:,294) * Y(:,24)) + (RC(:,297) * Y(:,27)) &
        + (RC(:,241) * Y(:,56) * Y(:,5)) + (RC(:,291) * Y(:,22)) &
        + (RC(:,238) * Y(:,119) * Y(:,5)) + (RC(:,240) * Y(:,54) * Y(:,5)) &
        + (RC(:,231) * Y(:,114) * Y(:,5)) + (RC(:,232) * Y(:,48) * Y(:,5)) &
        + (RC(:,226) * Y(:,112) * Y(:,5)) + (RC(:,230) * Y(:,45) * Y(:,5)) &
        + (RC(:,223) * Y(:,110) * Y(:,5)) + (RC(:,225) * Y(:,112) * Y(:,5)) &
        + (RC(:,220) * Y(:,44) * Y(:,5)) + (RC(:,221) * Y(:,44) * Y(:,5)) &
        + (RC(:,212) * Y(:,92) * Y(:,5)) + (RC(:,215) * Y(:,106) * Y(:,5)) &
        + (RC(:,210) * Y(:,103) * Y(:,5)) + (RC(:,211) * Y(:,90) * Y(:,5)) &
        + (RC(:,208) * Y(:,35) * Y(:,5)) + (RC(:,209) * Y(:,95) * Y(:,5)) &
        + (RC(:,206) * Y(:,31) * Y(:,5)) + (RC(:,207) * Y(:,33) * Y(:,5)) &
        + (RC(:,204) * Y(:,69) * Y(:,5)) + (RC(:,205) * Y(:,31) * Y(:,5)) &
        + (RC(:,202) * Y(:,65) * Y(:,5)) + (RC(:,203) * Y(:,68) * Y(:,5)) &
        + (RC(:,200) * Y(:,62) * Y(:,5)) + (RC(:,201) * Y(:,65) * Y(:,5)) &
        + (RC(:,193) * Y(:,26) * Y(:,5)) + (RC(:,195) * Y(:,29) * Y(:,5)) &
        + (RC(:,191) * Y(:,24) * Y(:,5)) + (RC(:,192) * Y(:,27) * Y(:,5)) &
        + (RC(:,161) * Y(:,56) * Y(:,8)) + (RC(:,190) * Y(:,22) * Y(:,5)) &
        + (RC(:,156) * Y(:,119) * Y(:,8)) + (RC(:,158) * Y(:,54) * Y(:,8)) &
        + (RC(:,148) * Y(:,114) * Y(:,8)) + (RC(:,149) * Y(:,48) * Y(:,8)) &
        + (RC(:,143) * Y(:,112) * Y(:,8)) + (RC(:,147) * Y(:,45) * Y(:,8)) &
        + (RC(:,140) * Y(:,110) * Y(:,8)) + (RC(:,142) * Y(:,112) * Y(:,8)) &
        + (RC(:,137) * Y(:,44) * Y(:,8)) + (RC(:,138) * Y(:,44) * Y(:,8)) &
        + (RC(:,129) * Y(:,92) * Y(:,8)) + (RC(:,132) * Y(:,106) * Y(:,8)) &
        + (RC(:,127) * Y(:,103) * Y(:,8)) + (RC(:,128) * Y(:,90) * Y(:,8)) &
        + (RC(:,125) * Y(:,35) * Y(:,8)) + (RC(:,126) * Y(:,95) * Y(:,8)) &
        + (RC(:,123) * Y(:,31) * Y(:,8)) + (RC(:,124) * Y(:,33) * Y(:,8)) &
        + (RC(:,121) * Y(:,69) * Y(:,8)) + (RC(:,122) * Y(:,31) * Y(:,8)) &
        + (RC(:,119) * Y(:,65) * Y(:,8)) + (RC(:,120) * Y(:,68) * Y(:,8)) &
        + (RC(:,117) * Y(:,62) * Y(:,8)) + (RC(:,118) * Y(:,65) * Y(:,8)) &
        + (RC(:,110) * Y(:,26) * Y(:,8)) + (RC(:,112) * Y(:,29) * Y(:,8)) &
        + (RC(:,108) * Y(:,24) * Y(:,8)) + (RC(:,109) * Y(:,27) * Y(:,8)) &
        + (RC(:,97) * Y(:,40) * Y(:,3)) + (RC(:,107) * Y(:,22) * Y(:,8)) &
        + (RC(:,93) * Y(:,78) * Y(:,3)) + (RC(:,95) * Y(:,3) * Y(:,79)) &
        + (RC(:,90) * Y(:,3) * Y(:,76)) + (RC(:,91) * Y(:,3) * Y(:,77)) &
        + (RC(:,82) * Y(:,3) * Y(:,39)) + (RC(:,85) * Y(:,5) * Y(:,39)) &
        + (RC(:,77) * Y(:,61) * Y(:,3)) + (RC(:,79) * Y(:,64) * Y(:,3)) &
        + (RC(:,61) * Y(:,6) * Y(:,43)) + (RC(:,74) * Y(:,59) * Y(:,3)) &
        + (RC(:,38) * Y(:,18)) + (RC(:,53) * Y(:,6) * Y(:,30)) &
        + (RC(:,28) * Y(:,3) * Y(:,5)) + (RC(:,31) * Y(:,15)) &
        + (RC(:,19) * Y(:,3) * Y(:,11)) + (RC(:,20) * Y(:,3) * Y(:,12)) &
        + (RC(:,17) * Y(:,3) * Y(:,6)) + (RC(:,18) * Y(:,3) * Y(:,10)) 
        L(:) = 0.0 &
        + (RC(:,289) * Y(:,122)) + (RC(:,290) * Y(:,55)) &
        + (RC(:,286) * Y(:,121)) + (RC(:,287) * Y(:,54)) + (RC(:,288) * Y(:,56)) &
        + (RC(:,283) * Y(:,117)) + (RC(:,284) * Y(:,118)) + (RC(:,285) * Y(:,119)) &
        + (RC(:,280) * Y(:,49)) + (RC(:,281) * Y(:,50)) + (RC(:,282) * Y(:,116)) &
        + (RC(:,277) * Y(:,45)) + (RC(:,278) * Y(:,114)) + (RC(:,279) * Y(:,48)) &
        + (RC(:,274) * Y(:,36)) + (RC(:,275) * Y(:,37)) + (RC(:,276) * Y(:,38)) &
        + (RC(:,271) * Y(:,44)) + (RC(:,272) * Y(:,110)) + (RC(:,273) * Y(:,112)) &
        + (RC(:,268) * Y(:,75)) + (RC(:,269) * Y(:,107)) + (RC(:,270) * Y(:,108)) &
        + (RC(:,265) * Y(:,72)) + (RC(:,266) * Y(:,106)) + (RC(:,267) * Y(:,74)) &
        + (RC(:,262) * Y(:,90)) + (RC(:,263) * Y(:,92)) + (RC(:,264) * Y(:,70)) &
        + (RC(:,259) * Y(:,35)) + (RC(:,260) * Y(:,95)) + (RC(:,261) * Y(:,103)) &
        + (RC(:,256) * Y(:,69)) + (RC(:,257) * Y(:,31)) + (RC(:,258) * Y(:,33)) &
        + (RC(:,253) * Y(:,62)) + (RC(:,254) * Y(:,65)) + (RC(:,255) * Y(:,68)) &
        + (RC(:,250) * Y(:,91)) + (RC(:,251) * Y(:,93)) + (RC(:,252) * Y(:,94)) &
        + (RC(:,247) * Y(:,26)) + (RC(:,248) * Y(:,29)) + (RC(:,249) * Y(:,89)) &
        + (RC(:,244) * Y(:,22)) + (RC(:,245) * Y(:,24)) + (RC(:,246) * Y(:,27)) &
        + (RC(:,29) * Y(:,8)) + (RC(:,30) * Y(:,4)) + (RC(:,33) * Y(:,5)) &
        + (RC(:,23) * Y(:,9)) + (RC(:,24) * Y(:,9)) + (RC(:,24) * Y(:,9)) &
        + (RC(:,21) * Y(:,6)) + (RC(:,22) * Y(:,3)) + (RC(:,23) * Y(:,9)) 
        Y(:,9) = (YP(:,9) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          H2               Y(:,10) 
        P(:) = (DJ(:,10) * Y(:,39)) 
        L(:) = 0.0 &
        + (RC(:,18) * Y(:,3)) 
        Y(:,10) = (YP(:,10) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CO               Y(:,11) 
        P(:) = (DJ(:,92) * Y(:,181)) &
        + (DJ(:,40) * Y(:,120)) + (DJ(:,73) * Y(:,173)) &
        + (DJ(:,30) * Y(:,113) * 2.00) &
        + (DJ(:,25) * Y(:,98)) + (DJ(:,29) * Y(:,109)) &
        + (DJ(:,24) * Y(:,60) * 2.00) &
        + (DJ(:,12) * Y(:,71)) + (DJ(:,22) * Y(:,102)) &
        + (DJ(:,10) * Y(:,39)) + (DJ(:,11) * Y(:,42)) &
        + (RC(:,480) * Y(:,3) * Y(:,202)) + (DJ(:,9) * Y(:,39)) &
        + (RC(:,474) * Y(:,3) * Y(:,199)) + (RC(:,475) * Y(:,3) * Y(:,200)) &
        + (RC(:,442) * Y(:,3) * Y(:,173)) + (RC(:,473) * Y(:,3) * Y(:,198)) &
        + (RC(:,367) * Y(:,3) * Y(:,98)) + (RC(:,374) * Y(:,6) * Y(:,109)) &
        + (RC(:,362) * Y(:,6) * Y(:,46)) + (RC(:,366) * Y(:,3) * Y(:,60) * 2.00) &
        + (RC(:,340) * Y(:,114)) + (RC(:,348) * Y(:,121)) &
        + (RC(:,231) * Y(:,114) * Y(:,5)) + (RC(:,239) * Y(:,121) * Y(:,5)) &
        + (RC(:,157) * Y(:,121) * Y(:,8)) + (RC(:,223) * Y(:,110) * Y(:,5)) &
        + (RC(:,140) * Y(:,110) * Y(:,8)) + (RC(:,148) * Y(:,114) * Y(:,8)) &
        + (RC(:,82) * Y(:,3) * Y(:,39)) + (RC(:,85) * Y(:,5) * Y(:,39)) &
        + (RC(:,73) * Y(:,53) * Y(:,6)) + (RC(:,74) * Y(:,59) * Y(:,3)) &
        + (RC(:,57) * Y(:,6) * Y(:,34)) + (RC(:,61) * Y(:,6) * Y(:,43)) &
        + (RC(:,53) * Y(:,6) * Y(:,30)) + (RC(:,55) * Y(:,6) * Y(:,32)) 
        L(:) = 0.0 &
        + (RC(:,19) * Y(:,3)) 
        Y(:,11) = (YP(:,11) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          H2O2             Y(:,12) 
        P(:) = (RC(:,363) * Y(:,6) * Y(:,46)) + (RC(:,375) * Y(:,6) * Y(:,109)) &
        + (RC(:,66) * Y(:,47) * Y(:,6)) + (RC(:,71) * Y(:,53) * Y(:,6)) &
        + (RC(:,23) * Y(:,9) * Y(:,9)) + (RC(:,24) * Y(:,9) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,20) * Y(:,3)) + (DJ(:,3)) 
        Y(:,12) = (YP(:,12) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HONO             Y(:,13) 
        P(:) = (RC(:,25) * Y(:,3) * Y(:,8)) + (RC(:,26) * Y(:,4)) 
        L(:) = 0.0 &
        + (RC(:,34) * Y(:,3)) + (DJ(:,7)) 
        Y(:,13) = (YP(:,13) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HNO3             Y(:,14) 
        P(:) = (RC(:,422) * Y(:,5) * Y(:,197)) &
        + (RC(:,417) * Y(:,5) * Y(:,195)) + (RC(:,419) * Y(:,5) * Y(:,66)) &
        + (RC(:,388) * Y(:,5) * Y(:,120)) + (RC(:,414) * Y(:,5) * Y(:,63)) &
        + (RC(:,373) * Y(:,5) * Y(:,109)) + (RC(:,387) * Y(:,5) * Y(:,51)) &
        + (RC(:,361) * Y(:,5) * Y(:,46)) + (RC(:,365) * Y(:,5) * Y(:,102)) &
        + (RC(:,86) * Y(:,5) * Y(:,42)) + (RC(:,87) * Y(:,5) * Y(:,71)) &
        + (RC(:,27) * Y(:,3) * Y(:,4)) + (RC(:,85) * Y(:,5) * Y(:,39)) 
        L(:) = 0.0 &
        + (RC(:,35) * Y(:,3)) + (RC(:,39)) + (DJ(:,8)) 
        Y(:,14) = (YP(:,14) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HO2NO2           Y(:,15) 
        P(:) = (RC(:,30) * Y(:,9) * Y(:,4)) 
        L(:) = 0.0 &
        + (RC(:,31)) + (RC(:,32) * Y(:,3)) 
        Y(:,15) = (YP(:,15) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          SO2              Y(:,16) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,36) * Y(:,2)) + (RC(:,37) * Y(:,3)) 
        Y(:,16) = (YP(:,16) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          SO3              Y(:,17) 
        P(:) = (RC(:,36) * Y(:,2) * Y(:,16)) + (RC(:,38) * Y(:,18)) 
        L(:) = 0.0 &
        + (RC(:,41)) 
        Y(:,17) = (YP(:,17) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HSO3             Y(:,18) 
        P(:) = (RC(:,37) * Y(:,3) * Y(:,16)) 
        L(:) = 0.0 &
        + (RC(:,38)) 
        Y(:,18) = (YP(:,18) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NA               Y(:,19) 
        P(:) = (RC(:,40) * Y(:,7)) &
        + (RC(:,39) * Y(:,14)) + (RC(:,40) * Y(:,7)) 
        L(:) = 0.0 
        Y(:,19) = (YP(:,19) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          SA               Y(:,20) 
        P(:) = (RC(:,41) * Y(:,17)) 
        L(:) = 0.0 
        Y(:,20) = (YP(:,20) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH4              Y(:,21) 
        P(:) = 0.0
        L(:) = 0.0 &
        + (RC(:,42) * Y(:,3)) 
        Y(:,21) = (YP(:,21) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3O2            Y(:,22) 
        P(:) = (DJ(:,61) * Y(:,159)) &
        + (DJ(:,13) * Y(:,73)) + (DJ(:,36) * Y(:,97)) &
        + (RC(:,423) * Y(:,3) * Y(:,144)) + (DJ(:,11) * Y(:,42)) &
        + (RC(:,322) * Y(:,70)) + (RC(:,381) * Y(:,3) * Y(:,97)) &
        + (RC(:,165) * Y(:,217)) + (RC(:,213) * Y(:,70) * Y(:,5)) &
        + (RC(:,101) * Y(:,3) * Y(:,82)) + (RC(:,130) * Y(:,70) * Y(:,8)) &
        + (RC(:,99) * Y(:,3) * Y(:,80)) + (RC(:,100) * Y(:,3) * Y(:,81)) &
        + (RC(:,57) * Y(:,6) * Y(:,34)) + (RC(:,98) * Y(:,41) * Y(:,3)) &
        + (RC(:,42) * Y(:,3) * Y(:,21)) + (RC(:,55) * Y(:,6) * Y(:,32)) 
        L(:) = 0.0 &
        + (RC(:,292)) + (RC(:,293)) &
        + (RC(:,190) * Y(:,5)) + (RC(:,244) * Y(:,9)) + (RC(:,291)) &
        + (RC(:,107) * Y(:,8)) + (RC(:,164) * Y(:,4)) + (RC(:,166) * Y(:,8)) 
        Y(:,22) = (YP(:,22) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H6             Y(:,23) 
        P(:) = 0.0
        L(:) = 0.0 &
        + (RC(:,43) * Y(:,3)) 
        Y(:,23) = (YP(:,23) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H5O2           Y(:,24) 
        P(:) = (DJ(:,64) * Y(:,162)) &
        + (DJ(:,57) * Y(:,148)) + (DJ(:,62) * Y(:,160)) &
        + (DJ(:,38) * Y(:,99)) + (DJ(:,45) * Y(:,127)) &
        + (DJ(:,17) * Y(:,88)) + (DJ(:,33) * Y(:,96)) &
        + (DJ(:,12) * Y(:,71)) + (DJ(:,14) * Y(:,101)) &
        + (RC(:,378) * Y(:,3) * Y(:,96)) + (RC(:,383) * Y(:,3) * Y(:,99)) &
        + (RC(:,303) * Y(:,29)) + (RC(:,323) * Y(:,72)) &
        + (RC(:,194) * Y(:,29) * Y(:,5)) + (RC(:,214) * Y(:,72) * Y(:,5)) &
        + (RC(:,111) * Y(:,29) * Y(:,8)) + (RC(:,131) * Y(:,72) * Y(:,8)) &
        + (RC(:,43) * Y(:,3) * Y(:,23)) + (RC(:,102) * Y(:,3) * Y(:,83)) 
        L(:) = 0.0 &
        + (RC(:,296)) &
        + (RC(:,245) * Y(:,9)) + (RC(:,294)) + (RC(:,295)) &
        + (RC(:,108) * Y(:,8)) + (RC(:,167) * Y(:,8)) + (RC(:,191) * Y(:,5)) 
        Y(:,24) = (YP(:,24) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C3H8             Y(:,25) 
        P(:) = 0.0
        L(:) = 0.0 &
        + (RC(:,44) * Y(:,3)) + (RC(:,45) * Y(:,3)) 
        Y(:,25) = (YP(:,25) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          IC3H7O2          Y(:,26) 
        P(:) = (RC(:,44) * Y(:,3) * Y(:,25)) 
        L(:) = 0.0 &
        + (RC(:,302)) &
        + (RC(:,247) * Y(:,9)) + (RC(:,300)) + (RC(:,301)) &
        + (RC(:,110) * Y(:,8)) + (RC(:,169) * Y(:,8)) + (RC(:,193) * Y(:,5)) 
        Y(:,26) = (YP(:,26) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN10O2           Y(:,27) 
        P(:) = (DJ(:,35) * Y(:,97)) + (DJ(:,65) * Y(:,163)) &
        + (DJ(:,15) * Y(:,186)) + (DJ(:,16) * Y(:,187)) &
        + (RC(:,45) * Y(:,3) * Y(:,25)) + (RC(:,380) * Y(:,3) * Y(:,97)) 
        L(:) = 0.0 &
        + (RC(:,299)) &
        + (RC(:,246) * Y(:,9)) + (RC(:,297)) + (RC(:,298)) &
        + (RC(:,109) * Y(:,8)) + (RC(:,168) * Y(:,8)) + (RC(:,192) * Y(:,5)) 
        Y(:,27) = (YP(:,27) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NC4H10           Y(:,28) 
        P(:) = 0.0
        L(:) = 0.0 &
        + (RC(:,46) * Y(:,3)) 
        Y(:,28) = (YP(:,28) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN13O2           Y(:,29) 
        P(:) = (DJ(:,95) * Y(:,184)) &
        + (DJ(:,37) * Y(:,99)) + (DJ(:,66) * Y(:,164)) &
        + (RC(:,358) * Y(:,3) * Y(:,104)) + (RC(:,382) * Y(:,3) * Y(:,99)) &
        + (RC(:,242) * Y(:,122) * Y(:,5)) + (RC(:,351) * Y(:,122)) &
        + (RC(:,46) * Y(:,3) * Y(:,28)) + (RC(:,163) * Y(:,122) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,303)) + (RC(:,304)) &
        + (RC(:,194) * Y(:,5)) + (RC(:,195) * Y(:,5)) + (RC(:,248) * Y(:,9)) &
        + (RC(:,111) * Y(:,8)) + (RC(:,112) * Y(:,8)) + (RC(:,170) * Y(:,8)) 
        Y(:,29) = (YP(:,29) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H4             Y(:,30) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,54) * Y(:,6)) &
        + (RC(:,47) * Y(:,3)) + (RC(:,50) * Y(:,5)) + (RC(:,53) * Y(:,6)) 
        Y(:,30) = (YP(:,30) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HOCH2CH2O       Y(:,31) 
        P(:) = (RC(:,466) * Y(:,3) * Y(:,192)) &
        + (RC(:,105) * Y(:,3) * Y(:,86)) + (RC(:,106) * Y(:,3) * Y(:,87)) &
        + (RC(:,103) * Y(:,3) * Y(:,84)) + (RC(:,104) * Y(:,3) * Y(:,85)) &
        + (RC(:,47) * Y(:,3) * Y(:,30)) + (RC(:,92) * Y(:,3) * Y(:,77)) 
        L(:) = 0.0 &
        + (RC(:,314)) + (RC(:,315)) &
        + (RC(:,205) * Y(:,5)) + (RC(:,206) * Y(:,5)) + (RC(:,257) * Y(:,9)) &
        + (RC(:,122) *  Y(:,8)) + (RC(:,123) * Y(:,8)) + (RC(:,173) * Y(:,8)) 
        Y(:,31) = (YP(:,31) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C3H6             Y(:,32) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,56) * Y(:,6)) &
        + (RC(:,48) * Y(:,3)) + (RC(:,51) * Y(:,5)) + (RC(:,55) * Y(:,6)) 
        Y(:,32) = (YP(:,32) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN9O2            Y(:,33) 
        P(:) = (RC(:,96) * Y(:,3) * Y(:,79)) + (RC(:,368) * Y(:,3) * Y(:,100)) &
        + (RC(:,48) * Y(:,3) * Y(:,32)) + (RC(:,94) * Y(:,78) * Y(:,3)) 
        L(:) = 0.0 &
        + (RC(:,258) * Y(:,9)) + (RC(:,316)) &
        + (RC(:,124) * Y(:,8)) + (RC(:,174) * Y(:,8)) + (RC(:,207) * Y(:,5)) 
        Y(:,33) = (YP(:,33) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TBUT2ENE         Y(:,34) 
        P(:) = 0.0
        L(:) = 0.0 &
        + (RC(:,58) * Y(:,6)) &
        + (RC(:,49) * Y(:,3)) + (RC(:,52) * Y(:,5)) + (RC(:,57) * Y(:,6)) 
        Y(:,34) = (YP(:,34) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN12O2           Y(:,35) 
        P(:) = (RC(:,369) * Y(:,3) * Y(:,189)) + (RC(:,371) * Y(:,3) * Y(:,191)) &
        + (RC(:,198) * Y(:,93) * Y(:,5)) + (RC(:,305) * Y(:,93)) &
        + (RC(:,49) * Y(:,3) * Y(:,34)) + (RC(:,115) * Y(:,93) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,259) * Y(:,9)) + (RC(:,317)) &
        + (RC(:,125) * Y(:,8)) + (RC(:,175) * Y(:,8)) + (RC(:,208) * Y(:,5)) 
        Y(:,35) = (YP(:,35) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRN6O2           Y(:,36) 
        P(:) = (RC(:,50) * Y(:,5) * Y(:,30)) 
        L(:) = 0.0 &
        + (RC(:,336)) &
        + (RC(:,144) * Y(:,8)) + (RC(:,227) * Y(:,5)) + (RC(:,274) * Y(:,9)) 
        Y(:,36) = (YP(:,36) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRN9O2           Y(:,37) 
        P(:) = (RC(:,51) * Y(:,5) * Y(:,32)) 
        L(:) = 0.0 &
        + (RC(:,337)) &
        + (RC(:,145) * Y(:,8)) + (RC(:,228) * Y(:,5)) + (RC(:,275) * Y(:,9)) 
        Y(:,37) = (YP(:,37) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRN12O2          Y(:,38) 
        P(:) = (RC(:,52) * Y(:,5) * Y(:,34)) 
        L(:) = 0.0 &
        + (RC(:,338)) &
        + (RC(:,146) * Y(:,8)) + (RC(:,229) * Y(:,5)) + (RC(:,276) * Y(:,9)) 
        Y(:,38) = (YP(:,38) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HCHO             Y(:,39) 
        P(:) = (DJ(:,93) * Y(:,182)) + (DJ(:,96) * Y(:,185)) &
        + (DJ(:,80) * Y(:,170)) + (DJ(:,91) * Y(:,180)) &
        + (DJ(:,75) * Y(:,155)) + (DJ(:,79) * Y(:,169) * 2.00) &
        + (DJ(:,74) * Y(:,154) * 2.00) &
        + (DJ(:,63) * Y(:,161)) + (DJ(:,69) * Y(:,166)) &
        + (DJ(:,41) * Y(:,123)) + (DJ(:,53) * Y(:,144)) &
        + (DJ(:,31) * Y(:,115)) + (DJ(:,32) * Y(:,115)) &
        + (DJ(:,22) * Y(:,102)) + (DJ(:,23) * Y(:,46)) &
        + (RC(:,475) * Y(:,3) * Y(:,200)) + (DJ(:,18) * Y(:,111)) &
        + (RC(:,449) * Y(:,3) * Y(:,170)) + (RC(:,473) * Y(:,3) * Y(:,198)) &
        + (RC(:,424) * Y(:,3) * Y(:,144)) &
        + (RC(:,448) * Y(:,3) * Y(:,169) * 2.00) &
        + (RC(:,392) * Y(:,3) * Y(:,123)) + (RC(:,410) * Y(:,3) * Y(:,141)) &
        + (RC(:,362) * Y(:,6) * Y(:,46)) + (RC(:,363) * Y(:,6) * Y(:,46)) &
        + (RC(:,349) * Y(:,54)) + (RC(:,352) * Y(:,55)) &
        + (RC(:,337) * Y(:,37)) + (RC(:,347) * Y(:,119)) &
        + (RC(:,336) * Y(:,36) * 2.00) &
        + (RC(:,334) * Y(:,112)) + (RC(:,335) * Y(:,112)) &
        + (RC(:,325) * Y(:,74)) + (RC(:,330) * Y(:,44)) &
        + (RC(:,316) * Y(:,33)) + (RC(:,324) * Y(:,106)) &
        + (RC(:,314) * Y(:,31) * 2.00) &
        + (RC(:,291) * Y(:,22)) + (RC(:,292) * Y(:,22)) &
        + (RC(:,240) * Y(:,54) * Y(:,5)) + (RC(:,243) * Y(:,55) * Y(:,5)) &
        + (RC(:,228) * Y(:,37) * Y(:,5)) + (RC(:,238) * Y(:,119) * Y(:,5)) &
        + (RC(:,227) * Y(:,36) * Y(:,5)) + (RC(:,227) * Y(:,36) * Y(:,5)) &
        + (RC(:,225) * Y(:,112) * Y(:,5)) + (RC(:,226) * Y(:,112) * Y(:,5)) &
        + (RC(:,216) * Y(:,74) * Y(:,5)) + (RC(:,221) * Y(:,44) * Y(:,5)) &
        + (RC(:,207) * Y(:,33) * Y(:,5)) + (RC(:,215) * Y(:,106) * Y(:,5)) &
        + (RC(:,205) * Y(:,31) * Y(:,5) * 2.00) &
        + (RC(:,162) * Y(:,56) * Y(:,8)) + (RC(:,190) * Y(:,22) * Y(:,5)) &
        + (RC(:,158) * Y(:,54) * Y(:,8)) + (RC(:,160) * Y(:,55) * Y(:,8)) &
        + (RC(:,145) * Y(:,37) * Y(:,8)) + (RC(:,156) * Y(:,119) * Y(:,8)) &
        + (RC(:,144) * Y(:,36) * Y(:,8)) + (RC(:,144) * Y(:,36) * Y(:,8)) &
        + (RC(:,142) * Y(:,112) * Y(:,8)) + (RC(:,143) * Y(:,112) * Y(:,8)) &
        + (RC(:,133) * Y(:,74) * Y(:,8)) + (RC(:,138) * Y(:,44) * Y(:,8)) &
        + (RC(:,124) * Y(:,33) * Y(:,8)) + (RC(:,132) * Y(:,106) * Y(:,8)) &
        + (RC(:,122) * Y(:,31) * Y(:,8) * 2.00) &
        + (RC(:,90) * Y(:,3) * Y(:,76)) + (RC(:,107) * Y(:,22) * Y(:,8)) &
        + (RC(:,71) * Y(:,53) * Y(:,6)) + (RC(:,72) * Y(:,53) * Y(:,6)) &
        + (RC(:,55) * Y(:,6) * Y(:,32)) + (RC(:,56) * Y(:,6) * Y(:,32)) &
        + (RC(:,53) * Y(:,6) * Y(:,30)) + (RC(:,54) * Y(:,6) * Y(:,30)) 
        L(:) = 0.0 &
        + (DJ(:,10)) &
        + (RC(:,82) * Y(:,3)) + (RC(:,85) * Y(:,5)) + (DJ(:,9)) 
        Y(:,39) = (YP(:,39) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HCOOH            Y(:,40) 
        P(:) = (RC(:,74) * Y(:,59) * Y(:,3)) &
        + (RC(:,54) * Y(:,6) * Y(:,30)) + (RC(:,62) * Y(:,6) * Y(:,43)) 
        L(:) = 0.0 &
        + (RC(:,97) * Y(:,3)) 
        Y(:,40) = (YP(:,40) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3CO2H          Y(:,41) 
        P(:) = (RC(:,56) * Y(:,6) * Y(:,32)) + (RC(:,58) * Y(:,6) * Y(:,34)) 
        L(:) = 0.0 &
        + (RC(:,98) * Y(:,3)) 
        Y(:,41) = (YP(:,41) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3CHO           Y(:,42) 
        P(:) = (DJ(:,81) * Y(:,171) * 2.00) &
        + (DJ(:,77) * Y(:,157)) + (DJ(:,80) * Y(:,170)) &
        + (DJ(:,76) * Y(:,156) * 2.00) &
        + (DJ(:,57) * Y(:,148)) + (DJ(:,75) * Y(:,155)) &
        + (DJ(:,45) * Y(:,127)) + (DJ(:,54) * Y(:,145)) &
        + (DJ(:,20) * Y(:,104)) + (DJ(:,42) * Y(:,124)) &
        + (RC(:,474) * Y(:,3) * Y(:,199)) + (DJ(:,19) * Y(:,188)) &
        + (RC(:,449) * Y(:,3) * Y(:,170)) &
        + (RC(:,450) * Y(:,3) * Y(:,171) * 2.00) &
        + (RC(:,393) * Y(:,3) * Y(:,124)) + (RC(:,425) * Y(:,3) * Y(:,145)) &
        + (RC(:,338) * Y(:,38) * 2.00) &
        + (RC(:,327) * Y(:,107)) + (RC(:,337) * Y(:,37)) &
        + (RC(:,318) * Y(:,95)) + (RC(:,326) * Y(:,75)) &
        + (RC(:,317) * Y(:,35) * 2.00) &
        + (RC(:,303) * Y(:,29)) + (RC(:,316) * Y(:,33)) &
        + (RC(:,294) * Y(:,24)) + (RC(:,295) * Y(:,24)) &
        + (RC(:,229) * Y(:,38) * Y(:,5) * 2.00) &
        + (RC(:,218) * Y(:,107) * Y(:,5)) + (RC(:,228) * Y(:,37) * Y(:,5)) &
        + (RC(:,209) * Y(:,95) * Y(:,5)) + (RC(:,217) * Y(:,75) * Y(:,5)) &
        + (RC(:,207) * Y(:,33) * Y(:,5)) + (RC(:,208) * Y(:,35) * Y(:,5) * 2.00) &
        + (RC(:,191) * Y(:,24) * Y(:,5)) + (RC(:,194) * Y(:,29) * Y(:,5)) &
        + (RC(:,146) * Y(:,38) * Y(:,8) * 2.00) &
        + (RC(:,135) * Y(:,107) * Y(:,8)) + (RC(:,145) * Y(:,37) * Y(:,8)) &
        + (RC(:,126) * Y(:,95) * Y(:,8)) + (RC(:,134) * Y(:,75) * Y(:,8)) &
        + (RC(:,125) * Y(:,35) * Y(:,8) * 2.00) &
        + (RC(:,111) * Y(:,29) * Y(:,8)) + (RC(:,124) * Y(:,33) * Y(:,8)) &
        + (RC(:,91) * Y(:,3) * Y(:,77)) + (RC(:,108) * Y(:,24) * Y(:,8)) &
        + (RC(:,57) * Y(:,6) * Y(:,34)) + (RC(:,58) * Y(:,6) * Y(:,34)) 
        L(:) = 0.0 &
        + (RC(:,83) * Y(:,3)) + (RC(:,86) * Y(:,5)) + (DJ(:,11)) 
        Y(:,42) = (YP(:,42) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C5H8             Y(:,43) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,62) * Y(:,6)) &
        + (RC(:,59) * Y(:,3)) + (RC(:,60) * Y(:,5)) + (RC(:,61) * Y(:,6)) 
        Y(:,43) = (YP(:,43) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RU14O2           Y(:,44) 
        P(:) = (RC(:,59) * Y(:,3) * Y(:,43)) 
        L(:) = 0.0 &
        + (RC(:,329)) + (RC(:,330)) &
        + (RC(:,220) * Y(:,5)) + (RC(:,221) * Y(:,5)) + (RC(:,271) * Y(:,9)) &
        + (RC(:,137) * Y(:,8)) + (RC(:,138) * Y(:,8)) + (RC(:,180) * Y(:,8)) 
        Y(:,44) = (YP(:,44) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRU14O2          Y(:,45) 
        P(:) = (RC(:,60) * Y(:,5) * Y(:,43)) 
        L(:) = 0.0 &
        + (RC(:,339)) &
        + (RC(:,147) * Y(:,8)) + (RC(:,230) * Y(:,5)) + (RC(:,277) * Y(:,9)) 
        Y(:,45) = (YP(:,45) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          UCARB10          Y(:,46) 
        P(:) = (DJ(:,69) * Y(:,166)) &
        + (RC(:,330) * Y(:,44)) + (RC(:,481) * Y(:,3) * Y(:,201)) &
        + (RC(:,138) * Y(:,44) * Y(:,8)) + (RC(:,221) * Y(:,44) * Y(:,5)) &
        + (RC(:,61) * Y(:,6) * Y(:,43)) + (RC(:,62) * Y(:,6) * Y(:,43)) 
        L(:) = 0.0 &
        + (RC(:,363) * Y(:,6)) + (DJ(:,23)) &
        + (RC(:,360) * Y(:,3)) + (RC(:,361) * Y(:,5)) + (RC(:,362) * Y(:,6)) 
        Y(:,46) = (YP(:,46) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          APINENE          Y(:,47) 
        P(:) = 0.0
        L(:) = 0.0 &
        + (RC(:,66) * Y(:,6)) + (RC(:,67) * Y(:,6)) &
        + (RC(:,63) * Y(:,3)) + (RC(:,64) * Y(:,5)) + (RC(:,65) * Y(:,6)) 
        Y(:,47) = (YP(:,47) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN28O2          Y(:,48) 
        P(:) = (RC(:,63) * Y(:,47) * Y(:,3)) 
        L(:) = 0.0 &
        + (RC(:,232) * Y(:,5)) + (RC(:,279) * Y(:,9)) + (RC(:,341)) &
        + (RC(:,149) * Y(:,8)) + (RC(:,150) * Y(:,8)) + (RC(:,185) * Y(:,8)) 
        Y(:,48) = (YP(:,48) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRTN28O2         Y(:,49) 
        P(:) = (RC(:,64) * Y(:,47) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,342)) &
        + (RC(:,151) * Y(:,8)) + (RC(:,233) * Y(:,5)) + (RC(:,280) * Y(:,9)) 
        Y(:,49) = (YP(:,49) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN26O2          Y(:,50) 
        P(:) = (RC(:,483) * Y(:,203)) + (DJ(:,39) * Y(:,51)) &
        + (RC(:,387) * Y(:,5) * Y(:,51)) + (RC(:,455) * Y(:,3) * Y(:,176)) &
        + (RC(:,65) * Y(:,47) * Y(:,6)) + (RC(:,384) * Y(:,3) * Y(:,51)) 
        L(:) = 0.0 &
        + (RC(:,343)) + (RC(:,482) * Y(:,4)) &
        + (RC(:,152) * Y(:,8)) + (RC(:,234) * Y(:,5)) + (RC(:,281) * Y(:,9)) 
        Y(:,50) = (YP(:,50) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TNCARB26         Y(:,51) 
        P(:) = (DJ(:,85) * Y(:,174)) + (DJ(:,86) * Y(:,175)) &
        + (RC(:,454) * Y(:,3) * Y(:,174)) + (RC(:,456) * Y(:,3) * Y(:,175)) &
        + (RC(:,342) * Y(:,49)) + (RC(:,408) * Y(:,3) * Y(:,139)) &
        + (RC(:,233) * Y(:,49) * Y(:,5)) + (RC(:,341) * Y(:,48)) &
        + (RC(:,151) * Y(:,49) * Y(:,8)) + (RC(:,232) * Y(:,48) * Y(:,5)) &
        + (RC(:,66) * Y(:,47) * Y(:,6)) + (RC(:,149) * Y(:,48) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,384) * Y(:,3)) + (RC(:,387) * Y(:,5)) + (DJ(:,39)) 
        Y(:,51) = (YP(:,51) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RCOOH25          Y(:,52) 
        P(:) = (RC(:,67) * Y(:,47) * Y(:,6)) + (RC(:,490) * Y(:,206)) 
        L(:) = 0.0 &
        + (RC(:,389) * Y(:,3)) + (RC(:,489)) 
        Y(:,52) = (YP(:,52) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          BPINENE          Y(:,53) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,71) * Y(:,6)) + (RC(:,72) * Y(:,6)) + (RC(:,73) * Y(:,6)) &
        + (RC(:,68) * Y(:,3)) + (RC(:,69) * Y(:,5)) + (RC(:,70) * Y(:,6)) 
        Y(:,53) = (YP(:,53) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTX28O2          Y(:,54) 
        P(:) = (RC(:,68) * Y(:,53) * Y(:,3)) + (RC(:,462) * Y(:,3) * Y(:,182)) 
        L(:) = 0.0 &
        + (RC(:,240) * Y(:,5)) + (RC(:,287) * Y(:,9)) + (RC(:,349)) &
        + (RC(:,158) * Y(:,8)) + (RC(:,159) * Y(:,8)) + (RC(:,187) * Y(:,8)) 
        Y(:,54) = (YP(:,54) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRTX28O2         Y(:,55) 
        P(:) = (RC(:,69) * Y(:,53) * Y(:,5)) + (RC(:,465) * Y(:,3) * Y(:,185)) 
        L(:) = 0.0 &
        + (RC(:,352)) &
        + (RC(:,160) * Y(:,8)) + (RC(:,243) * Y(:,5)) + (RC(:,290) * Y(:,9)) 
        Y(:,55) = (YP(:,55) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTX24O2          Y(:,56) 
        P(:) = (RC(:,70) * Y(:,53) * Y(:,6)) + (RC(:,390) * Y(:,3) * Y(:,57)) 
        L(:) = 0.0 &
        + (RC(:,241) * Y(:,5)) + (RC(:,288) * Y(:,9)) + (RC(:,350)) &
        + (RC(:,161) * Y(:,8)) + (RC(:,162) * Y(:,8)) + (RC(:,188) * Y(:,8)) 
        Y(:,56) = (YP(:,56) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TXCARB24         Y(:,57) 
        P(:) = (DJ(:,96) * Y(:,185)) &
        + (RC(:,410) * Y(:,3) * Y(:,141)) + (DJ(:,93) * Y(:,182)) &
        + (RC(:,349) * Y(:,54)) + (RC(:,352) * Y(:,55)) &
        + (RC(:,240) * Y(:,54) * Y(:,5)) + (RC(:,243) * Y(:,55) * Y(:,5)) &
        + (RC(:,158) * Y(:,54) * Y(:,8)) + (RC(:,160) * Y(:,55) * Y(:,8)) &
        + (RC(:,71) * Y(:,53) * Y(:,6)) + (RC(:,73) * Y(:,53) * Y(:,6)) 
        L(:) = 0.0 &
        + (RC(:,390) * Y(:,3)) 
        Y(:,57) = (YP(:,57) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TXCARB22         Y(:,58) 
        P(:) = (DJ(:,52) * Y(:,142)) + (DJ(:,94) * Y(:,183)) &
        + (RC(:,411) * Y(:,3) * Y(:,142)) + (RC(:,463) * Y(:,3) * Y(:,183)) &
        + (RC(:,241) * Y(:,56) * Y(:,5)) + (RC(:,350) * Y(:,56)) &
        + (RC(:,72) * Y(:,53) * Y(:,6)) + (RC(:,161) * Y(:,56) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,391) * Y(:,3)) 
        Y(:,58) = (YP(:,58) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H2             Y(:,59) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,74) * Y(:,3)) + (RC(:,75) * Y(:,3)) 
        Y(:,59) = (YP(:,59) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB3            Y(:,60) 
        P(:) = (DJ(:,83) * Y(:,152)) &
        + (DJ(:,50) * Y(:,137)) + (DJ(:,82) * Y(:,151)) &
        + (RC(:,452) * Y(:,3) * Y(:,152)) + (DJ(:,49) * Y(:,136)) &
        + (RC(:,406) * Y(:,3) * Y(:,137)) + (RC(:,451) * Y(:,3) * Y(:,151)) &
        + (RC(:,311) * Y(:,69)) + (RC(:,405) * Y(:,3) * Y(:,136)) &
        + (RC(:,308) * Y(:,65)) + (RC(:,310) * Y(:,68)) &
        + (RC(:,203) * Y(:,68) * Y(:,5)) + (RC(:,307) * Y(:,62)) &
        + (RC(:,200) * Y(:,62) * Y(:,5)) + (RC(:,201) * Y(:,65) * Y(:,5)) &
        + (RC(:,118) * Y(:,65) * Y(:,8)) + (RC(:,120) * Y(:,68) * Y(:,8)) &
        + (RC(:,75) * Y(:,59) * Y(:,3)) + (RC(:,117) * Y(:,62) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,366) * Y(:,3)) + (DJ(:,24)) 
        Y(:,60) = (YP(:,60) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          BENZENE          Y(:,61) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,76) * Y(:,3)) + (RC(:,77) * Y(:,3)) 
        Y(:,61) = (YP(:,61) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA13O2           Y(:,62) 
        P(:) = (RC(:,76) * Y(:,61) * Y(:,3)) 
        L(:) = 0.0 &
        + (RC(:,253) * Y(:,9)) + (RC(:,307)) &
        + (RC(:,117) * Y(:,8)) + (RC(:,181) * Y(:,8)) + (RC(:,200) * Y(:,5)) 
        Y(:,62) = (YP(:,62) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          AROH14           Y(:,63) 
        P(:) = (RC(:,77) * Y(:,61) * Y(:,3)) 
        L(:) = 0.0 &
        + (RC(:,413) * Y(:,3)) + (RC(:,414) * Y(:,5)) 
        Y(:,63) = (YP(:,63) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TOLUENE          Y(:,64) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,78) * Y(:,3)) + (RC(:,79) * Y(:,3)) 
        Y(:,64) = (YP(:,64) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA16O2           Y(:,65) 
        P(:) = (RC(:,78) * Y(:,64) * Y(:,3)) 
        L(:) = 0.0 &
        + (RC(:,308)) + (RC(:,309)) &
        + (RC(:,201) * Y(:,5)) + (RC(:,202) * Y(:,5)) + (RC(:,254) * Y(:,9)) &
        + (RC(:,118) * Y(:,8)) + (RC(:,119) * Y(:,8)) + (RC(:,182) * Y(:,8)) 
        Y(:,65) = (YP(:,65) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          AROH17           Y(:,66) 
        P(:) = (RC(:,79) * Y(:,64) * Y(:,3)) 
        L(:) = 0.0 &
        + (RC(:,418) * Y(:,3)) + (RC(:,419) * Y(:,5)) 
        Y(:,66) = (YP(:,66) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          OXYL             Y(:,67) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,80) * Y(:,3)) + (RC(:,81) * Y(:,3)) 
        Y(:,67) = (YP(:,67) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA19AO2          Y(:,68) 
        P(:) = (RC(:,80) * Y(:,67) * Y(:,3)) 
        L(:) = 0.0 &
        + (RC(:,255) * Y(:,9)) + (RC(:,310)) &
        + (RC(:,120) * Y(:,8)) + (RC(:,183) * Y(:,8)) + (RC(:,203) * Y(:,5)) 
        Y(:,68) = (YP(:,68) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA19CO2          Y(:,69) 
        P(:) = (RC(:,81) * Y(:,67) * Y(:,3)) 
        L(:) = 0.0 &
        + (RC(:,256) * Y(:,9)) + (RC(:,311)) &
        + (RC(:,121) * Y(:,8)) + (RC(:,184) * Y(:,8)) + (RC(:,204) * Y(:,5)) 
        Y(:,69) = (YP(:,69) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3CO3           Y(:,70) 
        P(:) = (DJ(:,71) * Y(:,168)) &
        + (DJ(:,40) * Y(:,120) * 2.00) &
        + (DJ(:,31) * Y(:,115)) + (DJ(:,32) * Y(:,115)) &
        + (DJ(:,27) * Y(:,189)) + (DJ(:,29) * Y(:,109)) &
        + (DJ(:,25) * Y(:,98)) + (DJ(:,26) * Y(:,100) * 2.00) &
        + (DJ(:,19) * Y(:,188)) + (DJ(:,23) * Y(:,46)) &
        + (DJ(:,17) * Y(:,88)) + (DJ(:,18) * Y(:,111)) &
        + (DJ(:,14) * Y(:,101)) + (DJ(:,15) * Y(:,186)) &
        + (RC(:,468) * Y(:,198)) + (DJ(:,13) * Y(:,73)) &
        + (RC(:,374) * Y(:,6) * Y(:,109)) + (RC(:,431) * Y(:,3) * Y(:,159)) &
        + (RC(:,362) * Y(:,6) * Y(:,46)) + (RC(:,367) * Y(:,3) * Y(:,98)) &
        + (RC(:,331) * Y(:,110)) + (RC(:,333) * Y(:,112)) &
        + (RC(:,325) * Y(:,74)) + (RC(:,326) * Y(:,75)) &
        + (RC(:,222) * Y(:,110) * Y(:,5)) + (RC(:,224) * Y(:,112) * Y(:,5)) &
        + (RC(:,216) * Y(:,74) * Y(:,5)) + (RC(:,217) * Y(:,75) * Y(:,5)) &
        + (RC(:,139) * Y(:,110) * Y(:,8)) + (RC(:,141) * Y(:,112) * Y(:,8)) &
        + (RC(:,133) * Y(:,74) * Y(:,8)) + (RC(:,134) * Y(:,75) * Y(:,8)) &
        + (RC(:,83) * Y(:,3) * Y(:,42)) + (RC(:,86) * Y(:,5) * Y(:,42)) 
        L(:) = 0.0 &
        + (RC(:,322)) + (RC(:,467) * Y(:,4)) &
        + (RC(:,130) * Y(:,8)) + (RC(:,213) * Y(:,5)) + (RC(:,264) * Y(:,9)) 
        Y(:,70) = (YP(:,70) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H5CHO          Y(:,71) 
        P(:) = (DJ(:,78) * Y(:,158) * 2.00) &
        + (DJ(:,55) * Y(:,146)) + (DJ(:,77) * Y(:,157)) &
        + (DJ(:,21) * Y(:,105)) + (DJ(:,43) * Y(:,125)) &
        + (RC(:,394) * Y(:,3) * Y(:,125)) + (RC(:,426) * Y(:,3) * Y(:,146)) &
        + (RC(:,318) * Y(:,95)) + (RC(:,319) * Y(:,103) * 2.00) &
        + (RC(:,297) * Y(:,27)) + (RC(:,298) * Y(:,27)) &
        + (RC(:,210) * Y(:,103) * Y(:,5) * 2.00) &
        + (RC(:,192) * Y(:,27) * Y(:,5)) + (RC(:,209) * Y(:,95) * Y(:,5)) &
        + (RC(:,126) * Y(:,95) * Y(:,8)) &
        + (RC(:,127) * Y(:,103) * Y(:,8) * 2.00) &
        + (RC(:,93) * Y(:,78) * Y(:,3)) + (RC(:,109) * Y(:,27) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,84) * Y(:,3)) + (RC(:,87) * Y(:,5)) + (DJ(:,12)) 
        Y(:,71) = (YP(:,71) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H5CO3          Y(:,72) 
        P(:) = (RC(:,470) * Y(:,199)) &
        + (RC(:,327) * Y(:,107)) + (RC(:,432) * Y(:,3) * Y(:,160)) &
        + (RC(:,135) * Y(:,107) * Y(:,8)) + (RC(:,218) * Y(:,107) * Y(:,5)) &
        + (RC(:,84) * Y(:,3) * Y(:,71)) + (RC(:,87) * Y(:,5) * Y(:,71)) 
        L(:) = 0.0 &
        + (RC(:,323)) + (RC(:,469) * Y(:,4)) &
        + (RC(:,131) * Y(:,8)) + (RC(:,214) * Y(:,5)) + (RC(:,265) * Y(:,9)) 
        Y(:,72) = (YP(:,72) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3COCH3         Y(:,73) 
        P(:) = (DJ(:,90) * Y(:,179)) + (DJ(:,95) * Y(:,184)) &
        + (DJ(:,44) * Y(:,126)) + (DJ(:,56) * Y(:,147)) &
        + (RC(:,464) * Y(:,3) * Y(:,184)) + (RC(:,484) * Y(:,3) * Y(:,203)) &
        + (RC(:,412) * Y(:,3) * Y(:,143)) + (RC(:,427) * Y(:,3) * Y(:,147)) &
        + (RC(:,395) * Y(:,3) * Y(:,126)) + (RC(:,409) * Y(:,3) * Y(:,140)) &
        + (RC(:,346) * Y(:,118)) + (RC(:,351) * Y(:,122)) &
        + (RC(:,300) * Y(:,26)) + (RC(:,301) * Y(:,26)) &
        + (RC(:,237) * Y(:,118) * Y(:,5)) + (RC(:,242) * Y(:,122) * Y(:,5)) &
        + (RC(:,163) * Y(:,122) * Y(:,8)) + (RC(:,193) * Y(:,26) * Y(:,5)) &
        + (RC(:,159) * Y(:,54) * Y(:,8)) + (RC(:,162) * Y(:,56) * Y(:,8)) &
        + (RC(:,150) * Y(:,48) * Y(:,8)) + (RC(:,155) * Y(:,118) * Y(:,8)) &
        + (RC(:,95) * Y(:,3) * Y(:,79)) + (RC(:,110) * Y(:,26) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,88) * Y(:,3)) + (DJ(:,13)) 
        Y(:,73) = (YP(:,73) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN8O2            Y(:,74) 
        P(:) = (DJ(:,92) * Y(:,181)) &
        + (DJ(:,28) * Y(:,190) * 2.00) &
        + (DJ(:,21) * Y(:,105)) + (DJ(:,27) * Y(:,189)) &
        + (DJ(:,16) * Y(:,187)) + (DJ(:,20) * Y(:,104)) &
        + (RC(:,239) * Y(:,121) * Y(:,5)) + (RC(:,348) * Y(:,121)) &
        + (RC(:,88) * Y(:,3) * Y(:,73)) + (RC(:,157) * Y(:,121) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,325)) &
        + (RC(:,133) * Y(:,8)) + (RC(:,216) * Y(:,5)) + (RC(:,267) * Y(:,9)) 
        Y(:,74) = (YP(:,74) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN11O2           Y(:,75) 
        P(:) = (RC(:,89) * Y(:,101) * Y(:,3)) + (RC(:,355) * Y(:,3) * Y(:,88)) 
        L(:) = 0.0 &
        + (RC(:,326)) &
        + (RC(:,134) * Y(:,8)) + (RC(:,217) * Y(:,5)) + (RC(:,268) * Y(:,9)) 
        Y(:,75) = (YP(:,75) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3OH            Y(:,76) 
        P(:) = (RC(:,293) * Y(:,22)) 
        L(:) = 0.0 &
        + (RC(:,90) * Y(:,3)) 
        Y(:,76) = (YP(:,76) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H5OH           Y(:,77) 
        P(:) = (RC(:,296) * Y(:,24)) 
        L(:) = 0.0 &
        + (RC(:,91) * Y(:,3)) + (RC(:,92) * Y(:,3)) 
        Y(:,77) = (YP(:,77) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NPROPOL          Y(:,78) 
        P(:) = (RC(:,299) * Y(:,27)) 
        L(:) = 0.0 &
        + (RC(:,93) * Y(:,3)) + (RC(:,94) * Y(:,3)) 
        Y(:,78) = (YP(:,78) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          IPROPOL          Y(:,79) 
        P(:) = (RC(:,302) * Y(:,26)) 
        L(:) = 0.0 &
        + (RC(:,95) * Y(:,3)) + (RC(:,96) * Y(:,3)) 
        Y(:,79) = (YP(:,79) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3CL            Y(:,80) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,99) * Y(:,3)) 
        Y(:,80) = (YP(:,80) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH2CL2           Y(:,81) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,100) * Y(:,3)) 
        Y(:,81) = (YP(:,81) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CHCL3            Y(:,82) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,101) * Y(:,3)) 
        Y(:,82) = (YP(:,82) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3CCL3          Y(:,83) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,102) * Y(:,3)) 
        Y(:,83) = (YP(:,83) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TCE              Y(:,84) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,103) * Y(:,3)) 
        Y(:,84) = (YP(:,84) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TRICLETH         Y(:,85) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,104) * Y(:,3)) 
        Y(:,85) = (YP(:,85) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CDICLETH         Y(:,86) 
        P(:) = 0.0
        L(:) = 0.0 &
        + (RC(:,105) * Y(:,3)) 
        Y(:,86) = (YP(:,86) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TDICLETH         Y(:,87) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,106) * Y(:,3)) 
        Y(:,87) = (YP(:,87) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB11A          Y(:,88) 
        P(:) = (DJ(:,58) * Y(:,148)) &
        + (RC(:,428) * Y(:,3) * Y(:,148)) + (DJ(:,46) * Y(:,127)) &
        + (RC(:,304) * Y(:,29)) + (RC(:,396) * Y(:,3) * Y(:,127)) &
        + (RC(:,112) * Y(:,29) * Y(:,8)) + (RC(:,195) * Y(:,29) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,355) * Y(:,3)) + (DJ(:,17)) 
        Y(:,88) = (YP(:,88) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN16O2           Y(:,89) 
        P(:) = (RC(:,359) * Y(:,3) * Y(:,105)) + (DJ(:,67) * Y(:,165)) 
        L(:) = 0.0 &
        + (RC(:,249) * Y(:,9)) + (RC(:,312)) &
        + (RC(:,113) * Y(:,8)) + (RC(:,171) * Y(:,8)) + (RC(:,196) * Y(:,5)) 
        Y(:,89) = (YP(:,89) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN15AO2          Y(:,90) 
        P(:) = (DJ(:,59) * Y(:,149)) &
        + (RC(:,312) * Y(:,89)) + (RC(:,385) * Y(:,3) * Y(:,193)) &
        + (RC(:,113) * Y(:,89) * Y(:,8)) + (RC(:,196) * Y(:,89) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,262) * Y(:,9)) + (RC(:,320)) &
        + (RC(:,128) * Y(:,8)) + (RC(:,178) * Y(:,8)) + (RC(:,211) * Y(:,5)) 
        Y(:,90) = (YP(:,90) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN19O2           Y(:,91) 
        P(:) = (RC(:,150) * Y(:,48) * Y(:,8)) + (RC(:,159) * Y(:,54) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,250) * Y(:,9)) + (RC(:,313)) &
        + (RC(:,114) * Y(:,8)) + (RC(:,172) * Y(:,8)) + (RC(:,197) * Y(:,5)) 
        Y(:,91) = (YP(:,91) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN18AO2          Y(:,92) 
        P(:) = (RC(:,313) * Y(:,91)) + (DJ(:,60) * Y(:,150)) &
        + (RC(:,114) * Y(:,91) * Y(:,8)) + (RC(:,197) * Y(:,91) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,263) * Y(:,9)) + (RC(:,321)) &
        + (RC(:,129) * Y(:,8)) + (RC(:,179) * Y(:,8)) + (RC(:,212) * Y(:,5)) 
        Y(:,92) = (YP(:,92) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN13AO2          Y(:,93) 
        P(:) = (RC(:,162) * Y(:,56) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,305)) &
        + (RC(:,115) * Y(:,8)) + (RC(:,198) * Y(:,5)) + (RC(:,251) * Y(:,9)) 
        Y(:,93) = (YP(:,93) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN16AO2          Y(:,94) 
        P(:) = (RC(:,328) * Y(:,108)) &
        + (RC(:,136) * Y(:,108) * Y(:,8)) + (RC(:,219) * Y(:,108) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,306)) &
        + (RC(:,116) * Y(:,8)) + (RC(:,199) * Y(:,5)) + (RC(:,252) * Y(:,9)) 
        Y(:,94) = (YP(:,94) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN15O2           Y(:,95) 
        P(:) = (DJ(:,47) * Y(:,128)) &
        + (RC(:,306) * Y(:,94)) + (RC(:,370) * Y(:,3) * Y(:,190)) &
        + (RC(:,116) * Y(:,94) * Y(:,8)) + (RC(:,199) * Y(:,94) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,260) * Y(:,9)) + (RC(:,318)) &
        + (RC(:,126) * Y(:,8)) + (RC(:,176) * Y(:,8)) + (RC(:,209) * Y(:,5)) 
        Y(:,95) = (YP(:,95) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          UDCARB8          Y(:,96) 
        P(:) = (DJ(:,49) * Y(:,136)) + (DJ(:,82) * Y(:,151)) &
        + (RC(:,405) * Y(:,3) * Y(:,136)) + (RC(:,451) * Y(:,3) * Y(:,151)) &
        + (RC(:,307) * Y(:,62)) + (RC(:,309) * Y(:,65)) &
        + (RC(:,202) * Y(:,65) * Y(:,5)) + (RC(:,204) * Y(:,69) * Y(:,5)) &
        + (RC(:,121) * Y(:,69) * Y(:,8)) + (RC(:,200) * Y(:,62) * Y(:,5)) &
        + (RC(:,117) * Y(:,62) * Y(:,8)) + (RC(:,119) * Y(:,65) * Y(:,8)) 
        L(:) = 0.0 &
        + (DJ(:,34)) &
        + (RC(:,378) * Y(:,3)) + (RC(:,379) * Y(:,3)) + (DJ(:,33)) 
        Y(:,96) = (YP(:,96) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          UDCARB11         Y(:,97) 
        P(:) = (DJ(:,84) * Y(:,153)) &
        + (DJ(:,51) * Y(:,138)) + (DJ(:,83) * Y(:,152)) &
        + (RC(:,453) * Y(:,3) * Y(:,153)) + (DJ(:,50) * Y(:,137)) &
        + (RC(:,407) * Y(:,3) * Y(:,138)) + (RC(:,452) * Y(:,3) * Y(:,152)) &
        + (RC(:,308) * Y(:,65)) + (RC(:,406) * Y(:,3) * Y(:,137)) &
        + (RC(:,118) * Y(:,65) * Y(:,8)) + (RC(:,201) * Y(:,65) * Y(:,5)) 
        L(:) = 0.0 &
        + (DJ(:,36)) &
        + (RC(:,380) * Y(:,3)) + (RC(:,381) * Y(:,3)) + (DJ(:,35)) 
        Y(:,97) = (YP(:,97) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB6            Y(:,98) 
        P(:) = (DJ(:,70) * Y(:,167)) + (DJ(:,84) * Y(:,153)) &
        + (RC(:,453) * Y(:,3) * Y(:,153)) + (DJ(:,51) * Y(:,138)) &
        + (RC(:,407) * Y(:,3) * Y(:,138)) + (RC(:,434) * Y(:,3) * Y(:,162)) &
        + (RC(:,375) * Y(:,6) * Y(:,109)) + (RC(:,377) * Y(:,3) * Y(:,115)) &
        + (RC(:,356) * Y(:,3) * Y(:,111)) + (RC(:,363) * Y(:,6) * Y(:,46)) &
        + (RC(:,309) * Y(:,65)) + (RC(:,334) * Y(:,112)) &
        + (RC(:,202) * Y(:,65) * Y(:,5)) + (RC(:,225) * Y(:,112) * Y(:,5)) &
        + (RC(:,119) * Y(:,65) * Y(:,8)) + (RC(:,142) * Y(:,112) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,367) * Y(:,3)) + (DJ(:,25)) 
        Y(:,98) = (YP(:,98) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          UDCARB14         Y(:,99) 
        P(:) = (RC(:,310) * Y(:,68)) + (RC(:,311) * Y(:,69)) &
        + (RC(:,120) * Y(:,68) * Y(:,8)) + (RC(:,203) * Y(:,68) * Y(:,5)) 
        L(:) = 0.0 &
        + (DJ(:,38)) &
        + (RC(:,382) * Y(:,3)) + (RC(:,383) * Y(:,3)) + (DJ(:,37)) 
        Y(:,99) = (YP(:,99) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB9            Y(:,100) 
        P(:) = (RC(:,357) * Y(:,3) * Y(:,188)) + (RC(:,435) * Y(:,3) * Y(:,163)) &
        + (RC(:,121) * Y(:,69) * Y(:,8)) + (RC(:,204) * Y(:,69) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,368) * Y(:,3)) + (DJ(:,26)) 
        Y(:,100) = (YP(:,100) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          MEK              Y(:,101) 
        P(:) = 0.0 
        L(:) = 0.0 &
        + (RC(:,89) * Y(:,3)) + (DJ(:,14)) 
        Y(:,101) = (YP(:,101) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HOCH2CHO         Y(:,102) 
        P(:) = (DJ(:,71) * Y(:,168)) &
        + (DJ(:,29) * Y(:,109)) + (DJ(:,70) * Y(:,167)) &
        + (RC(:,399) * Y(:,3) * Y(:,130)) + (RC(:,443) * Y(:,3) * Y(:,154)) &
        + (RC(:,374) * Y(:,6) * Y(:,109)) + (RC(:,375) * Y(:,6) * Y(:,109)) &
        + (RC(:,332) * Y(:,110)) + (RC(:,333) * Y(:,112)) &
        + (RC(:,315) * Y(:,31)) + (RC(:,331) * Y(:,110)) &
        + (RC(:,222) * Y(:,110) * Y(:,5)) + (RC(:,224) * Y(:,112) * Y(:,5)) &
        + (RC(:,141) * Y(:,112) * Y(:,8)) + (RC(:,206) * Y(:,31) * Y(:,5)) &
        + (RC(:,123) * Y(:,31) * Y(:,8)) + (RC(:,139) * Y(:,110) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,364) * Y(:,3)) + (RC(:,365) * Y(:,5)) + (DJ(:,22)) 
        Y(:,102) = (YP(:,102) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN18O2           Y(:,103) 
        P(:) = (DJ(:,48) * Y(:,129)) 
        L(:) = 0.0 &
        + (RC(:,261) * Y(:,9)) + (RC(:,319)) &
        + (RC(:,127) * Y(:,8)) + (RC(:,177) * Y(:,8)) + (RC(:,210) * Y(:,5)) 
        Y(:,103) = (YP(:,103) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB13           Y(:,104) 
        P(:) = (RC(:,446) * Y(:,3) * Y(:,157)) &
        + (RC(:,416) * Y(:,3) * Y(:,195)) + (RC(:,417) * Y(:,5) * Y(:,195)) &
        + (RC(:,320) * Y(:,90)) + (RC(:,402) * Y(:,3) * Y(:,133)) &
        + (RC(:,128) * Y(:,90) * Y(:,8)) + (RC(:,211) * Y(:,90) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,358) * Y(:,3)) + (DJ(:,20)) 
        Y(:,104) = (YP(:,104) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB16           Y(:,105) 
        P(:) = (RC(:,447) * Y(:,3) * Y(:,158)) + (RC(:,484) * Y(:,3) * Y(:,203)) &
        + (RC(:,421) * Y(:,3) * Y(:,197)) + (RC(:,422) * Y(:,5) * Y(:,197)) &
        + (RC(:,321) * Y(:,92)) + (RC(:,403) * Y(:,3) * Y(:,134)) &
        + (RC(:,129) * Y(:,92) * Y(:,8)) + (RC(:,212) * Y(:,92) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,359) * Y(:,3)) + (DJ(:,21)) 
        Y(:,105) = (YP(:,105) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HOCH2CO3         Y(:,106) 
        P(:) = (RC(:,433) * Y(:,3) * Y(:,161)) + (RC(:,472) * Y(:,200)) &
        + (RC(:,364) * Y(:,3) * Y(:,102)) + (RC(:,365) * Y(:,5) * Y(:,102)) 
        L(:) = 0.0 &
        + (RC(:,324)) + (RC(:,471) * Y(:,4)) &
        + (RC(:,132) * Y(:,8)) + (RC(:,215) * Y(:,5)) + (RC(:,266) * Y(:,9)) 
        Y(:,106) = (YP(:,106) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN14O2           Y(:,107) 
        P(:) = (RC(:,353) * Y(:,3) * Y(:,186)) 
        L(:) = 0.0 &
        + (RC(:,327)) &
        + (RC(:,135) * Y(:,8)) + (RC(:,218) * Y(:,5)) + (RC(:,269) * Y(:,9)) 
        Y(:,107) = (YP(:,107) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN17O2           Y(:,108) 
        P(:) = (RC(:,354) * Y(:,3) * Y(:,187)) 
        L(:) = 0.0 &
        + (RC(:,328)) &
        + (RC(:,136) * Y(:,8)) + (RC(:,219) * Y(:,5)) + (RC(:,270) * Y(:,9)) 
        Y(:,108) = (YP(:,108) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          UCARB12          Y(:,109) 
        P(:) = (RC(:,438) * Y(:,3) * Y(:,166)) + (DJ(:,68) * Y(:,166)) &
        + (RC(:,329) * Y(:,44)) + (RC(:,404) * Y(:,3) * Y(:,135)) &
        + (RC(:,137) * Y(:,44) * Y(:,8)) + (RC(:,220) * Y(:,44) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,375) * Y(:,6)) + (DJ(:,29)) &
        + (RC(:,372) * Y(:,3)) + (RC(:,373) * Y(:,5)) + (RC(:,374) * Y(:,6)) 
        Y(:,109) = (YP(:,109) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RU12O2           Y(:,110) 
        P(:) = (RC(:,439) * Y(:,3) * Y(:,167)) + (RC(:,477) * Y(:,201)) &
        + (RC(:,372) * Y(:,3) * Y(:,109)) + (RC(:,373) * Y(:,5) * Y(:,109)) 
        L(:) = 0.0 &
        + (RC(:,332)) + (RC(:,476) * Y(:,4)) &
        + (RC(:,223) * Y(:,5)) + (RC(:,272) * Y(:,9)) + (RC(:,331)) &
        + (RC(:,139) * Y(:,8)) + (RC(:,140) * Y(:,8)) + (RC(:,222) * Y(:,5)) 
        Y(:,110) = (YP(:,110) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB7            Y(:,111) 
        P(:) = (RC(:,480) * Y(:,3) * Y(:,202)) &
        + (RC(:,400) * Y(:,3) * Y(:,131)) + (RC(:,444) * Y(:,3) * Y(:,155)) &
        + (RC(:,332) * Y(:,110)) + (RC(:,335) * Y(:,112)) &
        + (RC(:,223) * Y(:,110) * Y(:,5)) + (RC(:,226) * Y(:,112) * Y(:,5)) &
        + (RC(:,140) * Y(:,110) * Y(:,8)) + (RC(:,143) * Y(:,112) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,356) * Y(:,3)) + (DJ(:,18)) 
        Y(:,111) = (YP(:,111) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RU10O2           Y(:,112) 
        P(:) = (RC(:,440) * Y(:,3) * Y(:,168)) + (RC(:,479) * Y(:,202)) &
        + (RC(:,360) * Y(:,3) * Y(:,46)) + (RC(:,361) * Y(:,5) * Y(:,46)) 
        L(:) = 0.0 &
        + (RC(:,335)) + (RC(:,478) * Y(:,4)) &
        + (RC(:,273) * Y(:,9)) + (RC(:,333)) + (RC(:,334)) &
        + (RC(:,224) * Y(:,5)) + (RC(:,225) * Y(:,5)) + (RC(:,226) * Y(:,5)) &
        + (RC(:,141) * Y(:,8)) + (RC(:,142) * Y(:,8)) + (RC(:,143) * Y(:,8)) 
        Y(:,112) = (YP(:,112) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NUCARB12         Y(:,113) 
        P(:) = (DJ(:,72) * Y(:,172)) &
        + (RC(:,339) * Y(:,45)) + (RC(:,441) * Y(:,3) * Y(:,172)) &
        + (RC(:,147) * Y(:,45) * Y(:,8)) + (RC(:,230) * Y(:,45) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,376) * Y(:,3)) + (DJ(:,30)) 
        Y(:,113) = (YP(:,113) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRU12O2          Y(:,114) 
        P(:) = (RC(:,376) * Y(:,3) * Y(:,113)) 
        L(:) = 0.0 &
        + (RC(:,340)) &
        + (RC(:,148) * Y(:,8)) + (RC(:,231) * Y(:,5)) + (RC(:,278) * Y(:,9)) 
        Y(:,114) = (YP(:,114) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NOA              Y(:,115) 
        P(:) = (DJ(:,30) * Y(:,113)) + (DJ(:,73) * Y(:,173)) &
        + (RC(:,340) * Y(:,114)) + (RC(:,442) * Y(:,3) * Y(:,173)) &
        + (RC(:,148) * Y(:,114) * Y(:,8)) + (RC(:,231) * Y(:,114) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,377) * Y(:,3)) + (DJ(:,31)) + (DJ(:,32)) 
        Y(:,115) = (YP(:,115) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN25O2          Y(:,116) 
        P(:) = (RC(:,457) * Y(:,3) * Y(:,177)) + (DJ(:,87) * Y(:,176)) &
        + (RC(:,343) * Y(:,50)) + (RC(:,389) * Y(:,3) * Y(:,52)) &
        + (RC(:,152) * Y(:,50) * Y(:,8)) + (RC(:,234) * Y(:,50) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,282) * Y(:,9)) + (RC(:,344)) &
        + (RC(:,153) * Y(:,8)) + (RC(:,186) * Y(:,8)) + (RC(:,235) * Y(:,5)) 
        Y(:,116) = (YP(:,116) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN24O2          Y(:,117) 
        P(:) = (DJ(:,88) * Y(:,177)) &
        + (RC(:,344) * Y(:,116)) + (RC(:,458) * Y(:,3) * Y(:,178)) &
        + (RC(:,153) * Y(:,116) * Y(:,8)) + (RC(:,235) * Y(:,116) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,345)) &
        + (RC(:,154) * Y(:,8)) + (RC(:,236) * Y(:,5)) + (RC(:,283) * Y(:,9)) 
        Y(:,117) = (YP(:,117) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN23O2          Y(:,118) 
        P(:) = (DJ(:,89) * Y(:,178)) &
        + (RC(:,345) * Y(:,117)) + (RC(:,459) * Y(:,3) * Y(:,179)) &
        + (RC(:,154) * Y(:,117) * Y(:,8)) + (RC(:,236) * Y(:,117) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,346)) &
        + (RC(:,155) * Y(:,8)) + (RC(:,237) * Y(:,5)) + (RC(:,284) * Y(:,9)) 
        Y(:,118) = (YP(:,118) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN14O2          Y(:,119) 
        P(:) = (DJ(:,90) * Y(:,179)) &
        + (RC(:,346) * Y(:,118)) + (RC(:,460) * Y(:,3) * Y(:,180)) &
        + (RC(:,155) * Y(:,118) * Y(:,8)) + (RC(:,237) * Y(:,118) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,347)) &
        + (RC(:,156) * Y(:,8)) + (RC(:,238) * Y(:,5)) + (RC(:,285) * Y(:,9)) 
        Y(:,119) = (YP(:,119) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TNCARB10         Y(:,120) 
        P(:) = (RC(:,347) * Y(:,119)) + (DJ(:,91) * Y(:,180)) &
        + (RC(:,156) * Y(:,119) * Y(:,8)) + (RC(:,238) * Y(:,119) * Y(:,5)) 
        L(:) = 0.0 &
        + (RC(:,386) * Y(:,3)) + (RC(:,388) * Y(:,5)) + (DJ(:,40)) 
        Y(:,120) = (YP(:,120) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN10O2          Y(:,121) 
        P(:) = (RC(:,461) * Y(:,3) * Y(:,181)) &
        + (RC(:,386) * Y(:,3) * Y(:,120)) + (RC(:,388) * Y(:,5) * Y(:,120)) 
        L(:) = 0.0 &
        + (RC(:,348)) &
        + (RC(:,157) * Y(:,8)) + (RC(:,239) * Y(:,5)) + (RC(:,286) * Y(:,9)) 
        Y(:,121) = (YP(:,121) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTX22O2          Y(:,122) 
        P(:) = (RC(:,391) * Y(:,3) * Y(:,58)) 
        L(:) = 0.0 &
        + (RC(:,289) * Y(:,9)) + (RC(:,351)) &
        + (RC(:,163) * Y(:,8)) + (RC(:,189) * Y(:,8)) + (RC(:,242) * Y(:,5)) 
        Y(:,122) = (YP(:,122) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3NO3           Y(:,123) 
        P(:) = (RC(:,166) * Y(:,22) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,392) * Y(:,3)) + (DJ(:,41)) 
        Y(:,123) = (YP(:,123) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H5NO3          Y(:,124) 
        P(:) = (RC(:,167) * Y(:,24) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,393) * Y(:,3)) + (DJ(:,42)) 
        Y(:,124) = (YP(:,124) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN10NO3          Y(:,125) 
        P(:) = (RC(:,168) * Y(:,27) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,394) * Y(:,3)) + (DJ(:,43)) 
        Y(:,125) = (YP(:,125) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          IC3H7NO3         Y(:,126) 
        P(:) = (RC(:,169) * Y(:,26) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,395) * Y(:,3)) + (DJ(:,44)) 
        Y(:,126) = (YP(:,126) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN13NO3          Y(:,127) 
        P(:) = (RC(:,170) * Y(:,29) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,396) * Y(:,3)) + (DJ(:,45)) + (DJ(:,46)) 
        Y(:,127) = (YP(:,127) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN16NO3          Y(:,128) 
        P(:) = (RC(:,171) * Y(:,89) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,397) * Y(:,3)) + (DJ(:,47)) 
        Y(:,128) = (YP(:,128) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN19NO3          Y(:,129) 
        P(:) = (RC(:,172) * Y(:,91) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,398) * Y(:,3)) + (DJ(:,48)) 
        Y(:,129) = (YP(:,129) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HOC2H4NO3        Y(:,130) 
        P(:) = (RC(:,173) * Y(:,31) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,399) * Y(:,3)) 
        Y(:,130) = (YP(:,130) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN9NO3           Y(:,131) 
        P(:) = (RC(:,174) * Y(:,33) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,400) * Y(:,3)) 
        Y(:,131) = (YP(:,131) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN12NO3          Y(:,132) 
        P(:) = (RC(:,175) * Y(:,35) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,401) * Y(:,3)) 
        Y(:,132) = (YP(:,132) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN15NO3          Y(:,133) 
        P(:) = (RC(:,176) * Y(:,95) * Y(:,8)) + (RC(:,178) * Y(:,90) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,402) * Y(:,3)) 
        Y(:,133) = (YP(:,133) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN18NO3          Y(:,134) 
        P(:) = (RC(:,177) * Y(:,103) * Y(:,8)) + (RC(:,179) * Y(:,92) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,403) * Y(:,3)) 
        Y(:,134) = (YP(:,134) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RU14NO3          Y(:,135) 
        P(:) = (RC(:,180) * Y(:,44) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,404) * Y(:,3)) 
        Y(:,135) = (YP(:,135) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA13NO3          Y(:,136) 
        P(:) = (RC(:,181) * Y(:,62) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,405) * Y(:,3)) + (DJ(:,49)) 
        Y(:,136) = (YP(:,136) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA16NO3          Y(:,137) 
        P(:) = (RC(:,182) * Y(:,65) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,406) * Y(:,3)) + (DJ(:,50)) 
        Y(:,137) = (YP(:,137) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA19NO3          Y(:,138) 
        P(:) = (RC(:,183) * Y(:,68) * Y(:,8)) + (RC(:,184) * Y(:,69) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,407) * Y(:,3)) + (DJ(:,51)) 
        Y(:,138) = (YP(:,138) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN28NO3         Y(:,139) 
        P(:) = (RC(:,185) * Y(:,48) * Y(:,8)) + (RC(:,486) * Y(:,204)) 
        L(:) = 0.0 &
        + (RC(:,408) * Y(:,3)) + (RC(:,485)) 
        Y(:,139) = (YP(:,139) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN25NO3         Y(:,140) 
        P(:) = (RC(:,186) * Y(:,116) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,409) * Y(:,3)) 
        Y(:,140) = (YP(:,140) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTX28NO3         Y(:,141) 
        P(:) = (RC(:,187) * Y(:,54) * Y(:,8)) + (RC(:,488) * Y(:,205)) 
        L(:) = 0.0 &
        + (RC(:,410) * Y(:,3)) + (RC(:,487)) 
        Y(:,141) = (YP(:,141) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTX24NO3         Y(:,142) 
        P(:) = (RC(:,188) * Y(:,56) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,411) * Y(:,3)) + (DJ(:,52)) 
        Y(:,142) = (YP(:,142) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTX22NO3         Y(:,143) 
        P(:) = (RC(:,189) * Y(:,122) * Y(:,8)) 
        L(:) = 0.0 &
        + (RC(:,412) * Y(:,3)) 
        Y(:,143) = (YP(:,143) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3OOH           Y(:,144) 
        P(:) = (RC(:,244) * Y(:,22) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,423) * Y(:,3)) + (RC(:,424) * Y(:,3)) + (DJ(:,53)) 
        Y(:,144) = (YP(:,144) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H5OOH          Y(:,145) 
        P(:) = (RC(:,245) * Y(:,24) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,425) * Y(:,3)) + (DJ(:,54)) 
        Y(:,145) = (YP(:,145) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN10OOH          Y(:,146) 
        P(:) = (RC(:,246) * Y(:,27) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,426) * Y(:,3)) + (DJ(:,55)) 
        Y(:,146) = (YP(:,146) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          IC3H7OOH         Y(:,147) 
        P(:) = (RC(:,247) * Y(:,26) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,427) * Y(:,3)) + (DJ(:,56)) 
        Y(:,147) = (YP(:,147) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN13OOH          Y(:,148) 
        P(:) = (RC(:,248) * Y(:,29) * Y(:,9)) + (RC(:,251) * Y(:,93) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,428) * Y(:,3)) + (DJ(:,57)) + (DJ(:,58)) 
        Y(:,148) = (YP(:,148) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN16OOH          Y(:,149) 
        P(:) = (RC(:,249) * Y(:,89) * Y(:,9)) + (RC(:,252) * Y(:,94) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,429) * Y(:,3)) + (DJ(:,59)) 
        Y(:,149) = (YP(:,149) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN19OOH          Y(:,150) 
        P(:) = (RC(:,250) * Y(:,91) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,430) * Y(:,3)) + (DJ(:,60)) 
        Y(:,150) = (YP(:,150) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA13OOH          Y(:,151) 
        P(:) = (RC(:,253) * Y(:,62) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,451) * Y(:,3)) + (DJ(:,82)) 
        Y(:,151) = (YP(:,151) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA16OOH          Y(:,152) 
        P(:) = (RC(:,254) * Y(:,65) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,452) * Y(:,3)) + (DJ(:,83)) 
        Y(:,152) = (YP(:,152) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RA19OOH          Y(:,153) 
        P(:) = (RC(:,255) * Y(:,68) * Y(:,9)) + (RC(:,256) * Y(:,69) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,453) * Y(:,3)) + (DJ(:,84)) 
        Y(:,153) = (YP(:,153) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HOC2H4OOH        Y(:,154) 
        P(:) = (RC(:,257) * Y(:,31) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,443) * Y(:,3)) + (DJ(:,74)) 
        Y(:,154) = (YP(:,154) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN9OOH           Y(:,155) 
        P(:) = (RC(:,258) * Y(:,33) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,444) * Y(:,3)) + (DJ(:,75)) 
        Y(:,155) = (YP(:,155) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN12OOH          Y(:,156) 
        P(:) = (RC(:,259) * Y(:,35) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,445) * Y(:,3)) + (DJ(:,76)) 
        Y(:,156) = (YP(:,156) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN15OOH          Y(:,157) 
        P(:) = (RC(:,260) * Y(:,95) * Y(:,9)) + (RC(:,262) * Y(:,90) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,446) * Y(:,3)) + (DJ(:,77)) 
        Y(:,157) = (YP(:,157) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN18OOH          Y(:,158) 
        P(:) = (RC(:,261) * Y(:,103) * Y(:,9)) + (RC(:,263) * Y(:,92) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,447) * Y(:,3)) + (DJ(:,78)) 
        Y(:,158) = (YP(:,158) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3CO3H          Y(:,159) 
        P(:) = (RC(:,264) * Y(:,70) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,431) * Y(:,3)) + (DJ(:,61)) 
        Y(:,159) = (YP(:,159) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          C2H5CO3H         Y(:,160) 
        P(:) = (RC(:,265) * Y(:,72) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,432) * Y(:,3)) + (DJ(:,62)) 
        Y(:,160) = (YP(:,160) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          HOCH2CO3H        Y(:,161) 
        P(:) = (RC(:,266) * Y(:,106) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,433) * Y(:,3)) + (DJ(:,63)) 
        Y(:,161) = (YP(:,161) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN8OOH           Y(:,162) 
        P(:) = (RC(:,267) * Y(:,74) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,434) * Y(:,3)) + (DJ(:,64)) 
        Y(:,162) = (YP(:,162) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN11OOH          Y(:,163) 
        P(:) = (RC(:,268) * Y(:,75) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,435) * Y(:,3)) + (DJ(:,65)) 
        Y(:,163) = (YP(:,163) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN14OOH          Y(:,164) 
        P(:) = (RC(:,269) * Y(:,107) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,436) * Y(:,3)) + (DJ(:,66)) 
        Y(:,164) = (YP(:,164) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RN17OOH          Y(:,165) 
        P(:) = (RC(:,270) * Y(:,108) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,437) * Y(:,3)) + (DJ(:,67)) 
        Y(:,165) = (YP(:,165) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RU14OOH          Y(:,166) 
        P(:) = (RC(:,271) * Y(:,44) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,438) * Y(:,3)) + (DJ(:,68)) + (DJ(:,69)) 
        Y(:,166) = (YP(:,166) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RU12OOH          Y(:,167) 
        P(:) = (RC(:,272) * Y(:,110) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,439) * Y(:,3)) + (DJ(:,70)) 
        Y(:,167) = (YP(:,167) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RU10OOH          Y(:,168) 
        P(:) = (RC(:,273) * Y(:,112) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,440) * Y(:,3)) + (DJ(:,71)) 
        Y(:,168) = (YP(:,168) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRN6OOH          Y(:,169) 
        P(:) = (RC(:,274) * Y(:,36) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,448) * Y(:,3)) + (DJ(:,79)) 
        Y(:,169) = (YP(:,169) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRN9OOH          Y(:,170) 
        P(:) = (RC(:,275) * Y(:,37) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,449) * Y(:,3)) + (DJ(:,80)) 
        Y(:,170) = (YP(:,170) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRN12OOH         Y(:,171) 
        P(:) = (RC(:,276) * Y(:,38) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,450) * Y(:,3)) + (DJ(:,81)) 
        Y(:,171) = (YP(:,171) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRU14OOH         Y(:,172) 
        P(:) = (RC(:,277) * Y(:,45) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,441) * Y(:,3)) + (DJ(:,72)) 
        Y(:,172) = (YP(:,172) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRU12OOH         Y(:,173) 
        P(:) = (RC(:,278) * Y(:,114) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,442) * Y(:,3)) + (DJ(:,73)) 
        Y(:,173) = (YP(:,173) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN28OOH         Y(:,174) 
        P(:) = (RC(:,279) * Y(:,48) * Y(:,9)) + (RC(:,496) * Y(:,209)) 
        L(:) = 0.0 &
        + (RC(:,454) * Y(:,3)) + (RC(:,495)) + (DJ(:,85)) 
        Y(:,174) = (YP(:,174) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRTN28OOH        Y(:,175) 
        P(:) = (RC(:,280) * Y(:,49) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,456) * Y(:,3)) + (DJ(:,86)) 
        Y(:,175) = (YP(:,175) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN26OOH         Y(:,176) 
        P(:) = (RC(:,281) * Y(:,50) * Y(:,9)) + (RC(:,498) * Y(:,210)) 
        L(:) = 0.0 &
        + (RC(:,455) * Y(:,3)) + (RC(:,497)) + (DJ(:,87)) 
        Y(:,176) = (YP(:,176) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN25OOH         Y(:,177) 
        P(:) = (RC(:,282) * Y(:,116) * Y(:,9)) + (RC(:,502) * Y(:,212)) 
        L(:) = 0.0 &
        + (RC(:,457) * Y(:,3)) + (RC(:,501)) + (DJ(:,88)) 
        Y(:,177) = (YP(:,177) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN24OOH         Y(:,178) 
        P(:) = (RC(:,283) * Y(:,117) * Y(:,9)) + (RC(:,492) * Y(:,207)) 
        L(:) = 0.0 &
        + (RC(:,458) * Y(:,3)) + (RC(:,491)) + (DJ(:,89)) 
        Y(:,178) = (YP(:,178) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN23OOH         Y(:,179) 
        P(:) = (RC(:,284) * Y(:,118) * Y(:,9)) + (RC(:,504) * Y(:,213)) 
        L(:) = 0.0 &
        + (RC(:,459) * Y(:,3)) + (RC(:,503)) + (DJ(:,90)) 
        Y(:,179) = (YP(:,179) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN14OOH         Y(:,180) 
        P(:) = (RC(:,285) * Y(:,119) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,460) * Y(:,3)) + (DJ(:,91)) 
        Y(:,180) = (YP(:,180) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN10OOH         Y(:,181) 
        P(:) = (RC(:,286) * Y(:,121) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,461) * Y(:,3)) + (DJ(:,92)) 
        Y(:,181) = (YP(:,181) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTX28OOH         Y(:,182) 
        P(:) = (RC(:,287) * Y(:,54) * Y(:,9)) + (RC(:,494) * Y(:,208)) 
        L(:) = 0.0 &
        + (RC(:,462) * Y(:,3)) + (RC(:,493)) + (DJ(:,93)) 
        Y(:,182) = (YP(:,182) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTX24OOH         Y(:,183) 
        P(:) = (RC(:,288) * Y(:,56) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,463) * Y(:,3)) + (DJ(:,94)) 
        Y(:,183) = (YP(:,183) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTX22OOH         Y(:,184) 
        P(:) = (RC(:,289) * Y(:,122) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,464) * Y(:,3)) + (DJ(:,95)) 
        Y(:,184) = (YP(:,184) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          NRTX28OOH        Y(:,185) 
        P(:) = (RC(:,290) * Y(:,55) * Y(:,9)) 
        L(:) = 0.0 &
        + (RC(:,465) * Y(:,3)) + (DJ(:,96)) 
        Y(:,185) = (YP(:,185) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB14           Y(:,186) 
        P(:) = (RC(:,397) * Y(:,3) * Y(:,128)) + (RC(:,429) * Y(:,3) * Y(:,149)) 
        L(:) = 0.0 &
        + (RC(:,353) * Y(:,3)) + (DJ(:,15)) 
        Y(:,186) = (YP(:,186) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB17           Y(:,187) 
        P(:) = (RC(:,398) * Y(:,3) * Y(:,129)) + (RC(:,430) * Y(:,3) * Y(:,150)) 
        L(:) = 0.0 &
        + (RC(:,354) * Y(:,3)) + (DJ(:,16)) 
        Y(:,187) = (YP(:,187) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB10           Y(:,188) 
        P(:) = (RC(:,401) * Y(:,3) * Y(:,132)) + (RC(:,445) * Y(:,3) * Y(:,156)) 
        L(:) = 0.0 &
        + (RC(:,357) * Y(:,3)) + (DJ(:,19)) 
        Y(:,188) = (YP(:,188) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB12           Y(:,189) 
        P(:) = (RC(:,436) * Y(:,3) * Y(:,164)) 
        L(:) = 0.0 &
        + (RC(:,369) * Y(:,3)) + (DJ(:,27)) 
        Y(:,189) = (YP(:,189) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CARB15           Y(:,190) 
        P(:) = (RC(:,437) * Y(:,3) * Y(:,165)) 
        L(:) = 0.0 &
        + (RC(:,370) * Y(:,3)) + (DJ(:,28)) 
        Y(:,190) = (YP(:,190) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CCARB12          Y(:,191) 
        P(:) = (RC(:,412) * Y(:,3) * Y(:,143)) + (RC(:,464) * Y(:,3) * Y(:,184)) 
        L(:) = 0.0 &
        + (RC(:,371) * Y(:,3)) 
        Y(:,191) = (YP(:,191) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          ANHY             Y(:,192) 
        P(:) = (DJ(:,38) * Y(:,99)) &
        + (DJ(:,34) * Y(:,96)) + (DJ(:,36) * Y(:,97)) &
        + (RC(:,383) * Y(:,3) * Y(:,99)) + (RC(:,510) * Y(:,216)) &
        + (RC(:,379) * Y(:,3) * Y(:,96)) + (RC(:,381) * Y(:,3) * Y(:,97)) 
        L(:) = 0.0 &
        + (RC(:,466) * Y(:,3)) + (RC(:,509)) 
        Y(:,192) = (YP(:,192) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          TNCARB15         Y(:,193) 
        P(:) = (RC(:,409) * Y(:,3) * Y(:,140)) 
        L(:) = 0.0 &
        + (RC(:,385) * Y(:,3)) 
        Y(:,193) = (YP(:,193) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RAROH14          Y(:,194) 
        P(:) = (RC(:,413) * Y(:,3) * Y(:,63)) + (RC(:,414) * Y(:,5) * Y(:,63)) 
        L(:) = 0.0 &
        + (RC(:,415) * Y(:,4)) 
        Y(:,194) = (YP(:,194) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          ARNOH14          Y(:,195) 
        P(:) = (RC(:,415) * Y(:,194) * Y(:,4)) + (RC(:,506) * Y(:,214)) 
        L(:) = 0.0 &
        + (RC(:,416) * Y(:,3)) + (RC(:,417) * Y(:,5)) + (RC(:,505)) 
        Y(:,195) = (YP(:,195) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RAROH17          Y(:,196) 
        P(:) = (RC(:,418) * Y(:,3) * Y(:,66)) + (RC(:,419) * Y(:,5) * Y(:,66)) 
        L(:) = 0.0 &
        + (RC(:,420) * Y(:,4)) 
        Y(:,196) = (YP(:,196) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          ARNOH17          Y(:,197) 
        P(:) = (RC(:,420) * Y(:,196) * Y(:,4)) + (RC(:,508) * Y(:,215)) 
        L(:) = 0.0 &
        + (RC(:,421) * Y(:,3)) + (RC(:,422) * Y(:,5)) + (RC(:,507)) 
        Y(:,197) = (YP(:,197) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          PAN              Y(:,198) 
        P(:) = (RC(:,467) * Y(:,70) * Y(:,4)) 
        L(:) = 0.0 &
        + (RC(:,468)) + (RC(:,473) * Y(:,3)) 
        Y(:,198) = (YP(:,198) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          PPN              Y(:,199) 
        P(:) = (RC(:,469) * Y(:,72) * Y(:,4)) 
        L(:) = 0.0 &
        + (RC(:,470)) + (RC(:,474) * Y(:,3)) 
        Y(:,199) = (YP(:,199) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          PHAN             Y(:,200) 
        P(:) = (RC(:,471) * Y(:,106) * Y(:,4)) 
        L(:) = 0.0 &
        + (RC(:,472)) + (RC(:,475) * Y(:,3)) 
        Y(:,200) = (YP(:,200) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RU12PAN          Y(:,201) 
        P(:) = (RC(:,476) * Y(:,110) * Y(:,4)) 
        L(:) = 0.0 &
        + (RC(:,477)) + (RC(:,481) * Y(:,3)) 
        Y(:,201) = (YP(:,201) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          MPAN             Y(:,202) 
        P(:) = (RC(:,478) * Y(:,112) * Y(:,4)) 
        L(:) = 0.0 &
        + (RC(:,479)) + (RC(:,480) * Y(:,3)) 
        Y(:,202) = (YP(:,202) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          RTN26PAN         Y(:,203) 
        P(:) = (RC(:,482) * Y(:,50) * Y(:,4)) + (RC(:,500) * Y(:,211)) 
        L(:) = 0.0 &
        + (RC(:,483)) + (RC(:,484) * Y(:,3)) + (RC(:,499)) 
        Y(:,203) = (YP(:,203) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P2604            Y(:,204) 
        P(:) = (RC(:,485) * Y(:,139)) 
        L(:) = 0.0 &
        + (RC(:,486)) 
        Y(:,204) = (YP(:,204) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P4608            Y(:,205) 
        P(:) = (RC(:,487) * Y(:,141)) 
        L(:) = 0.0 &
        + (RC(:,488)) 
        Y(:,205) = (YP(:,205) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P2631            Y(:,206) 
        P(:) = (RC(:,489) * Y(:,52)) 
        L(:) = 0.0 &
        + (RC(:,490)) 
        Y(:,206) = (YP(:,206) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P2635            Y(:,207) 
        P(:) = (RC(:,491) * Y(:,178)) 
        L(:) = 0.0 &
        + (RC(:,492)) 
        Y(:,207) = (YP(:,207) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P4610            Y(:,208) 
        P(:) = (RC(:,493) * Y(:,182)) 
        L(:) = 0.0 &
        + (RC(:,494)) 
        Y(:,208) = (YP(:,208) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P2605            Y(:,209) 
        P(:) = (RC(:,495) * Y(:,174)) 
        L(:) = 0.0 &
        + (RC(:,496)) 
        Y(:,209) = (YP(:,209) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P2630            Y(:,210) 
        P(:) = (RC(:,497) * Y(:,176)) 
        L(:) = 0.0 &
        + (RC(:,498)) 
        Y(:,210) = (YP(:,210) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P2629            Y(:,211) 
        P(:) = (RC(:,499) * Y(:,203)) 
        L(:) = 0.0 &
        + (RC(:,500)) 
        Y(:,211) = (YP(:,211) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P2632            Y(:,212) 
        P(:) = (RC(:,501) * Y(:,177)) 
        L(:) = 0.0 &
        + (RC(:,502)) 
        Y(:,212) = (YP(:,212) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P2637            Y(:,213) 
        P(:) = (RC(:,503) * Y(:,179)) 
        L(:) = 0.0 &
        + (RC(:,504)) 
        Y(:,213) = (YP(:,213) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P3612            Y(:,214) 
        P(:) = (RC(:,505) * Y(:,195)) 
        L(:) = 0.0 &
        + (RC(:,506)) 
        Y(:,214) = (YP(:,214) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P3613            Y(:,215) 
        P(:) = (RC(:,507) * Y(:,197)) 
        L(:) = 0.0 &
        + (RC(:,508)) 
        Y(:,215) = (YP(:,215) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P3442            Y(:,216) 
        P(:) = (RC(:,509) * Y(:,192)) 
        L(:) = 0.0 &
        + (RC(:,510)) 
        Y(:,216) = (YP(:,216) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          CH3O2NO2         Y(:,217) 
        P(:) = (RC(:,164) * Y(:,22) * Y(:,4)) 
        L(:) = 0.0 &
        + (RC(:,165)) 
        Y(:,217) = (YP(:,217) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          EMPOA            Y(:,218) 
        P(:) = 0.0 
        !     dry deposition of EMPOA
        L(:) = 0.0
        Y(:,218) = (YP(:,218) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !          P2007            Y(:,219) 
        ! P(:) = EM(:,219) &
        ! + (RC(:,511) * Y(:,167)) 
        ! L(:) = 0.0 &
        ! + (RC(:,512)) 
        ! Y(:,219) = (YP(:,219) + DTS * P(:)) /(1.0+ DTS * L(:)) 

        !  DEFINE TOTAL CONCENTRATION OF PEROXY RADICALS

        RO2(:) =  Y(:,22) + Y(:,24)+ Y(:,27) + Y(:,26) + Y(:,29) + &
        Y(:,93) + Y(:,94) + Y(:,89) + Y(:,91) + Y(:,31) + Y(:,33) + &
        Y(:,35) + Y(:,95) + Y(:,103) + Y(:,90) + Y(:,92) + Y(:,70) + &
        Y(:,72) + Y(:,75) + Y(:,107) + Y(:,108) + Y(:,106) + Y(:,44) + &
        Y(:,110) + Y(:,112) + Y(:,36) + Y(:,37) + Y(:,38) + Y(:,48) + &
        Y(:,45) + Y(:,114) + Y(:,62) + Y(:,65) + Y(:,68) + Y(:,69) + &
        Y(:,74) + Y(:,50) + Y(:,49) + Y(:,116) + Y(:,117) + Y(:,118) + &
        Y(:,119) + Y(:,121) + Y(:,54) + Y(:,56) + Y(:,122) + Y(:,55) 

        ! CALCULATE FLUX TERMS

        !      O + O2 + M = O3 + M
        FL(:,1)=FL(:,1)+RC(:,1)*Y(:,2)*DTS
        !      O + N2 + M = O3 + M
        FL(:,2)=FL(:,2)+RC(:,2)*Y(:,2)*DTS
        !      O + O3 = 
        FL(:,3)=FL(:,3)+RC(:,3)*Y(:,2)*Y(:,6)*DTS
        !       O + NO = NO2
        FL(:,4)=FL(:,4)+RC(:,4)*Y(:,2)*Y(:,8)*DTS
        !      O + NO2 = NO
        FL(:,5)=FL(:,5)+RC(:,5)*Y(:,2)*Y(:,4)*DTS
        !      O + NO2 = NO3
        FL(:,6)=FL(:,6)+RC(:,6)*Y(:,2)*Y(:,4)*DTS
        !      O1D + O2 + M = O + M
        FL(:,7)=FL(:,7)+RC(:,7)*Y(:,1)*DTS
        !      O1D + N2 + M = O + M
        FL(:,8)=FL(:,8)+RC(:,8)*Y(:,1)*DTS
        !      NO + O3 = NO2
        FL(:,9)=FL(:,9) + RC(:,9)*Y(:,8)*Y(:,6)*DTS
        !      NO2 + O3 = NO3
        FL(:,10)=FL(:,10)+RC(:,10)*Y(:,4)*Y(:,6)*DTS
        !      NO + NO = NO2 + NO2
        FL(:,11)=FL(:,11)+RC(:,11)*Y(:,4)*Y(:,4)*DTS
        !      NO + NO3 = NO2 + NO2
        FL(:,12)=FL(:,12)+RC(:,12)*Y(:,4)*Y(:,5)*DTS
        !      NO2 + NO3 = NO + NO2
        FL(:,13)=FL(:,13)+RC(:,13)*Y(:,4)*Y(:,5)*DTS
        !      NO2 + NO3 = N2O5
        FL(:,14)=FL(:,14)+RC(:,14)*Y(:,4)*Y(:,5)*DTS
        !      N2O5 = NO2 + NO3
        FL(:,15)=FL(:,15)+RC(:,15)*Y(:,7)*DTS
        !      O1D = OH + OH
        FL(:,16)=FL(:,16)+RC(:,16)*Y(:,1)*DTS
        !      OH + O3 = HO2
        FL(:,17)=FL(:,17)+RC(:,17)*Y(:,3)*Y(:,6)*DTS
        !      OH + H2 = HO2
        FL(:,18)=FL(:,18)+RC(:,18)*Y(:,3)*Y(:,10)*DTS
        !       OH + CO = HO2
        FL(:,19)=FL(:,19)+RC(:,19)*Y(:,3)*Y(:,11)*DTS
        !      OH + H2O2 = HO2
        FL(:,20)=FL(:,20)+RC(:,20)*Y(:,3)*Y(:,12)*DTS
        !      HO2 + O3 = OH
        FL(:,21)=FL(:,21)+RC(:,21)*Y(:,9)*Y(:,6)*DTS
        !      OH + HO2 = 
        FL(:,22)=FL(:,22)+RC(:,22)*Y(:,3)*Y(:,9)*DTS
        !      HO2 + HO2 = H2O2
        FL(:,23)=FL(:,23)+RC(:,23)*Y(:,9)*Y(:,9)*DTS
        !      HO2 + HO2 = H2O2
        FL(:,24)=FL(:,24)+RC(:,24)*Y(:,9)*Y(:,9)*DTS
        !      OH + NO = HONO
        FL(:,25)=FL(:,25)+RC(:,25)*Y(:,3)*Y(:,8)*DTS
        !      NO2 = HONO
        FL(:,26)=FL(:,26)+RC(:,26)*Y(:,4)*DTS
        !      OH + NO2 = HNO3
        FL(:,27)=FL(:,27)+RC(:,27)*Y(:,3)*Y(:,4)*DTS
        !      OH + NO3 = HO2 + NO2
        FL(:,28)=FL(:,28)+RC(:,28)*Y(:,3)*Y(:,5)*DTS
        !      HO2 + NO = OH + NO2
        FL(:,29)=FL(:,29)+RC(:,29)*Y(:,9)*Y(:,8)*DTS
        !      HO2 + NO2 = HO2NO2
        FL(:,30)=FL(:,30)+RC(:,30)*Y(:,9)*Y(:,4)*DTS
        !      HO2NO2 = HO2 + NO2
        FL(:,31)=FL(:,31)+RC(:,31)*Y(:,15)*DTS
        !      OH + HO2NO2 = NO2 
        FL(:,32)=FL(:,32)+RC(:,32)*Y(:,3)*Y(:,15)*DTS
        !      HO2 + NO3 = OH + NO2
        FL(:,33)=FL(:,33)+RC(:,33)*Y(:,9)*Y(:,5)*DTS
        !      OH + HONO = NO2
        FL(:,34)=FL(:,34)+RC(:,34)*Y(:,3)*Y(:,13)*DTS
        !      OH + HNO3 = NO3
        FL(:,35)=FL(:,35)+RC(:,35)*Y(:,3)*Y(:,14)*DTS
        !      O + SO2 = SO3
        FL(:,36)=FL(:,36)+RC(:,36)*Y(:,2)*Y(:,16)*DTS
        !      OH + SO2 = HSO3 
        FL(:,37)=FL(:,37)+RC(:,37)*Y(:,3)*Y(:,16)*DTS
        !      HSO3 = HO2 + SO3
        FL(:,38)=FL(:,38)+RC(:,38)*Y(:,18)*DTS
        !      HNO3 = NA
        FL(:,39)=FL(:,39)+RC(:,39)*Y(:,14)*DTS
        !      N2O5 = NA + NA
        FL(:,40)=FL(:,40)+RC(:,40)*Y(:,7)*DTS
        !      SO3 = SA
        FL(:,41)=FL(:,41)+RC(:,41)*Y(:,17)*DTS
        !      OH + CH4 = CH3O2
        FL(:,42)=FL(:,42)+RC(:,42)*Y(:,3)*Y(:,21)*DTS
        !      OH + C2H6 = C2H5O2
        FL(:,43)=FL(:,43)+RC(:,43)*Y(:,3)*Y(:,26)*DTS
        !      OH + C3H8 = IC3H7O2
        FL(:,44)=FL(:,44)+RC(:,44)*Y(:,3)*Y(:,25)*DTS
        !      OH + C3H8 = RN10O2 
        FL(:,45)=FL(:,45)+RC(:,45)*Y(:,3)*Y(:,25)*DTS
        !      OH + NC4H10 = RN13O2
        FL(:,46)=FL(:,46)+RC(:,46)*Y(:,3)*Y(:,28)*DTS
        !      OH + C2H4 = HOCH2CH2O2
        FL(:,47)=FL(:,47)+RC(:,47)*Y(:,3)*Y(:,30)*DTS
        !      OH + C3H6 = RN9O2
        FL(:,48)=FL(:,48)+RC(:,48)*Y(:,3)*Y(:,32)*DTS
        !      OH + TBUT2ENE = RN12O2
        FL(:,49)=FL(:,49)+RC(:,49)*Y(:,3)*Y(:,34)*DTS
        !      NO3 + C2H4 = NRN6O2
        FL(:,50)=FL(:,50)+RC(:,50)*Y(:,3)*Y(:,30)*DTS
        !      NO3 + C3H6 = NRN9O2
        FL(:,51)=FL(:,51)+RC(:,51)*Y(:,5)*Y(:,32)*DTS
        !      NO3 + TBUT2ENE = NRN12O2
        FL(:,52)=FL(:,52)+RC(:,52)*Y(:,5)*Y(:,34)*DTS
        !      O3 + C2H4 = HCHO + CO + HO2 + OH
        FL(:,53)=FL(:,53)+RC(:,53)*Y(:,6)*Y(:,30)*DTS
        !      O3 + C2H4 = HCHO + HCOOH
        FL(:,54)=FL(:,54)+RC(:,54)*Y(:,6)*Y(:,30)*DTS
        !      O3 + C3H6 = HCHO + CO + CH3O2 + OH
        FL(:,55)=FL(:,55)+RC(:,55)*Y(:,6)*Y(:,32)*DTS
        !      O3 + C3H6 = HCHO + CH3CO2H
        FL(:,56)=FL(:,56)+RC(:,56)*Y(:,6)*Y(:,32)*DTS
        !      O3 + TBUT2ENE = CH3CHO + CO + CH3O2 + OH
        FL(:,57)=FL(:,57)+RC(:,57)*Y(:,6)*Y(:,34)*DTS
        !      O3 + TBUT2ENE = CH3CHO + CH3CO2H
        FL(:,58)=FL(:,58)+RC(:,58)*Y(:,6)*Y(:,34)*DTS
        !      OH + C5H8 = RU14O2
        FL(:,59)=FL(:,59)+RC(:,59)*Y(:,3)*Y(:,43)*DTS
        !      NO3 + C5H8 = NRU14O2
        FL(:,60)=FL(:,60)+RC(:,60)*Y(:,5)*Y(:,43)*DTS
        !      O3 + C5H8 = UCARB10 + CO + HO2 + OH
        FL(:,61)=FL(:,61)+RC(:,61)*Y(:,6)*Y(:,43)*DTS
        !      O3 + C5H8 = UCARB10 + HCOOH
        FL(:,62)=FL(:,62)+RC(:,62)*Y(:,6)*Y(:,43)*DTS
        !      APINENE + OH = RTN28O2
        FL(:,63)=FL(:,63)+RC(:,63)*Y(:,47)*Y(:,3)*DTS
        !      APINENE + NO3 = NRTN28O2
        FL(:,64)=FL(:,64)+RC(:,64)*Y(:,47)*Y(:,5)*DTS
        !      APINENE + O3 = OH + RTN26O2 
        FL(:,65)=FL(:,65)+RC(:,65)*Y(:,47)*Y(:,6)*DTS
        !      APINENE + O3 = TNCARB26 + H2O2
        FL(:,66)=FL(:,66)+RC(:,66)*Y(:,47)*Y(:,6)*DTS
        !      APINENE + O3 = RCOOH25 
        FL(:,67)=FL(:,67)+RC(:,67)*Y(:,47)*Y(:,6)*DTS
        !      BPINENE + OH = RTX28O2
        FL(:,68)=FL(:,68)+RC(:,68)*Y(:,53)*Y(:,3)*DTS
        !      BPINENE + NO3 = NRTX28O2
        FL(:,69)=FL(:,69)+RC(:,69)*Y(:,53)*Y(:,5)*DTS
        !      BPINENE + O3 =  RTX24O2 + OH
        FL(:,70)=FL(:,70)+RC(:,70)*Y(:,53)*Y(:,6)*DTS
        !      BPINENE + O3 =  HCHO + TXCARB24 + H2O2
        FL(:,71)=FL(:,71)+RC(:,71)*Y(:,53)*Y(:,6)*DTS
        !      BPINENE + O3 =  HCHO + TXCARB22
        FL(:,72)=FL(:,72)+RC(:,72)*Y(:,53)*Y(:,6)*DTS
        !      BPINENE + O3 =  TXCARB24 + CO 
        FL(:,73)=FL(:,73)+RC(:,73)*Y(:,53)*Y(:,6)*DTS
        !      C2H2 + OH = HCOOH + CO + HO2
        FL(:,74)=FL(:,74)+RC(:,74)*Y(:,59)*Y(:,3)*DTS
        !      C2H2 + OH = CARB3 + OH
        FL(:,75)=FL(:,75)+RC(:,75)*Y(:,59)*Y(:,3)*DTS
        !      BENZENE + OH = RA13O2
        FL(:,76)=FL(:,76)+RC(:,76)*Y(:,61)*Y(:,3)*DTS
        !      BENZENE + OH = AROH14 + HO2
        FL(:,77)=FL(:,77)+RC(:,77)*Y(:,61)*Y(:,3)*DTS
        !      TOLUENE + OH = RA16O2
        FL(:,78)=FL(:,78)+RC(:,78)*Y(:,64)*Y(:,3)*DTS
        !      TOLUENE + OH = AROH17 + HO2
        FL(:,79)=FL(:,79)+RC(:,79)*Y(:,64)*Y(:,3)*DTS
        !      OXYL + OH = RA19AO2
        FL(:,80)=FL(:,80)+RC(:,80)*Y(:,67)*Y(:,3)*DTS
        !      OXYL + OH = RA19CO2
        FL(:,81)=FL(:,81)+RC(:,81)*Y(:,67)*Y(:,3)*DTS
        !      OH + HCHO = HO2 + CO
        FL(:,82)=FL(:,82)+RC(:,82)*Y(:,3)*Y(:,39)*DTS
        !      OH + CH3CHO = CH3CO3
        FL(:,83)=FL(:,83)+RC(:,83)*Y(:,3)*Y(:,42)*DTS
        !      OH + C2H5CHO = C2H5CO3
        FL(:,84)=FL(:,84)+RC(:,84)*Y(:,3)*Y(:,71)*DTS
        !      NO3 + HCHO = HO2 + CO + HNO3
        FL(:,85)=FL(:,85)+RC(:,85)*Y(:,5)*Y(:,39)*DTS
        !      NO3 + CH3CHO = CH3CO3 + HNO3
        FL(:,86)=FL(:,86)+RC(:,86)*Y(:,5)*Y(:,42)*DTS
        !      NO3 + C2H5CHO = C2H5CO3 + HNO3
        FL(:,87)=FL(:,87)+RC(:,87)*Y(:,5)*Y(:,71)*DTS
        !      OH + CH3COCH3 = RN8O2
        FL(:,88)=FL(:,88)+RC(:,88)*Y(:,3)*Y(:,73)*DTS
        !      MEK + OH = RN11O2
        FL(:,89)=FL(:,89)+RC(:,89)*Y(:,101)*Y(:,3)*DTS
        !      OH + CH3OH = HO2 + HCHO
        FL(:,90)=FL(:,90)+RC(:,90)*Y(:,3)*Y(:,76)*DTS
        !      OH + C2H5OH = CH3CHO + HO2
        FL(:,91)=FL(:,91)+RC(:,91)*Y(:,3)*Y(:,76)*DTS
        !      OH + C2H5OH = HOCH2CH2O2 
        FL(:,92)=FL(:,92)+RC(:,92)*Y(:,3)*Y(:,77)*DTS
        !     NPROPOL + OH = C2H5CHO + HO2 
        FL(:,93)=FL(:,93)+RC(:,93)*Y(:,3)*Y(:,78)*DTS
        !      NPROPOL + OH = RN9O2
        FL(:,94)=FL(:,94)+RC(:,94)*Y(:,3)*Y(:,78)*DTS
        !      OH + IPROPOL = CH3COCH3 + HO2
        FL(:,95)=FL(:,95)+RC(:,95)*Y(:,3)*Y(:,79)*DTS
        !      OH + IPROPOL = RN9O2
        FL(:,96)=FL(:,96)+RC(:,96)*Y(:,3)*Y(:,79)*DTS
        !      HCOOH + OH = HO2
        FL(:,97)=FL(:,97)+RC(:,97)*Y(:,3)*Y(:,40)*DTS
        !      CH3CO2H + OH = CH3O2
        FL(:,98)=FL(:,98)+RC(:,98)*Y(:,3)*Y(:,41)*DTS
        !      OH + CH3CL = CH3O2 
        FL(:,99)=FL(:,99)+RC(:,99)*Y(:,3)*Y(:,80)*DTS
        !      OH + CH2CL2 = CH3O2
        FL(:,100)=FL(:,100)+RC(:,100)*Y(:,3)*Y(:,81)*DTS
        !      OH + CHCL3 = CH3O2
        FL(:,101)=FL(:,101)+RC(:,101)*Y(:,3)*Y(:,80)*DTS
        !      OH + CH3CCL3 = C2H5O2
        FL(:,102)=FL(:,102)+RC(:,102)*Y(:,3)*Y(:,83)*DTS
        !      OH + TCE = HOCH2CH2O2 
        FL(:,103)=FL(:,103)+RC(:,103)*Y(:,3)*Y(:,84)*DTS
        !      OH + TRICLETH = HOCH2CH2O2
        FL(:,104)=FL(:,104)+RC(:,104)*Y(:,3)*Y(:,85)*DTS
        !      OH + CDICLETH = HOCH2CH2O2
        FL(:,105)=FL(:,105)+RC(:,105)*Y(:,3)*Y(:,86)*DTS
        !      OH + TDICLETH = HOCH2CH2O2
        FL(:,106)=FL(:,106)+RC(:,106)*Y(:,3)*Y(:,87)*DTS
        !      CH3O2 + NO = HCHO + HO2 + NO2
        FL(:,107)=FL(:,107)+RC(:,107)*Y(:,8)*Y(:,22)*DTS
        !      C2H5O2 + NO = CH3CHO + HO2 + NO2
        FL(:,108)=FL(:,108)+RC(:,108)*Y(:,8)*Y(:,24)*DTS
        !      RN10O2 + NO = C2H5CHO + HO2 + NO2
        FL(:,109)=FL(:,109)+RC(:,109)*Y(:,8)*Y(:,27)*DTS
        !      IC3H7O2 + NO = CH3COCH3 + HO2 + NO2
        FL(:,110)=FL(:,110)+RC(:,110)*Y(:,8)*Y(:,26)*DTS
        !      RN13O2 + NO = CH3CHO + C2H5O2 + NO2 
        FL(:,111)=FL(:,111)+RC(:,111)*Y(:,8)*Y(:,29)*DTS
        !      RN13O2 + NO = CARB11A + HO2 + NO2
        FL(:,112)=FL(:,112)+RC(:,112)*Y(:,8)*Y(:,29)*DTS
        !      RN16O2 + NO = RN15AO2 + NO2 
        FL(:,113)=FL(:,113)+RC(:,113)*Y(:,8)*Y(:,89)*DTS
        !      RN19O2 + NO = RN18AO2 + NO2
        FL(:,114)=FL(:,114)+RC(:,114)*Y(:,8)*Y(:,91)*DTS
        !      RN13AO2 + NO = RN12O2 + NO2
        FL(:,115)=FL(:,115)+RC(:,115)*Y(:,8)*Y(:,93)*DTS
        !      RN16AO2 + NO = RN15O2 + NO2 
        FL(:,116)=FL(:,116)+RC(:,116)*Y(:,8)*Y(:,94)*DTS
        !      RA13O2 + NO = CARB3 + UDCARB8 + HO2 + NO2
        FL(:,117)=FL(:,117)+RC(:,117)*Y(:,8)*Y(:,62)*DTS
        !      RA16O2 + NO = CARB3 + UDCARB11 + HO2 + NO2 
        FL(:,118)=FL(:,118)+RC(:,118)*Y(:,8)*Y(:,65)*DTS
        !      RA16O2 + NO = CARB6 + UDCARB8 + HO2 + NO2  
        FL(:,119)=FL(:,119)+RC(:,119)*Y(:,8)*Y(:,65)*DTS
        !      RA19AO2 + NO = CARB3 + UDCARB14 + HO2 + NO2
        FL(:,120)=FL(:,120)+RC(:,120)*Y(:,8)*Y(:,68)*DTS
        !      RA19CO2 + NO = CARB9 + UDCARB8 + HO2 + NO2
        FL(:,121)=FL(:,121)+RC(:,121)*Y(:,8)*Y(:,69)*DTS
        !      HOCH2CH2O2 + NO = HCHO + HCHO + HO2 + NO2
        FL(:,122)=FL(:,122)+RC(:,122)*Y(:,8)*Y(:,31)*DTS
        !      HOCH2CH2O2 + NO = HOCH2CHO + HO2 + NO2
        FL(:,123)=FL(:,123)+RC(:,123)*Y(:,8)*Y(:,31)*DTS
        !       RN9O2 + NO = CH3CHO + HCHO + HO2 + NO2 
        FL(:,124)=FL(:,124)+RC(:,124)*Y(:,8)*Y(:,33)*DTS
        !      RN12O2 + NO = CH3CHO + CH3CHO + HO2 + NO2
        FL(:,125)=FL(:,125)+RC(:,125)*Y(:,8)*Y(:,35)*DTS
        !      RN15O2 + NO = C2H5CHO + CH3CHO + HO2 + NO2
        FL(:,126)=FL(:,126)+RC(:,126)*Y(:,8)*Y(:,95)*DTS
        !      RN18O2 + NO = C2H5CHO + C2H5CHO + HO2 + NO2 
        FL(:,127)=FL(:,127)+RC(:,127)*Y(:,8)*Y(:,103)*DTS
        !      RN15AO2 + NO = CARB13 + HO2 + NO2 
        FL(:,128)=FL(:,128)+RC(:,128)*Y(:,8)*Y(:,90)*DTS
        !      RN18AO2 + NO = CARB16 + HO2 + NO2
        FL(:,129)=FL(:,129)+RC(:,129)*Y(:,8)*Y(:,92)*DTS
        !      CH3CO3 + NO = CH3O2 + NO2 
        FL(:,130)=FL(:,130)+RC(:,130)*Y(:,8)*Y(:,70)*DTS


    END SUBROUTINE DERIV

    SUBROUTINE WRITE(DTS)
        IMPLICIT NONE
        REAL :: DTS
        ! INTEGER :: TS, NCELL, PP, PC, TC, S

        PI = 4.0*ATAN(1.0)
        DO S = 1,NS
            Y_PPB(:,S) = Y(:,S) * 1.0E+9 / M(:) 
        END DO

        DO S = 1,NEMI
            EMI_PPB(:,S) = EMI(:,S) * 1.0E+9 / M(:)
        END DO

        IERR = NF90_PUT_VAR(NCID, VARID_LEVEL, LEVEL(:), START=[1], COUNT=[NCELL])
        CALL CHECK(IERR, 'LEVEL WRITE FAILED')
        IERR = NF90_PUT_VAR(NCID, VARID_LON, LON(:), START=[1], COUNT=[NCELL])
        CALL CHECK(IERR, 'LON WRITE FAILED')
        IERR = NF90_PUT_VAR(NCID, VARID_LAT, LAT(:), START=[1], COUNT=[NCELL])
        CALL CHECK(IERR, 'LAT WRITE FAILED')

        IERR = NF90_PUT_VAR(NCID, VARID_Y, Y_PPB(:, :), START=[1, 1, TS], COUNT=[NCELL, NS, 1])
        CALL CHECK(IERR, 'Y WRITE FAILED')
        IERR = NF90_PUT_VAR(NCID, VARID_SZA, SZA(:), START=[1, TS], COUNT=[NCELL, 1])
        CALL CHECK(IERR, 'SZA WRITE FAILED')
        IERR = NF90_PUT_VAR(NCID, VARID_J, J(:, :), START=[1, 1, TS], COUNT=[NCELL, 5, 1])
        CALL CHECK(IERR, 'J WRITE FAILED')
        IERR = NF90_PUT_VAR(NCID, VARID_DJ, DJ(:, :), START=[1, 1, TS], COUNT=[NCELL, 5, 1])
        CALL CHECK(IERR, 'DJ WRITE FAILED')
        IERR = NF90_PUT_VAR(NCID, VARID_RC, RC(:, :), START=[1, 1, TS], COUNT=[NCELL, 5, 1])
        CALL CHECK(IERR, 'RC WRITE FAILED')
        IERR = NF90_PUT_VAR(NCID, VARID_EMI, EMI_PPB(:, :), START=[1, 1, TS], COUNT=[NCELL, NEMI, 1])
        CALL CHECK(IERR, 'EMI WRITE FAILED')

        TIME1 = TIME1 + DTS

    END SUBROUTINE WRITE
    
    ! POST INTEGRATION
    SUBROUTINE DEALLOCATE
        IMPLICIT NONE
        
        IERR = NF90_CLOSE(NCID)
        CALL CHECK(IERR, 'NCID CLOSE FAILED')


    END SUBROUTINE DEALLOCATE

END MODULE BOXM

PROGRAM BOXM_RUN
    USE BOXM
    IMPLICIT NONE
    

    ! PRE-INTEGRATION
    CALL OPEN_NC()
    CALL GET_DIMS()
    CALL INIT_VARS()
    CALL GET_PRE_INT_VARIDS()
    CALL GET_INT_VARIDS()
    CALL GET_PRE_INT_VARS()

    ! ! INTEGRATION
    DTS = 20.0
    DO TS = 1, NTS
        IF (MOD(REAL(TS), 3600 / DTS) == 0) THEN
            PRINT *, 'TS = ', TS, ' (Hour = ', TS * DTS / 3600, ')'
        END IF
        CALL GET_INT_VARS(TS)
        CALL CALC_AEROSOL()
        CALL CHEMCO()
        CALL CALC_J()
        CALL PHOTOL()
        IF (TS /= 1) THEN
            CALL DERIV(DTS)
        END IF
        CALL WRITE(DTS)
    END DO

    ! CALL DEALLOCATE()

END PROGRAM BOXM_RUN
