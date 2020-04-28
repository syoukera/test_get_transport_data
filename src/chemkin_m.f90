module chemkin_params
    ! unit nuber for input
    integer, parameter :: unit_tplink = 20
    integer, parameter :: unit_cklink = 21
    integer, parameter :: unit_stdout = 6

    ! lenth of work array obtained from cklen
    integer len_logi_cklen
    integer len_int_cklen
    integer len_real_cklen
    integer len_char_cklen
    
    ! length of work array
    integer len_logi_ckwk
    integer len_int_ckwk
    integer len_real_ckwk
    integer len_char_ckwk

    ! chemkin work array
    integer,      allocatable :: int_ckwk(:)
    real(8),      allocatable :: real_ckwk(:)
    character(6), allocatable :: char_ckwk(:)*16
    
    ! transport work array
    integer,      allocatable :: int_tpwk(:)
    real(8),      allocatable :: real_tpwk(:)

    ! index of chemistry
    integer mm ! number of elements
    integer kk ! number of species
    integer ii ! number of reactions
    integer nift ! number of fitting coefficients
    
    ! array pointers
    integer nsys, neq, lsdas, lidas, lrdas
    integer iprd, iprck, ipwt, ipwdot, ipu, iptot
    integer nrpar, nrdas, nsdas, natol, nrtol, nxmol, nz, nzp, lrtot
    integer nipar, nidas, litot, ipcck, nksym, lctot

    ! constant physical values
    real(8) pressure              ! Dyne/cm**2
    real(8), allocatable :: wt(:) ! Molecular weights gm/mole

    contains

    subroutine initialize_chemkin_workarray(pressure_Pa)
        real(8), intent(in) :: pressure_Pa
        ! to allign with SENKIN
        integer, parameter  :: ipick       = 20

        ! open input unit
        open(unit_cklink, form='unformatted', file='./link/cklink')
        open(unit_tplink, form='unformatted', file='./link/tplink')

        !   ------- get chemkin index ---------

        call cklen(unit_cklink, unit_stdout, len_int_cklen, len_real_cklen, len_char_cklen)
        call get_ckindex()

        !   ------- calculate pointer index ---------

        call calc_pointer()

        !   ------- initialize chemkin work array ---------

        len_int_ckwk  = 7000
        len_real_ckwk = 10000
        len_char_ckwk = 200

        allocate(int_ckwk(len_int_ckwk+ipick))
        allocate(real_ckwk(len_real_ckwk))
        allocate(char_ckwk(len_char_ckwk))

        call ckinit(len_int_ckwk, len_real_ckwk, len_char_ckwk, unit_cklink, &
                    unit_stdout, int_ckwk(ipick), real_ckwk, char_ckwk)

        !   ------- initialize transport work array ---------

        call mclen(unit_tplink, unit_stdout, len_int_tpwk, len_real_tpwk)
        
        allocate(int_tpwk(len_int_tpwk))
        allocate(real_tpwk(len_real_tpwk))

        call mcinit(unit_tplink, unit_stdout, len_int_tpwk, len_real_tpwk, &
                      int_tpwk, real_tpwk)

        !   ------- initialize pointer integer ---------

        call ckindx(int_ckwk, real_ckwk, mm, kk, ii, nfit)
    
        !   ------- initialize pointer integer ---------
        pressure = pressure_Pa*10d0
        allocate(wt(kk))
        call ckwt(int_ckwk, real_ckwk, wt)

    end subroutine initialize_chemkin_workarray

    
    subroutine get_ckindex()

        ! temporary work array
        integer      :: int_ckwk_temp(len_int_cklen)
        real(8)      :: real_ckwk_temp(len_real_cklen)
        character(6) :: char_ckwk_temp(len_char_cklen)*16

        call ckinit(len_int_cklen, len_real_cklen, len_char_cklen, unit_cklink, &
                    unit_stdout, int_ckwk_temp, real_ckwk_temp, char_ckwk_temp)
        call ckindx(int_ckwk_temp, real_ckwk_temp, mm, kk, ii, nfit)

    end subroutine get_ckindex

    subroutine calc_pointer()

        NSYS  = KK + 1
        NEQ = NSYS
        LSDAS  = 1
        LIDAS  = 20 + NEQ
        LRDAS  = 40 + NEQ*9 + NSYS*NSYS
! C  
! C           SET RPAR POINTERS
! C       NOTE: MUST HAVE IPRD 1ST IN RPAR
! C  

        IPRD   = 1
        IPRCK  = IPRD  + II
        IPWT   = IPRCK + LENR
        IPWDOT = IPWT  + KK
        IPU    = IPWDOT+ KK
        IPTOT  = IPU   + KK
! C  
! C           APPORTION REAL WORK SPACE
! C  
        NRPAR  = 1
        NRDAS  = NRPAR + IPTOT
        NSDAS  = NRDAS + LRDAS
        NATOL  = NSDAS + LSDAS
        NRTOL  = NATOL + NEQ
        NXMOL  = NRTOL + NEQ
        NZ     = NXMOL + KK
        NZP    = NZ    + NEQ
        LRTOT  = NZP   + NEQ
! C  
! C          APPORTION INTEGER WORK SPACE
! C  
        NIPAR  = 1
        NIDAS  = NIPAR + LENI + IPICK
        LITOT  = NIDAS + LIDAS
! C  
! C          APPORTION CHARACTER WORK SPACE
! C  
        IPCCK = 1
        NKSYM = IPCCK + LENC
        LCTOT = NKSYM + KK

    end subroutine calc_pointer


end module chemkin_params