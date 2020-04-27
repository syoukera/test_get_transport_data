!     general parameter indempendent from each cells 
module chemkin_params
    ! unit nuber for input
    integer, parameter :: unit_tplink = 20
    integer, parameter :: unit_cklink = 21

    ! chemkin work array
    integer,      allocatable :: int_ckwk(:)
    real(8),      allocatable :: real_ckwk(:)
    character(6), allocatable :: char_ckwk(:)*16
    
    ! transport work array
    integer,      allocatable :: int_tpwk(:)
    real(8),      allocatable :: real_tpwk(:)
end module chemkin_params

!   ----------------------------------------

subroutine initialize_chemkin_workarray()
    use chemkin_params

    integer, parameter :: unit_stdout = 6

    ! length of work array
    integer len_logi_ckwk
    integer len_int_ckwk
    integer len_real_ckwk
    integer len_char_ckwk

    ! open input unit
    open(unit_cklink, form='unformatted', file='./link/cklink')
    open(unit_tplink, form='unformatted', file='./link/tplink')

    !   ------- initialize chemkin work array ---------

    CALL cklen(unit_cklink, unit_stdout, len_int_ckwk, len_real_ckwk, len_char_ckwk)

    allocate(int_ckwk(len_int_ckwk))
    allocate(real_ckwk(len_real_ckwk))
    allocate(char_ckwk(len_char_ckwk))

    call ckinit(len_int_ckwk, len_real_ckwk, len_char_ckwk, unit_cklink, &
                unit_stdout, int_ckwk, real_ckwk, char_ckwk)

    !   ------- initialize transport work array ---------

    call mclen(unit_tplink, unit_stdout, len_int_tpwk, len_real_tpwk)
    
    allocate(int_tpwk(len_int_tpwk))
    allocate(real_tpwk(len_real_tpwk))

    call mcinit(unit_tplink, unit_stdout, len_int_tpwk, len_real_tpwk, &
                int_tpwk, real_tpwk)

end subroutine initialize_chemkin_workarray