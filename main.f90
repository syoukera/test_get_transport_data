!     test code for obtaining chemkin transport data
!
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

program test_main
      integer,parameter :: num_species = 53
      
      ! phisycal values from CFD calculation
      real(8) :: T_CFD = 1000d0         ! K
      real(8) :: P_CFD = 1.01325d5*40   ! Pa
      real(8) Y_CFD(num_species)        ! Mass fractions
      real(8) :: delta_t_CFD = 1.0d-12  ! s
      ! real(8) :: TOLS_CFD(4)            ! Tolerances

      ! logical :: MAKE_OUTPUT = .false.

      ! transport data obtained from this program
      ! real(8) :: D_diff(num_species)
      ! real(8) lambda
      
      ! Assurme Y has a same secuence as species in ckout
      data Y_CFD /0.00d+00,0.00d+00,0.00d+00,2.20d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,5.51d-02,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00, &
                  0.00d+00,0.00d+00,7.24d-01,0.00d+00,0.00d+00,0.00d+00,0.00d+00,0.00d+00/

      call initialize_chemkin_library

      call get_tranport_data()
end program test_main

!   ----------------------------------------

subroutine initialize_chemkin_library()
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

end subroutine initialize_chemkin_library

!   ----------------------------------------

subroutine get_tranport_data()
      ! real(8), intent(out) :: D_diff
      ! real(8)

      ! call mcmdif(P, T, X, KDIM, IMCWRK, RMCWRK, D)
      ! call mcacon
end subroutine get_tranport_data