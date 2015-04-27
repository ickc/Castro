module align_array_module

  implicit none

#ifdef BL_ALIGN_64_BYTE
  integer, parameter :: alignbyte = 64
#elif BL_ALIGN_32_BYTE
  integer, parameter :: alignbyte = 32
#elif BL_ALIGN_16_BYTE
  integer, parameter :: alignbyte = 16
#else
  integer, parameter :: alignbyte = 8
#endif

  integer, parameter :: nalign_double = alignbyte/8

contains

  ! Add padding to the first dimension so that loops like below have good alignment
  ! for array a(loa:hia,:,...)
  !
  ! do j = ...
  !   do i = lo, hi
  !     a(i,j,...) = ...
  !
  subroutine align_padding(start, lo, hi, loa, hia)
    integer, intent(in) :: start, lo, hi
    integer, intent(out) :: loa, hia
    integer :: rem, nx
    if (start .eq. lo) then
       loa = lo
    else
       loa = min(lo, start-nalign_double)
    end if
    nx = hi - loa + 1  ! We want this to be a multiple of nalign_double
    rem = modulo(nx, nalign_double)
    if (rem .eq. 0) then
       hia = hi
    else
       hia = hi + (nalign_double-rem)
    end if
  end subroutine align_padding

  ! If a(lo:) is aligned, what's the first index after 'start' that is aligned?
  function align_next_index(lo, start) result (start_align)
    integer, intent(in) :: lo, start
    integer :: start_align
    integer :: rem
    rem = modulo(start-lo, nalign_double)
    if (rem .eq. 0) then
       start_align = start
    else
       start_align = start + (nalign_double - rem)
    end if
  end function align_next_index

end module align_array_module
