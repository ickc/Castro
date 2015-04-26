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

end module align_array_module
