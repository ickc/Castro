! problem-specific Fortran stuff goes here

subroutine problem_checkpoint(int_dir_name, len) bind(C, name="problem_checkpoint")

  ! called by the IO processor during checkpoint

  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo



end subroutine problem_checkpoint


subroutine problem_restart(int_dir_name, len) bind(C, name="problem_restart")

  ! called by ALL processors during restart 

  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

end subroutine problem_restart
