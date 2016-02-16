subroutine stockInitialState
  use parallel
  use hamiltonian

  implicit none

  call startCount("stockState")

  if (hami%flag_spectrum) then
     invariant_matrix_spectrum_initial=invariant_matrix_spectrum
  end if

  call stopCount
end subroutine stockInitialState
