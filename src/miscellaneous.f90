module miscellaneous
  implicit none
  private

  public :: helpf
contains
  subroutine helpf()
   write(*,*) "Possible commands:"
   write(*,*) "-mpi <int>"
   write(*,*) "-defgrid <int>"
   write(*,*) "-guess <guess options>"
   write(*,*) "-polar"
   write(*,*) "-hyppol"
   write(*,*) "-polgrad"
   write(*,*) "-dipgrad"
   write(*,*) "-geoopt"
   write(*,*) "-nocosx"
   write(*,*) "-tightscf"
   write(*,*) "-strongscf"
   write(*,*) "-v"
   write(*,*) "-suborca # (make 'inp.inp'-Output)"
   write(*,*) "-nouseshark # (use different integral library)"
   write(*,*) "-smallauxbasis # (use auxbasis from ~/.auxbasis_vDZP instead of def2/J)"
   write(*,*) "-plot # (Plot the electron density with the following settings)"
  end subroutine helpf
end module miscellaneous