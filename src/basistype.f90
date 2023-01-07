module basistype
   use mctc_io
   use mctc_env
   implicit none
   private

   public :: basis_type

   type :: basis_type
      integer               :: atmax ! atom index
      integer, allocatable  :: npr(:,:) ! number of primitive functions per atom
      integer, allocatable  :: nbf(:) ! total number of primitive functions
      integer, allocatable  :: angmom(:,:) ! angular momentum of each primitive function
      real(wp), allocatable :: exp(:,:,:) ! exponent of each primitive function
      real(wp), allocatable :: coeff(:,:,:) ! contraction coefficient of each primitive function
   end type basis_type

end module basistype
