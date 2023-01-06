module ioroutines
   use mctc_io
   use mctc_env
   implicit none
   private

   public :: rdfile
contains
   subroutine rdfile(filename,imol)

      type(structure_type),intent(out)  :: imol
      type(error_type), allocatable     :: error
      character(len=*), intent(in)      :: filename

      call read_structure(imol, filename, error)
      if (allocated(error)) then
         print '(a)', error%message
         error stop
      end if

      write(*,*) imol%nat

      ! test input
      call write_structure(imol,"coord.rev",error,2)
      if (allocated(error)) then
         print '(a)', error%message
         error stop
      end if
      print *, "Hello, revo4wb3c!"
   end subroutine rdfile
end module ioroutines
