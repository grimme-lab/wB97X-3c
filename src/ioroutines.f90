module ioroutines
   use mctc_io
   use mctc_env
   use basistype, only: basis_type
   implicit none
   private

   public :: rdfile,rdbas
contains
   subroutine rdfile(filename,imol,nel)

      type(structure_type),intent(out)  :: imol
      type(error_type), allocatable     :: error
      character(len=*), intent(in)      :: filename
      integer,intent(out)               :: nel
      integer                           :: i

      call read_structure(imol, filename, error)
      if (allocated(error)) then
         print '(a)', error%message
         error stop
      end if

      nel = 0
      do i=1,imol%nat
         nel = nel + imol%num(imol%id(i))
      end do

   end subroutine rdfile

   subroutine rdbas(ibasis,verb,basisfilename)

      type(basis_type),intent(out)           :: ibasis
      logical, intent(in), optional          :: verb
      character(len=*), intent(in),optional  :: basisfilename
      type(error_type), allocatable          :: error

      character(len=:), allocatable          :: fname,homedir
      character(len=120)                     :: atmp
      character(len=1)                       :: ltmp

      logical                                :: da

      integer, allocatable                   :: nbf(:),npr(:,:)
      character(len=1), allocatable          :: angmom(:,:)

      integer                                :: myunit,char_length,ierr,iread
      integer                                :: iat,tmpnpr,imax
      integer                                :: i,j,l,k

      allocate(nbf(100),npr(100,20),angmom(100,20))

      if (.not. present(basisfilename)) then
         call get_environment_variable("HOME", length=char_length)
         allocate(character(len=char_length) :: homedir)
         CALL get_environment_variable("HOME", value=homedir, status=ierr)
         select case (ierr)
          case (0)
            ! do nothing
          case (1)
            call fatal_error(error, "HOME environment variable is not set!")
            if (allocated(error)) then
               print '(a)', error%message
               error stop "I/O error stop."
            end if
          case (2)
            call fatal_error(error, "this compiler does not support environment variables")
            if (allocated(error)) then
               print '(a)', error%message
               error stop "I/O error stop."
            end if
         end select
         fname=trim(homedir)//'/.basis_vDZP'
      else
         fname=basisfilename
      end if

      write(*,'(a,a)') "Used basis set file: ",fname

      inquire(file=fname,exist=da)
      if(da)then
         open(newunit=myunit,file=fname,action='read',status='old',iostat=ierr)
      else
         call fatal_error(error, "basis set file cannot be opened!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      endif
      iread =  0
      l     =  0
      npr   =  0
      nbf   =  0
      imax  =  0
      do while (iread >= 0)
         read(myunit,'(a)',iostat=iread) atmp
         l = l + 1
         if (iread < 0) then
            exit
         end if
         if (iread > 0) then
            write(*,*) "Current line number: ",l
            error stop "I/O error in basis set read in."
         end if
         if(index(atmp,'*').ne.0) then
            cycle
         else
            read(atmp,*,iostat=iread) iat
            if (iread /= 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in basis set read in."
            end if
            nbf(iat) = 0
            npr(iat,:) = 0

            if (iat > imax) then
               imax = iat
            end if

            read(myunit,'(a)',iostat=iread) atmp
            l = l + 1
            if (iread /= 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in basis set read in."
            end if
            k = 0
            do while (index(atmp,'*').eq.0)
               read(atmp,*) tmpnpr,ltmp
               k = k + 1
               nbf(iat) = k
               npr(iat,nbf(iat)) = tmpnpr
               angmom(iat,nbf(iat)) = ltmp
               do j=1,tmpnpr
                  read(myunit,*,iostat=iread) atmp
               end do
               read(myunit,'(a)',iostat=iread) atmp
               l = l + 1
               if (iread /= 0) then
                  write(*,*) "Current line number: ",l
                  error stop "I/O error in basis set read in."
               end if
            enddo
         endif
      enddo

      ibasis%atmax = imax
      allocate(ibasis%npr(imax,20),ibasis%angmom(imax,20),ibasis%nbf(imax))
      allocate(ibasis%exp(imax,20,20),ibasis%coeff(imax,20,20))
      close(myunit)

      ibasis%npr = 0
      ibasis%nbf = 0
      ibasis%angmom = 0
      ibasis%exp = 0.0_wp
      ibasis%coeff = 0.0_wp

      do i=1,ibasis%atmax
         ibasis%nbf(i) = nbf(i)
         do j=1,ibasis%nbf(i)
            ibasis%npr(i,j) = npr(i,j)
            if (angmom(i,j) == "s") ibasis%angmom(i,j) = 0
            if (angmom(i,j) == "p") ibasis%angmom(i,j) = 1
            if (angmom(i,j) == "d") ibasis%angmom(i,j) = 2
            if (angmom(i,j) == "f") ibasis%angmom(i,j) = 3
            if (angmom(i,j) == "g") ibasis%angmom(i,j) = 4
            if (angmom(i,j) == "h") ibasis%angmom(i,j) = 5
            if (angmom(i,j) == "i") ibasis%angmom(i,j) = 6
            if (angmom(i,j) == "j") error stop "too high primitives"
         ibasis%npr(iat,:) = npr(iat,:)
         enddo
      enddo


      open(newunit=myunit,file=fname,action='read',status='old',iostat=ierr)
      iread =  0
      l     =  0
      if (verb) then
         write(*,*) "Z, # basis function, # primitive, exponent, coefficient"
      endif
      do while (iread >= 0)
         read(myunit,'(a)',iostat=iread) atmp
         l = l + 1
         if (iread < 0) then
            write(*,'(a,1x,i3,/)') "End of file reached after reading basis set for element:", iat
            exit
         end if
         if (iread > 0) then
            write(*,*) "Current line number: ",l
            error stop "I/O error in basis set read in."
         end if
         if(index(atmp,'*').ne.0) then
            cycle
         else
            read(atmp,*,iostat=iread) iat
            if (iread /= 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in basis set read in."
            end if
            do i=1,ibasis%nbf(iat)
               read(myunit,*,iostat=iread) atmp
               l = l + 1
               do j=1,ibasis%npr(iat,i)
                  read(myunit,*,iostat=iread) ibasis%exp(iat,i,j),ibasis%coeff(iat,i,j)
                  if (verb) then
                     write(*,*) iat, i, j, ibasis%exp(iat,i,j),ibasis%coeff(iat,i,j)
                  end if
                  l = l + 1
               end do
            enddo
         endif
      enddo

   end subroutine rdbas
end module ioroutines
