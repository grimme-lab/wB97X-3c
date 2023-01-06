program main
   use mctc_io
   use mctc_env
   use ioroutines, only: rdfile
   use miscellaneous, only: helpf
   implicit none
   integer              :: narg,i,myunit
   integer              :: mpi,defgrid,charge,nopen

   character(len=80)    :: atmp,guess,filen,outn

   real(wp)             :: coremem

   logical              :: indguess, polar, beta, polgrad, dipgrad, geoopt, nocosx
   logical              :: tightscf, strongscf, verbose, suborca, nouseshark, sauxbas
   logical              :: largeaux, ploteldens, help, uhfgiven, da

   type(structure_type) :: mol
   type(error_type), allocatable   :: error

   filen = 'coord' ! input  filename
   outn  = 'wb97x3c.inp'   ! output filename
   mpi   =  4      ! #procs
   coremem   = 5000      ! #procs
   polar = .false. ! polarizability calc
   beta = .false. ! hyperpolarizabilities
   polgrad = .false. ! polarizability derivatives
   dipgrad = .false. ! dipole moment gradients
   geoopt = .false. ! turn on !OPT keyword
   nocosx = .false. ! turns on RIJCOSX, seminumerical exchange
   verbose = .false. ! verbose output
   nouseshark = .false. ! use different integral library
   sauxbas = .false. ! use small auxbasis from ~/.auxbasis_vDZP
   tightscf = .false. ! SCF conv criterium
   strongscf = .false. ! SCF conv criterium
   indguess = .false. ! SCF conv criterium
   uhfgiven = .false. ! SCF conv criterium
   help = .false. ! SCF conv criterium
   largeaux = .false.
   ploteldens = .false.

   ! get number of arguments
   narg = command_argument_count()

   do i=1,narg
      call get_command_argument(i,atmp)
      if(index(atmp,'-mpi').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) mpi
      endif
      call get_command_argument(i,atmp)
      if(index(atmp,'-struc').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) filen
         filen=trim(filen)
      endif
      if(index(atmp,'-memory').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) coremem
      endif
      if(index(atmp,'-defgrid').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) defgrid
      endif
      if(index(atmp,'-guess').ne.0) then
         indguess=.true.
         call get_command_argument(i+1,atmp)
         guess=trim(atmp)
      endif
      if(index(atmp,'-polar').ne.0) polar=.true.
      if(index(atmp,'-hyppol').ne.0) beta=.true.
      if(index(atmp,'-polgrad').ne.0) polgrad=.true.
      if(index(atmp,'-dipgrad').ne.0) dipgrad=.true.
      if(index(atmp,'-geoopt').ne.0) geoopt=.true.
      if(index(atmp,'-nocosx').ne.0) nocosx=.true.
      if(index(atmp,'-tightscf').ne.0) tightscf=.true.
      if(index(atmp,'-strongscf').ne.0) strongscf=.true.
      if(index(atmp,'-v').ne.0) verbose=.true.
      if(index(atmp,'-suborca').ne.0) suborca=.true.
      if(index(atmp,'-nouseshark').ne.0) nouseshark=.true.
      if(index(atmp,'-smallauxbasis').ne.0) sauxbas=.true.
      if(index(atmp,'-help').ne.0) help=.true.
      if(index(atmp,'-largeaux').ne.0) largeaux=.true.
      if(index(atmp,'-plot').ne.0) ploteldens=.true.
   enddo
   inquire(file='.UHF',exist=da)
   if(da)then
      open(newunit=myunit,file='.UHF')
      read(21,*) nopen
      uhfgiven=.true.
   endif

   inquire(file='.CHRG',exist=da)
   if(da)then
      open(newunit=myunit,file='.CHRG')
      read(21,*) charge
   endif

   if (help) then
      call helpf()
      stop
   endif

   call rdfile(filen,mol)

end program main
