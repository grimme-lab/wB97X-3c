      implicit none
      integer i,j,k,l,nat,nel,charge,zmax,nopen,nn,line,mpi,jat,ecpp
      integer elem(100)           
      integer iat(10000)          
      real*8  xyz(3,10000)
      real*8 bohr
      real*8 xx(10)
      character*1 ang(5)                           
      DATA ang/'S','P','D','F','G'/
      
      character*1 ang2(6)                           
      DATA ang2/'s','p','d','f','g','h'/
      
      character*2 asym                             
      CHARACTER*2 ELEMNT(107)
      CHARACTER*80 filen, outn 
      CHARACTER*128 atmp,guess        
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/
      logical da,polar,beta,dipgrad,polgrad,verbose,geoopt,nocosx
      logical tightscf,strongscf,indguess,help,uhfgiven,suborca
      logical nouseshark,sauxbas, largeaux

      integer lao, nao, npr, carg,defgrid,coremem
      integer myunit
      integer npp, n_core, lmax, lpp
      integer nao_aux(94)
      real*8  ec, ec2
      real*8  auxec(3,6,10,94)
      integer lao_aux(10,94),npr_aux(10,94)
      common /basis/ ec(3,6,10,94), lao(10,94), npr(10,94), nao(94) 
      common /basis/ lmax(94), n_core(94), lpp(100,94), ec2(3,6,10,94)
      common /basis/ npp(94) 

      charge=0
      nopen =0
      bohr=0.529177210903
      defgrid=2

      filen ='coord' ! input  filename
      outn  ='inp'   ! output filename
      mpi   = 4      ! #procs
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

      inquire(file='.UHF',exist=da)
      if(da)then
        open(unit=21,file='.UHF')
        read(21,*) nopen
        uhfgiven=.true.
      endif

      inquire(file='.CHRG',exist=da)
      if(da)then
        open(unit=21,file='.CHRG')
        read(21,*) charge
      endif

      guess="hueckel"
      carg = command_argument_count()
      do i=1,carg
          call getarg(i,atmp)
          if(index(atmp,'-mpi').ne.0) then
              call getarg(i+1,atmp)
              read(atmp,*) mpi
          endif
          if(index(atmp,'-memory').ne.0) then
              call getarg(i+1,atmp)
              read(atmp,*) coremem
          endif
          if(index(atmp,'-defgrid').ne.0) then
              call getarg(i+1,atmp)
              read(atmp,*) defgrid
          endif
          if(index(atmp,'-guess').ne.0) then
              indguess=.true.
              call getarg(i+1,atmp)
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
      enddo

      if (help) then
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
          write(*,*) "-smallauxbasis # (use auxbasis from 
     . ~/.auxbasis_vDZP instead of def2/J)"
          stop
      endif
      if (suborca) outn='inp.inp'   ! output filename

c read coords
      inquire(file=filen,exist=da)
      if(.not.da) stop 'coord file not found'
      open(unit=142,file=filen)
      if(index(filen,'.xyz').ne.0)then
      call rd0xyz(filen,nat)
      call rdxyz(filen,nat,xyz,iat)
      else
      call rd0(filen,nat)
      call rd(filen,nat,xyz,iat)
      endif

      elem=0
      nel=0
      j=0
      do i=1,nat
         nel=nel+iat(i)
         elem(iat(i))=elem(iat(i))+1
         if(iat(i).gt.j) j=iat(i)
      enddo
      zmax=j

      nel=nel-charge

      if(.not.uhfgiven.and.(nel.eq.1.or.mod(nel,2).ne.0))nopen=1

c start writing
      open(unit=7,file=outn)
      if (largeaux) then
          write(7,'(''! RKS WB97X-D4 def2/J def2-TZVP'')')
      else
          write(7,'(''! RKS WB97X-D4 def2/J'')')
      endif
      if (verbose) then
          write(7,'(''! PRINTBASIS LARGEPRINT'')')
      endif

      if (tightscf) then
          write(7,'(''! TightSCF'',2x,a,i1)') "DEFGRID", defgrid
      elseif (strongscf) then
          write(7,'(''! StrongSCF'',2x,a,i1)') "DEFGRID", defgrid
      else
          write(7,'(''! NormalSCF'',2x,a,i1)') "DEFGRID", defgrid
      endif
      if(geoopt) write(7,'(''! Opt'')')
      if(nocosx) write(7,'(''! NOCOSX'')')
      if(dipgrad) write(7,'(''! Freq'')')
      if(polgrad) write(7,'(''! NumFreq'')')
      if(nouseshark) write(7,'(''! NoUseShark'')')
      write(7,'( )') 

      if(mpi.gt.0) write(7,'(''%pal'',/,''nprocs'',i4,/,''end'',/)') mpi
      !if(mpi.gt.0) write(7,'(''%pal'',/,''nprocs '',i2,/,''end'',/)') mpi
      write(7,'(''%MaxCore '',i6,/)') coremem
c method input
!    .''! rks wb97x-d3 SDD Printbasis Largeprint nososcf slowconv'')') 
    
      write(7,'(a)')    "%method"
      write(7,'(a)')    "  D4A1    0.2464"
      write(7,'(a)')    "  D4A2    4.737"
      write(7,'(a)')    "  D4S6    1.00"
      write(7,'(a)')    "  D4S8    0.00"
      write(7,'(a)')    "  D4S9    1.00"
      write(7,'(a,/)')  "end"

      if(.not.indguess) then
          do i=37,45
              if (elem(i).ge.1) then
                  indguess=.true.
                  guess='hueckel'
              endif
          enddo
          do i=48,54
              if (elem(i).ge.1) then
                  indguess=.true.
                  guess='hueckel'
              endif
          enddo
          if (elem(21).ge.1.) then
              indguess=.true.
              guess='hcore'
          endif
          if (elem(47).ge.1) then
              indguess=.true.
              guess='hcore'
          endif
          if (elem(74).ge.1) then
              indguess=.true.
              guess='hcore'
          endif
          if (elem(82).ge.1) then
              indguess=.true.
              guess='hcore'
          endif
          if (elem(83).ge.1) then
              indguess=.true.
              guess='hcore'
          endif
      endif
      if(indguess)then
          write(7,'(''%scf'')')
          write(7,'(''  guess '',a20)') guess
          write(7,'(''end'',/)')
      endif

      if (beta) then
          write(7, '(a)') "%scf"
          write(7, '(a)') "efield XXX,YYY,ZZZ"
          write(7, '(a,/)') "end"
      endif
      if(polar.or.polgrad.or.beta) 
     & write(7,'(''%elprop polar 1'',/,''end'',/)')
      if (verbose) then
          write(7,'(''%output'')')
          write(7,'(''       print[P_Hirshfeld] 1'')')
          write(7,'(''       print[P_BondOrder_M] 1'')')
          write(7,'(''       print[P_basis] 2'')')
          write(7,'(''end'')')
      endif
c read file ~/.basis and convert to ORCA style
      call rdbas
      write(7,'(''%basis'')')
      do jat=1,94
         if(elem(jat).gt.0.and.nao(jat).gt.0)then
         write(7,'(''  NewGTO '',a2)') asym(jat)
         do i=1,nao(jat)           
            write(7,'(a2,2x,i3)') ang(lao(i,jat)),npr(i,jat)
            do j=1,npr(i,jat)
               write(7,'(i3,2x,2f14.8)') j,ec(1:2,j,i,jat)
            enddo
         enddo
         write(7,'(''  end'')')
         endif
      enddo
      if (sauxbas) then
      call rdauxbas(nao_aux,auxec,npr_aux,lao_aux)
      do jat=1,94
         if(elem(jat).gt.0.and.nao_aux(jat).gt.0)then
         write(7,'(''  NewAuxJGTO '',a2)') asym(jat)
         do i=1,nao_aux(jat)           
            write(7,'(a2,2x,i3)') ang(lao_aux(i,jat)),npr_aux(i,jat)
            do j=1,npr_aux(i,jat)
               write(7,'(i3,2x,2f14.8)') j,auxec(1:2,j,i,jat)
            enddo
         enddo
         write(7,'(''  end'')')
         endif
      enddo
      endif

c read file ~/.ecp and convert to ORCA style
      call rdecp
      do jat=1,94
         if(elem(jat).gt.0.and.npp(jat).gt.0)then
         write(7,'(''  NewECP '',a2)') asym(jat)
         write(7,'(''  N_core '',i2)') n_core(jat)
         write(7,'(''  lmax '',a2)') ang2(lmax(jat)+1)
         do i=1,lpp(1,jat)           
            write(7,'(a2,2x,i3)') ang2(i),npr(i,jat)
            do j=1,npr(i,jat)
               ecpp=idint(ec2(3,i,j,jat))
               write(7,'(i3,2x,2f14.8,2x,i3)') j,ec2(1:2,i,j,jat),ecpp 
            enddo
         enddo
         write(7,'(''  end'')')
         endif
      enddo
      write(7,'(''end'')')

      write(7,'(''* xyz'',2i4)')charge,nopen+1
      do i=1,nat
         write(7,'(a2,2x,3F22.14)') elemnt(iat(i)),xyz(1:3,i)
      enddo
      write(7,'(''*'')')

      close(7)

      end
c*************************************************************************
      subroutine rdecp
      implicit none
      character*80 atmp,ecpfilename,homedir
      integer nn,i,iat,np,l
      real*8 xx(10)
      logical da, da2
      integer lao, nao, npr
      real*8  ec, ecp, ec2
      character null0
      common /basis/ ec(3,6,10,94), lao(10,94), npr(10,94), nao(94)
      integer npp, n_core, lmax, lpp
      common /basis/  lmax(94), n_core(94), lpp(100,94), ec2(3,6,10,94)
      common /basis/ npp(94) 
      character*1 ang2(6)                           
      DATA ang2/'s','p','d','f','g','h'/
      


      lpp = 0
      npr = 0
      npp = 0

      CALL get_environment_variable("HOME", homedir)
      ecpfilename=trim(homedir)//'/.ecp'

      inquire(file='~/.ecp',exist=da)
      if(da)then
          open(unit=45,file='~/.ecp')
      else
          inquire(file=ecpfilename,exist=da2)
          if(da2)then
              open(unit=45,file=ecpfilename)
          else
              write(*,*) "~/.ecp file in $HOME is missing!"
              error stop "I/O error stop."
          endif
      endif
 10   read(45,'(a)',end=20) atmp
      i=0
      if(index(atmp,'*').ne.0) then        
         read(45,*) iat             
         read(45,*) null0, null0, n_core(iat), null0, null0, lmax(iat)
 12      read(45,'(a)',end=20) atmp
!          no space in first column allowed in ecp file
         if(index(atmp,'s').eq.1) l=1
         if(index(atmp,'p').eq.1) l=2
         if(index(atmp,'d').eq.1) l=3
         if(index(atmp,'f').eq.1) l=4
         if(index(atmp,'g').eq.1) l=5
         if(index(atmp,'h').eq.1) l=6  
         if(index(atmp,'*').ne.0) then
            goto 10     
         endif
         if(index(atmp,ang2(lmax(iat)+1)).eq.0) then
         backspace 45
         lpp(npp(iat)+1,iat)=l
         npp(iat)=npp(iat)+1
         npr(l,iat)=npr(l,iat) + 1
 
         read(45,*) ec2(2,l,npr(l,iat),iat)
     1   ,ec2(3,l,npr(l,iat),iat),ec2(1,l,npr(l,iat),iat)
         ENDIF
         GOTO 12
      ENDIF

      goto 10
 20   close(45)
      
      return
      end



C     *****************************************************************         

      subroutine rdbas
      implicit none
      character*80 atmp,homedir,basisfilename
      integer nn,i,iat,np,l
      real*8 xx(10)
      logical da, da2
      integer lao, nao, npr
      real*8  ec, ec2
      common /basis/ ec(3,6,10,94), lao(10,94), npr(10,94), nao(94) 
      integer npp, n_core, lmax, lpp
      common /basis/ lmax(94), n_core(94), lpp(10,94), ec2(3,6,10,94)
      common /basis/ npp(94) 



      lao = 0
      nao = 0
      npr = 0
      ec  = 0

      CALL get_environment_variable("HOME", homedir)
      basisfilename=trim(homedir)//'/.basis_vDZP'

      inquire(file='~/.basis_vDZP',exist=da)
      if(da)then
          open(unit=44,file='~/.basis_vDZP')
      else
          inquire(file=basisfilename,exist=da2)
          if(da2)then
              open(unit=44,file=basisfilename)
          else
              write(*,*) "~/.basis_vDZP file in $HOME is missing!"
              error stop "I/O error stop."
          endif
      endif
 10   read(44,'(a)',end=20) atmp
      if(index(atmp,'*').ne.0) then        
         read(44,*) iat             
 12      read(44,'(a)',end=20) atmp
         if(index(atmp,'*').ne.0) goto 10     
         call readl(atmp,xx,nn)
         if(nn.eq.1)then
            np = idint(xx(1))
            if(index(atmp,'s').ne.0) l=1         
            if(index(atmp,'p').ne.0) l=2         
            if(index(atmp,'d').ne.0) l=3         
            if(index(atmp,'f').ne.0) l=4         
            nao(iat)=nao(iat)+1
            lao(nao(iat),iat)=l           
            npr(nao(iat),iat)=np          
            do i=1,np
               read(44,*) ec(1,i,nao(iat),iat),ec(2,i,nao(iat),iat)
            enddo
         endif
         goto 12
      endif
      goto 10
 20   close(44)
      
      return
      end

C     *****************************************************************         

      subroutine rdauxbas(nao_aux,auxec,npr_aux,lao_aux)
      implicit none
      character*80 atmp,homedir,auxbasisfilename
      integer nn,i,iat,np,l
      real*8 xx(10)
      logical da, da2
      integer myunit
      integer lao_aux(10,94),npr_aux(10,94)
      integer nao_aux(94) 
      real*8  auxec(3,6,10,94)


      lao_aux = 0
      nao_aux = 0
      npr_aux = 0
      auxec  = 0

      CALL get_environment_variable("HOME", homedir)
      auxbasisfilename=trim(homedir)//'/.auxbasis_vDZP'

      inquire(file='~/.auxbasis_vDZP',exist=da)
      if(da)then
          open(newunit=myunit,file='~/.auxbasis_vDZP')
      else
          inquire(file=auxbasisfilename,exist=da2)
          if(da2)then
              open(newunit=myunit,file=auxbasisfilename)
          endif
      endif
100   read(myunit,'(a)',end=202) atmp
      if(index(atmp,'*').ne.0) then        
         read(myunit,*) iat             
120      read(myunit,'(a)',end=202) atmp
         if(index(atmp,'*').ne.0) goto 100
         call readl(atmp,xx,nn)
         if(nn.eq.1)then
         np = idint(xx(1))
         if(index(atmp,'s').ne.0) l=1         
         if(index(atmp,'p').ne.0) l=2         
         if(index(atmp,'d').ne.0) l=3         
         if(index(atmp,'f').ne.0) l=4         
         nao_aux(iat)=nao_aux(iat)+1
         lao_aux(nao_aux(iat),iat)=l           
         npr_aux(nao_aux(iat),iat)=np          
         do i=1,np
         read(myunit,*) auxec(1,i,nao_aux(iat),iat),
     1 auxec(2,i,nao_aux(iat),iat)
         enddo
         endif
         goto 120
      endif
      goto 100
202   close(myunit)
      
      return
      end


C     *****************************************************************         

      SUBROUTINE backstring(A1,A2,lena2)
      CHARACTER*(*) A1
      CHARACTER*(*) A2
      integer n,lena2
      n=0
      DO J=1,len(a1)
         if(a1(j:j).ne.' ')then
            n=n+1
            a2(n:n)=a1(j:j)
         endif
      enddo
      DO J=1,lena2   
         a2(j:j)=' '
      enddo
      a1=a2
      a2='                                                            '
      n=0
      DO J=1,len(a1)
         if(a1(j:j).ne.' ')then
            n=n+1
            a2(n:n)=a1(j:j)
         endif
      enddo

      END


C     *****************************************************************         
                                                                                
      SUBROUTINE READL(A1,X,N)                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      CHARACTER*(*) A1                                                      
      DIMENSION X(*)                                                            
      I=0                                                                       
      IS=1                                                                      
  10  I=I+1                                                                     
      X(I)=READAA(A1,IS,IB,IE)                                               
      IF(IB.GT.0 .AND. IE.GT.0) THEN                                            
                                IS=IE                                           
                                GOTO 10                                         
      ENDIF                                                                     
      N=I-1                                                                     
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      FUNCTION READAA(A,ISTART,IEND,IEND2)                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 READAA                                                             
      CHARACTER*(*) A                                                      
      NINE=ICHAR('9')                                                           
      IZERO=ICHAR('0')                                                          
      MINUS=ICHAR('-')                                                          
      IDOT=ICHAR('.')                                                           
      ND=ICHAR('D')                                                             
      NE=ICHAR('E')                                                             
      IBL=ICHAR(' ')                                                            
      IEND=0                                                                    
      IEND2=0                                                                   
      IDIG=0                                                                    
      C1=0                                                                      
      C2=0                                                                      
      ONE=1.D0                                                                  
      X = 1.D0                                                                  
      NL=LEN(A) 
      DO 10 J=ISTART,NL-1                                                       
         N=ICHAR(A(J:J))                                                          
         M=ICHAR(A(J+1:J+1)) 
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20                      
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO                            
     1 .OR. M.EQ.IDOT)) GOTO 20                                                 
   10 CONTINUE                                                                  
      READAA=0.D0                                                               
      RETURN                                                                    
   20 CONTINUE                                                                  
      IEND=J                                                                    
      DO 30 I=J,NL                                                              
         N=ICHAR(A(I:I))                                                          
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C1=C1*10+N-IZERO                                                    
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN                                     
            ONE=-1.D0                                                           
         ELSEIF(N.EQ.IDOT) THEN                                                 
            GOTO 40                                                             
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   30 CONTINUE                                                                  
   40 CONTINUE                                                                  
      IDIG=0                                                                    
      DO 50 II=I+1,NL                                                           
         N=ICHAR(A(II:II))                                                         
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C2=C2*10+N-IZERO                                                    
            X = X /10                                                           
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN                                    
            X=-X                                                                
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   50 CONTINUE                                                                  
C                                                                               
C PUT THE PIECES TOGETHER                                                       
C                                                                               
   60 CONTINUE                                                                  
      READAA= ONE * ( C1 + C2 * X)                                              
      DO 55 J=IEND,NL                                                           
         N=ICHAR(A(J:J))                                                          
         IEND2=J                                                                
         IF(N.EQ.IBL)RETURN                                                     
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57                                           
      RETURN                                                                    
                                                                                
   57 C1=0.0D0                                                                  
      ONE=1.0D0                                                                 
      DO 31 I=J+1,NL                                                            
         N=ICHAR(A(I:I))                                                          
         IEND2=I                                                                
         IF(N.EQ.IBL)GOTO 70                                                    
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO                      
         IF(N.EQ.MINUS)ONE=-1.0D0                                               
   31 CONTINUE                                                                  
   61 CONTINUE                                                                  
   70 READAA=READAA*10**(ONE*C1)                                                
      RETURN                                                                    
      END                                                                       

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read coordinates

      subroutine rdxyz(fname,n,xyz,iat)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,n),iat(n),xx(10)
      character*128 line
      character*(*) fname

      ich=142
      open(unit=ich,file=fname)
      read(ich,'(a)')line
      read(ich,'(a)')line
      do i=1,n
         read(ich,'(a)')line
         call READL(line,xx,nn)    
         if(nn.ne.3) stop 'read error'
         xyz(1:3,i)=xx(1:3)
         call elem(line,j)
         iat(i)=j
      enddo

      close(ich)
      end

      subroutine rd(fname,n,xyz,iat)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,n),iat(n),xx(10)
      character*128 line
      character*(*) fname
      bohr=0.52917726

      ich=142
      open(unit=ich,file=fname)
      read(ich,'(a)')line
      do i=1,n
         read(ich,*) xyz(1:3,i)
      enddo

      rewind ich
      read(ich,'(a)')line
      do i=1,n
         read(ich,'(a)')line
         call elem(line,j)
         iat(i)=j
      enddo

      xyz=xyz*bohr
      close(ich)
      end

      subroutine rd0xyz(fname,n)
      character*(*) fname
      open(unit=142,file=fname)
      read(142,*) n
      close(142)
      end

      subroutine rd0(fname,n)
      implicit real*8 (a-h,o-z)
      dimension xx(10)
      character*128 line
      character*(*) fname

      ich=142
      open(unit=ich,file=fname)
      n=0
 100  read(ich,'(a)',end=200)line
         if(index(line,'$set').ne.0)goto 200
         if(index(line,'$redu').ne.0)goto 200
         if(index(line,'$user').ne.0)goto 200
         call readl(line,xx,nn)
         if(nn.ne.3) goto 100
         n=n+1
      goto 100
 200  continue

      close(ich)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ELEM(KEY1, NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) KEY1
      CHARACTER*2 ELEMNT(107),E

      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/

      nat=0
      e='  '
      k=1
      DO J=1,len(key1)
         if (k.gt.2)exit
         N=ICHAR(key1(J:J))
         if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
            e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
            k=k+1
         endif
         if(n.ge.ichar('a') .and. n.le.ichar('z') )then
            e(k:k)=key1(j:j)
            k=k+1
         endif
      enddo

      DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

C     *****************************************************************         

      FUNCTION ASYM(I)
      CHARACTER*2 ASYM
      CHARACTER*2 ELEMNT(107), AS
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/
      AS=ELEMNT(I)
      CALL UPPER(AS)
      ASYM=AS
      if(i.eq.103) asym='XX'
      RETURN
      END

      SUBROUTINE UPPER(AS)
      CHARACTER*2 AS
      NSP=ICHAR(' ')
      ND=ICHAR('A')-ICHAR('a')
      DO 10 I=1,2
         J=ICHAR(AS(I:I))
         IF(J.NE.NSP)J=J+ND
         AS(I:I)=CHAR(J)
  10  CONTINUE
      RETURN
      END

  
