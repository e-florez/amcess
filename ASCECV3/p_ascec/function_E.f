******************************************************************************
*                            Modificacion "Version 3"
*     Por: Andy Danian Zapata Escobar
*     04/01/2015
*
*     Para usar el nwchen y gamess, en el calculo de la energia del sistema
*     de estudio.
*     Cuando se usa el gamess solo se pude usar gaussianas sto o de pople,
*     con el nwchem cualquier base, que este implementada en el
*     programa nwchen.
*     Con el gamess se puede usar hf, semi-empiricos, dft y solo mp2; con
*     el nwchen solo se puede usar hf y dft.

      function funkGAU(natom,r,jo) !Call system
      implicit real*8 (a-h,o-z)
      include 'param.inc'     
      character*10 stat1
      character*20 names
	character*3 stat2   
      character*150 stat3 !DALTON
      character*2 sym, symbol, element 
      character*8 hamil
      character*15 alias	
      character*20 bde
      character*2  cc
      character*5  ngva_char
      character*1  acc,ccc,ccn,cccc
      character*5  sp 
      integer proc,j,thamil,tbase,ng,ngv,sistat1,tb,nli
      integer cout, k
      character*8 mspin
      common/elemente/element(max_atom),mspin,namolec(max_atom)
      common/names/iznu(max_atom),sym(120),xbox,nnu(max_atom) !iznu: tipo de nucleo, nnu: cantidad tn
      common/ruta/ia,iQ,iS2,hamil,alias !ia:1,2,3,4,5,6 (g09,MSINDO.1.1.opt,nwchem,gamess,adf,dal) 
      common/proce/proc,thamil,ielem    !thamil:0,1,2,3 (HF,PM3,DFT,MP2) ielem: # tipo átomos
      common/db/bde(max_atom),tbase(max_atom),iz(max_atom)
      common/molecu/rcm(max_mole,3),imolec(0:max_mole),nmo,ds,dphi
      character*1000 nproc, aliasg
      dimension r(max_atom,3) 
      data B2/3.166D-6/		!Const. Boltzmann en Hartree/K
CDALTON
      CHARACTER(len=*),PARAMETER :: FMT1="(A)"
      CHARACTER(len=*),PARAMETER :: FMT2="(I2)"
      CHARACTER(len=*),PARAMETER :: FMT3="(A,I2,A)"
      CHARACTER(len=*),PARAMETER :: FMT4="(A,I2,A,I2,A)"
      CHARACTER(len=*),PARAMETER :: FMT5="(A,I3,A,A)"
      CHARACTER(len=*),PARAMETER ::FMT6="(A,3X,F20.6,3X,F20.6,3X,F20.6)"
      CHARACTER(len=*),PARAMETER ::FMT7="(A,3X,E20.6,3X,E20.6,3X,E20.6)"
      CHARACTER(len=150) :: commanddalton
C======================================================================       
C Contrucci'on input del GAMESS y calculo
C
      IF(ia.eq.4)then !Call Gamess

       open(34,file='anneal.inp')
       vert=xbox/2.d0
       write(34,*)'$CONTRL RUNTYP=Energy $END'
       SELECT CASE (thamil)
         CASE (0)
          if((iS2.eq.1).AND.(iQ.eq.0))then
           write(34,*)'$CONTRL SCFTYP=RHF $END'
          else
           write(34,*)'$CONTRL SCFTYP=UHF $END'
          endif
         CASE (1)  
          if((iS2.eq.1).AND.(iQ.eq.0))then
            write(34,*)'$BASIS GBASIS=',hamil,' $END'
          else
            write(34,*)'$CONTRL SCFTYP=UHF $END'
            write(34,*)'$BASIS GBASIS=',hamil,' $END'
          endif
         CASE (2)
          if((iS2.eq.1).AND.(iQ.eq.0))then
           write(34,*)'$CONTRL SCFTYP=RHF $END'
           write(34,*)'$DFT  DFTTYP=',hamil,' $END'
          else
           write(34,*)'$CONTRL SCFTYP=UHF $END'
           write(34,*)'$DFT  DFTTYP=',hamil,' $END'
          endif
         CASE (3)
          if((iS2.eq.1).AND.(iQ.eq.0))then
           write(34,*)'$CONTRL SCFTYP=RHF $END'
           write(34,*)'$CONTRL MPLEVL=2 $END'
          else
           write(34,*)'$CONTRL SCFTYP=UHF $END'
           write(34,*)'$CONTRL MPLEVL=2 $END'
          endif
         CASE DEFAULT
         write(*,*)"EL TIPO HAMILTONIANO NO ESTA IMPLEMENTADO EN 
     $  EL ASCEC Y GAMESS"
        END SELECT

       if(( thamil .NE. 1).OR.( thamil .NE. 1))then
          
          call system ("awk 'NR==11{print $4}' *.in > tb")
          open(4,file='tb')
          read(4,*)tb
          close(4)
          call system ('rm -f tb')

        SELECT CASE (tb)

          CASE (0)
          call system ("awk 'NR==11{print $3}' *.in | cut -c5 > ng")
          open(4,file='ng')
          read(4,*)ng
          close(4)
          call system ('rm -f ng')
          write(34,*)'$BASIS GBASIS=sto NGAUSS=',ng,' $END'

          CASE (1)
          call system ("awk 'NR==11{print $3}' *.in| cut -c1 > ng")
          open(4,file='ng')
          read(4,*)ng
          close(4)
          call system ('rm -f ng')
          call system ("awk 'NR==11{print $3}' *.in | cut -c3-5 | 
     $    sed 's/g//' | sed 's/+//' | sed 's/*//' > ngv")
          open(4,file='ngv')
          read(4,*)ngv
          close(4)
          call system ('rm -f ngv')
          write(ngva_char,'(i4)')ngv
          write(34,*)'$BASIS GBASIS=N',adjustl(ngva_char),' NGAUSS=',ng
     $  ,' $END'
         call system ("awk 'NR==11{print $3}' *.in | sed 's/[0-9]//g' | 
     $  sed 's/[a-z]//g' | sed 's/-//' | sed 's/*//' | sed 's/*//' 
     $  | wc -c > fd")
         open(4,file='fd')
         read(4,*)fd
         close(4)
         call system ('rm -r fd')
         call system ("awk 'NR==11{print $3}' *.in | sed 's/[0-9]//g' | 
     $ sed 's/[a-z]//g' | sed 's/-//' | sed 's/+//' | sed 's/+//' 
     $ | wc -c > fp") 
         open(4,file='fp')
         read(4,*)fp
         close(4)
         call system ('rm -r fp')
         if (( fd.eq.2 ).AND.( fp.eq.2 ))then
         write(34,*)'$BASIS DIFFSP=.TRUE. $END' 
         write(34,*)'$BASIS NDFUNC=1 $END' 
         elseif (( fd.eq.3 ).AND.( fp.eq.3 ))then
         write(34,*)'$BASIS DIFFSP=.TRUE. DIFFS=.TRUE. $END'
         write(34,*)'$BASIS NDFUNC=1 NPFUNC=1 $END'
         elseif (( fd.eq.1 ).AND.( fp.eq.2 ))then
         write(34,*)'$BASIS NDFUNC=1 $END'
         elseif (( fd.eq.2 ).AND.( fp.eq.1 ))then
         write(34,*)'$BASIS DIFFSP=.TRUE. $END' 
         elseif (( fd.eq.1 ).AND.( fp.eq.3 ))then
         write(34,*)'$BASIS NDFUNC=1 NPFUNC=1 $END' 
         elseif (( fd.eq.3 ).AND.( fp.eq.1 ))then
         write(34,*)'$BASIS DIFFSP=.TRUE. DIFFS=.TRUE. $END'
         elseif (( fd.eq.2 ).AND.( fp.eq.3 ))then
         write(34,*)'$BASIS DIFFSP=.TRUE. $END'
         write(34,*)'$BASIS NDFUNC=1 NPFUNC=1 $END'
         elseif (( fd.eq.3).AND.( fp.eq.2 ))then
         write(34,*)'$BASIS DIFFSP=.TRUE. DIFFS=.TRUE. $END'
         write(34,*)'$BASIS NDFUNC=1 $END'
         endif

          CASE(2)
          
         call system ("awk 'NR==11{print $3}' *.in | cut -c1 | 
     $   tr '[A-Z]' '[a-z]' > acc")
          open(4,file='acc')
          read(4,*)acc
          close(4)
          call system ('rm -f acc')
         

         if ( acc == "a" ) then
         call system ("awk 'NR==11{print $3}' *.in | cut -c9 
     $   | tr '[A-Z]' '[a-z]' > ccc")
          open(4,file='ccc')
          read(4,*)ccc
          close(4)
          call system ('rm -f ccc')
          if ( ccc == "c" ) then


          call system ("awk 'NR==11{print $3}' *.in | cut -c5-6 |
     $    tr '[A-Z]' '[a-z]' > cc")
          open(4,file='cc')
          read(4,*)cc
          close(4)
          call system ('rm -f cc')

          call system ("awk 'NR==11{print $3}' *.in | cut -c11 |
     $    tr '[A-Z]' '[a-z]' > ccn")
          open(4,file='ccn')
          read(4,*)ccn
          close(4)
          call system ('rm -f ccn')
          write(34,*)'$BASIS GBASIS=',acc,cc
     $    ,ccn,ccc,' $END'

          endif
         endif

         if ( acc == "a" ) then

         call system ("awk 'NR==11{print $3}' *.in | cut -c9 
     $   | tr '[A-Z]' '[a-z]' > ccc")
          open(4,file='ccc')
          read(4,*)ccc
          close(4)
          call system ('rm -f ccc')

          if ( ccc /= "c" ) then

          call system ("awk 'NR==11{print $3}' *.in | cut -c5-6 |
     $    tr '[A-Z]' '[a-z]' > cc")
          open(4,file='cc')
          read(4,*)cc
          close(4)
          call system ('rm -f cc')

          call system ("awk 'NR==11{print $3}' *.in | cut -c10 |
     $    tr '[A-Z]' '[a-z]' > ccn")
          open(4,file='ccn')
          read(4,*)ccn
          close(4)
          call system ('rm -f ccn')
          write(34,*)'$BASIS GBASIS=',acc,cc,ccn,' $END'
         endif 
        endif

         call system ("awk 'NR==11{print $3}' *.in | cut -c5 |
     $    tr '[A-Z]' '[a-z]' > cccc")
          open(4,file='cccc')
          read(4,*)cccc
          close(4)
          call system ('rm -f cccc')
  
         if (( acc /= "a" ).AND.( cccc == "c" ))then

          call system ("awk 'NR==11{print $3}' *.in | cut -c1-2 |
     $    tr '[A-Z]' '[a-z]' > cc")
          open(4,file='cc')
          read(4,*)cc
          close(4)
          call system ('rm -f cc')

          call system ("awk 'NR==11{print $3}' *.in | cut -c7 |
     $    tr '[A-Z]' '[a-z]' > ccn")
          open(4,file='ccn')
          read(4,*)ccn
          close(4)
          call system ('rm -f ccn')
          write(34,*)'$BASIS GBASIS=',cc,ccn,cccc,' $END'
        
        endif


         if (( acc /= "a" ).AND.( cccc /= "c" ))then

          call system ("awk 'NR==11{print $3}' *.in | cut -c1-2 |
     $    tr '[A-Z]' '[a-z]' > cc")
          open(4,file='cc')
          read(4,*)cc
          close(4)
          call system ('rm -f cc')

          call system ("awk 'NR==11{print $3}' *.in | cut -c6 |
     $    tr '[A-Z]' '[a-z]' > ccn")
          open(4,file='ccn')
          read(4,*)ccn
          close(4)
          call system ('rm -f ccn')
          write(34,*)'$BASIS GBASIS=',cc,ccn,' $END'
        
        endif
        write(34,*)"$CONTRL ISPHER=1 $END"

          CASE (3)
!           call system ("awk 'NR==11{print $3}' *.in | 
!     $    tr '[A-Z]' '[a-z]' > sp")
!          open(4,file='sp')
!          read(4,*)sp
!          close(4)
!          call system ('rm -f sp')
!           if(sp == "sbkjc")then           
            write(34,*)'$BASIS GBASIS=',bde(1),' $END'
            write(34,*)'$CONTRL ECP=',bde(1),' $END'
!          else
!            write(34,*)'$BASIS GBASIS=HW $END'
!            write(34,*)'$CONTRL ECP=HW $END'
!           endif

          CASE DEFAULT
         write(*,*)"EL TIPO BASE NO ESTA IMPLEMENTADO EN 
     $  EL ASCEC Y GAMESS",bde(1)
         END SELECT
       endif
       write(34,*)'$CONTRL ICHARG=',iQ,' MULT=',iS2,' $END'
       write(34,*)'$DATA'
       write(34,*)'Molecule specification'
       write(34,*)'C1'
        do i=1,natom
         write(34,100)element(i),iznu(i),r(i,1),r(i,2),r(i,3)
        enddo
  100  Format (A,I3,3X,F20.6,3X,F20.6,3X,F20.6)
        write(34,200)-vert,-vert,-vert
        write(34,200)vert,-vert,-vert
        write(34,200)vert, vert,-vert
        write(34,200)vert, vert, vert
        write(34,200)-vert, vert,-vert
        write(34,200)-vert, vert, vert
        write(34,200)-vert,-vert, vert
        write(34,200)vert,-vert, vert
  200   Format ('Xx  0',3X,F20.6,3X,F20.6,3X,F20.6)
       write(34,*)'$END'
       close(34)
C
C Termina construcci'on del INPUT
C

C
C**** Se llama al gammes para el single-point
C
      write(nproc,'(i5)')proc
      nproc=adjustl(nproc)
      call system ('/usr/src/gamess/rungms anneal.inp 01 
     $ '//trim(nproc)//' > anneal.log ')
      call system ("grep 'EXECUTION OF GAMESS TERMINATED' 
     $ anneal.log | awk '{print $5}' | sed 's/\-//g' > statu")
      open(4,file='statu')
      read(4,*)stat1
      close(4)
      call system ('rm -f statu') 
C
C*** Si el calculo del SP acaba mal
C
       IF(stat1.EQ.'ABNORMALLY')then
        jo=0
        goto 49	 
       ENDIF

       SELECT CASE (thamil)
         CASE (0)
          
          call system ("grep 'TOTAL ENERGY =' anneal.log |
     $ awk 'NR==1 {print $4}' > energy")
          open(5,file='energy')
          read(5,*)eval
          close(5)
          call system ('rm -f energy')

        CASE (1)
          call system ("grep -e 'FINAL.*ENERGY IS' anneal.log |
     $ awk '{print $5}' > energy")
          open(5,file='energy')
          read(5,*)eval
          close(5)
          call system ('rm -f energy')
        
        CASE (2)
          call system ("grep -e 'FINAL.*ENERGY IS.*AFTER' anneal.log |
     $ awk '{print $5}' > energy")
          open(5,file='energy')
          read(5,*)eval
          close(5)
          call system ('rm -f energy')


        CASE (3)
          call system ("grep 'E(MP2)=' anneal.log |
     $ awk '{print $2}' > energy")
          open(5,file='energy')
          read(5,*)eval
          close(5)
          call system ('rm -f energy')

         CASE DEFAULT
         write(*,*)"EL TIPO HAMILTONIANO NO ESTA IMPLEMENTADO EN 
     $  EL ASCEC Y GAMESS"
      END SELECT

      funkGAU=eval/B2
      jo=1	

       ELSEIF(ia.eq.3)then !Call Nwchem

        open(33,file='anneal.nw')
        vert=xbox/2.d0
        write(33,*)'echo'
        write(33,*)''
        write(33,*)'start anneal'
        write(33,*)''
        call system ("echo $USER > name")
        open(4,file='name')
        read(4,*)names
        close(4)
        write(33,*)'scratch_dir /home/scratch/nwchem/',names
        call system ('rm -rf name')
        write(33,*)''
        write(33,*)'memory stack 500 mb heap 100 mb global 1000 mb'
        write(33,*)''
        write(33,*)'title "Annealing" '
        write(33,*)'charge ',iQ
        write(33,*)''
        write(33,*)'geometry units angstroms print xyz autosym'
        do i=1,natom
         write(33,*)element(i),' ',(r(i,j),j=1,3)
        enddo
        write(33,*)'Xx  ',-vert,-vert,-vert
        write(33,*)'Xx  ', vert,-vert,-vert
        write(33,*)'Xx  ', vert, vert,-vert
        write(33,*)'Xx  ', vert, vert, vert
        write(33,*)'Xx  ',-vert, vert,-vert
        write(33,*)'Xx  ',-vert, vert, vert
        write(33,*)'Xx  ',-vert,-vert, vert
        write(33,*)'Xx  ', vert,-vert, vert
        write(33,*)'end'
        write(33,*)''

         write(33,*)'basis'
          do i=1,ielem
           write(33,*)sym(iz(i)),' library ',bde(i)  
          enddo 
         write(33,*)'end'
         write(33,*)''
         j=0
          do i=1,ielem
           if(tbase(i).eq.3)then
             j = j + 1
           endif
          enddo
          if(j > 0)then
           write(33,*)'ecp'
            do i=1,ielem
             if(tbase(i)==3)then
              write(33,*)sym(iz(i)),' library ',bde(i)
             endif
            enddo
           write(33,*)'end'
           write(33,*)''
          endif

        SELECT CASE (thamil)
         CASE (0)
         write(33,*)'scf'
         write(33,*)'',mspin
         write(33,*)'end'
         write(33,*)''
         write(33,*)'task scf energy'
         CASE (2)
         write(33,*)'dft'
         write(33,*)'xc ', hamil 
         write(33,*)'mult', iS2
         write(33,*)'end'
         write(33,*)''
         write(33,*)'task dft energy'
         
         CASE(3)
          if((hamil == "mp2").OR.(hamil == "MP2"))then
           write(33,*)'mp2'
           write(33,*)'freeze core atomic'
           write(33,*)'end'
           write(33,*)''          
           write(33,*)'task mp2 energy'
          else          
           write(33,*)'tce'
           write(33,*)'',hamil
           write(33,*)'freeze core atomic'
           write(33,*)'end'
           write(33,*)''          
           write(33,*)'task tce energy'
         endif
  
         CASE DEFAULT
         write(*,*)"EL TIPO HAMILTONIANO NO ESTA IMPLEMENTADO EN 
     $  EL ASCEC Y NWCHEM"
        END SELECT
        close(33)
C**************************************************************
C SP con NWCHEM
C**************************************************************
      write(nproc,'(i5)')proc
      nproc=adjustl(nproc)
      call system ('mpirun -np '//trim(nproc)//' nwchem anneal.nw 
     $ > anneal.log')
C 
C ** Termino mal
      call system ("grep ' There is an error ' anneal.log | 
     $ awk '{print $4}' > statu")
      call system ("wc statu | awk '{print $3}' > sistatu")
      open(4,file='sistatu')
      read(4,*)sistat1
      close(4)
      call system ('rm -f sistatu')
      if(sistat1.eq.6)then
       open(4,file='statu')
       read(4,*)stat1
       close(4)
       call system ('rm -f statu')
       if(stat1.EQ.'error')then
        jo=0
        goto 49	 
       endif
      endif
      call system ('rm -f statu')

      SELECT CASE (thamil)
    
      CASE (0)      
      call system ("grep 'Total SCF energy' anneal.log |
     $ awk 'NR==1 {print $5}' > energy")

      CASE (2)
      call system ("grep 'Total DFT energy' anneal.log |
     $ awk 'NR==1 {print $5}' > energy")

      CASE (3)
      call system ("grep 'Total MP2 energy' anneal.log |
     $ awk '{print $4}' > energy")
 
      CASE DEFAULT
         write(*,*)" " 
      END SELECT

      open(5,file='energy')
      read(5,*)eval
      close(5)
      call system ('rm -f energy')
      funkGAU=eval/B2
      jo=1	
C=============================================================================
C INPUT y calculo con el G09
C
	ELSEIF(ia.eq.1)then !Call Gaussian	
 
      open(3,file='anneal.com')      
      vert=xbox/2.d0
      write(3,*)'%nproc=',proc
      if(iS2.eq.1)then
      write(3,*)'# ',hamil,'/',bde(1)
      else
      write(3,*)'# ',hamil,'/',bde(1),'scf=qc'
      endif
      write(3,*)''
      write(3,*)'Annealing'
      write(3,*)''       
      write(3,*)iQ,iS2
      do i=1,natom
      write(3,*)iznu(i),' 0 ',(r(i,j),j=1,3)
      enddo
      write(3,*)'x  0 ',-vert,-vert,-vert
      write(3,*)'x  0 ', vert,-vert,-vert
      write(3,*)'x  0 ', vert, vert,-vert
      write(3,*)'x  0 ', vert, vert, vert
      write(3,*)'x  0 ',-vert, vert,-vert
      write(3,*)'x  0 ',-vert, vert, vert
      write(3,*)'x  0 ',-vert,-vert, vert
      write(3,*)'x  0 ', vert,-vert, vert
      write(3,*)''           
      close(3) 
C***************************************************
C SP con G09     
      if ((alias.eq."g03").OR.(alias.eq."G03")) then
       call system ('g03 anneal.com')
      else
       call system ('g09 anneal.com')
      endif
C Se busca si t'ermino bien el calculo
C Se guarda la palabra clave en el archivo statu
      call system ('grep termination anneal.log | 
     $ sed s/termination.*// > statu')
C Se lee la palabra clave en statu
      open(4,file='statu')
      read(4,*)stat1
      close(4)
      call system ('rm -f statu') 
C Si sale ERROR se va al final de este archivo
      if(stat1.eq.'Error')then
      jo=0
      goto 49	 
      endif         	     
C Se guarda el valor de la energ'ia
      call system ('grep Energy= anneal.log | 
     $ cut -c10-24 > energy')     
      call system ('grep Done: anneal.log |
     $ sed s/.*=// | sed s/A.U.*// >> energy')
C Se lee el valor de energ'ia
      open(5,file='energy')
      read(5,*)eval
      close(5)
      call system ('rm -f energy')
C 
      funkGAU=eval/B2
      jo=1	

	ELSEIF(ia.eq.2)then !Call MSINDO.1.1.opt
	open(32,file='anneal.inp')
	write(32,*)'Annealing :Neu'
	write(32,*)hamil,' CHARGE=',iQ,' MULTIP=',iS2
	write(32,*)'MAXCYC=100'
	write(32,*)'CARTES'
	write(32,*)'NATOMS=',natom
	write(32,*)':End'
      do i=1,natom
      write(32,*)iznu(i),(r(i,j),j=1,3)
      enddo	
	write(32,'(a4)')'ENDE'
      write(32,*)''	
	close(32)	
	call system ('MSINDO.1.1.opt < anneal.inp > anneal.out')
	call system ('grep -a ERROR anneal.out > statu2')	
	call system ('grep -a CONVERGED anneal.out >> statu2')
      open(42,file='statu2')
      read(42,*)stat2
      close(42)
      call system ('rm -f statu2')	
      if(stat2.eq.'***')then
      jo=0
      goto 49	 
      endif	
	call system ('grep -a TOTAL anneal.out | grep -a ENERGY 
     $ | cut -c22-36 > energy2')
      open(52,file='energy2')
      read(52,*)eval
      close(52)
      call system ('rm -f energy2')
      funkGAU=eval/B2
      jo=1	
C=======================================================================
CINPUT del ADF y Calculo con el ADF
C
       ELSEIF(ia.eq.5)then !Call Nwchem
C Abro el anneal.run, en el cual escribo la informaci'on del input del ADF
        open(39,file='anneal.run')
C
        vert=xbox/2.d0 !calculo para la caja
C
C        write(39,*)'rm -r logfile'
C        write(39,*)'rm -r TAPE21'
C        write(39,*)' rm -r TAPE13'
C        write(39,*)''
C Adf en mi laptop
c        write(39,'(A,I2,A)')"faketime '2013-01-01' /home/danian/Document
c     &s/programas/adf2013.01/bin/adf -n",proc," <<eor>> anneal.out"
C ADF en el cluster-IMIT
        write(39,'(A,I2,A)')"/home/andyzapa/adf/dfaketime/bin/faketime 
     &'2013-01-01' /home/andyzapa/adf2013.01/bin/adf -n ",proc," <<eor>>
     & anneal.out"
        write(39,*)''
        write(39,*)'TITLE Annealing ASCEC'
        write(39,*)''
C Se escribe la carga y multiplicidad
        write(39,'(A,I1,A)')'CHARGE {',iQ,',0}'
        write(39,*)''
C Se adiciona la geometr'ia con la coordenadas de los 'atomos fantasmas
        write(39,*)'ATOMS'
        do i=1,natom
         write(39,*)element(i),' ',(r(i,j),j=1,3)
        enddo
        write(39,*)'Xx  ',-vert,-vert,-vert
        write(39,*)'Xx  ', vert,-vert,-vert
        write(39,*)'Xx  ', vert, vert,-vert
        write(39,*)'Xx  ', vert, vert, vert
        write(39,*)'Xx  ',-vert, vert,-vert
        write(39,*)'Xx  ',-vert, vert, vert
        write(39,*)'Xx  ',-vert,-vert, vert
        write(39,*)'Xx  ', vert,-vert, vert
        write(39,*)'END'
        write(39,*)''
C Se pone la base
        write(39,*)'BASIS'
        write(39,*)'type ',bde(1)
        write(39,*)'core None'
        write(39,*)'END'
        write(39,*)''
C
        write(39,*)'symmetry auto'
        write(39,*)''
C 
C Hamiltoniano
C
        SELECT CASE (thamil)
         CASE (0)
          write(39,*)'RELATIVISTIC Scalar ZORA'
          write(39,*)''
          write(39,*)'XC'
          write(39,*)'GGA ',hamil
          write(39,*)'END'
          write(39,*)''

         CASE (2)
          write(39,*)'RELATIVISTIC Spin-Orbit ZORA'
          write(39,*)''
          write(39,*)'XC'
          write(39,*)'GGA ',hamil
          write(39,*)'END'
          write(39,*)''
         
         CASE DEFAULT
          write(*,*)"EL TIPO HAMILTONIANO NO ESTA IMPLEMENTADO EN 
     $               EL ASCEC CON ADF"
        END SELECT
C
c          write(39,*)'SAVE TAPE21 TAPE13'
          write(39,*)'SCF'
          write(39,*)'CONVERGE 1e-6 1e-6'
          write(39,*)'ITERATION 200'
          write(39,*)'END'
          write(39,*)''
C
          write(39,*)'INTEGRATION 6.0'
          write(39,*)''
C
c          write(39,*)'eor'
c          write(39,*)''
C
C Cierro archivo del input del ADF
         close(39)
C********************************************
C Calculos SP con ADF
         call system ('bash anneal.run')
C********************************************
C Se busca si término bien o no, se guarda la palabra clave 
C en el archivo statu
         call system ("tail -n2 logfile | awk 'NR==1' | 
     & awk '{print $3}' > statu")
C Se lee la palabra clave en el archivo statu
       open(4,file='statu')
       read(4,*)stat1
       close(4)
       call system ('rm -f statu') 
C
       IF(stat1.ne.'NORMAL')THEN
        jo=0
        GOTO 49 !Se va al final de este archivo
       ENDIF
C     
        call system("grep 'Bond Energy' logfile | tail -n3 | 
     & awk 'NR==1' | awk '{print $5}' > energy")    
C
       call system("wc -l energy | awk '{print $1}' > lineas")

      open(8,file='lineas')
      read(8,*)nli
      close(8)
      call system ('rm -f lineas')

      IF(nli.eq.0) THEN
       jo=0
       call system ('rm -f energy')
       GOTO 49
      ELSE
       open(5,file='energy')
       read(5,*)eval
       close(5)
       call system ('rm -f energy')
      ENDIF
C
      IF (eval.EQ.0.0D0) THEN
      call system ('rm -f anneal.out TAPE21 KidOut* logfile t21.*')
       GOTO 49 
      ENDIF
      IF (eval.GT.0.0D0) eval = -1.0D0/eval
C        
      funkGAU=eval/B2
      jo=1


      call system ('rm -f anneal.out TAPE21 KidOut* logfile t21.*')
c===========================================================
********DALTON
c===========================================================
      ELSEIF(ia.eq.6)THEN
!-Input Programs
       OPEN(200,FILE="anneald.dal",STATUS='UNKNOWN')
       CALL SYSTEM('rm -f anneal*.out')
!- start: .dal
       WRITE(200,FMT1)'**DALTON INPUT'
       WRITE(200,FMT1)'.RUN WAVE FUNCTIONS'
       WRITE(200,FMT1)'**WAVE FUNCTIONS'

       IF (thamil==0) THEN
        WRITE(200,FMT1)'.HF'
       ELSEIF(thamil==3) THEN
        WRITE(200,FMT1)'.HF'
        WRITE(200,FMT1)'.MP2'
       ELSEIF(thamil==2) THEN
        WRITE(200,FMT1)'.DFT'
        WRITE(200,FMT1)hamil
       ELSE
        WRITE(*,*)'No está implementado en el DALTON o ASCEC'
       ENDIF

c       WRITE(200,FMT1)'*CONFIGURATION INPUT'
c       WRITE(200,FMT1)'.SPIN MULTIPLICITY'
c       WRITE(200,FMT2)Multiplicity
       WRITE(200,FMT1)'*SCF INPUT'
       WRITE(200,FMT1)'.THRESH'
       WRITE(200,FMT1)' 1.0D-7'
       WRITE(200,FMT1)'**END OF DALTON INPUT'

       CLOSE(UNIT=200)
!- end: .dal
!- start: .mol
       OPEN(300,FILE="annealm.mol",STATUS='UNKNOWN')
       WRITE(300,FMT1)'ATOMBASIS'
       WRITE(300,FMT1)'Configuration'
       WRITE(300,FMT1)'Coordenadas'
       
       l = 0
       DO i = 1, nmo
        DO j=1,namolec(i)
         l = l + 1
        ENDDO
       ENDDO
       l = l !+ 8
       IF(iQ==0) THEN
        WRITE(300,FMT3)'Angstrom Atomtypes=',l,
     &                 ' Generators=0'
       ELSE
        WRITE(300,FMT4)'Angstrom Atomtypes=',l,' Charge=',Charge, 
     &      ' Generators=0'
       ENDIF

      l = 1
      DO i=1,nmo
       DO j=1,namolec(i)
        cout=0
        DO k = 1, ielem !# de tipo de elementos
         IF (iz(k).eq.iznu(l).and.cout.eq.0) THEN
          WRITE(300,FMT5)'Charge=',iznu(l),' Atoms=1 Basis=',bde(k)
          WRITE(300,FMT6)element(l),r(l,1),r(l,2),r(l,3)
          l = l + 1
          cout = 1
         ENDIF
        ENDDO
       ENDDO
      ENDDO
       !WRITE(300,FMT5)'Charge=0.0 Atoms=8 Basis=STO-2G'
       !write(300,FMT7)'X ',-vert,-vert,-vert
       !write(300,FMT7)'X ', vert,-vert,-vert
       !write(300,FMT7)'X ', vert, vert,-vert
       !write(300,FMT7)'X ', vert, vert, vert
       !write(300,FMT7)'X ',-vert, vert,-vert
       !write(300,FMT7)'X ',-vert, vert, vert
       !write(300,FMT7)'X ',-vert,-vert, vert
       !write(300,FMT7)'X ', vert,-vert, vert

       CLOSE(300)
!---------------------------------------------------------------------------------------------
!Calculate
       write(nproc,'(i5)')proc
       nproc=adjustl(nproc)
       CALL SYSTEM("echo $PATH_DALTON > statu")
       OPEN(4,FILE='statu')
       READ(4,'(A)')stat3
       CLOSE(4)
       CALL SYSTEM('rm -f statu')
 
       commanddalton=""//trim(stat3)//" -noarch -t /home/scratch/dalton
     & -N "//trim(nproc)//" anneald annealm"

       CALL EXECUTE_COMMAND_LINE(commanddalton)

       l = 0
       CALL SYSTEM("grep ABORTED anneald_annealm.out | wc -l > statu")
       OPEN(4,FILE='statu')
       READ(4,*)l
       CLOSE(4)
       CALL SYSTEM('rm -f statu')
       IF (l.GT.0) THEN
        jo = 0
        GOTO 49
       ENDIF

       CALL SYSTEM("grep '@    Electronic energy:' anneald_annealm.out |
     &          awk '{print $2}' > statu")
       OPEN(4,FILE='statu')
       READ(4,*)stat1
       CLOSE(4)
       CALL SYSTEM('rm -f statu')

       IF(stat1=='Electronic') THEN
        CALL SYSTEM("grep '@    Electronic energy:' anneald_annealm.out|
     &      awk '{print $4}' >> energy")
         OPEN(5,FILE='energy')
         READ(5,*)eval
         CLOSE(5)
         funkGAU=eval/B2
         jo=1
         CALL SYSTEM('rm -f energy')
       ELSE
        jo = 0
        GOTO 49
       ENDIF


      ELSE
C ===============Programa no articulado======================
        write(6,*)'WARNING: No program to evalue energy'
        STOP
      ENDIF
!       call system ('unset NPROC')	
  49  continue      	 

      return
      end
*******************************************************************************
