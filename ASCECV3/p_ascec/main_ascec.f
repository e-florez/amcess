*******************************************************************************
*                               MODIFICACION "Version 3"
*     Por: Andy Danian Zapata Escobar
*     06/18/2020
*
*     Adecuacion para poder articular el ascec con el nwchem y gamess, ademas 
*     de adicionar el alias del program, el tipo de hamiltoniano y la base,
*     se debe adicionar los simbolos de los elementos en la geometria y
*     tambian se cambio el numero de la semilla, ahora es la hora con el año
*     (hminsegaño).
*     Adici'on del ADF
*

      program main_ascec
      implicit real*8 (a-h,o-z)
      integer h, proc, thamil,tbase
      include 'param.inc'     
      dimension rp(max_atom,3),rf(max_atom,3)
      !dimension nnu(10)
      common/elemente/element(max_atom),mspin,namolec(max_atom)
      common/names/iznu(max_atom),sym(120),xbox,nnu(max_atom)
      common/rando/idum
      common/label/jdum
      common/proce/proc,thamil,ielem
      common/db/bde(max_atom),tbase(max_atom),iz(max_atom)
      common/size/Ratom(20)
      common/ruta/ia,iQ,iS2,hamil,alias
      common/molecu/rcm(max_mole,3),imolec(0:max_mole),nmo,ds,dphi	
      common/weigth/wt(max_atom)            
      character*2 sym
      character*2 element
      character*10 jdum, idum
      character*8 hamil
      character*8 tbde
      character*15 alias
      character*20 bde
      character*20 molecula(max_mole)
      character*8 mspin     
      data B1/8.6173D-5/	!Const. Boltzmann en eV/K
      data B2/3.166D-6/		!Const. Boltzmann en Hartree/K
      data Ratom/0.373,1.280,1.520,1.133,0.830,0.772,0.710,0.604,0.709,
     $ 1.500,1.537,1.600,1.431,1.170,0.093,1.040,0.994,1.740,2.27,1.97/
      open(6,file='ascec.out')
!      open(14,file='tvse.dat')
      
!      call system('rm -f mto.xyz')
!      call system('rm -f result.xyz')
!      call system('rm -f anneal.com')
!      call system('rm -f anneal.log')
!      call system('rm -f rless.out')      	        
      call system('echo ""')
      call system('echo "             ASCEC V. 03.     "')
      call system('echo ""')
      call system('echo "Autores: Jhon Fredy Perez"')
      call system('echo "         Andy Danian Zapata Escobar"')
      call system('echo "         Albeiro Restrepo"')
      call system('echo ""')
      call system('echo "Institucion: Universidad de Antioquia."')
      call system('echo "Grupo: Quimica-Fisica Teorica."')
      call system('echo "Pagina: quimica.udea.edu.co/~qfteorica/"')
      call system('echo ""')
      call system('echo "Institucion: Universidad del Nordeste de la"') 
      call system('echo "                   Argentina."')
      call system('echo "Grupo: IMIT-CONICET."')
      call system('echo ""')
      call system('echo "Este programa se puede usar para realizar"')
      call system('echo "un muestreo estocastico de las posibles"')
      call system('echo "conformaciones del sistema atomico o "')
      call system('echo "molecular de estudio, en base a una"')
      call system('echo "modificacion del criterio de metropolis"')
      call system('echo "(mayor informacion, en la carpeta de"')
      call system('echo "papers y examples). Comando ejecucion,"')
      call system('echo "ascecv03 < abc.in > abc.out."')
      call system('echo ""')
      call system('echo "Nota: No poner como extension inp, porque"')
      call system('echo "esta extension es asignado a los archivos"')
      call system('echo "del gamess. Una plantilla del archivo abc.in"')
      call system('echo "la puede encontrar en examples".')
      call system('echo ""')
      call system('echo ""')
      call system('echo "GAMESS"')
      call system('echo ""')
      call system('echo "El ASCEC tiene implementados los formatos"')
      call system('echo "del gamess de los metodos hf, dft, mp2 (el"') 
      call system('echo "gamess solo tiene implementado mp2) y"')
      call system('echo "semi-empiricos, en capa cerrada y abierta,"')
      call system('echo "ademas, de los formatos de las bases sto,"')
      call system('echo "pople, duning corre. cons. y algunos "')
      call system('echo "pseudopotenciales (SBKJL valence"')
      call system('echo "{alias=sbkjl} y Hay/Wadt Valence {alias=hw})"')
      call system('echo ""')
      call system('echo ""')
      call system('echo "NWCHEM"')
      call system('echo ""')
      call system('echo "El ASCEC tiene implementados los formatos "')
      call system('echo "del nwchem de los metodos hf, dft mp2 y mpn,"')
      call system('echo "para n mayores a 2, en capa cerrada y"')
      call system('echo "abierta, ademas, de los formatos de las "')
      call system('echo "bases sto, pople, dunings corre. cons.,"')
      call system('echo "pseudopotenciales (Hay-Wadt MB (n+1) ECP,"')
      call system('echo "LANL2DZ ECP, SBKJC VDZ ECP, CRENBL ECP,"')
      call system('echo "Stuttgart RLC ECP, entre otras. El ASCEC "')
      call system('echo "tambien tiene implementado el formato del "')
      call system('echo "nwchem para usar un tipo base por clase"')
      call system('echo "atomo."')
      call system('echo ""')
      call system('echo ""')
      call system('echo "ADF"')
      call system('echo ""')
      call system('echo "El ASCEC tiene implementados los formatos "')
      call system('echo "del adf SR y SR-SO, con capa cerrada."')
      call system('echo ""')
      call system('echo ""')
      call system('echo "GAUSSIAN Y MSINDO1.1.opt"')
      call system('echo ""')
      call system('echo "Ademas el ASCEC esta programado para"')     
      call system('echo "ejecutarse con el gaussian y el MSINDO"')
      call system('echo ""')     

      call symelem			!Simbolos de los elementos quimicos
      call wtelem			!Masas atomicas de los elementos qcos.
! Lectura de datos ---------------------------------------------------------
      read(*,*)ivalE,mto                !Evalua E? (0:No -> mto movimientos)
      read(*,*)ia,alias			!1,2,3,4,5 (g09,MSINDO.1.1.opt,nwchem,gamess,adf)
      read(*,*)xbox			!Longitud del cubbo
      read(*,*)iRT			!Ruta de enfriamiento
      read(*,*)To1,dT,nT1      	        !T_inicial / delta_T / núm_T
      read(*,*)To2,pdis,nT		!T_inicial / %dis_T / núm_T
      read(*,*)maxstep			!Numero maximo de pasos (inicialmente)
      read(*,*)proc                     !Numero procesadores
      read(*,*)hamil,thamil
      read(*,*)ielem			!Numero de elementos     
      k=0	
      if(thamil==1)then
       do i=1,ielem
         read(*,*)iz(i),nnu(i)  !Tipo de nucleo / Numero de nucleos 
	  do j=1,nnu(i)
	   k=k+1
	   iznu(k)=iz(i)			!Numero atomico de cada nucleo
	  enddo
       enddo
      natom=k				!Numero total de nucleos
      else
       do i=1,ielem
         read(*,*)iz(i),nnu(i),bde(i),tbase(i)  !Tipo de nucleo / Numero de nucleos / Base para cada nucleo / Tipo base cada elemento
	  do j=1,nnu(i)
	   k=k+1
	   iznu(k)=iz(i)			!Numero atomico de cada nucleo
	  enddo
       enddo
      natom=k				!Numero total de nucleos
      endif
      read(*,*)iQ,iS2			!carga / espinín

      IF(ia==3)then
       SELECT CASE (int(iS2)) 
         CASE (1)
         mspin = "SINGLET"
         CASE (2)
         mspin = "DOUBLET"
         CASE (3)
         mspin = "TRIPLET"
         CASE (4)
         mspin = "QUARTET"
         CASE (5)
         mspin = "QUINTET"
         CASE (6)
         mspin = "SEXTET"
         CASE (7)
         mspin = "SEPTET"
         CASE (8)
         mspin = "OCTET"
         CASE (9)
         mspin = "NOPEN"
         CASE DEFAULT
         write(*,*)"LA MULTIPLICIDAD NO TIENE TRADUCCION LITERAL 
     $  PARA NWCHEM"
        END SELECT
      END IF

      read(*,*)
      read(*,*)nmo			!numero de moleculas
      read(*,*)
      inu=0
      imolec(0)=0
      do imo=1,nmo
       read(*,*)molecula(imo),jb
       namolec(imo) = jb
       imolec(imo)=imolec(imo-1)+jb	
       do h=1,jb
        inu=inu+1
        read(*,*)iznu(inu),element(inu),rp(inu,1),rp(inu,2)  !Z, Simbolo, Coordenadas
     $  ,rp(inu,3)   
       enddo
      enddo
      read(*,*)
      read(*,*)ds,dphi !Max. displacement (A); Max. rotation angle (radians)

c=========Empieza Esciruta en Pantalla
      write(6,*)
      write(6,*)'***********************************************'
!      call system('date')
      write(6,*)'Annealing Simulado Con Energia Cuantica **ASCEC**'
	write(6,*)'ASCEC-V02: February-2010'
      write(6,*)
      write(6,*)'Elemental composition of the system:'
      do i=1,ielem
       write(6,*)sym(iz(i)),' ',nnu(i)
!       write(6,*)iz(i),' ',nnu(i)
      enddo
      write(6,*)'There are a total of',natom,' nuclei'
      write(6,*)"Cube's length =",xbox,' A'
      call masscen(natom,rp)		!Halla el centro de masa de cada mole.	
      write(6,*)
      write(6,*)'Number of molecules:',nmo
      write(6,*)'Molecular composition'
      do imo=1,nmo
       write(6,*)molecula(imo),
     $ (sym(iznu(j)),j=imolec(imo-1)+1,imolec(imo))
!       write(6,*)molecula(imo),
!     $ (iznu(j),j=imolec(imo-1)+1,imolec(imo))
      enddo
      write(6,*)	
      write(6,*)'Maximum displacement of each mass center =',ds,' A'
      write(6,*)'Maximum rotation angle =',
     $ dphi,' radians' 

C
      IF(ivalE.eq.0)THEN
C
!      call system (' date | cut -c15-16,18-19,25-28 > dum')
       call system (" date | awk '{print $4, $6}' | sed 's/://g' | 
     $ sed 's/ //g' > dum")
       open(11,file='dum',status='unknown')
       read(11,*)idum
       close(11)
       open(11,file='dum',status='unknown')
       read(11,*)jdum
       close(11)
       write(6,*)'Seed =',idum
       call system('rm -f dum')	
       write(6,*)
       write(6,*)'** Energy will not be evaluated **'	
       DO i=1,mto
        call config(natom,rp)
        call drawmol(natom,rp,i,0,0.d0)  	
       ENDDO
       write(6,*)'Will produce ',mto,' random configurations'
       write(6,*)'Coordenates stored in mto_'//jdum//'.xyz'
       GOTO 200
C
      ENDIF
C
      IF(iRT.eq.1)THEN
       write(6,*)'Linear quenching route. To ='
     $           ,To1,' Tf =',Tf,' dT =',dT
       Tm1=To1
      ENDIF
C 
      IF(iRT.eq.2)THEN
      write(6,*)'Geometrical quenching route. To ='
     $           ,To2,' %dism =',pdis,' nT =',nT
      Tm1=To2
      ENDIF    	

      write(6,*)
      write(6,*)'Energy calculated with ',alias
      write(6,*)'Hamiltonian: ',hamil
       write(6,*)'Basis set: '
       do i=1,ielem
        write(6,*)sym(iz(i)),'  ',bde(i)
       enddo
      write(6,*)'Charge =',iQ,' Multiplicity =',iS2
      write(6,*)'Processors: ',proc 	 
!      call system ('rm -f anneal.com')
!      call system ('rm -f anneal.log')      
!      call system ('touch anneal.com')      
	
!**  Inicio de annealing	   
!      call system (' date | cut -c15-16,18-19,25-28 > dum')
      call system ("date | awk '{print $4, $6}' | sed 's/://g' | 
     $      sed 's/ //g' > dum")
      open(2,file='dum',status='unknown')
      read(2,*)idum
      close(2)
      open(2,file='dum',status='unknown')
      read(2,*)jdum
      close(2)
	write(6,*)'Seed =',idum
      call system('rm -f dum')   
C
! Se abre los archivos tvse_semilla.dat
!                     rless_semilla.out
C
       open(14,file='tvse_'//jdum//'.dat',status='unknown')
       open(17,file='rless_'//jdum//'.out',status='unknown')
C
C Se escribe en rless_semilla.out, geometr'ia
C
       do imo=1,nmo
        write(17,'(a20,20x,I3)')molecula(imo),imolec(imo)-imolec(imo-1)
        do inu=imolec(imo-1)+1,imolec(imo)
	 k=iznu(inu)
         write(17,'(I2,2x,3F12.6)')k,(rp(inu,j),j=1,3)
        enddo
       enddo
C Se cierra el archivo rless_semilla.out
       close(17) 
C
      do i=1,natom
       do j=1,3
        rf(i,j)=rp(i,j)
       enddo
      enddo
C
C La funci'on funkGAU, se encarga de crear el archivo de entrada
C del programa seleccionada
C
      Ep=funkGAU(natom,rp,jo)    
!      Ep=0.
      nstep=1
      icount=0
      iboltz=0
      oldE=Ep
      write(6,*)''
      if(jo.eq.0)then
       write(6,*)'No wavefunction convergency for initial configuration'
       else
       write(6,*)'Energy of initial configuration =',Ep*B2,' u.a.'
      endif
      write(14,'(I8,F10.2,F20.6)')nstep,Tm1,Ep*B2       
      call drawmol(natom,rp,nstep,1,Ep*B2)      

      write(6,*)
      write(6,*)'History: (T(K),E(u.a.),n-eval)'      
      write(6,*)                
      if(iRT.eq.1)then 
      goto 11
      endif     
      if(iRT.eq.2)then
      goto 12
      endif       
! Ruta de enfriamiento lineal ------------------------------------------------
  11  continue     	
      do iT=1,nT       
  21   call config(natom,rf)
       Ef=funkGAU(natom,rf,jo)	
       if(jo.eq.0)goto 21
       nstep=nstep+1
!       write(14,'(I8,F10.2,F20.6)')nstep,Ti,Ef*B2      
       icount=icount+1
       if(icount.eq.maxstep)then
       write(6,'(F10.2,I10,a20)')Ti,icount,'steps'
       maxstep=maxstep*3/4
       if(maxstep.lt.10)maxstep=10       
       goto 41
       endif
       del=Ef-oldE
       if(del.lt.0.d0)then
       write(14,'(I8,F10.2,F20.6)')nstep,Ti,Ef*B2
       Ep=Ef              
       oldE=Ep
       do i=1,natom
        do j=1,3
         rp(i,j)=rf(i,j)
        enddo
       enddo       
       open(17,file='rless_'//jdum//'.out',status='unknown')
       do imo=1,nmo
        write(17,'(a20,20x,I3)')molecula(imo),imolec(imo)-imolec(imo-1)
        do inu=imolec(imo-1)+1,imolec(imo)
	 k=iznu(inu)
         write(17,'(I2,2x,3F12.6)')k,(rp(inu,j),j=1,3)
        enddo
       enddo
       close(17)        
       call drawmol(natom,rp,nstep,1,Ep*B2)  
       goto 31
       else
       pE=dexp(-(del)/Ti)
       crt=del/dabs(Ef)
       if(crt.lt.pE)then
       iboltz=iboltz+1	  
       write(14,'(I8,F10.2,F20.6)')nstep,Ti,Ef*B2
       Ep=Ef      
       do i=1,natom
        do j=1,3
         rp(i,j)=rf(i,j)
        enddo
       enddo
       call drawmol(natom,rp,nstep,1,Ep*B2)              
       write(6,*)'Max.-Boltz.',crt,pE	    
       endif	    
       goto 21
       endif
  31   continue
       write(6,'(G10.2,F20.6,I8)')Ti,Ep*B2,icount
  41   continue
       icount=0
       Ti=To1-dT        	   
      enddo
      goto 100
!-----------------------------------------------------------------------------

! Ruta de enfriamiento por progresion ----------------------------------------
  12  continue
      Ti=To2 	      
      do iT=1,nT              
  22   call config(natom,rf)    
       Ef=funkGAU(natom,rf,jo)	       
!	Ef=ran0(idum)
       nstep=nstep+1
!       write(14,'(I8,F10.2,F20.6)')nstep,Ti,Ef*B2      
       icount=icount+1
       if(icount.eq.maxstep)then
       write(6,'(F10.2,I10,a20)')Ti,icount,'steps'
       maxstep=maxstep*3/4
       if(maxstep.lt.10)maxstep=10
       goto 42
       endif
       del=Ef-oldE
       if(del.lt.0.d0)then
       write(14,'(I8,F10.2,F20.6)')nstep,Ti,Ef*B2
       Ep=Ef
       oldE=Ep
       do i=1,natom
        do j=1,3
         rp(i,j)=rf(i,j)
        enddo
       enddo       
       open(17,file='rless_'//jdum//'.out',status='unknown')
       do imo=1,nmo
        write(17,'(a20,20x,I3)')molecula(imo),imolec(imo)-imolec(imo-1)
        do inu=imolec(imo-1)+1,imolec(imo)
	 k=iznu(inu)
         write(17,'(I2,2x,3F12.6)')k,(rp(inu,j),j=1,3)
        enddo
       enddo
       close(17)        
       call drawmol(natom,rp,nstep,1,Ep*B2)         
       goto 32
       else
       pE=dexp(-(del)/Ti)
       crt=del/dabs(Ef)
       if(crt.lt.pE)then
       write(14,'(I8,F10.2,F20.6)')nstep,Ti,Ef*B2
       Ep=Ef
       iboltz=iboltz+1	  
       do i=1,natom
	do j=1,3
         rp(i,j)=rf(i,j)
	enddo
       enddo
       call drawmol(natom,rp,nstep,1,Ep*B2)              
       write(6,*)'Max.-Boltz.',crt,pE
       endif	    
       goto 22
       endif
  32   continue
       write(6,'(G10.2,F20.6,I8)')Ti,Ep*B2,icount	
  42   continue	
       icount=0

       Ti=Ti-Ti*pdis/100.d0 !Cambio de T    
      enddo
      goto 100
!-----------------------------------------------------------------------------

 100  write(6,*)
      write(6,*)'** Normal annealing termination **'
      write(6,*)'Energy was evaluated',nstep,' times'
      write(6,*)'Energy evolution in tvse_'//jdum//'.dat'
      write(6,*)'Configurations accepted by Max.-Boltz. statistics'
     $, iboltz
      write(6,*)'Accepted configurations in result_'//jdum//'.xyz'
      write(6,*)'Lowest energy configuration in rless_'//jdum//
     $'.out'
      write(6,*)'Lowest energy =',oldE*B2,' u.a.'
!      write(6,*)('Configuracion de minima energia:')
!      call system('cat historial.out rless.out') 
      write(6,*)''
!      call system('date')
 200  continue
      write(6,*)'*************************************************'
      write(6,*)

      close(14)	
      close(6)
      stop
      end
      	
*******************************************************************************
