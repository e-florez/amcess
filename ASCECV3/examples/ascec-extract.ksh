#!/bin/ksh
########################################################
#            MODIFICACION DEL ASCEC-EXTRACT
# 
#POR: Andy Danian Zapata Escobar
#       04/01/2015
#
#

echo " "
echo "               ASCEC EXTRACT"
echo " "
echo "Autor:  Jhon Fredy Perez"
echo "        Andy Danian Zapata Escobar"
echo "        Albeiro Restrepo"
echo " "
echo "Institucion: Universidad de Antioquia."
echo "Grupo: Quimica-Fisica Teorica."
echo "Pagina: quimica.udea.edu.co/~qfteorica/"
echo " "
echo "Script para crear un archivos con el formato de los inputs
(.nw y .inp) del NWChem o el Gamess, por cada geometria obtenida 
por el ascecv03. Ejemplo de ejecucion:

 ksh ascec-extract.ksh result...xyz"
echo " "

echo "Nota: Solo esta implementado el formato del nwchem, para asignar
una base a cada clase de elemento, del sistema de estudio."
echo ""

echo "Programa a usar (gaussian:0/nwchem:1/gamess:2)"
read programa 
echo ""
if [[ $programa -eq 2 ]]; then
  echo "Este script tiene implemetado los formatos de las bases sto, 
pople y Dunning's correlation consistent, para el gamess, además, 
de los métodos RHF, DFT, MP2 (gamess solo tiene implemetado mp2) y 
coupled cluster, para capa cerrada; en el caso de capa abierta se 
tiene UHF, DFT, MP2 y ROHF con CCSD (gamess en capa abierta solo 
soporta para coupled cluster CCSD, mayor informacion 
http://www.msg.ameslab.gov/gamess/documentation.html)."
  echo " "
fi
if [[ $programa -eq 1 ]]; then
  echo "El nwchem tiene implemetado las bases sto, pople, Dunning's
correlation consistent, Jensen, iglo-II y iglo-iii, además, de 
los métodos hf, dft, mp2, coupled cluster y mpn, para n mayores
de 2; para los dos últimos métodos se hace uso de la teoría 
Tensor Contraction Engine, para solucionar la ecuación Schrödinger 
(mayor informarcion http://www.nwchem-sw.org/index.php/Release62:NWChem_Documentation)."
  echo ''
fi

echo "Si la base que necesita no se encuentra implementada, ya sea 
por formato o otro inconveniente, la puede descargar de la siguiente
pagina web en el formato del programa deseado: 
https://bse.pnl.gov/bse/portal."
echo ' '

echo "Tipo Método (hf=0/dft=1/mp(mp2)=2/cc(ccsdt)=3/sm(pm3)=4)"
read tm
if [[ $tm -ne 0 ]]; then
echo "Método"
read m
fi
if [[ $programa -eq 1 ]]; then
echo "Desea poner una base para cada atomo (si:0/no:1)"
read sino
 if [[ $sino -eq 0 ]]; then
  echo 'Numero atomos diferentes' 
  read nt
  echo 'Atomos (Separados por espacio)'
  read  -A at
  nat=`echo "1-$nt" | bc -l`
 echo "Tipo de bases en el orden de los atomos (separadas por espacio,
sto=0/pople=1/cc(cc-pvdz)=2/sp(sbkjl)=3)"
 read -A tb
  echo 'Bases en el orden de los atomos (Separados por espacio)'
  read  -A ba
 else 
 echo "Tipo base (sto=0/pople=1/cc(cc-pvdz)=2/sp(sbkjl)=3)"
 read tb
  echo 'Base'
  read b
 fi
fi

if [[ $programa -eq 2 ]]; then
 echo "Tipo base (sto=0/pople=1/cc(cc-pvdz)=2/sp(sbkjl)=3)"
 read tb
 echo "Base"
 read b
fi

if [[ $programa -eq 0 ]]; then
 echo "Base"
 read b
fi

echo "Carga"
read carga
echo "Multiplicidad"
read spin
echo "Cantidad procesadores para el cálculo"
read np

if [[ $programa -eq 2  ]]; then #GAMESS
rm -rf launch_g.sh
integer Nconf Natom i j l1 l2 k
Nconf=`grep "Configuration" $1 | wc -l | awk '{print $1-1}'`
Natom=`sed -n '1p' $1`
i=1
k=1
j=4+Natom 
rm -f rungauss
while (( i <= Nconf ))
 do
 titulo=`sed -n ''$j'p' $1`
 l1=j+1
 l2=l1+Natom-1
 j=j+Natom+2 
 print ' $SYSTEM MWORDS=16 $END' > s$k.inp
 print ' $CONTRL RUNTYP=OPTIMIZE $END' >> s$k.inp
 print ' $STATPT HSSEND=.TRUE. $END' >> s$k.inp
 print " \$CONTRL MAXIT=100 ICHARG=$carga MULT=$spin \$END" >> s$k.inp

 if [[ $tm -eq "hf" ]]; then
  if [[ $carga = 0 && $spin = 1 ]]; then
   print ' $CONTRL SCFTYP=RHF $END' >> s$k.inp
  else
   print ' $CONTRL SCFTYP=UHF $END' >> s$k.inp 
  fi
 fi

 if [[ $tm -eq 1 ]]; then
  if [[ $carga = 0 && $spin = 1 ]]; then
   print ' $CONTRL SCFTYP=RHF $END' >> s$k.inp
   print " \$DFT  DFTTYP=$m \$END" >> s$k.inp 
  else
   print ' $CONTRL SCFTYP=UHF $END' >> s$k.inp 
   print " \$DFT  DFTTYP=$m \$END" >> s$k.inp
  fi
 fi

 
 if [[ $tm -eq 2 ]]; then
  if [[ $carga = 0 && $spin = 1 ]]; then
   print ' $CONTRL SCFTYP=RHF $END' >> s$k.inp
   print ' $CONTRL MPLEVL=2 $END' >> s$k.inp
  else
   print ' $CONTRL SCFTYP=UHF $END' >> s$k.inp 
   print ' $CONTRL MPLEVL=2 $END' >> s$k.inp
  fi
 fi

 if [[ $tm -eq 3 ]]; then
  if [[ $carga = 0 && $spin = 1 ]]; then
   print ' $CONTRL SCFTYP=RHF $END' >> s$k.inp
   print " \$CONTRL CCTYP=$m \$END" >> s$k.inp
  else
   print ' $CONTRL SCFTYP=ROHF $END' >> s$k.inp 
   print " \$CONTRL CCTYP=CCSD \$END" >> s$k.inp
  fi
 fi

 if [[ $tm -eq 4 ]]; then
   print " \$BASIS GBASIS=$m \$END" >> s$k.inp
 fi

 if [[ $tm -ne 4 ]]; then
  if [[ $tb -eq 0 ]]; then  
   ng=`echo $b | cut -c5-6 | sed 's/g//'`
   print " \$BASIS GBASIS=sto NGAUSS=$ng \$END" >> s$k.inp
  fi

  if [[ $tb -eq 1 ]]; then
   ngv=`echo $b | cut -c3-5 | sed 's/g//' | sed 's/+//' | sed 's/*//'`
   ng=`echo $b | cut -c1`
   print " \$BASIS GBASIS=N$ngv NGAUSS=$ng \$END" >> s$k.inp
   fd=`echo $b | sed 's/[0-9]//g' | sed 's/[a-z]//g' | sed 's/-//' | sed 's/*//' | sed 's/*//' | wc -c`
   fp=`echo $b | sed 's/[0-9]//g' | sed 's/[a-z]//g' | sed 's/-//' | sed 's/+//' | sed 's/+//' | wc -c` 

    if [[ $fd == 2 && $fp == 2 ]]; then
     print ' $BASIS DIFFSP=.TRUE. $END' >> s$k.inp
     print ' $BASIS NDFUNC=1 $END' >> s$k.inp
    fi 
    if [[ $fd == 3 && $fp == 3 ]]; then
     print ' $BASIS DIFFSP=.TRUE. DIFFS=.TRUE. $END' >> s$k.inp
     print ' $BASIS NDFUNC=1 NPFUNC=1 $END' >> s$k.inp
    fi 
    if [[ $fd == 1 && $fp == 2 ]]; then
     print ' $BASIS NDFUNC=1 $END' >> s$k.inp
    fi 
    if [[ $fd == 2 && $fp == 1 ]]; then
     print ' $BASIS DIFFSP=.TRUE. $END' >> s$k.inp
    fi 
    if [[ $fd == 1 && $fp == 3 ]]; then
     print ' $BASIS NDFUNC=1 NPFUNC=1 $END' >> s$k.inp
    fi 
    if [[ $fd == 3 && $fp == 1 ]]; then
     print ' $BASIS DIFFSP=.TRUE. DIFFS=.TRUE. $END' >> s$k.inp
    fi 
    if [[ $fd == 2 && $fp == 3 ]]; then
     print ' $BASIS DIFFSP=.TRUE. $END' >> s$k.inp
     print ' $BASIS NDFUNC=1 NPFUNC=1 $END' >> s$k.inp
    fi 
    if [[ $fd == 3 && $fp == 2 ]]; then
     print ' $BASIS DIFFSP=.TRUE. DIFFS=.TRUE. $END' >> s$k.inp
     print ' $BASIS NDFUNC=1 $END' >> s$k.inp
    fi 

  fi 

  if [[ $tb -eq 2 ]]; then
   acc=`echo $b | cut -c1 | tr '[A-Z]' '[a-z]'`
   ccc=`echo $b | cut -c9 | tr '[A-Z]' '[a-z]'`
    if [[ $acc == "a" && $ccc == "c" ]]; then 
      cc=`echo $b | cut -c5-6 | tr '[A-Z]' '[a-z]'`
      ccn=`echo $b | cut -c11 | tr '[A-Z]' '[a-z]'`
      print " \$BASIS GBASIS=$acc$cc$ccn$ccc \$END" >> s$k.inp
    fi

     if [[ $acc == "a" && $ccc != "c" ]]; then 
      cc=`echo $b  | cut -c5-6 | tr '[A-Z]' '[a-z]'`
      ccn=`echo $b | cut -c10 | tr '[A-Z]' '[a-z]'`
      print " \$BASIS GBASIS=$acc$cc$ccn \$END" >> s$k.inp
     fi
    ccc=`echo $b | cut -c5 | tr '[A-Z]' '[a-z]'`
    if [[ $acc != "a" && $ccc == "c" ]]; then 
      cc=`echo $b | cut -c1-2 | tr '[A-Z]' '[a-z]'`
      ccn=`echo $b | cut -c7 | tr '[A-Z]' '[a-z]'`
      print " \$BASIS GBASIS=$cc$ccn$ccc \$END" >> s$k.inp
    fi

    if [[ $acc != 'a' && $ccc != 'c' ]]; then 
      cc=`echo $b | cut -c1-2 | tr '[A-Z]' '[a-z]'`
      ccn=`echo $b | cut -c6 | tr '[A-Z]' '[a-z]'`
      print " \$BASIS GBASIS=$cc$ccn \$END" >> s$k.inp
    fi
     print " \$CONTRL ISPHER=1 \$END" >> s$k.inp
   fi 

  if [[ $tb -eq 3 ]]; then
    print " \$BASIS GBASIS=$b \$END" >> s$k.inp
    print " \$CONTRL ECP=$b \$END" >> s$k.inp
  fi

fi
 
 print ' $DATA' >> s$k.inp
 print 'Molecule specification' >> s$k.inp
 echo 'C1' >> s$k.inp
 sed -n ''$l1','$l2'p' $1 >> s$k.inp
 print ' $END' >> s$k.inp

 i=i+1
 echo "gamess 01 $np s$k.nw > s$k.log" >> launch_g.sh
 k=k+1
done
fi


if [[ $programa -eq 1 ]]; then #NWCHEM
rm -rf launch_g.sh
integer Nconf Natom i j l1 l2 k
Nconf=`grep "Configuration" $1 | wc -l | awk '{print $1-1}'`
Natom=`sed -n '1p' $1`
i=1
k=1
j=4+Natom 
rm -f rungauss
while (( i <= Nconf ))
 do
 titulo=`sed -n ''$j'p' $1`
 l1=j+1
 l2=l1+Natom-1
 j=j+Natom+2
 
 echo -e "echo \n\nstart s$k \n" > s$k.nw
 echo -e "scratch_dir /home/scratch/nwchem/$USER \n" >> s$k.nw
 echo -e "memory stack 500 mb heap 100 mb global 1000 mb\n" >> s$k.nw
 echo -e "title \"$titulo\" \ncharge $carga \n" >> s$k.nw
 echo -e "geometry units angstroms print xyz autosym" >> s$k.nw 
 sed -n ''$l1','$l2'p' $1 >> s$k.nw
 echo -e "end\n" >> s$k.nw

 if [[ $sino -eq 0 ]]; then
 echo "basis" >> s$k.nw

  for iii in {0..$nat}
  do
   echo " ${at[$iii]}  library ${ba[$iii]}" >> s$k.nw
  done

  echo -e "end \n" >> s$k.nw
  jjj=0

  for iii in {0..$nat}
  do
   if [[ ${tb[$iii]} -eq 3 ]]; then
    jjjj=`echo $jjj+1 | bc -l` 
   fi
  done

  if [[ $jjjj -gt 0 ]]; then
    echo "ecp" >> s$k.nw
 
  for iii in {0..$nat}
   do
    if [[ ${tb[$iii]} -eq 3 ]]; then  
      echo " ${at[$iii]}  library ${ba[$iii]}" >> s$k.nw
    fi
   done 
  
  echo -e "end \n" >> s$k.nw
  fi

 else
  if [[ $tb -eq 3 ]]; then
   echo -e "basis \n * library $b \nend \n" >> s$k.nw
   echo -e "ecp \n * library $b \nend \n" >> s$k.nw
  else
   echo -e "basis \n * library $b \nend \n" >> s$k.nw
  fi 
 fi

#MODULO DE DFT 
 if [[ $tm -eq 1 ]]; 
 then
 print "dft \n xc $m \n mult $spin \nend \n" >> s$k.nw
 echo -e "task dft optimize" >> s$k.nw
 echo -e "task dft freq\n" >> s$k.nw 
 fi

#MODULO DE MULTIPLICIDAD
 if [[ $tm -ne 1 ]]; then
  if [[ $spin = 1 ]]; then
  print "scf\n SINGLET" >> s$k.nw
  fi
  if [[ $spin = 2 ]]; then
  print "scf\n DOUBLET" >> s$k.nw
  fi
  if [[ $spin = 3 ]]; then
  print "scf\n TRIPLET" >> s$k.nw 
  fi
  if [[ $spin = 4 ]]; then
  print "scf\n QUARTET" >> s$k.nw
  fi
  if [[ $spin = 5 ]]; then
  print "scf\n QUINTET" >> s$k.nw
  fi
  if [[ $spin = 6 ]]; then
  print "scf\n SEXTET" >> s$k.nw
  fi
  if [[ $spin = 7 ]]; then
  print "scf\n SEPTET" >> s$k.nw
  fi
  if [[ $spin = 8 ]]; then
  print "scf\n OCTET" >> s$k.nw
  fi
  if [[ $spin = 9 ]]; then
  print "scf\n NOPEN" >> s$k.nw
  fi

#MODULO DE HF
 if [[ $tm -eq 0 ]]; then 
 print "end\n" >> s$k.nw
 echo -e "task scf optimize" >> s$k.nw
 echo -e "task scf freq\n" >> s$k.nw 
 fi

#MODULO DE MP 
  if [[ $tm -eq 2 ]]; then
  if [[ $spin != 1 ]]; then 
  print " uhf\nend \n" >> s$k.nw
  else
  print "end\n" >> s$k.nw
  fi
   if [[ $m == "mp2" || $m == "MP2" ]]; then
    print "mp2\n # Exclude core electrons from MP treatment \n freeze core atomic \nend \n" >> s$k.nw
    echo -e "task mp2 optimize" >> s$k.nw
    echo -e "task mp2 freq\n" >> s$k.nw 
   else
    print "tce\n $m \n # Exclude core electrons from MP treatment \n freeze core atomic \nend \n" >> s$k.nw
    echo -e "task tce optimize" >> s$k.nw
    echo -e "task tce freq\n" >> s$k.nw 
   fi
  fi

#MODULO DE CC 
 if [[ $tm -eq 3 ]]; then
  if [[ $spin != 1 ]]; then 
  print " uhf\nend \n" >> s$k.nw
  else
  print "end\n" >> s$k.nw
  fi
 print "tce\n $m \n # Exclude core electrons from CC treatment \n freeze core atomic \nend \n" >> s$k.nw
 echo -e "task tce optimize" >> s$k.nw
 echo -e "task tce freq\n" >> s$k.nw 
 fi

 fi
 i=i+1
 echo "mpirum -np $np nwchem s$k.nw" >> launch_g.sh
 k=k+1
done
 
fi

if [[ $programa -eq 0 ]]; then #GAUSSIAN

integer Nconf Natom i j l1 l2 k
Nconf=`grep "Configuration" $1 | wc -l | awk '{print $1-1}'`
Natom=`sed -n '1p' $1`
i=1
k=1
j=4+Natom 
rm -f rungauss
while (( i <= Nconf ))
 do
 titulo=`sed -n ''$j'p' $1`
 l1=j+1
 l2=l1+Natom-1
 j=j+Natom+2
 if [[ $tm -eq 0 ]]; then 
  echo -e "%chk=s$k.chk \n %nproc=$np \n %mem=200mw \n #p hf/$b opt freq\n\n $titulo\n\n$carga $spin" > s$k.com
 else 
  echo -e "%chk=s$k.chk \n %nproc=$np \n %mem=200mw \n #p $m/$b opt freq\n\n $titulo\n\n$carga $spin" > s$k.com 
 fi
 sed -n ''$l1','$l2'p' $1 >> s$k.com
 echo -e "\n\n" >> s$k.com
 i=i+1
 echo "g03 s$k.com" >> launch_g.sh
 k=k+1
done

fi
