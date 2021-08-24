      subroutine config(natom,r)
      implicit real*8 (a-h,o-z)
      include 'param.inc'          
      dimension rflo(max_atom,3)			
      dimension r(max_atom,3)
      dimension r2(max_atom,2:max_atom)
      dimension rel(max_atom,3)
      dimension rcm2(max_mole,3)      	
      character*2 sym,element
      !common/names/iznu(max_atom),element(max_atom),sym(120),xbox,mspin
      common/elemente/element(max_atom),mspin
      common/names/iznu(max_atom),sym(120),xbox
      common/rando/idum
      common/size/Ratom(20)
      common/molecu/rcm(max_mole,3),imolec(0:max_mole),nmo,ds,dphi	
      data pi/3.141592653589793238462643d0/

!****	
      DO imo=1,nmo
       DO k=1,3
       rcm2(imo,k)=rcm(imo,k)
        DO inu=imolec(imo-1)+1,imolec(imo)
	 rflo(inu,k)=r(inu,k)
	ENDDO
       ENDDO	
      ENDDO 
      DO imo=1,nmo		!Mto. cada molecula
       icount=0
  50   continue
       DO k=1,3
	rcm(imo,k)=rcm2(imo,k)
        DO inu=imolec(imo-1)+1,imolec(imo)
	 rflo(inu,k)=r(inu,k)
	ENDDO	  
       ENDDO	 	 	
       call trans(imo,rflo,rel)
       call rotac(imo,rflo,rel)
       if(imo.eq.1)goto 51
       DO jmo=1,imo-1
  	DO inu=imolec(imo-1)+1,imolec(imo)
	 DO jnu=imolec(jmo-1)+1,imolec(jmo)
	   r2(inu,jnu)=dsqrt(
     $    (rflo(inu,1)-rflo(jnu,1))**2 +
     $    (rflo(inu,2)-rflo(jnu,2))**2 +
     $    (rflo(inu,3)-rflo(jnu,3))**2 )
          rmin=(Ratom(iznu(inu))+Ratom(iznu(jnu)))*0.7
	  if(r2(inu,jnu).lt.rmin)then
	  icount=icount+1
	  if(icount.gt.100000)then
	   write(6,*)''
	   write(6,*)'** Warning **'
	   write(6,*)'More than ',icount-1,' movements.'
           write(6,*)'Molecula',imo,' keeps violating the space of
     $molecule',jmo
           write(6,*)'Involved atoms:',inu,' and',jnu
           icount=0
	   write(6,*)'Maximum displacement parameter will be enlarged 
     $by 20% '
	   ds=ds*1.2d0	   
	   if(ds.gt.xbox)then
           ds=xbox	   
           write(6,*)'Maximum value for ds = xbox reached'
	   endif	   
     	   write(6,*)'New value for ds =',ds
	   write(6,*)'*****************'
	   write(6,*)''	   
!	   stop
	  endif
	  goto 50
	  endif	    
	 ENDDO
	ENDDO 
       ENDDO
  51   continue	 
       DO k=1,3
	rcm2(imo,k)=rcm(imo,k)
        DO inu=imolec(imo-1)+1,imolec(imo)
	 r(inu,k)=rflo(inu,k)
	ENDDO
       ENDDO
      ENDDO
!****      
      return
      end	
!********************

!********************
!**** Traslada moléculas aleatoriamente      
      subroutine trans(imo,r,rel)
      implicit real*8 (a-h,o-z)
      include 'param.inc'               
      dimension r(max_atom,3)	
      dimension rel(max_atom,3)      	
      character*2 sym
      common/names/iznu(max_atom),sym(120),xbox
      common/rando/idum
      common/molecu/rcm(max_mole,3),imolec(0:max_mole),nmo,ds,dphi	
      data pi/3.141592653589793238462643d0/
	
      DO i=1,3
       z=ran0(idum)
       if(z.lt.0.5)then
       is=-1
       else
       is=1
       endif
       s=is*ran0(idum)*ds
       rcm(imo,i)=rcm(imo,i)+s    		!Mto CM.
       if(abs(rcm(imo,i)).gt.xbox/2.d0)then
       rcm(imo,i)=rcm(imo,i)-2.d0*s
       s=-s
       endif
       DO inu=imolec(imo-1)+1,imolec(imo)
        r(inu,i)=r(inu,i)+s
	rel(inu,i)=r(inu,i)-rcm(imo,i)	   
       ENDDO
      ENDDO

      return
      end 
!********************

!********************
!**** Rota moleculas aleatoriamente (se mantiene inalterado el CM.)
      subroutine rotac(imo,r,rel)
      implicit real*8 (a-h,o-z)
      include 'param.inc'               
      dimension r(max_atom,3)	
      dimension rel(max_atom,3)
      character*2 sym
      common/names/iznu(max_atom),sym(120),xbox
      common/rando/idum
      common/molecu/rcm(max_mole,3),imolec(0:max_mole),nmo,ds,dphi	
      data pi/3.141592653589793238462643d0/

      if(ran0(idum).lt.0.5)then
      ix=1
      else
      ix=-1
      endif
      if(ran0(idum).lt.0.5)then
      iy=1
      else
      iy=-1
      endif
      if(ran0(idum).lt.0.5)then
      iz=1
      else
      iz=-1
      endif
      xa=ix*ran0(idum)*dphi
      ya=iy*ran0(idum)*dphi
      za=iz*ran0(idum)*dphi	  
      DO inu=imolec(imo-1)+1,imolec(imo)
       X=rel(inu,1)
       Y=rel(inu,2)
       Z=rel(inu,3)      
       rel(inu,1)=
     $  X*(cos(za)*cos(ya))+
     $  Y*(-cos(za)*sin(ya)*sin(xa)+sin(za)*cos(xa))+
     $  Z*(cos(za)*sin(ya)*cos(xa)+sin(za)*sin(xa))
       rel(inu,2)=
     $  X*(-sin(za)*cos(ya))+
     $  Y*(sin(za)*sin(ya)*sin(xa)+cos(za)*cos(xa))+
     $  Z*(-sin(za)*sin(ya)*cos(xa)+cos(za)*sin(xa)) 
       rel(inu,3)=
     $  X*(-sin(ya))+
     $  Y*(-cos(ya)*sin(xa))+
     $  Z*(cos(ya)*cos(xa))
       DO k=1,3
        r(inu,k)=rel(inu,k)+rcm(imo,k)
       ENDDO
      ENDDO
		  
      return
      end
!********************
