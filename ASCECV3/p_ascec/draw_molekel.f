*******************************************************************************
!! Draw the cluster
      subroutine drawmol(natom,r,nstep,iw,E)
      implicit real*8 (a-h,o-z)
      include 'param.inc'     
      dimension r(max_atom,3)
      character*2 sym, element
      character*8 jdum, hamil
      character*15 alias
      common/label/jdum
      common/elemente/element(max_atom),mspin
      common/names/iznu(max_atom),sym(120),xbox
      common/ruta/ia,iQ,iS2,hamil,alias

      if(ia.ne.4)then

       if(iw.eq.0)then	          	
        open(7,file='mto_'//jdum//'.xyz',status='unknown')
        write(7,*)natom       
        write(7,*)'Configuration:',nstep
        do i=1,natom
         write(7,'(1x,a2,1x,3f12.6)')sym(iznu(i)),(r(i,j),j=1,3)    
        enddo
       endif  

       if(iw.eq.1)then
        open(10,file='result_'//jdum//'.xyz',status='unknown')      
        write(10,*)natom       
        write(10,'(A15,I5,A15,F16.8)')'Configuration:',nstep,
     #' Energy(u.a.)=',E      
        do i=1,natom
         write(10,'(1x,a3,1x,3f12.6)')sym(iznu(i)),(r(i,j),j=1,3)    
        enddo
       endif        

      else      

        if(iw.eq.0)then	          	
         open(7,file='mto_'//jdum//'.xyz',status='unknown')
         write(7,*)natom       
         write(7,*)'Configuration:',nstep
         do i=1,natom
          write(7,'(1x,a2,1x,i3,1x,3f12.6)')sym(iznu(i)),iznu(i),
     $   (r(i,j),j=1,3)    
         enddo
        endif  

        if(iw.eq.1)then
         open(10,file='result_'//jdum//'.xyz',status='unknown')      
         write(10,*)natom       
         write(10,'(A15,I5,A15,F16.8)')'Configuration:',nstep,
     #' Energy(u.a.)=',E      
         do i=1,natom
          write(10,'(1x,a3,1x,i3,1x,3f12.6)')sym(iznu(i)),iznu(i),
     $   (r(i,j),j=1,3)    
         enddo
        endif
        
      endif
      return
      end
*******************************************************************************
