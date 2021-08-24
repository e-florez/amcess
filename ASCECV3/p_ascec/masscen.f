      subroutine masscen(natom,r)
      implicit real*8 (a-h,o-z)
      include 'param.inc'     
      dimension r(max_atom,3)	
      common/molecu/rcm(max_mole,3),imolec(0:max_mole),nmo,ds,dphi	
      !common/names/iznu(max_atom),element(max_atom),sym(120),xbox,mspin	
      common/elemente/element(max_atom),mspin	
      common/names/iznu(max_atom),sym(120),xbox	
      common/weigth/wt(max_atom)
	
      do imo=1,nmo
       do j=1,3
	rcm(imo,j)=0.
	TW=0.
	 do inu=imolec(imo-1)+1,imolec(imo)
	  rcm(imo,j)=rcm(imo,j)+r(inu,j)*wt(iznu(inu))
	  TW=TW+wt(iznu(inu))
	 enddo
	 rcm(imo,j)=rcm(imo,j)/TW
!	write(6,*)rcm(imo,j),TW
       enddo
      enddo
	
      return
      end
