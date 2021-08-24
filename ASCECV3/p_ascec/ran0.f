*******************************************************************************
      function ran0(idum)
      integer idum,ia,im,iq,ir,mask
      real*8 ran0,am
      parameter(ia=16807,im=2147483647,am=1.d0/im,
     $          iq=127773,ir=2836,mask=123459876)
     
      integer k
      idum=ieor(idum,mask)
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      ran0=am*idum
      idum=ieor(idum,mask)      
      return
      end
*******************************************************************************
