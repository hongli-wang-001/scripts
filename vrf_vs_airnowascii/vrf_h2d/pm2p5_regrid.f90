       implicit none
       integer, parameter :: mix=2000,mjy=350

       integer :: iflaga(mix,mjy),iflagb(mix,mjy)
       real    :: vb(mix,mjy,2)

       real    :: wk(5)
       INTEGER :: Reason,i,j,ii,jj,nsam

       real    :: slat,slon
       real    :: corr,cx

       slat=20.0
       slon=100.0

       iflaga=0
       iflagb=0
       vb=0.0


        open(12,file="vobs.fhhh.txt")
DO
   READ(12,*,IOSTAT=Reason) wk
   !print*,wk
   IF (Reason > 0)  THEN
      print*,"check if your file exist"
   ELSE IF (Reason < 0) THEN
      print*," end of file reached ..."
      exit
   ELSE
      if(wk(5).lt.1.0)cycle
      if(wk(1).gt.slat)then
        if(wk(2).lt.0)wk(2)=wk(2)+360.0
        jj=int((wk(1)-slat)*10.0+0.5)
        ii=int((wk(2)-slon)*10.0+0.5)
        !print*,"ii,jj= ",ii,jj
        if(ii.gt.mix)cycle
        if(jj.gt.mjy)cycle
        iflagb(ii,jj)=iflagb(ii,jj)+1
        vb(ii,jj,1)=vb(ii,jj,1)+wk(3)
        vb(ii,jj,2)=vb(ii,jj,2)+wk(4)
      endif
   END IF
END DO
   print*,"Reading 1 Finished"

   do i=1,mix
   do j=1,mjy
      nsam = iflagb(i,j)
      if(iflagb(i,j).gt.0)then
      nsam=iflagb(i,j)
      !print*,i,j,iflagb(i,j),vb(i,j,: )/iflagb(i,j)
      !write(11,'(I6,1X,I6,1x,I6,1F10.5)')i,j,iflagb(i,j),vb(i,j,: )/iflagb(i,j)
      !print*,i,j,nsam
      write(11,'(2F10.5,1x,F6.1,1x,2F10.5,F7.2)')slon+0.1*i,slat+0.1*j,1.0*iflagb(i,j),vb(i,j,: )/iflagb(i,j)!,cx 
      end if
   end do 
   end do

 end 

subroutine corr_1d(x,y,n,corr)
integer n
real x(n),y(n),xb(n),yb(n)
real xa,ya
   print*,n
   print*,x
   print*,y
   if(n >= 1)then
     xa=sum(x)/n
     ya=sum(y)/n
     xb=x-xa
     yb=y-ya
     corr=sum(xb*yb)/sqrt((sum(xb*xb)*sum(yb*yb)))
   else
     corr=-9999

   end if
   print*,corr
   return
end subroutine corr_1d




