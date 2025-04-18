       implicit none
       integer, parameter :: mix=2000,mjy=350

       integer :: iflaga(mix,mjy),iflagb(mix,mjy)
       real    :: vb(mix,mjy,2)

       real    :: wk(6),wk1(200),wk2(200)
       INTEGER :: Reason,i,j,k,ii,jj,nsam,nn

       real, allocatable :: xobs(:,:)
       real    :: slat,slon
       real    :: corr,cx,xa,bias,rms

       slat=20.0
       slon=100.0

       iflaga=0
       iflagb=0
       vb=0.0
       nsam=0

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
      nsam=nsam+1
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
   close(12)
   print*,"Reading 1 Finished: ",nsam
   allocate(xobs(nsam,6))
        open(12,file="vobs.fhhh.txt")
do i=1,nsam
   READ(12,*,IOSTAT=Reason) xobs(i,:)  
END DO
   close(12)
   print*,"Reading 2 Finished: ",nsam

   do i=1,mix
   do j=1,mjy

      nn = 0
      do k=1,nsam
      wk=xobs(k,:)
      !print*,wk
      if(wk(5).lt.1.0)cycle
      if(wk(1).gt.slat)then
        if(wk(2).lt.0)wk(2)=wk(2)+360.0
        jj=int((wk(1)-slat)*10.0+0.5)
        ii=int((wk(2)-slon)*10.0+0.5)
        !print*,"ii,jj= ",ii,jj
        if(ii.gt.mix)cycle
        if(jj.gt.mjy)cycle
        if(ii.eq.i .and. jj.eq.j)then
        nn=nn+1
        wk1(nn)=wk(3)
        wk2(nn)=wk(4)
        xobs(k,5)=-99
        end if
      endif

      end do 
      if(nn.gt.0)then
      cx=-99.0
      call stat_1d(wk1,wk2,nn,cx,xa,bias,rms)
      print*,i,j,nn,bias,rms,cx
      if(abs(bias.gt.10))then
      write(12,*)slon+0.1*i,slat+0.1*j,1.0*nn,bias
      do k=1,nn
      write(12,*)wk1(k),wk2(k)
      end do
      end if
      if(nn >= 1)then
      write(11,'(2F10.5,1x,F6.1,1x,2F10.5,F10.2,F7.2)')slon+0.1*i,slat+0.1*j,1.0*nn,xa,bias,rms,cx 
      end if
      end if
   end do 
   end do

 end 

subroutine stat_1d(x,y,n,corr,xa,bias,rmse)
integer n
real x(n),y(n),xb(n),yb(n)
real xa,ya
real corr,bias,rmse
   !print*,n
   !print*,x
   !print*,y
     corr=-99.0
   if(n >= 1)then
     xa=sum(x)/n
     ya=sum(y)/n
     xb=x-xa
     yb=y-ya
     if(n >= 2) corr=sum(xb*yb)/sqrt((sum(xb*xb)*sum(yb*yb)))
     bias=ya-xa
     rmse=sqrt(sum((x-y)*(x-y))/n)

   end if
   !print*,corr
   return
end subroutine stat_1d




