! ifort -r8 -O3 dct2d.f90
! WRAPIT -in dct2df90.stub dct2d.f90 ! in: use intel compiler 
!  Bertrand DENIS MWR 2002 130:1812-1829
      subroutine dct2d(mxx,nyy,nspt,kd,dx,fin,spt1d,length1d)
!      integer, parameter :: mxx=201,nyy=101  ! nyy should be <= mxx
!     nspt has to be min(mxx,nyy)-1
      integer :: mxx,nyy,nspt,kd
      real    :: fin(mxx,nyy),dx
      ! min(mxx,nyy) to define spt2d
      real    :: spt1d(nspt),length1d(nspt)

      real    :: bm(mxx),bn(nyy)
      real    :: fin2(mxx,nyy),fout(mxx,nyy),spt2d(mxx,nyy)
      real    :: pi
!      integer :: kd   ! interval for accumulated variance spetra, 1,2,...
      integer :: i,j

!      kd=4

      pi=3.1415926

      call betam_betan(bm,bn,mxx,nyy)
!       print*,bm(1:2)  !,bm(1:2)**2.0*mxx !values of 1,2
!       print*,bn(1:2)  !,bn(1:2)**2.0*nyy
!case 1
!      fin = 2.0
!case 2 see WMR 2002 paper, maximum spetra is at k=p, here p is 4
!      do i=1,mxx
!      do j=1,nyy
!         fin(i,j)=cos(20*pi*(i-1.0)/mxx) 
!      end do
!      end do
!case 3
!      do i=1,mxx
!      do j=1,nyy
!         fin(i,j)=cos(20*pi*(i-1.0)/mxx)+sin(40.0*pi*(j-1.0)/nyy)
!      end do
!      end do
!
!      print*,fin(:,1)
      call dst(bm,bn,mxx,nyy, fin,fout)
!      print*,spt2d
! 2D variance
      spt2d= fout*fout/(1.0*mxx*nyy)
!      print*,"sum(spt2d)",sum(spt2d),spt2d(1,1)

! Recover oringinal input fin
!      fin2=fout
!      call dst_inv(bm,bn,mxx,nyy, fin2, fout)
!      print*,fout
!      print*,sum(fout-fin)/mxx/nyy

      call  two2one_dct(mxx,nyy,spt2d,nspt,kd,spt1d,length1d)
      open(12,file="spt1d.txt")
      do i=1,nspt,kd
         write(12,*)i,length1d(i)*dx,spt1d(i)
      end do
      end subroutine dct2d 

       subroutine dst(bm,bn,mx,ny,fin, fout)
       real bm(0:mx-1), bn(0:ny-1)
       real fin(0:mx-1,0:ny-1),fout(0:mx-1,0:ny-1),spt2d(0:mx-1,0:ny-1)

       real pi

       integer i,j,im,jn

       pi=3.1415926

       fout=0.0
       spt2d=0.0
       do im=0,mx-1
       do jn=0,ny-1
          do i=0,mx-1
          do j=0,ny-1
             fout(im,jn)= fout(im,jn) + fin(i,j)*cos(im*pi*(i+0.5)/mx)*cos(jn*pi*(j+0.5)/ny)
          end do
          end do
             fout(im,jn)= bm(im)*bn(jn)*fout(im,jn)
       end do
       end do
       
       end  subroutine dst 
       
       subroutine dst_inv(bm,bn,mx,ny,fin, fout)
       real bm(0:mx-1), bn(0:ny-1)
       real fin(0:mx-1,0:ny-1),fout(0:mx-1,0:ny-1)

       real pi

       integer i,j,im,jn

       pi=3.1415926

       fout=0.0
       do i=0,mx-1
       do j=0,ny-1
          do im=0,mx-1
          do jn=0,ny-1
            fout(i,j)= fout(i,j) + &
                     bm(im)*bn(jn)*fin(im,jn)*cos(im*pi*(i+0.5)/mx)*cos(jn*pi*(j+0.5)/ny)
          end do
          end do
       end do
       end do

       end  subroutine dst_inv


       subroutine betam_betan(bm,bn,mx,ny)
       real bm(0:mx-1), bn(0:ny-1)

       do i=0,mx-1
          if (i.eq.0) bm(i)=sqrt(1.0/mx)
          if (i.ge.1) bm(i)=sqrt(2.0/mx)
       end do

       do i=0,ny-1
          if (i.eq.0) bn(i)=sqrt(1.0/ny)
          if (i.ge.1) bn(i)=sqrt(2.0/ny)
       end do

       end subroutine betam_betan


       subroutine two2one_dct(mx,ny,fin,md,kd,spt1d,length1d)

       integer mx,ny,md
       real  fin(0:mx-1,0:ny-1), spt1d(md),length1d(md)
     
       integer im,jn,k,kd
!      md: minimum dimension-1
!      kd: interval for accumulated variace
!       kd=5 
!       print*,fin(0:5,0:1) 
       spt1d=0.0
       do k=1,md,kd
          alpha1=1.0*k/(md+1.0)
          alpha2=1.0*(k+kd)/(md+1.0)
          do im=0,mx-1
          do jn=0,ny-1
!             if(im.le.2*(k+1).or.jn.le.2*(k+1))then
               alpha=sqrt(1.0*im*im/mx/mx+1.0*jn*jn/ny/ny)
                if(alpha.ge.alpha1 .and. alpha.lt.alpha2)then
                   if(im.eq.0.or.jn.eq.0)then
                   spt1d(k)= spt1d(k)+fin(im,jn)/2.0
                   else
                   spt1d(k)= spt1d(k)+fin(im,jn)
                   end if
!                  print*,k,im,jn 
                end if
!              end if
          end do
          end do
!          length1d(k)=2.0/alpha1
          length1d(k)=2.0/(0.5*(alpha1+alpha2))
!          if(k.le.100)then
!          print*,k,"length and spetra= ",length1d(k), spt1d(k)
!          end if
        end do 

       end  subroutine two2one_dct 
