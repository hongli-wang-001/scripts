program pmbufr
!
! read all observations out from prepbufr. 
! read bufr table from prepbufr file
!
 implicit none

 integer, parameter :: mxmn=35, mxlv=255, maxcnt=5000, pm_limit=-15
 integer :: ireadmg,ireadsb,idate,nmsg,ntb,nsubset,nlvl
 integer :: i,j,k,NARG,iret_hdr,iret_ob,interval,cyc_year,cyc_mon,cyc_day,cyc_hr
 integer :: dhr,dhr_new,hr_new,idate_new,cyc_new_day,iret_hdr1,iret_ob1
 integer :: unit_table,unit_out,unit_site,unit_var,cnt_site,cnt_rec 

 character(80):: hdstr='SID XOB YOB DHR TYP T29 SQN PROCN RPT CAT TYPO TSIG'
 character(80):: obstr='TPHR QCIND COPOPM'
 character(8) :: subset,cyc_str,c_sid
 character(10):: idate_str,idate_new_str,file_var,cnt_str,site_loc,cyc,state
 character(2) :: hr_new_str,day_new_str,f_str 
 character(20),allocatable,dimension(:)::sta_id,sta_id1,var,var1,stat,mon,day,year,hour
 character(20) :: s_id

 real*8,parameter:: bmiss=10e10
 real(8),allocatable,dimension(:):: lat, lon, elv,val
 real(8) :: hdr(mxmn),obs(mxmn,mxlv)
 real(8) :: rstation_id

 ! use sta_id*20 to find out station lat/lon from station file, 
 ! but only c_sid*8 is encoded to bufr file since SID is 64 bit in bufr table
 equivalence(rstation_id,c_sid)


! read arguments

 NARG=IARGC()
 if (NARG==5) then   

    call getarg(1,site_loc)
    call getarg(2,cnt_str)
    read(cnt_str,*) cnt_site

    call getarg(3,file_var)
    call getarg(4,cnt_str)
    read(cnt_str,*) cnt_rec

    call getarg(5,cyc)
    write(*,'(a55,a10,3x,a10,3x,i7,3x,a10,i7)') 'cycle, file_var, cnt_record, site_loc, cnt_site:  ', cyc,file_var,cnt_rec, site_loc, cnt_site
    
 else
    write(*,*)
    write(*,*) 'Wrong arguments'
    write(*,*)
    call exit(2)
 endif


! open and read station file
 unit_site=10
 allocate(sta_id(cnt_site))
 allocate(var(cnt_site))
 allocate(stat(cnt_site))
 allocate(lat(cnt_site))
 allocate(lon(cnt_site))
 allocate(elv(cnt_site))
 open (unit = unit_site, file = site_loc,status = 'old')
 do i = 1,cnt_site  
      read(unit_site,*) sta_id(i), var(i), stat(i),lat(i),lon(i) 
      !write(*,'(a10,i5,4x,3a20,3f10.2)')  'site:  ', i,sta_id(i),var(i),stat(i),lat(i),lon(i)
 end do 


! open and read pm2.5 ascii file
 unit_var=11
 allocate(mon(cnt_rec))
 allocate(day(cnt_rec))
 allocate(year(cnt_rec))
 allocate(hour(cnt_rec))
 allocate(sta_id1(cnt_rec))
 allocate(var1(cnt_rec))
 allocate(val(cnt_rec))
 open (unit = unit_var, file = file_var,status = 'old')
 do i = 1,cnt_rec 
      read(unit_var,*) mon(i),day(i),year(i),hour(i), sta_id1(i),var1(1),val(i)
      !write(*,'(a10,i5,4x,4a7,a20,a7,f15.10)') 'obs  ', i, mon(i),day(i), year(i), hour(i), sta_id1(i),var1(1),val(i)
 end do




! open output files and write bufr table to output files
 
 unit_out=12
 unit_table=13
 open(unit_table,file='pm.table')
 open(unit_out,file='HourlyData_'//cyc//'.bufr',action='write',form='unformatted')
 call openbf(unit_out,'OUT',unit_table)

 !call maxout(20000) ! 20k byte any message size
 call datelen(10)

 ! unchaged variables for each profile at current time
 subset="ANOWPM"
 read(cyc,*)  idate    ! string to integer



 ! hdstr='SID XOB YOB DHR TYP T29 SQN PROCN RPT CAT TYPO TSIG'
 ! obstr='TPHR QCIND COPOPM'

! encode each record to bufr
 do j=1, cnt_rec

  hdr=bmiss
  obs=bmiss; 

  do k=1, cnt_site

    if (trim(sta_id1(j)) == trim(sta_id(k))) then
      hdr(2)=lon(k)
      hdr(3)=lat(k)
      state=stat(k)
      exit 
    endif

  enddo

  c_sid=sta_id1(j)(1:8); hdr(1)=rstation_id; hdr(4)=0.0; hdr(5)=102; hdr(10)=6.0! Single level report
  write(*,'(a10,2i6,i12,2x,a20,a10,2x,a10,2x,2f10.2,f25.12)') 'j,k= ',j,k,idate,sta_id1(j),c_sid,state,hdr(3),hdr(2),val(j)

  nlvl=1
  obs(1,nlvl)=-1.0; 

  if (val(j) > 0.0) then
    obs(2,nlvl)=0.0
    obs(3,nlvl)=val(j)*1e-9  ! "UG/M3" to KG/(M3)
    write(*,'(a25,f5.2,f25.12)') 'goodObs',obs(2,nlvl),obs(3,nlvl)
  else if (val(j) > pm_limit) then
    obs(2,nlvl)=0.0
    obs(3,nlvl)=0.0  
    write(*,'(a25,f5.2,f25.12)') 'limitObs',obs(2,nlvl),obs(3,nlvl)
  else   
    obs(2,nlvl)=bmiss
    obs(3,nlvl)=bmiss
    write(*,'(a25,f5.2,f25.12)') 'badObs',obs(2,nlvl),obs(3,nlvl)
  endif

  call openmb(unit_out,subset,idate)
  call ufbint(unit_out,hdr,mxmn,1,iret_hdr,hdstr)
  call ufbint(unit_out,obs,mxmn,nlvl,iret_ob,obstr)
  call writsb(unit_out)

 enddo

 call closbf(unit_out)
 close(unit_site)
 close(unit_var)



end program

