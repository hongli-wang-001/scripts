program pmbufr
!
! read all observations out from prepbufr. 
! read bufr table from prepbufr file
!
 implicit none

 integer, parameter :: mxmn=35, mxlv=255
 integer :: ireadmg,ireadsb,idate,nmsg,ntb,nsubset
 integer :: i,j,k,NARG,iret_hdr,iret_ob,interval,cyc_year,cyc_mon,cyc_day,cyc_hr
 integer :: dhr,dhr_new,hr_new,idate_new,cyc_new_day,iret_hdr1,iret_ob1
 integer :: unit_table,unitOut,unit_out(24),unit_in 

 character(80):: hdstr='SID XOB YOB DHR TYP T29 SQN PROCN RPT CAT TYPO TSIG'
 character(80):: obstr='TPHR QCIND COPOPM'
 character(8) :: subset,c_sid,cyc_interv
 character(10):: idate_str,idate_new_str
 character(2) :: hr_new_str,day_new_str,f_str 

 real(8) :: hdr(mxmn),obs(mxmn,mxlv)
 real(8) :: rstation_id,typ
 equivalence(rstation_id,c_sid)


! read arguments

 NARG=IARGC()
 if (NARG==1) then ! get cycle 
    call getarg(1,cyc_interv)
    read(cyc_interv,*) interval
    write(*,'(a20,i4)') 'cycle interval=',interval
 else
    write(*,*)
    write(*,*) 'Usage: wrtrad_latency.x  yyyymmddhh beginhr endhr'
    write(*,*)
    call exit(2)
 endif

 unit_in=10
 unit_table=24
 do i=1,24,interval
   unit_out(i)=50+i
 enddo

! dump bufr table
 !open(unit_table,file='pm.table')
 !open(unit_in,file='pm.bufr',form='unformatted',status='old')
 !call openbf(unit_in,'IN',unit_in)
 !call dxdump(unit_in,unit_table)


! open input and output files and write bufr table to output files

 open(unit_in,file='pm.bufr',form='unformatted',status='old')
 call openbf(unit_in,'IN',unit_in)

 open(unit_table,file='pm.table')
 do i=1,24,interval
   write(f_str,'(i2.2)') i
   open(unit_out(i),file='pm25.bufr.'//f_str,action='write',form='unformatted')
   call openbf(unit_out(i),'OUT',unit_table)
 enddo

 !call maxout(20000) ! 20k byte any message size
 call datelen(10)


! write to hourly data

 nmsg=0
 nsubset=0
 msg_report: do while (ireadmg(unit_in,subset,idate) == 0)
   nmsg=nmsg+1
   ntb = 0
   write(*,*)
   write(*,'(3a,i10)') 'msgtype=',subset,' cycle time =',idate

   write (idate_str,'(i10)')  idate    ! integer to string
   read(idate_str(1:4),*)  cyc_year    ! string to integer
   read(idate_str(5:6),*)  cyc_mon     ! string to integer
   read(idate_str(7:8),*)  cyc_day     ! string to integer
   read(idate_str(9:10),*) cyc_hr      ! string to integer

   sb_report: do while (ireadsb(unit_in) == 0)
     ntb = ntb+1
     nsubset=nsubset+1
     call ufbint(unit_in,hdr,mxmn,1   ,iret_hdr,hdstr)
     call ufbint(unit_in,obs,mxmn,mxlv,iret_ob,obstr)      !single level report, ireb_ob=1 

     dhr=int(hdr(4))
     hr_new=cyc_hr+dhr

     if (hr_new==24) then
        cyc_new_day=cyc_day+1
        write(day_new_str,'(i2.2)') cyc_new_day
        hr_new_str="00"
        idate_new_str=idate_str(1:6)//day_new_str//hr_new_str
     else if (hr_new <24 .and. hr_new>0) then
        write(hr_new_str,'(i2.2)') hr_new   ! integer to str, width 2 with zeros at the left
        idate_new_str=idate_str(1:8)//hr_new_str
     else
        write(*,*) 'Subset hour is wrong'
        call exit(2)
     endif

     read(idate_new_str,*) idate_new
     dhr_new=0
     hdr(4)=dhr_new
     unitOut=unit_out(hr_new)
     write(*,*) "subset_new ",dhr,idate_new,unitOut
   
     call openmb(unitOut,subset,idate_new) 
     call ufbint(unitOut,hdr,mxmn,1,iret_hdr1,hdstr)
     call ufbint(unitOut,obs,mxmn,iret_ob,iret_ob1,obstr)
     call writsb(unitOut)

            
   enddo sb_report

   write(*,*)  'message ',nmsg, '  total subset ',ntb
   enddo msg_report

   write(*,*)  'total message ',nmsg, 'total subset ',nsubset

 call closbf(unit_in)
 do i=1,24,interval
   call closbf(unit_out(i))
 enddo

end program

