program pmbufr_decode
!
! read all observations out from prepbufr. 
! read bufr table from prepbufr file
!
 implicit none

 integer, parameter :: mxmn=35, mxlv=255
 integer :: ireadmg,ireadsb,idate,nmsg,ntb,nsubset
 integer :: i,j,k,NARG,iret_hdr,iret_ob,interval,cyc_year,cyc_mon,cyc_day,cyc_hr
 integer :: dhr,dhr_new,hr_new,idate_new,cyc_new_day,iret_hdr1,iret_ob1
 integer :: unit_table,unitOut,unit_in

 character(80):: hdstr='SID XOB YOB DHR TYP T29 SQN PROCN RPT CAT TYPO TSIG'
 character(80):: obstr='TPHR QCIND COPOPM'
 character(8) :: subset,cyc_interv,c_sid
 character(10):: idate_str,idate_new_str
 character(2) :: hr_new_str,day_new_str,f_str 
 
 real(8) :: hdr(mxmn),obs(mxmn,mxlv)
 real(8) :: rstation_id,typ
 equivalence(rstation_id,c_sid)

 unit_in=10
 unit_table=24

! dump bufr table
 open(unit_table,file='pm.table.decode')
 open(unit_in,file='pm.bufr',form='unformatted',status='old')
 call openbf(unit_in,'IN',unit_in)
 call dxdump(unit_in,unit_table)


! open input and output files and write bufr table to output files

   !call maxout(20000)
   call datelen(10) 
   nmsg=0
   nsubset=0

    write(*,'(a130)') "-----------------------------------------------------------------------------------------------------------------------------------"
    write(*,'(2a10,15a10)') "Subset=","SID","Lat","Lon","DhrTime","RepType","DmpRepTyp","RepSeqNum","ProcNum","RepTime","Category","TypePoll","TimeSigf"
    write(*,*)
    write(*,'(a10,3a30)') "Obs=","TimePeriodDisplacement","QualityControlIndication","ConcentrationOfPollutant"
    write(*,'(a130)') "-----------------------------------------------------------------------------------------------------------------------------------"


   msg_report: do while (ireadmg(unit_in,subset,idate) == 0)
    nmsg=nmsg+1
    ntb = 0
    write(*,*)
    write(*,'(3a,i10)') 'msgtype=',subset,' cycle time =',idate

    sb_report: do while (ireadsb(unit_in) == 0)
     ntb = ntb+1
     nsubset=nsubset+1
     call ufbint(unit_in,hdr,mxmn,1   ,iret_hdr,hdstr)
     call ufbint(unit_in,obs,mxmn,mxlv,iret_ob,obstr)      !single level report, ireb_ob=1 

     rstation_id=hdr(1)                      ! SID XOB YOB DHR TYP T29 SQN PROCN RPT CAT TYPO TSIG
     write(*,'(2a10,4f10.2,4f20.2,f7.2,2f20.2)')  "Subset=",c_sid,hdr(3),hdr(2),hdr(4),hdr(5),hdr(6),hdr(7),hdr(8),hdr(9),hdr(10),hdr(11),hdr(12)
     do k=1,iret_ob                          ! single level report, ireb_ob=1
       write(*,'(i3,a10,f10.2,2f25.12)') k,'Obs=',(obs(i,k),i=1,3)
     enddo
            
   enddo sb_report

   write(*,*)  'message ',nmsg, '  total subset ',ntb
   enddo msg_report

   write(*,*)  'total message ',nmsg, 'total subset ',nsubset
 call closbf(unit_in)

end program
