#!/bin/ksh
#set -euax

log="log.20190807.tm00"
grep subset= ${log} | awk '$5==-1{print $0}' |wc
grep subset= ${log} | awk '$5==-2{print $0}' |wc
grep subset= ${log} | awk '$5==-3{print $0}' |wc
grep subset= ${log} | awk '$5==-4{print $0}' |wc
grep subset= ${log} | awk '$5==-5{print $0}' |wc
grep subset= ${log} | awk '$5==-6{print $0}' |wc
grep subset= ${log} | awk '$5==-7{print $0}' |wc
grep subset= ${log} | awk '$5==-8{print $0}' |wc
grep subset= ${log} | awk '$5==-9{print $0}' |wc
grep subset= ${log} | awk '$5==-10{print $0}' |wc
grep subset= ${log} | awk '$5==-11{print $0}' |wc
grep subset= ${log} | awk '$5==0{print $0}' |wc
grep subset= ${log} | awk '$5==1{print $0}' |wc
grep subset= ${log} | awk '$5==2{print $0}' |wc
grep subset= ${log} | awk '$5==3{print $0}' |wc
grep subset= ${log} | awk '$5==4{print $0}' |wc
grep subset= ${log} | awk '$5==5{print $0}' |wc
grep subset= ${log} | awk '$5==6{print $0}' |wc
grep subset= ${log} | awk '$5==7{print $0}' |wc
grep subset= ${log} | awk '$5==8{print $0}' |wc
grep subset= ${log} | awk '$5==9{print $0}' |wc
grep subset= ${log} | awk '$5==10{print $0}' |wc
grep subset= ${log} | awk '$5==11{print $0}' |wc
grep subset= ${log} | awk '$5==12{print $0}' |wc

