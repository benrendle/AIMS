
!   generateList.f95   2018-10-09

!===============================================================================
!   USE: generate the list of models for AIMS based on given criteria (Teff, Xc, mHe, etc)
!   based on a list of summary files, use in conjunction with the script "readgrid" to loop
!   over summary files in a given repertory (tracks in the current grid)
!===============================================================================

! Disclaimer: the code is unoptimized and poorly written.

program generateList

implicit none
integer :: io,nlinesSum,i,nmod,narg,stat,iarg,lsummaryfile
! Solar Constants (change depending on references used)
real(8), parameter :: Lsun=3.828d33
real(8), parameter :: Rsun=6.957d10
real(8), parameter :: Msun=1.988475415d33
real(8), parameter :: Teffsun=5772d0
real(8),dimension(:),pointer :: Age,logTeff,logL,Xc,Yc,Lumg,Rhoc,Tc,R,Mconv,logg,ZXsurfstar,Diff,AlphaConv,AlphaOver,&
 &Rconv,Z0,X0,Mass,P0sum,P0recalc,Mcomix,mHe, Numax
character(256) :: nameseq,summaryfile
integer, dimension(:),pointer :: modNumber
character(256), dimension(:),pointer :: Modelname
character(4),dimension(:),pointer :: mNumber
character(60),dimension(:),pointer :: model
nlinesSum = 0

narg=command_argument_count()

if(narg.gt.1) write(*,*) 'Too many arguments!'
do iarg=1,narg
  if(iarg.eq.1) call get_command_argument(iarg,summaryfile,status=stat)
end do
lsummaryfile=len_trim(summaryfile)
nameseq=summaryfile(1:(lsummaryfile-8))
write(*,*) trim(nameseq)
write(*,*) trim(summaryfile)

! writing format from summary.f in CLES
901 format(1x,i4,22es16.5)

open (1, file = summaryfile)
do
  read(1,*,iostat=io)
  if (io/=0) exit
  nlinesSum = nlinesSum + 1
end do
close(1)
nmod =nlinesSum-3
allocate(Age(nmod),logTeff(nmod),logL(nmod),Xc(nmod),Yc(nmod),Lumg(nmod),Rhoc(nmod),Tc(nmod),&
 &R(nmod),Mconv(nmod),logg(nmod),ZXsurfstar(nmod),Diff(nmod),AlphaConv(nmod),AlphaOver(nmod),&
 &Rconv(nmod),Z0(nmod),X0(nmod),Mass(nmod),P0sum(nmod),P0recalc(nmod),Mcomix(nmod),mHe(nmod),&
 &modNumber(nmod),Modelname(nmod),mnumber(nmod), Numax(nmod),model(nmod))

open (1, file = summaryfile)
do i=1,3
 read(1,*)
end do
do i=1,nmod
 read(1,901) modNumber(i),Age(i),logTeff(i),logL(i),Xc(i),Yc(i),Lumg(i),Rhoc(i),Tc(i),&
 &R(i),Mconv(i),logg(i),ZXsurfstar(i),Diff(i),AlphaConv(i),AlphaOver(i),&
 &Rconv(i),Z0(i),X0(i),Mass(i),P0sum(i),Mcomix(i),mHe(i)
end do
close(1)

do i=1,nmod
  Numax(i)=3104d0*(Mass(i)/MSun)*((1d0/R(i))**2)*(Teffsun/(10**logTeff(i)))**0.5
  write (mnumber(i),'(I4.4)') modNumber(i)
  Modelname(i)=trim(nameseq)//'/AIMS/'//trim(nameseq)//'-atm-'//mnumber(i)
  !write (*,*) LEN(Modelname(i))
  !write (*,*) LEN(trim(Modelname(i)))
  !Modelname(i)=trim(nameseq)//'/AIMS/'//trim(nameseq)//'-atm-'//mnumber(i)//'.freq'
end do


! Generation of the list, adapt name and criteria accordingly
open (3, file ='/home/bmr135/bison/Sim2/AIMS_Gael/NGC6791_in',position='append')
!write (3,*) '/home/ADF/bmr135/Sim3/cles-19.1-Up/scripts/	.mod'
do i=1,nmod
model(i) = trim(Modelname(i))
!if((Z0(i).eq.(0.0032))) then
!  if((mHe(i).gt.(0.08)).and. ((10.0**logTeff(i)).lt. (5300d0)) .and. (Numax(i).gt.(60d0))) then
! write(3,*) trim(Modelname(i)),Mass(i),(R(i)*Rsun),(10.0**logL(i))*Lsun,Z0(i),X0(i),(Age(i)*1d-6),(10.0**logTeff(i)),(mHe(i)*Msun)&
!&,Xc(i),P0sum(i)
!  endif
!elseif((Z0(i).eq.(0.0057))) then
!  if((mHe(i).gt.(0.08)).and. ((10.0**logTeff(i)).lt. (5200d0)) .and. (Numax(i).gt.(60d0))) then
! write(3,*) trim(Modelname(i)),Mass(i),(R(i)*Rsun),(10.0**logL(i))*Lsun,Z0(i),X0(i),(Age(i)*1d-6),(10.0**logTeff(i)),(mHe(i)*Msun)&
!&,Xc(i),P0sum(i)
!  endif
!else
!  if((mHe(i).gt.(0.08)).and. ((10.0**logTeff(i)).lt. (5200d0)) .and. (Numax(i).gt.(60d0))) then
! write(3,*) trim(Modelname(i)),Mass(i),(R(i)*Rsun),(10.0**logL(i))*Lsun,Z0(i),X0(i),(Age(i)*1d-6),(10.0**logTeff(i)),(mHe(i)*Msun)&
!&,Xc(i),P0sum(i)
!  endif
!if((Xc(i).gt.(0.1)).and. (abs(X0(i)-Xc(i)).gt. (0.05d0))) then
if (mHe(i).gt.0.05) then 
 !.and.(mHe(i).lt.0.265)) then
 !write(*,*) Numax(i), 10.0**logTeff(i)
 write(3,*) model(i),Mass(i),(R(i)*Rsun),(10.0**logL(i))*Lsun,Z0(i),X0(i),(Age(i)*1d-6),(10.0**logTeff(i))&
 &,(mHe(i)*Msun),Xc(i),P0sum(i)
 endif
end do
close(3)
end program generateList
