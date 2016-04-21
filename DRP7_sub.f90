real(kind=8) FUNCTION DRP7(u0,delta_x,flag)
IMPLICIT NONE
!coeff(1,:):0-6
!coeff(2,:):-1~5
!%coeff(3,:):-2~4
!%coeff(4,:):-3~3
!%coeff(5,:):-4~2
!%coeff(6,:):-5~1
!%coeff(7,:):-6~0
!%flag: the value changes from 0 to -6
 integer::i,j,flag
 real(kind=8)::u0(7),delta_x
 real(kind=8)::coeff(7,7)
 real(kind=8):: temp_coeff1(7)=(/ -2.19228033900d0,4.74861440100d0,-5.10885191500d0,&
		&4.46156710400d0,-2.83349874100d0,1.12832886100d0,-0.20387637100d0 /)
 real(kind=8):: temp_coeff2(7)=(/ -0.20933762200d0,-1.08487567600d0,2.14777605000d0,&
		&-1.38892832200d0,0.76894976600d0,-0.28181465000d0,0.048230454000d0 /)
 real(kind=8):: temp_coeff3(7)=(/ 0.049041958000d0,-0.4688403500d0,-0.47476091400d0,&
		&1.27327473700d0,-0.51848452600d0,0.16613853300d0,-0.026369431000d0 /)
 real(kind=8):: temp_coeff4(7)=(/ -0.02084314277031176d0,0.166705904414580469d0,&
				&-0.77088238051822552d0,0.0d0,0.77088238051822552d0,&
				&-0.166705904414580469d0,0.02084314277031176d0 /)
 !character(len=80)::filename='DRP7_coeff.dat'
 real(kind=8)::start,finish
!
!call cpu_time(start) 
 do i=1,7
 	coeff(1,i)=temp_coeff1(i)
	coeff(2,i)=temp_coeff2(i)
	coeff(3,i)=temp_coeff3(i)
	coeff(4,i)=temp_coeff4(i)
 end do
! open(unit=10,file=filename) 
! do i=1,4
! 	read(unit=10,fmt=*)coeff(i,:)
! end do
! close(unit=10)
!call cpu_time(finish)
!write(*,*)'time6:',finish-start
!
 do i=1,7
	coeff(5,i)=-coeff(3,8-i)
	coeff(6,i)=-coeff(2,8-i)
	coeff(7,i)=-coeff(1,8-i)
 end do
 DRP7=0.0;
 do i=1,7
     DRP7=DRP7+coeff(1-flag,i)*u0(i);
 end do
 DRP7=DRP7/delta_x;
END FUNCTION DRP7
