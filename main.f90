	PROGRAM MAIN
		USE constants_Model
		IMPLICIT NONE
		INTEGER::I,J,KK,lag_num,time_step,flag
		INTEGER::index_XY(3,2)
		REAL(KIND=8)::X_Coor(M),Y_Coor(N)
		REAL(KIND=8),allocatable::Body_x(:),Body_y(:),AA(:,:)
		REAL(KIND=8)::U(M,N),V(M,N),P(M,N),Qu(M,N),Qv(M,N),Qp(M,N)
		REAL(KIND=8)::dU(M,N),dV(M,N)
		REAL(KIND=8)::Delta_s, Delta_t, Total_Time, t
		REAL(KIND=8),allocatable:: P_Observe(:,:)
		CHARACTER(len=80)::g,gg,str1,filename
		REAL(KIND=8),external::Numerical_Probe
		REAL(KIND=8)::START,FINISH
		!
		CALL CPU_TIME(START)
		!Get the Cartesian grid
		CALL Grid(X_Coor,Y_Coor,M,N,Length,Width)
		!
		!model
		lag_num=int(PI*2.0*R/Delta_x)
		write(*,*)lag_num
		allocate(Body_x(lag_num),Body_y(lag_num),AA(lag_num,lag_num))
		CALL Model(Body_X,Body_Y,lag_num,Delta_s)
		!
		!initializaition
		CALL Initialization(U, V, P, Qu, Qv, Qp, X_Coor,Y_Coor, M,N)
		!
		Delta_t=0.5d0*Delta_x
		Total_Time=15d0
		t=0.0d0
		KK=0
		g='.dat'
		gg='p_'
		OPEN(UNIT=11,FILE='Probe_ABC.dat')
		CALL Search_Index(X_Coor,Y_Coor, M,N, index_xy)
		allocate(P_Observe(int(Total_Time/Delta_t),3))
		!
		DO time_step=1, int(Total_Time/Delta_t)
			t=t+Delta_t
			!step-1: calculate the acoustic field without regard the body
			CALL LDDRK(U,V,P,Qu,Qv,Qp,X_coor,Y_coor,M,N,Body_x,Body_y,lag_num,Delta_x,Delta_t)
			!
			!step-2: calculate the correctional quantity for the acoustic field
			flag=time_step
			CALL Get_dUdV(dU,dV,U,V,X_Coor,Y_Coor,Body_x,Body_y,M,N,lag_num,Delta_x,Delta_s,Delta_t,flag,AA)
			!
			!STEP-3: correct the acoustic field for u, v, p
			DO I=1,M
				DO J=1,N
					U(I,J)=U(I,J)+dU(I,J)
					V(I,J)=V(I,J)+dV(I,J)
				END DO
			END DO
			!
			!Finally: output the results
			!the transient pressure
			IF((time_step==1) .or. (mod(time_step,int(0.5d0/Delta_t))==0))THEN
				write(*,*)time_step
				KK=KK+1
				if(kk/10>=1)then
     	 			write(str1,'(I2)')kk
       			else
            		write(str1,'(I1)')kk
        		end if
				write(filename,*)(trim(gg)//trim(str1)//trim(g))
				OPEN(UNIT=10,FILE=trim(filename))
				do I=1,M
					do J=1,N
						if(J==N)then
							write(unit=10,fmt='(f18.8)',advance='yes')p(I,J)
						else
	 						write(unit=10,fmt='(f18.8)',advance='no')p(I,J)
						end if
					end do
				end do
				close(unit=10)
			END IF
			!
			!the observing point: A, B, C(numerical probe)
			DO flag=1,3
				P_Observe(time_step,flag)=P(index_xy(flag,1),index_xy(flag,2))
			END DO
			!write(*,fmt="(3f18.8)")P_Observe(time_step,:)
		END DO
		DO I=1,int(Total_Time/Delta_t)
			WRITE(UNIT=11,FMT="(3F18.13)")P_Observe(I,:)
		END DO
		CLOSE(UNIT=11)
		CALL CPU_TIME(FINISH)
		write(*,*)FINISH-START
	END PROGRAM MAIN
