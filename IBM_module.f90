 
    SUBROUTINE Normal_Vector(fnx,fny,Body_x,Body_y,lag_num)
		!-----------------introduction---------------
		!get the outer normal vector for body
		!fnx: the x component for the unit normal 
		!vector
		!fny: the y component for the unit normal
		!vector
		!--------------------------------------------
		IMPLICIT NONE
		INTEGER::I,J,lag_num
        REAL(KIND=8),parameter::PI=3.141592653d0
		REAL(KIND=8)::fnx(lag_num),fny(lag_num),Body_x(lag_num),Body_y(lag_num)
		!
		DO I=1,lag_num
			fnx(I)=dcos((I-1)*2.0d0*PI/lag_num)
			fny(I)=dsin((I-1)*2.0d0*PI/lag_num)
		END DO
	END SUBROUTINE Normal_Vector
	!
	!
	SUBROUTINE Get_Un(Un,U,V,X_Coor,Y_Coor,M,N,Delta_x,Body_x,Body_y,lag_num)
		IMPLICIT NONE
		INTEGER::I,J,K,M,N,lag_num
		REAL(KIND=8)::Un(lag_num),U(M,N),V(M,N)
		REAL(KIND=8)::X_Coor(M),Y_Coor(N),Body_x(lag_num),Body_y(lag_num)
		REAL(KIND=8)::fnx(lag_num),fny(lag_num)
		REAL(KIND=8)::coordinate(1,2),points(1,2)
		REAL(KIND=8)::Delta_x,summ
        REAL(KIND=8),parameter::PI=3.141592653d0
        REAL(KIND=8),external::Distribution
		!
		CALL Normal_Vector(fnx,fny,Body_x,Body_y,lag_num)
		!
		DO K=1,lag_num
			summ=0.0d0
			points(1,1)=Body_x(k)
			points(1,2)=Body_y(k)
			DO I=1,M
				DO J=1,N
					coordinate(1,1)=X_Coor(I)
					coordinate(1,2)=Y_Coor(J)
					summ=summ+fnx(k)*U(I,J)*Distribution(points,coordinate,Delta_x)*Delta_x**2 + &
						&fny(k)*V(I,J)*Distribution(points,coordinate,Delta_x)*Delta_x**2
				END DO
			END DO
			Un(k)=summ
		END DO
	END SUBROUTINE 
	!
	!	
	SUBROUTINE Get_dUndt(dUndt,U,V,X_Coor,Y_Coor,Body_x,Body_y,lag_num,M,N,Delta_x,Delta_t)
		IMPLICIT NONE
		INTEGER::I,J,M,N,lag_num
		REAL(KIND=8)::dUndt(lag_num),Un(lag_num)
		REAL(KIND=8)::U(M,N),V(M,N),X_Coor(M),Y_Coor(N),Body_x(lag_num),Body_y(lag_num)
		REAL(KIND=8)::Delta_x,Delta_t,Delta_s
        REAL(KIND=8),parameter::PI=3.141592653d0
		!
		CALL Get_Un(Un,U,V,X_Coor,Y_Coor,M,N,Delta_x,Body_x,Body_y,lag_num)
		!
		DO I=1,lag_num
			dUndt(I)=-Un(I)/Delta_t
		END DO
	END SUBROUTINE Get_dUndt
	!
	!
	SUBROUTINE Get_A(A,X_Coor,Y_Coor,Body_x,Body_y,M,N,lag_num,Delta_x,Delta_s,Delta_t)
		IMPLICIT NONE
		INTEGER::I,J,K,L,M,N,lag_num
		REAL(KIND=8)::A(lag_num,lag_num)
		REAL(KIND=8)::X_Coor(M),Y_Coor(N),Body_x(lag_num),Body_y(lag_num)
		REAL(KIND=8)::fnx(lag_num),fny(lag_num)
		REAL(KIND=8)::coordinate1(1,2),coordinate2(1,2),points(1,2)
		REAL(KIND=8)::delta_s,Delta_x,Delta_t,summ
        REAL(KIND=8),parameter::PI=3.141592653d0
        REAL(KIND=8),external::Distribution
		!
		CALL Normal_Vector(fnx,fny,Body_x,Body_y,lag_num)
		write(*,*)'here we are'
		DO K=1,lag_num
			coordinate1(1,1)=Body_x(K)
			coordinate1(1,2)=Body_y(k)
			DO L=1,lag_num
				coordinate2(1,1)=Body_x(L)
				coordinate2(1,2)=Body_y(L)
				summ=0.0d0
				DO I=1,M
					DO J=1,N
						points(1,1)=X_Coor(I)
						points(1,2)=Y_Coor(J)
						summ=summ+Distribution(coordinate1,points,Delta_x)* &
							& Distribution(coordinate2,points,Delta_x)*Delta_x**2*Delta_s
					END DO
				END DO
				A(K,L)=(fnx(L)*fnx(K)+fny(L)*fny(K))*summ
			END DO
		END DO
		write(*,*)'here we are again!'
	END SUBROUTINE Get_A
	!
	!
	SUBROUTINE Get_Fn(Fn,U,V,X_Coor,Y_Coor,Body_x,Body_y,M,N,lag_num,Delta_x,Delta_s,Delta_t,flag,AA)
		IMPLICIT NONE
		INTEGER::I,J,M,N,lag_num,flag,INFO,RANK,LWORK
		REAL(KIND=8)::AA(lag_num,lag_num),A(lag_num,lag_num),At(lag_num,lag_num),B(lag_num,1)
		REAL(KIND=8)::Fn(lag_num),U(M,N),V(M,N),X_Coor(M),Y_Coor(N)
		REAL(KIND=8)::Body_x(lag_num),Body_y(lag_num)
		REAL(KIND=8)::dUndt(lag_num)
		REAL(KIND=8)::Delta_x,Delta_s,Delta_t
		REAL(KIND=8)::RCOND
		REAL(KIND=8),allocatable::S(:),WORK(:)
        REAL(KIND=8),parameter::PI=3.141592653d0
		!REAL(KIND=8),external:: dgelss
		!
		IF(flag==1)THEN
			!CALL Get_A(A,X_Coor,Y_Coor,Body_x,Body_y,M,N,lag_num,Delta_x,Delta_s,Delta_t)
			CALL Get_A(AA,X_Coor,Y_Coor,Body_x,Body_y,M,N,lag_num,Delta_x,Delta_s,Delta_t)
			!DO I=1,lag_num
			!	DO J=1,lag_num
			!		AA(I,J)=A(I,J)
			!	END DO
			!END DO
			write(*,*)'flag',flag
			open(unit=12,file="A.dat")
			do I=1,lag_num
				do J=1,lag_num
					if(J==lag_num)then
						write(unit=12,fmt='(f18.8)',advance='yes')AA(I,J)
					else
	 					write(unit=12,fmt='(f18.8)',advance='no')AA(I,J)
					end if
				end do
			end do
			close(unit=12)
		END IF
		!
		CALL Get_dUndt(dUndt,U,V,X_Coor,Y_Coor,Body_x,Body_y,lag_num,M,N,Delta_x,Delta_t)
		!
		DO I=1,lag_num
			DO J=1,lag_num
				At(I,J)=AA(I,J)
			END DO
		END DO
		!
		DO I=1,lag_num
			B(I,1)=dUndt(I)
		END DO
		!
		LWORK=64*lag_num
		allocate(S(lag_num),WORK(LWORK))
		RCOND=0.01d0
		CALL DGELSS(lag_num, lag_num, 1, At, lag_num, B, lag_num, S, RCOND, RANK, WORK, LWORK, INFO)
		!
		DO I=1,lag_num
			Fn(I)=B(I,1)
		END DO
		deallocate(S,WORK)
	END SUBROUTINE 
	!Get_dUdV
	!
	!
	!
