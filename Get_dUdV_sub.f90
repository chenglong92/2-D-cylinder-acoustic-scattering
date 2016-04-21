	SUBROUTINE Get_dUdV(dU,dV,U,V,X_Coor,Y_Coor,Body_x,Body_y,M,N,lag_num,Delta_x,Delta_s,Delta_t,flag,AA)
		IMPLICIT NONE
		INTEGER::I,J,K,M,N,lag_num,flag
		REAL(KIND=8)::dU(M,N),dV(M,N),U(M,N),V(M,N)
		REAL(KIND=8)::X_Coor(M),Y_Coor(N),Body_x(lag_num),Body_y(lag_num)
		REAL(KIND=8)::AA(lag_num,lag_num)
		REAL(KIND=8)::Fn(lag_num),sumFx,sumFy
		REAL(KIND=8)::fnx(lag_num),fny(lag_num)
		REAL(KIND=8)::Delta_x,Delta_s,Delta_t
		REAL(KIND=8)::coordinate(1,2),points(1,2)
        REAL(KIND=8),external:: Distribution
		!
		CALL Get_Fn(Fn,U,V,X_Coor,Y_Coor,Body_x,Body_y,M,N,lag_num,Delta_x,Delta_s,Delta_t,flag,AA)
		!write(*,*)'here'
		!
		CALL Normal_Vector(fnx,fny,Body_x,Body_y,lag_num)
		!
		DO I=1,M
			DO J=1,N
				sumFx=0.0d0
				sumFy=0.0d0
				points(1,1)=X_Coor(I)
				points(1,2)=Y_Coor(J)
				DO K=1,lag_num
					coordinate(1,1)=Body_x(K)
					coordinate(1,2)=Body_y(K)
					sumFx=sumFx+Fn(K)*fnx(K)*Delta_s*Distribution(coordinate,points,Delta_x)
					sumFy=sumFy+Fn(K)*fny(K)*Delta_s*Distribution(coordinate,points,Delta_x)
				END DO
				dU(I,J)=sumFx*Delta_t
				dV(I,J)=sumFy*Delta_t
			END DO
		END DO
	END SUBROUTINE Get_dUdV
