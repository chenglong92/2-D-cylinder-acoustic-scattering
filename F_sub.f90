	SUBROUTINE F(beita,k1,k2,k3,k4,k5,k6,U,V,P,Qu,Qv,Qp,X_coor,Y_coor,Delta_t,Body_x,Body_y,lag_num)
		USE constants_Model
		IMPLICIT NONE
		INTEGER::I,J,flag,lag_num
		REAL(KIND=8)::X_coor(M),Y_coor(N),Body_x(lag_num),Body_y(lag_num)
		REAL(KIND=8)::U(M,N),V(M,N),P(M,N),Qu(M,N),Qv(M,N),Qp(M,N)
		REAL(KIND=8)::Delta_t,beita
		REAL(KIND=8)::K1(M,N),K2(M,N),K3(M,N),K4(M,N),K5(M,N),K6(M,N)
		REAL(KIND=8)::K11(M,N),K12(M,N),K13(M,N),K14(M,N),K15(M,N),K16(M,N)
		REAL(KIND=8)::absorb_x,absorb_y,absorb_max,absorb_beita
		parameter(absorb_max=2.0d0/Delta_x, absorb_beita=2.0d0)
		REAL(KIND=8)::temp(7)
		REAL(KIND=8),external::DRP7
		!
		!------------------------------------------------------------
		!for boundary condition, we use the periodic boundary condition
		!for the PML boundaries since the numerical solution decays 
		!exponentially toward all the boundaries
		!------------------------------------------------------------
		!initialize
		k11=k1
		k12=k2
		k13=k3
		k14=k4
		k15=k5
		k16=k6
		k1=0.0d0
		k2=0.0d0
		k3=0.0d0
		k4=0.0d0
		k5=0.0d0
		k6=0.0d0
		!
		!solve the Euler zone, here absorbing coefficient is 0
		!
		DO I=1,M
			DO J=1,N
				CALL Absorbing_Coeff(absorb_x,absorb_y,I,J,M,N,PML_M,PML_N,absorb_max,absorb_beita,Delta_x,Body_x,Body_y,lag_num)
				!derivative for x direction
				!write(*,*)'absorbing:',absorb_x,absorb_y
				flag=0
				CALL Periodic_Boundary(temp,I,J,P,K13,beita,M,N,flag)
				k1(I,J)=k1(I,J)+DRP7(temp,Delta_x,-3)
					!
				CALL Periodic_Boundary(temp,I,J,Qp,K16,beita,M,N,flag)
				k1(I,J)=k1(I,J)+absorb_y*DRP7(temp,Delta_x,-3)
					!
				CALL Periodic_Boundary(temp,I,J,U,K11,beita,M,N,flag)
				k3(I,J)=k3(I,J)+DRP7(temp,Delta_x,-3)
					!
				CALL Periodic_Boundary(temp,I,J,Qu,K14,beita,M,N,flag)
				k3(I,J)=k3(I,J)+absorb_y*DRP7(temp,Delta_x,-3)
					!
				!derivative for y direction
				flag=1
				CALL Periodic_Boundary(temp,I,J,P,K13,beita,M,N,flag)
				k2(I,J)=k2(I,J)+DRP7(temp,Delta_x,-3)
					!
				CALL Periodic_Boundary(temp,I,J,Qp,K16,beita,M,N,flag)
				k2(I,J)=k2(I,J)+absorb_x*DRP7(temp,Delta_x,-3)
					!
				CALL Periodic_Boundary(temp,I,J,V,K12,beita,M,N,flag)
				k3(I,J)=k3(I,J)+DRP7(temp,Delta_x,-3)
					!
				CALL Periodic_Boundary(temp,I,J,Qv,K15,beita,M,N,flag)
				k3(I,J)=k3(I,J)+absorb_x*DRP7(temp,Delta_x,-3)
					!
				!for polynomial terms
				k1(I,J)=k1(I,J)+(absorb_x+absorb_y)*(U(I,J)+beita*K11(I,J))&
							&+absorb_x*absorb_y*(Qu(I,J)+beita*K14(I,J))
				K2(I,J)=K2(I,J)+(absorb_x+absorb_y)*(V(I,J)+beita*K12(I,J))+&
							&absorb_x*absorb_y*(Qv(I,J)+beita*k15(I,J))
				K3(I,J)=K3(I,J)+(absorb_x+absorb_y)*(P(I,J)+beita*K13(I,J))+&
							&absorb_x*absorb_y*(Qp(I,J)+beita*k16(I,J))
				k4(I,J)=K4(I,J)+(U(I,J)+beita*K11(I,J))
				K5(I,J)=K5(I,J)+(V(I,J)+beita*K12(I,J))
				K6(I,J)=K6(I,J)+(P(I,J)+beita*K13(I,J))
			END DO
		END DO
		!
		K1=-K1*Delta_t
		k2=-k2*Delta_t
		k3=-k3*Delta_t
		k4=k4*Delta_t
		k5=k5*Delta_t
		k6=k6*Delta_t
	END SUBROUTINE F
