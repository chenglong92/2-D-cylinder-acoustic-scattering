	SUBROUTINE LDDRK(U,V,P,Qu,Qv,Qp,X_coor,Y_coor,M,N,Body_x,Body_y,lag_num,Delta_x,Delta_t)
		USE constants_LDDRK
		!----------------------introduction-----------------
		!this subroutine considers the case which of no flow
		!input: U,V,P,Qu,Qv,Qp,m,n,Delta_x,Delta_t
		!output: U,V,P,Qu,Qv,Qp
		!---------------------------------------------------
		!the governing equation is:
		!	
		!
		!---------------------------------------------------
		IMPLICIT NONE
		INTEGER::I,J,M,N,lag_num
		REAL(KIND=8)::X_coor(M),Y_coor(N),Body_x(lag_num),Body_y(lag_num)
		REAL(KIND=8)::U(M,N),V(M,N),P(M,N),Qu(M,N),Qv(M,N),Qp(M,N)
		REAL(KIND=8)::U1(M,N),V1(M,N),P1(M,N),Qu1(M,N),Qv1(M,N),Qp1(M,N)
		REAL(KIND=8)::Delta_x,Delta_t
		REAL(KIND=8)::K1(M,N),K2(M,N),K3(M,N),K4(M,N),K5(M,N),K6(M,N)
		!
		!initialization
		k1=0.0d0
		k2=0.0d0
		k3=0.0d0
		k4=0.0d0
		k5=0.0d0
		k6=0.0d0
		U1=U
		V1=V
		P1=P
		Qu1=Qu
		Qv1=Qv
		Qp1=Qp
		!
		DO I=1,4
			CALL F(beita(I),k1,k2,k3,k4,k5,k6,U1,V1,P1,Qu1,Qv1,Qp1,X_coor,Y_coor,Delta_t,Body_x,Body_y,lag_num)
			U=U+w(I)*k1
			V=V+w(I)*k2
			P=P+w(I)*k3
			Qu=Qu+w(I)*k4
			Qv=Qv+w(I)*k5
			Qp=Qp+w(I)*k6
		END DO
	END SUBROUTINE LDDRK
