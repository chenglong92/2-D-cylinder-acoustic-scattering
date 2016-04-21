	SUBROUTINE Initialization(U, V, P, Qu, Qv, Qp, X_Coor,Y_Coor, M,N)
	IMPLICIT NONE
	INTEGER::I,J,M,N
	REAL(KIND=8)::X_Coor(M),Y_Coor(N)
	REAL(KIND=8)::U(M,N),V(M,N),P(M,N),Qu(M,N),Qv(M,N),Qp(M,N)
	REAL(KIND=8)::kexi=1.0d-3, x0=4.0d0, y0=0.0d0, b=0.2d0
	!
	U=0.0d0
	V=0.0d0
	Qu=0.0d0
	Qv=0.0d0
	Qp=0.0d0
	!
	DO I=1,M
		DO J=1,N
			P(I,J)=kexi*dexp(-log(2.0d0)*((X_Coor(I)-4.0d0)**2+Y_Coor(J)**2)/b**2)
		END DO
	END DO
	END SUBROUTINE Initialization
