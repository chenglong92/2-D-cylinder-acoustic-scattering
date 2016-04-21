	SUBROUTINE Grid(X_Coor,Y_Coor,M,N,Length,Width)
	!------------------------introduction--------------
	!input: M,N,Length,Width,Delta_X,Delta_Y
	!output: X_Coor,Y_Coor
	!--------------------------------------------------
	IMPLICIT NONE
	INTEGER::I,M,N
	REAL(KIND=8)::X_Coor(M),Y_Coor(N)
	REAL(KIND=8)::Length,Width
	!
	DO I=1,M
		X_Coor(I)=-Length/2.0d0+Length*(I-1)*1.0d0/(M-1)
	END DO
	!
	DO I=1,N
		Y_Coor(I)=-Width/2.0d0+Width*(I-1)*1.0d0/(N-1)
	END DO
	END SUBROUTINE Grid
