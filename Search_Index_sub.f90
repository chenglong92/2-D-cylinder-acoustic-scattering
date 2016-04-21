	SUBROUTINE Search_Index(X_Coor,Y_Coor, M,N, index_x)
		IMPLICIT NONE
		INTEGER::I,J,M,N
		INTEGER::index_x(3,2)
		REAL(KIND=8)::X_Coor(M),Y_Coor(N)
		REAL(KIND=8)::A_Coor(1,2),B_Coor(1,2),C_Coor(1,2)
		!
		A_Coor(1,1)=0.0d0
		A_Coor(1,2)=5.0d0
		B_Coor(1,1)=-5.0d0/dsqrt(2.0d0)
		B_Coor(1,2)=5.0d0/dsqrt(2.0d0)
		C_Coor(1,1)=-5.0d0
		C_Coor(1,2)=0.0d0
		DO I=1,M-1
			IF((X_Coor(I)-A_Coor(1,1))*(X_Coor(I+1)-A_Coor(1,1))<=0)THEN
				index_x(1,1)=I
			END IF
			IF((X_Coor(I)-B_Coor(1,1))*(X_Coor(I+1)-B_Coor(1,1))<=0)THEN
				index_x(2,1)=I
			END IF
			IF((X_Coor(I)-C_Coor(1,1))*(X_Coor(I+1)-C_Coor(1,1))<=0)THEN
				index_x(3,1)=I
			END IF
		END DO
		!y
		DO I=1,N-1
			IF((Y_Coor(I)-A_Coor(1,2))*(Y_Coor(I+1)-A_Coor(1,2))<=0)THEN
				index_x(1,2)=I
			END IF
			IF((Y_Coor(I)-B_Coor(1,2))*(Y_Coor(I+1)-B_Coor(1,2))<=0)THEN
				index_x(2,2)=I
			END IF
			IF((Y_Coor(I)-C_Coor(1,2))*(Y_Coor(I+1)-C_Coor(1,2))<=0)THEN
				index_x(3,2)=I
			END IF
		END DO
	END SUBROUTINE Search_Index
