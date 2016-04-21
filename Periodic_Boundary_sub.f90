	SUBROUTINE Periodic_Boundary(temp,I,J,UVPQ,K,beita,M,N,flag)
	!-------------------------introduction--------------------
	!here use the periodic boundary condition on the boundary
	!flag=0: x direction
	!flag=1: y direction
	!---------------------------------------------------------
	IMPLICIT NONE
	INTEGER::I,J,M,N,flag
	REAL(KIND=8)::UVPQ(M,N),K(M,N)
	REAL(KIND=8)::temp(7),beita
	IF(flag==0)THEN
		IF(I<4)THEN
			temp(1:(4-I))=UVPQ((M-3+I):M,J)+beita*k((M-3+I):M,J)
			temp((4-I+1):7)=UVPQ(1:(I+3),J)+beita*k(1:(I+3),J)
		ELSEIF(I>M-3)THEN
			temp((M+5-I):7)=UVPQ(1:(I+3-M),J)+beita*k(1:(I+3-M),J)
			temp(1:(M+4-I))=UVPQ((I-3):M,J)+beita*k((I-3):M,J)
		ELSE
			temp=UVPQ((I-3):(I+3),J)+beita*k((I-3):(I+3),J)
		END IF
	ELSEIF(flag==1)THEN
		IF(J<4)THEN
			temp(1:(4-J))=UVPQ(I,(N-3+J):N)+beita*k(I,(N-3+J):N)
			temp((4-J+1):7)=UVPQ(I,1:(J+3))+beita*k(I,1:(J+3))
		ELSEIF(J>N-3)THEN
			temp((N+5-J):7)=UVPQ(I,1:(J+3-N))+beita*k(I,1:(J+3-N))
			temp(1:(N+4-J))=UVPQ(I,(J-3):N)+beita*k(I,(J-3):N)
		ELSE
			temp=UVPQ(I,(J-3):(J+3))+beita*k(I,(J-3):(J+3))
		END IF
	END IF
	END SUBROUTINE Periodic_Boundary
