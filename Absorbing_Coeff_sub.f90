	SUBROUTINE Absorbing_Coeff(absorb_x,absorb_y,I,J,M,N,PML_M,PML_N,absorb_max,absorb_beita,Delta_x,Body_x,Body_y,lag_num)
		!-------------------------introduction-----------------------------------------------------------------
		!define the absorbing coefficience for all points including the boundary and the body
		!input: I,J,M,N,PML_M,PML_N,absorb_max,absorb_beita,Body_x,Body_y,lag_num,Delta_x
		!output: absorb_x,absorb_y
		!I,J: the index of calculating point
		!M,N: the x-direction and y-direction grid number
		!PML_M,PML_N: the PML zone grid number
		!absorb_max: the maximum absorbing coefficient
		!absorb_beita: the absorbing index
		!Body_x, Body_y: the body coordinates
		!Delta_x: the cartesian grid gap
		!-----------------------------------------------------------------------------------------------------
		IMPLICIT NONE
		INTEGER::II,JJ
		INTEGER::I,J,M,N,PML_M,PML_N,lag_num
		REAL(KIND=8)::absorb_x,absorb_y
		REAL(KIND=8)::absorb_max,absorb_beita
		REAL(KIND=8)::Body_x(lag_num),Body_y(lag_num)
		REAL(KIND=8)::Delta_x,R,A(1,2),IB(lag_num,2),centr(2)
		!get R,A,IB
		R=(maxval(Body_x)-minval(Body_x))/2.0d0
		A(1,1)=-(M-1)*Delta_x/2.0d0+(I-1)*Delta_x
		A(1,2)=-(N-1)*Delta_x/2.0d0+(J-1)*Delta_x
		DO II=1,lag_num
			IB(II,1)=Body_x(II)
			IB(II,2)=Body_y(II)
		END DO
		centr=0.0d0
		!
		IF((I>PML_M .and. I<(M-PML_M+1)) .and. (J>PML_N .and. J<(N-PML_N+1)))THEN
		!judge whether grid(I,J) is in the body, if it is, then (I,J) belongs to PML zone 
			IF(((A(1,1)-centr(1))**2+(A(1,2)-centr(2))**2)<=R**2)THEN
				absorb_x=absorb_max*((R-dsqrt((A(1,1)-centr(1))**2+(A(1,2)-centr(2))**2))/R)**absorb_beita
				absorb_y=absorb_x
			ELSE
				absorb_x=0.0d0
				absorb_y=0.0d0
			END IF
		ELSEIF((I<=PML_M) .and. (J>PML_N .and. J<(N-PML_N+1)))THEN
			absorb_x=absorb_max*(abs(I-PML_M)*1.0d0/(PML_M-1))**absorb_beita
			absorb_y=0.0d0
		ELSEIF((I>M-PML_M) .and. (J>PML_N .and. J<(N-PML_N+1)))THEN
			absorb_x=absorb_max*(abs(I-(M-PML_M+1))*1.0d0/(PML_M-1))**absorb_beita
			absorb_y=0.0d0
		ELSEIF((I<=PML_M) .and. (J<=PML_N))THEN
			absorb_x=absorb_max*(abs(I-PML_M)*1.0d0/(PML_M-1))**absorb_beita
			absorb_y=absorb_max*(abs(J-PML_N)*1.0d0/(PML_N-1))**absorb_beita
		ELSEIF((I<=PML_M) .and. J>(N-PML_N))THEN
			absorb_x=absorb_max*(abs(I-PML_M)*1.0d0/(PML_M-1))**absorb_beita
			absorb_y=absorb_max*(abs(J-(N-PML_N+1))*1.0d0/(PML_N-1))**absorb_beita
		ELSEIF((I>M-PML_M) .and. (J<=PML_N))THEN
			absorb_x=absorb_max*(abs(I-(M-PML_M+1))*1.0d0/(PML_M-1))**absorb_beita
			absorb_y=absorb_max*(abs(J-PML_N)*1.0d0/(PML_N-1))**absorb_beita
		ELSEIF((I>M-PML_M) .and. (J>(N-PML_N)))THEN
			absorb_x=absorb_max*(abs(I-(M-PML_M+1))*1.0d0/(PML_M-1))**absorb_beita
			absorb_y=absorb_max*(abs(J-(N-PML_N+1))*1.0d0/(PML_N-1))**absorb_beita
		ELSEIF((I>PML_M .and. I<(M-PML_M+1)) .and. (J>(N-PML_N)))THEN
			absorb_x=0.0d0
			absorb_y=absorb_max*(abs(J-(N-PML_N+1))*1.0d0/(PML_N-1))**absorb_beita
		ELSEIF((I>PML_M .and. I<(M-PML_M+1)) .and.  (J<=PML_N))THEN
			absorb_x=0.0d0
			absorb_y=absorb_max*(abs(J-PML_N)*1.0d0/(PML_N-1))**absorb_beita
		END IF
	END SUBROUTINE Absorbing_Coeff
