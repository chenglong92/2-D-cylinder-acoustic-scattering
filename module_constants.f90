	MODULE constants_Model
		IMPLICIT NONE
		REAL(KIND=8)::Length,Width,PML_Length,PML_Width
		parameter(Length=12,Width=12,PML_Length=0.5d0,PML_Width=0.5d0)
		INTEGER::M,N
		parameter(M=553,N=553)
		REAL(KIND=8)::Delta_X,Delta_Y
		parameter(Delta_X=Length/(M-1),Delta_Y=Width/(N-1))
		INTEGER::PML_M,PML_N
		parameter(PML_M=int(PML_Length/Delta_X),PML_N=int(PML_Width/Delta_Y))
		REAL(KIND=8)::PI, R
		parameter(PI=3.141592653d0, R=0.5d0)
	END MODULE constants_Model
	!
	MODULE constants_LDDRK
		IMPLICIT NONE
		REAL(KIND=8)::beita(4),w(4)
		DATA beita /0.0d0, 0.5d0, 0.5d0, 1.0d0/
		DATA w /0.1630296d0, 0.348012d0, 0.3259288d0, 0.1630296d0/
	END MODULE constants_LDDRK
