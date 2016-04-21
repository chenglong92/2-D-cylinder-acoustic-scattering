	SUBROUTINE Model(Body_X,Body_Y,M,Delta_s)
	!-------------------introduction------------------
	!input: Delta_x
	!output: Body_X,Body_Y,M,Delta_s
	!here R is the radius of cylinder,centre is the 
	!centre of cylinder
	!-------------------------------------------------
	IMPLICIT NONE
	REAL(KIND=8)::PI
	parameter(PI=3.141592653d0)
	INTEGER::I,M
	REAL(KIND=8)::Delta_s,Delta_thita
	REAL(KIND=8)::R=0.5d0
	REAL(KIND=8)::Centre(2)=(/0.0d0,0.0d0/)
	REAL(KIND=8)::Body_x(M),Body_y(M)
	!
	!M=int(PI*2.0*R/Delta_x)
	!allocate(Body_x(M),Body_Y(M))
	!
	Delta_thita=2.0*PI/M
	DO I=1,M
		Body_x(I)=Centre(1)+R*dcos(Delta_thita*(I-1))
		Body_y(I)=Centre(2)+R*dsin(Delta_thita*(I-1))
	END DO
	Delta_s=dcos(Delta_thita)*R
	END SUBROUTINE Model
