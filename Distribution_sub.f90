REAL(KIND=8) FUNCTION Distribution(coordinate,points,h)
 IMPLICIT NONE
 INTEGER::I,J
 REAL(KIND=8)::h,r1,r2,f1,f2,coordinate(1,2),points(1,2)
!
 r1=dabs((points(1,1)-coordinate(1,1))/h);
 r2=dabs((points(1,2)-coordinate(1,2))/h);
 IF(r1<1.0d0)THEN
     f1=(3.0d0-2.0d0*r1+dsqrt(1+4.0d0*r1-4.0d0*r1**2))/8.0d0
 ELSE IF(r1>=1.0d0 .and. r1<2.0d0)THEN
     r1=2-r1
     f1=0.5-(3.0d0-2.0d0*r1+dsqrt(1+4.0d0*r1-4.0d0*r1**2))/8.0d0
 ELSE
    f1=0
 END IF
!
!
 IF(r2<1.0d0)THEN
     f2=(3.0d0-2.0d0*r2+dsqrt(1+4.0d0*r2-4.0d0*r2**2))/8.0d0
 ELSE IF(r2>=1.0d0 .and. r2<2.0d0)THEN
     r2=2-r2
     f2=0.5-(3.0d0-2.0d0*r2+dsqrt(1+4.0d0*r2-4.0d0*r2**2))/8.0d0
 ELSE
     f2=0
 END IF
!
!
 Distribution=f1*f2/h**2
END FUNCTION Distribution
