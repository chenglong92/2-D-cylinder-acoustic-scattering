******************************************************************************************
*This program solves the acoustic scattering problems 					**
*using Immersed boundary method and PML boundary condtion,				**
*if you are interested in the IB, You can get the papers for				**
*Peskin(C S Peskin. The immersed boundary method[J]. Acta Numerica, 2002, 11:479-517.), **
*if you are interested in IB with the application of					**
*acoustics, you can see the papers for Liang An and Jiang yongsong			**
*(梁岸, 钟国华, 孙晓峰. 基于浸入式边界方法的运动物体声辐射数值模拟[J]. 航空动力学报, 	**
*2011, 26(3):512-517.)									**
*（Sun X, Jiang Y, Liang A, et al. An immersed boundary computational model for acoustic**
* scattering problems with complex geometries[J]. Journal of the Acoustical Society of 	**
*America, 2012, 132(5):3190-3199.）							**
*----------------------------Algorithm------------------------------------------------- **
*LDDRK+DRP+PML(using periodic boundary condtion)+IBM					**
*---------------------------usage-------------------------------------------------------**
*you can revise the parameters like grid number which in the subroutine <module_constants.f90> 
* or the total calculation time which in the <main.f90>					**
*and you can run this program in the platform Linux using the <makefile> and just input "make"**
*in the current directory or in the platform Windows using the Visual Studio & Intel Fortran 90/95
******************************************************************************************

Note: This program is developed by Cheng Long from BUAA, and referred to the work by WuLong.
	So if you want to use this code, please note its' source. No copyright protected!
----------------------------------20150607---BUAA-FAEL407-------------------------------
