!  HarmonicanalysisTideLoad.f90 
!
!  FUNCTIONS:
!  HarmonicanalysisTideLoad - Entry point of console application.
!
!****************************************************************************
       program HarmonicanalysisTideLoad
       implicit none
	 character*800::dtmfl,tidefl
	 integer knd
	 real*8::dk,itd
!---------------------------------------------------------------------
      !���볱ϫ��������knd
      !knd=0 surface atmosphere tide, =1 ocean tide
      knd=1!knd=0����ѹ����=1����ϫ
      !����½����������������ļ��������ڷ���½������knd=0 or =1
      !The spatial resolution of the land-sea terrain grid should not be lower than that of the surface load grid.
      write(dtmfl,*)'sphetopo30m.dat'
      !����ϵ�зֳ����ͳ�������������ļ����ļ�,���������ͷ�ļ����߸���Ϊ������ETideLoadϵͳ��ʽʱ��
      !The seventh numerical value of spherical coordinate load ewh grid file header is the long integer time agreed by ETideLoad
      write(tidefl,*)'sphgridtidefls.txt'
      !���������ֹ�������в��׼��Ϊԭ��ֵ��׼��dk��
      !input the iteration condition parameters
      !Iteration termination condition: The standard deviation of the residual grid value is less than dk of the standard deviation
      !of the original grid value, or the difference of the residual standard deviation of the previous step iteration relative to 
      !the current step iteration is less than itd of the standard deviation of the original grid values.
      dk=1.d-3; itd=1.d-4
      write(*, *)"    Begin compulation......"
      write(*, *)"    name doodson  C00(+) relative error(+)%  C00(-) relative error(-)%"
      call tideharmanalysis(dtmfl,tidefl,knd,dk,itd)
      pause
      end
