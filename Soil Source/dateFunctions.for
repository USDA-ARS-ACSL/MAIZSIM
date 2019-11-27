c these functions come from the numerical recipes book. I modified the julian day base to be
c 0 for 3/1/1900 to be compatible with Microsoft Excel's method of calculating Julian dates
C for astronomy (not the julian calendar)
C Caldat returns a calendar date given a julian day
C Julian is a function that returns a julian day given a calendar day


      SUBROUTINE caldat(julian,mm,id,iyyy)
! passes julian and returns mm, dd, yyyy
      INTEGER id,iyyy,julian,mm,IGREG
      PARAMETER (IGREG=2299161)
!(IGREG=2299161)
      INTEGER ja,jalpha,jb,jc,jd,je
      !Julian=julian+2415079   !provide a reference of 3/1/1900
      Julian=julian+2415019   ! now 12/30/1899
      if(julian.ge.IGREG)then
        jalpha=int(((julian-1867216)-0.25d0)/36524.25d0)
        ja=julian+1+jalpha-int(0.25d0*jalpha)
      else if(julian.lt.0)then
        ja=julian+36525*(1-julian/36525)
      else
        ja=julian
      endif
      jb=ja+1524
      jc=int(6680.0d0+((jb-2439870)-122.1d0)/365.25d0)
      jd=365*jc+int(0.25d0*jc)
      je=int((jb-jd)/30.6001d0)
      id=jb-jd-int(30.6001d0*je)
      mm=je-1
      if(mm.gt.12)mm=mm-12
      iyyy=jc-4715
      if(mm.gt.2)iyyy=iyyy-1
      if(iyyy.le.0)iyyy=iyyy-1
      if(julian.lt.0)iyyy=iyyy-100*(1-julian/36525)
      return
      END
      
      FUNCTION julday(Date)
      INTEGER julday,id,iyyy,mm,IGREG
      PARAMETER (IGREG=15+31*(10+12*1582))
      INTEGER ja,jm,jy
      Character*10 date
      
      Call FromCalendar(date,mm,id,iyyy)
      
      
      jy=iyyy
      if (jy.eq.0) pause 'julday: there is no year zero'
      if (jy.lt.0) jy=jy+1
      if (mm.gt.2) then
        jm=mm+1
      else
        jy=jy-1
        jm=mm+13
      endif
      julday=365*jy+int(0.25d0*jy+2000.d0)+int(30.6001d0*jm)+id+1718995
      if (id+31*(mm+12*iyyy).ge.IGREG) then
        ja=int(0.01d0*jy)
        julday=julday+2-ja+int(0.25d0*ja)
        !julday=julday-2415019   ! 1/1/1900 JD zero
        Julday=Julday-2415019 ! 12/30/1899 JD is 0
        !julday=julday-2415079;      ! on 3/1/1900 jday = 1
      endif
      return
      END
C**********************************************************************************
C Subroutine fromCalendar
C**********************************************************************************      
      Subroutine FromCalendar(date,month,day,year)
C takes a calendar date string as mm/dd/yyyy and returns the month, day and year
      
      character*10 date
      Character*10 date2
      integer month, day, year,pos
      
       pos=index(date,'/')
    
      read(date(1:pos-1),*) month
      date2=date(pos+1:)
      pos=index(date2,'/')
      if (pos.gt.0) then
         read(date2(1:pos-1),*) day
         read(date2(pos+1:),*) year
        else                               ! case for no year
         read(date2,*) day
        endif
      return
       END      
       
       Function toCalendarDate(month,day,year)
       
       character*10 toCalendarDate, temp
       integer month, day, year
              
       write (temp,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy 
       toCalendarDate=temp
       return
       END
             
      
        integer Function CurYear(date)
      
          character*10 date
          integer month, day, iyear
          call FromCalendar(Date,month,day, CurYear)
          
            
          return 
       END
       
       
       