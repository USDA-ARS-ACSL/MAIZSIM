c==== Error messages
      Subroutine ErrMes(N,M)
      Include 'public.ins'
      If(N.ne.0) then
        Write(4,20) N,M
        Stop
      Endif
20    format('Error # ',i4,',   Line # ',i4)
      Return
      End
       
