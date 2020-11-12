
subroutine UpperCase(string,ilen)


!------------------------------------------------------------------------------
!UpperCase  Converts a string to all upper case characters.
!   UPPERCASE takes the string "string" of length ilen and converts
!   every character of it to upper case.
!   
!   See also ReadParams
!
!   todo: 
!
!==============================================================================
!  DIPLOMA THESIS WS01/02 ICOS                                    ETH-ZUERICH  
!------------------------------------------------------------------------------
!                                                                             
!              PROTEIN DIFFUSION INSIDE THE ENDOPLASMIC RETICULUM             
!                                                                             
!============================= ivo f. sbalzarini ==============================
!

!------------------------------------------------------------------------------
! Input/Output arguments
!------------------------------------------------------------------------------

! the string to be converted
CHARACTER(LEN=*), INTENT(INOUT)             :: string
! length of the string
INTEGER, INTENT(IN)                         :: ilen

!------------------------------------------------------------------------------
! Declaration of local variables
!------------------------------------------------------------------------------

! loop counter
INTEGER                                     :: i, j
! alphabet boundaries and shift lower->upper
INTEGER                                     :: i1,i2,i3,iadd

!------------------------------------------------------------------------------
! uppercase
!------------------------------------------------------------------------------

! determine alphabet boundaries
i1   = IACHAR('a') - 1
i2   = IACHAR('z') + 1
i3   = IACHAR('A')
! shift to upper case
iadd = i3 - i1 - 1
! shift all lower case characters to upper case
do i=1,ilen
   j = IACHAR(string(i:i))
   if (j.GT.i1.AND.j.LT.i2) then
      string(i:i) = CHAR(j+iadd)
   end if
end do

!------------------------------------------------------------------------------
! return
!------------------------------------------------------------------------------

9999  CONTINUE
return
     
END subroutine UpperCase
