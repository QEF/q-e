  function capital (character)
    !-----------------------------------------------------------------------
    !
    !   converts character to capital if lowercase
    !   copy character to output in all other cases
    !
    implicit none
    character (len=1) :: capital, character
    !
    character(len=26) :: minuscole='abcdefghijklmnopqrstuvwxyz', &
                         maiuscole='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer :: i
    !
    do i=1,26
       if (character.eq.minuscole(i:i)) then
          capital=maiuscole(i:i)
          return
       end if
    end do
    capital = character
    !
    return
    end
