subroutine add_numbers(a, b, result)
    implicit none
    real, intent(in) :: a, b
    real, intent(out) :: result
    result = a + b
end subroutine add_numbers