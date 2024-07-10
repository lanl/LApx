module averages_mod
    use kinds
    implicit none

contains 

function harmonicMean(values) result(mean)
    implicit none
    real(k_real), intent(in) :: values(:)
    real(k_real) :: mean
    
    mean = 1._k_real/(sum(1._k_real/values)/int2real(size(values)))
    
end function

function unweightedHarmonicMean(values) result(mean)
    implicit none
    real(k_real), intent(in) :: values(:)
    real(k_real) :: mean
    
    mean = 1._k_real/(sum(1._k_real/values))
    
end function
end module