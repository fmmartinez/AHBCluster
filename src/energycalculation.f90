module energycalculation
use definitions
implicit none

contains

function get_total_energy(pairs) result(e)
implicit none

      real(8) :: e
      real(8) :: rab,rah,rbh
      type(AtomPairData),dimension(:,:),intent(in) :: pairs
   
      rab = pairs(1,2)%rij
      rah = pairs(1,3)%rij
      rbh = pairs(2,3)%rij
      e = b*exp(-a*rab)+d*(1-exp(-(na*(rah-da)**2)/(2d0*rah)))+c*d*(1-exp(-nb*(rbh-db)**2/(2d0*rbh)))

end function get_total_energy
end module energycalculation
