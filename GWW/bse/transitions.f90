! this module contains the variables related to the solution of the BSE 
! in transition space

MODULE transitions

USE kinds, ONLY: DP
INTEGER, ALLOCATABLE :: ttab(:,:) !maps the (iv,ic) couple into the it index
INTEGER, ALLOCATABLE :: itiv(:)!for the it-th transition gives the
INTEGER, ALLOCATABLE :: itic(:)!corresponding valence and conduction band index respectively 
REAL(KIND=DP), ALLOCATABLE :: exch(:,:) !excitonic Hamiltonian in transition space


END MODULE transitions
