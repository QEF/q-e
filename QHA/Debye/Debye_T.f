        real*8 function Debye_T(x)
	implicit real*8(a-h,o-z)

	if(x.eq.0.d0) then 
	d=0.333333333333333
	else 
	d=-x/(dexp(x)-1) + 4*debye3(x)/3
	endif

	Debye_T=d
	return		
	end

