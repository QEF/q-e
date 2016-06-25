include ../make.inc

NAME = Bi

default : all

all : check_results

check_results: $(NAME).pw-out $(NAME).tddfpt-eels-out $(NAME).tddfpt_pp-out

clean :
	- /bin/rm -rf $(NAME).pw-out $(NAME).tddfpt-eels-out $(NAME).tddfpt_pp-out $(NAME).tddfpt-st-out *.plot_chi.dat *.plot_eps.dat out/* 
