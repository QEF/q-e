include ../make.inc

NAME = CH4

default : all

all : check_results

check_results: $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out

small_test: $(NAME).pw-out $(NAME).tddfpt-st-out

clean :
	- /bin/rm -rf $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out *.plot_chi.dat out/*
