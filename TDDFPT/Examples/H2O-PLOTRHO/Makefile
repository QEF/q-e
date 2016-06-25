include ../make.inc

NAME = H2O

default : all

all : check_results

check_results: $(NAME).pw-out $(NAME).dav-out $(NAME).tddfpt_pp-out $(NAME).pp-out

clean :
	- /bin/rm -rf $(NAME).pw-out $(NAME).dav-out $(NAME).tddfpt_pp-out $(NAME).eigen *.plot out/* *.cube $(NAME).pp-out
