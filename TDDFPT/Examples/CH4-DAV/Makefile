include ../make.inc

NAME = CH4

default : all

all : check_results

check_results: $(NAME).pw-out $(NAME).dav-out $(NAME).tddfpt_pp-out

clean :
	- /bin/rm -rf $(NAME).pw-out $(NAME).dav-out $(NAME).tddfpt_pp-out $(NAME).eigen *.plot out/*
