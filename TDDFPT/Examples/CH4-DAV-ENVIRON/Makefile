include ../make.inc

NAME = CH4

default : all

all : check_results

check_results: $(NAME).pw-out $(NAME).dav-out

clean :
	- /bin/rm -rf $(NAME).pw-out $(NAME).dav-out $(NAME).eigen *.plot out/*
