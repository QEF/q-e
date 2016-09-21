include ../make.inc

NAME = Benzene

default : all

all : check_results

check_results: $(NAME).pw-out $(NAME).dav-out

clean :
	- /bin/rm -rf $(NAME).pw-out $(NAME).dav-out $(NAME).eigen $(NAME).plot out/*
