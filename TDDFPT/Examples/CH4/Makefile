include ../make.sys

NAME = CH4

default : all

all : check_results

check_results: $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out
	$(check_pw) $(NAME).pw-out $(NAME).pw-ref
	$(check_tddfpt) $(NAME).tddfpt-out $(NAME).tddfpt-ref
	$(check_pp) $(NAME).tddfpt_pp-out $(NAME).tddfpt_pp-ref

clean :
	- /bin/rm -f $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out *.plot
