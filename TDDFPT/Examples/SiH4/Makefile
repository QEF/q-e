include ../make.sys

NAME = SiH4

default : all

all : check_results

check_results: $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out
	$(check_pw) $(NAME).pw-out $(NAME).pw-ref
	$(check_tddfpt) $(NAME).tddfpt-out $(NAME).tddfpt-ref

small_test: $(NAME).pw-out $(NAME).tddfpt-st-out
	@$(check_pw) $(NAME).pw-out $(NAME).pw-ref
	@$(check_tddfpt) $(NAME).tddfpt-st-out $(NAME).tddfpt-st-ref

clean :
	- /bin/rm -rf $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out $(NAME).tddfpt-st-out *.plot out/*
