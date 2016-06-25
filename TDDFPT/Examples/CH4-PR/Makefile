include ../make.inc

NAME = CH4-PR

default : all

all : stage1 stage2

stage1: $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out
	$(check_pw) $(NAME).pw-out $(NAME).pw-ref
	$(check_tddfpt) $(NAME).tddfpt-out $(NAME).tddfpt-ref

stage2: $(NAME).tddfpt-out-s2
	$(check_tddfpt) $(NAME).tddfpt-out-s2 $(NAME).tddfpt-ref-s2

clean :
	- /bin/rm -rf $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out $(NAME).tddfpt-st-out *.plot out/* $(NAME).tddfpt-out-s2 *.cube 
