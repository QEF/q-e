include ../make.inc

NAME = CH4-BOND

default : all

all : stage1 stage2 stage2-2

stage1: $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out
	gunzip $(NAME).pw-ref.gz
	$(check_pw) $(NAME).pw-out $(NAME).pw-ref
	gzip $(NAME).pw-ref
	gunzip $(NAME).tddfpt-ref.gz
	$(check_tddfpt) $(NAME).tddfpt-out $(NAME).tddfpt-ref
	gzip $(NAME).tddfpt-ref

stage2: $(NAME).tddfpt-out-s2
	$(check_move_cube) s2
	gunzip $(NAME).tddfpt-ref-s2.gz
	$(check_tddfpt) $(NAME).tddfpt-out-s2 $(NAME).tddfpt-ref-s2
	gzip $(NAME).tddfpt-ref-s2

stage2-2: $(NAME).tddfpt-out-s2-2
	$(check_move_cube) s2-2
	gunzip $(NAME).tddfpt-ref-s2-2.gz
	$(check_tddfpt) $(NAME).tddfpt-out-s2-2 $(NAME).tddfpt-ref-s2-2
	gzip $(NAME).tddfpt-ref-s2-2

clean :
	- /bin/rm -rf $(NAME).pw-out $(NAME).tddfpt-out $(NAME).tddfpt_pp-out $(NAME).tddfpt-st-out *.plot out/* $(NAME).tddfpt-out-s2 *.cube *.cube.gz s2 s2-2
