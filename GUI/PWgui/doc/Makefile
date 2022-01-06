TOPDIR        = ../../..

all: input_html make_user_guide

input_html:
	./htmlise.sh $(TOPDIR)/PW/Doc/INPUT_PW.txt      > INPUT_PW.html
	./htmlise.sh $(TOPDIR)/NEB/Doc/INPUT_NEB.txt    > INPUT_NEB.html
	./htmlise.sh $(TOPDIR)/PHonon/Doc/INPUT_PH.txt  > INPUT_PH.html
	./htmlise.sh $(TOPDIR)/PP/Doc/INPUT_PP.txt      > INPUT_PP.html
	./htmlise.sh $(TOPDIR)/PP/Doc/INPUT_PROJWFC.txt > INPUT_PROJWFC.html
	./htmlise.sh $(TOPDIR)/PP/Doc/INPUT_BANDS.txt   > INPUT_BANDS.html
	./htmlise.sh $(TOPDIR)/PP/Doc/INPUT_DOS.txt     > INPUT_DOS.html
	./htmlise.sh $(TOPDIR)/atomic/Doc/INPUT_LD1.txt > INPUT_LD1.html
	./htmlise.sh $(TOPDIR)/HP/Doc/INPUT_HP.txt      > INPUT_HP.html

qe_user_guide:
	- ( cd $(TOPDIR)/Doc; $(MAKE) user_guide )

make_user_guide: qe_user_guide
	rm -rf user_guide/
	cp -a $(TOPDIR)/Doc/user_guide .
	cd user_guide/; \
	for file in *.html; do \
		cp $$file /tmp/$$file; \
		cat /tmp/$$file | sed 's/&nbsp;/ /g' - > $$file; \
		rm -f /tmp/$$file; \
	done


clean: 
	-rm -f  *~
	-rm -rf user_guide/	


distclean: clean
veryclean: clean
