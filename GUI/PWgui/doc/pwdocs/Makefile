TOPDIR        = ../../../..
QE_DOC_DIR    = $(TOPDIR)/Doc

INPUT_FILES_HTML = \
        INPUT_PW.html \
	INPUT_NEB.html \
        INPUT_PH.html \
        INPUT_PP.html \
        INPUT_PROJWFC.html \
        INPUT_BANDS.html \
        INPUT_DOS.html \
        INPUT_LD1.html \
        INPUT_D3.html

INPUT_FILES_TXT = $(INPUT_FILES_HTML:.html=.txt)
LINK_FILES      = $(INPUT_FILES_TXT) 


all: input_html make_user_guide


links:
	( cd $(TOPDIR); $(MAKE) doc )
	for file in $(LINK_FILES); do \
		if test ! -f $$file; then ln -sf $(QE_DOC_DIR)/$$file . ; fi; \
	done

clean_links:
	for file in $(LINK_FILES); do \
		if test -f $$file ; then rm -f $$file ; fi ; \
	done


input_html: links 
	./htmlise.sh INPUT_PW.txt      > INPUT_PW.html
	./htmlise.sh INPUT_NEB.txt     > INPUT_NEB.html
	./htmlise.sh INPUT_PH.txt      > INPUT_PH.html
	./htmlise.sh INPUT_PP.txt      > INPUT_PP.html
	./htmlise.sh INPUT_PROJWFC.txt > INPUT_PROJWFC.html
	./htmlise.sh INPUT_BANDS.txt   > INPUT_BANDS.html
	./htmlise.sh INPUT_DOS.txt     > INPUT_DOS.html
	./htmlise.sh INPUT_D3.txt      > INPUT_D3.html
	./htmlise.sh INPUT_LD1.txt     > INPUT_LD1.html

make_user_guide: links
	rm -rf user_guide/
	cp -a $(QE_DOC_DIR)/user_guide .
	cd user_guide/; \
	for file in *.png; do \
		convert $$file $${file%.png}.gif; \
		rm -f $$file; \
	done; \
	for file in *.html; do \
		cp $$file /tmp/$$file; \
		cat /tmp/$$file | sed 's/\.png/\.gif/g' - | sed 's/&nbsp;/ /g' - > $$file; \
		rm -f /tmp/$$file; \
	done


clean: clean_links
	-rm -f INPUT_*.html *~
	-rm -f *.tex *.aux *.log *.out *.toc *.gif *.pdf *.png
	-rm -rf user_guide/	


distclean: clean
veryclean: clean