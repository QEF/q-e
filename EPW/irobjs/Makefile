# Makefile for IR object files

default: irobjs

all: irobjs

irobjs:
	@if test "`which curl`" = ""; then \
		if test "`which wget`" = ""; then \
			echo "### wget or curl not found: will not be able to download IR object files ###"; \
		else \
			wget -O irobjs.tar.gz https://docs.epw-code.org/_downloads/de79081ec6c4c58e05779d2c878e8689/q-e-develop-EPW-irobjs.tar.gz; \
		fi \
	else \
		curl -o irobjs.tar.gz https://docs.epw-code.org/_downloads/de79081ec6c4c58e05779d2c878e8689/q-e-develop-EPW-irobjs.tar.gz; \
	fi
	tar -zxvf irobjs.tar.gz
	mv ./q-e-develop-EPW-irobjs/EPW/irobjs/ir_nlambda*_ndigit*.dat ./
	rm -r q-e-develop-EPW-irobjs ; rm irobjs.tar.gz