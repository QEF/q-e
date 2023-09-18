
default: epw

all: epw

epw: 
	(cd src ; make all )
	- (mkdir -p bin)
	(cd bin ; ln -fs ../src/epw.x . )
	(cd ZG/src ; make )

clean:
	cd src ; rm -f *.o *.a *.mod  *~ *.F90
	cd ZG/src ; rm -f *.o *.a *.mod  *~

distclean : clean
	rm -f src/epw.x bin/epw.x
	rm -f ZG/src/*.x

release:
	cd ../ ; cp -r EPW EPW-release; cd EPW-release ; \
	rm -f src/*.o src/*.mod src/*.F90 src/*~ ; \
	rm bin/*.x ; \
	rm -rf examples/*/epw/out/* examples/*/epw/tmp/* \
	 examples/*/phonons/out/* examples/*/phonons/tmp/* \
	 examples/*/phonons/save/* ; \
	cd .. ; tar cfz EPW/EPW-release.tgz EPW-release ; \
	rm -rf EPW-release ; cd EPW
	
