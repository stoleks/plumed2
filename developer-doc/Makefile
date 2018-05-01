# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean coverage

all:
ifeq ($(make_doc),yes)
	./extract
	./go-doxygen
	cp *png html/
endif

clean:
	rm -fr automatic html latex *~ links-to-user-manual.md plumedcheck.md coverage coverage.info

coverage:
	lcov -c -d ../src/ --output-file coverage.info
# remove system headers:
	lcov -r coverage.info "/opt/*" "/usr/*" -o coverage-full.info
# remove other modules
	lcov -o coverage.info -r coverage-full.info "*/blas/*" "*/lapack/*" "*/molfile/*" "*/lepton/*"
# generate coverage
	genhtml -o coverage -t "plumed test coverage" coverage.info
# keep only other modules
	lcov -o coverage-others.info -e coverage-full.info "*/blas/*" "*/lapack/*" "*/molfile/*" "*/lepton/*"
# generate coverage
	genhtml -o coverage-others -t "plumed test coverage (other modules)" coverage-others.info
