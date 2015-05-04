#
# This makefile compiles all tests and runs them via 
# a diff compare
#

.PHONY: all demo heat web test clean

all: test

demo: dkrdem.c dkrdem.out.txt
	@make -f make_cddem
	
heat: dheat.c dheat.out.txt
	@make -f make_cdh
	
heatilu: dheatilu.c dheatilu.out.txt
	@make -f make_cdhilu

web: dweb.c dweb.out.txt
	@make -f make_cdw

webilu: dwebilu.c dwebilu.out.txt
	@make -f make_cdwilu
	
test: demo heat heatilu web webilu
	@./demo > dkrdem.new.out.txt
	@cmp -s dkrdem.out.txt dkrdem.new.out.txt; \
	RETVAL=$$?; \
	echo -n "Test Demo: "; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "OK"; \
	else \
		echo "WRONG"; \
	fi
	@./heat > dheat.new.out.txt
	@cmp -s dheat.out.txt dheat.new.out.txt; \
	RETVAL=$$?; \
	echo -n "Test HEAT: "; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "OK"; \
	else \
		echo "WRONG"; \
	fi
	@./heatilu > dheatilu.new.out.txt
	@cmp -s dheatilu.out.txt dheatilu.new.out.txt; \
	RETVAL=$$?; \
	echo -n "Test HEATILU: "; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "OK"; \
	else \
		echo "WRONG"; \
	fi
	@./web > dweb.new.out.txt
	@cat wdout wccout >> dweb.new.out.txt
	@cmp -s dweb.out.txt dweb.new.out.txt; \
	RETVAL=$$?; \
	echo -n "Test WEB: "; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "OK"; \
	else \
		echo "WRONG"; \
	fi
	@./webilu > dwebilu.new.out.txt
	@cat wdout wccout >> dwebilu.new.out.txt
	@cmp -s dwebilu.out.txt dwebilu.new.out.txt; \
	RETVAL=$$?; \
	echo -n "Test WEBILU: "; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "OK"; \
	else \
		echo "WRONG"; \
	fi
	
clean:
	@make -f make_cddem clean
	@make -f make_cdh clean
	@make -f make_cdhilu clean
	@make -f make_cdw clean
	@make -f make_cdwilu clean
	@rm -f dkrdem.new.out.txt \
	      dweb.new.out.txt \
		  dwebilu.new.out.txt \
		  dheat.new.out.txt \
		  dheatilu.new.out.txt \
		  wccout \
		  wdout