SRCCF = curvefitting.c
SRCCOR = correction.c
SRCPLT = myplot.c
SRCDAR = DynamicArrays.c

PROG = curvefit
PROGC = correction 
PROGP = myplot
RDEV = realDevice
ERR = detError

BST  = base.txt
CT = correct.txt
CRT  = correction.txt
RDT = realdevice.txt
ALT = alldata.txt
PNG = alldata.png
MEM = mem.txt
DATA = data.txt

CFLG = gcc -Wall -std=c99 -pedantic -g -O0
VFLG = valgrind --tool=memcheck --leak-check=yes

.SILENT:
all: $(PROG) $(PROGC) $(PROGP)
$(PROG) : $(SRCCF) $(SRCDAR)
	@echo "Compiling curvefitting"
	$(CFLG) $(SRCCF) $(SRCDAR) -o $(PROG) -DHW8 -lgsl -lm

$(PROGC): $(SRCCOR)
	@echo "Compiling correction"
	$(CFLG) $(SRCCOR) -o $(PROGC)

$(PROGP): $(SRCPLT)
	@echo "Compiling myplot"
	$(CFLG) $(SRCPLT) -o $(PROGP)

.PHONY: base
base: all
	@echo "Calculating Error from realDevice"
	./$(RDEV) | ./$(ERR) &> $(BST)
	@cat $(BST) 

.PHONY: correct
correct: all
	@echo "Calculating Error from correction"
	./$(RDEV) | ./$(PROGC) | ./$(ERR) &> $(CT)
	@cat $(CT)

.PHONY: plot
plot: all
	@echo "redirecting realDevice into realdevice.txt" 
	./$(RDEV) &> $(RDT) 
	
	@echo "redirecting correction into correction.txt" 
	./$(RDEV) | ./$(PROGC) &> $(CRT)
	paste $(RDT) $(CRT) > $(ALT) 
	
	@echo "creating alldata.png from all data"
	./$(PROGP) -i $(ALT) -o $(PNG)
	@display $(PNG) &

.PHONY: mem
mem: all
	@echo "Running Valgrind "
	$(VFLG) ./$(PROG) -n -o 3 -p $(DATA) > $(MEM) 2>&1
	
.PHONY: help
help:
	@echo "All makefile targets: all, base, correct, plot, mem, clean, help"

.PHONY: clean 
clean:
	@echo "Removing all temp file"
	rm -f $(PROG)
	rm -f $(PROGC)
	rm -f $(PROGP)
	rm -f $(BST)
	rm -f $(CRT)
	rm -f $(RDT)
	rm -f $(ALT)
	rm -f $(CT)
	rm -f $(PNG)
	rm -f $(MEM)
