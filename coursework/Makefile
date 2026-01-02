help:
	@echo "Complete your assignment, including creating files:"
	@echo "	 part1.c,"
	@echo "	 part2.c,"
	@echo "	 part1.pdf, and"
	@echo "	 part2.pdf,"
	@echo "in this directory."
	@echo "Then run 'make check' to ensure all files are present and produce correct outputs."
	@echo ""
	@echo "Once happy, run 'make submit' to create a zip file of your submission, "
	@echo "download the <studentID>.zip file with "
	@echo "		rsync -uzhr <username>@<server>:/path/to/<studentID>.zip .			"
	@echo "and upload <studentID>.zip to the submission point on Blackboard."
	@echo ""
	@echo "This Makefile is a simple check of the submission and is not a guarantee of marks nor can it submit the <studentID>.zip file for you."

all: part1.dat part2.dat serial.dat

submission:
	@echo "Creating $(USER).zip from part1.c part2.c part1.pdf part2.pdf:" && zip -r ${USER}.zip part1.c part2.c part1.pdf part2.pdf;

%.c:
	touch $@

%.dat: %
	@echo "Generating $@ from $<"

%: %.sbatch
	@echo "Submitting $<" && sbatch $<

%.sbatch: %.c params.h

serial: serial.c
	@echo "Running serial code:" && sbatch serial.sbatch

part1: part1.c
	@echo "Running part 1 code:" && sbatch part1.sbatch

part2: part1.c
	@echo "Running part 2 code:" && sbatch part2.sbatch

check: part1.dat part2.dat serial.dat
	@echo "Checking part1.dat:" && \
	diff -s part1.dat serial.dat; \
	echo "Checking part2.dat:" && \
	diff -s part2.dat serial.dat

clean_outputs:
	for file in part1 part2 serial; do \
		rm -f $$file.dat; \
	done \

clean_executables:
	for file in part1 part2 serial; do \
		rm -f $$file; \
	done \

clean: clean_outputs clean_executables
	rm -f resit.zip slurm*.out

.PHONY: help
