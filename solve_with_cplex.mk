# CPLEX has a limit of 560 characters per line, so we need to use fmt to wrap the lines!
# Luckily the GnuWin32 version exists on my Windows machines.

# We can't actually use the $$ shell variable sensibly because make starts a new shell for each recipe line!
# We could dig out make's own $$, but that would break if we tried to use make -j.  Using $@ in a recursively-expanded variable is simple and 100% safe :)
#CPLEXWORKDIR := cplex_work_dir.`hostname`.$$$$
CPLEXWORKDIR = $@.cplex_work_dir
CPLEXTIMELIMIT := 7200
# That's 1 hour
#CPLEXTIMELIMIT := 3600
# 4Gb work limit; may increase it later
#CPLEXWORKMEM := 4096
#CPLEXWORKMEM := 8192
# Some instances need almost 11GB (and some need even more, but I'll special-case those ones).
CPLEXWORKMEM := 11264

# Other quirks:
# - Can't seem to turn off the need to prompt for overwrite, hence we delete the target at the outset
# - Produces an XML file
# - Can't seem to control the type of output that "write" does -- specifying e.g. "sol" (which is apparently the default anyway) causes an error.
# - Writes temporary files to cur dir by default, causing problems for parallel usage
# - Let's use the $(GRBTHREADS) variable to control the number of threads.

# We now use Stephan's program, with -l O, for this as well.
# $(RECORDCMD) will run the command line in $(CMD) and save info about it to $@.run.
%.orig_cplex.sol: CMD=$(TIME) tools/solve_cplex -D "$(CPLEXWORKDIR)" -T $(CPLEXTIMELIMIT) -M $(CPLEXWORKMEM) -l o -o $@ $< > $(@:.sol=.stdout) 2> $(@:.sol=.stderr)
%.orig_cplex.sol: %.txt
	$(SAFEMKDIR) $(CPLEXWORKDIR)
	$(DELCMD) $@
#	$(TIME) cplex -c "set workdir $(CPLEXWORKDIR)" "set timelimit $(CPLEXTIMELIMIT)" "set workmem $(CPLEXWORKMEM)" "set threads $(GRBTHREADS)" "read $< lp" mipopt "write $@" > $(@:.sol=.stdout) 2> $(@:.sol=.log)
	$(RECORDCMD)
	-rmdir $(CPLEXWORKDIR)

%.ilp.wrapped: %.ilp
	$(TIME) fmt -s < $< > $@

# Convert from the crazy XML file to a version that should compare byte-for-byte identical with other canonicalised solutions.
%cplex.sol.canon: %cplex.sol
	$(TIME) canonicalise_cplex_solution.pl < $< > $@

# Stephan's strengthened version of CPLEX.
# $(RECORDCMD) will run the command line in $(CMD) and save info about it to $@.run.
# This doesn't actually produce a list of edges in the solution -- just the optimal score!  But that should be enough.
#%.cplex_stephan.sol: %.txt
#	$(TIME) solve_cplex -l C $< > $@ 2> $@.log
%.cut_cplex.sol: CMD=$(TIME) tools/solve_cplex -D "$(CPLEXWORKDIR)" -T $(CPLEXTIMELIMIT) -M $(CPLEXWORKMEM) -l c -o $@ $< > $(@:.sol=.stdout) 2> $(@:.sol=.stderr)
%.cut_cplex.sol: %.txt
	$(SAFEMKDIR) $(CPLEXWORKDIR)
	$(DELCMD) $@
#	$(TIME) cplex -c "set workdir $(CPLEXWORKDIR)" "set timelimit $(CPLEXTIMELIMIT)" "set workmem $(CPLEXWORKMEM)" "set threads $(GRBTHREADS)" "read $< lp" mipopt "write $@" > $(@:.sol=.stdout) 2> $(@:.sol=.log)
	$(RECORDCMD)
	-rmdir $(CPLEXWORKDIR)
