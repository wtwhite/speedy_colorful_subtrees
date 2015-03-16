# Makefile rules for solving maximum colourful subtree graphs with both the Gurobi ILP solver and lp_solve's ILP solver, using Kai's ilp.jar to convert a graph file to a .ilp input file for them.
# Can be run as e.g.
#    make PROBLEMSPAT="some_subdir\*.txt" all_gurobi
# to work with the problems in some subdirectory, or matching just some wildcard pattern, instead of the default pattern (pos44247*.txt in the current directory).
#
# To specify that reduced versions of each problem in PROBLEMSPAT should be created using ft_reduce, also specify REDUCE=1.

# Set $S to forward or backward slash based on detected OS.  This is necessary because $(wildcard ...) doesn't expand paths with backslashes in them,
# so we need to convert to and from forward slashes around calls to it!
# We also need $ES ("escaped slash"), since in certain places (e.g. target names) the backslash needs to be doubled (while in other places it doesn't)...
ifeq ($(OS:Win%=YES),YES)
# The following is the only way I found of getting just a single backslash (a backslash at end-of-line continues the line... 2 backslashes stay as 2 backslashes (sometimes)...)
S := \$(NONEXISTENT_VARIABLE)
ES := $S$S
Q := "
PERCENT := %%
TIME := stopwatch -v -s c
SAFEMKDIR := cmd /c mkdir
MOVECMD := cmd /c move
# Although del always seems to print an error message, even if /F or /Q or both are supplied, its errorlevel == 1, so maybe the script will just keep working as I generally hope...
DELCMD := del
else
S := /
ES := $S
Q := '
PERCENT := %
TIME := /usr/bin/time
SAFEMKDIR := mkdir -p
MOVECMD := mv
DELCMD := rm -f
endif

# If not specified, use the default ft_reduce executable (which should be in the $PATH, or on Windows, in the current dir).
# We need to give the user an option here because sometimes different machines have different OSes (or different versions of the OS) and need different executables!
FTREDUCE ?= ft_reduce

# If not specified, tell Gurobi to use just 1 thread.  You can specify GRBTHREADS=0 to use the default, which is the number of CPU cores.
GRBTHREADS ?= 1

# Used with the move_results dummy target.
MOVEDIR ?= MOVED

$(warning Detected OS uses <$S> for slashes in pathnames.)

# This recursively-expanded variable is used to run a command and then output information about it (including the command line itself) to a file called $@.run.
# It expects the command to be in $(CMD).  This should be set as a per-target or per-pattern variable assignment.
#RECORDCMD = $(CMD) && echo command_line=$(subst <,\<,$(subst >,\>,$(CMD))) > $@.run && echo job_id=$$JOB_ID >> $@.run && echo hostname=`hostname` >> $@.run && echo completed_time=`date` >> $@.run && echo cwd=`pwd` >> $@.run && env >> $@.run
#RECORDCMD = $(CMD) && echo command_line=$(subst <,\<,$(subst >,\>,$(CMD))) > $@.run && echo job_id=$$JOB_ID >> $@.run && echo hostname=`hostname` >> $@.run && echo completed_time=`date` >> $@.run && echo cwd=`pwd` >> $@.run && env >> $@.run
RECORDCMD = $(CMD) && echo 'command_line=$(CMD:'='"'"')' > $@.run && echo job_id=$$JOB_ID >> $@.run && echo hostname=`hostname` >> $@.run && echo completed_time=`date` >> $@.run && echo cwd=`pwd` >> $@.run && env >> $@.run

PROBLEMSPAT := pos44247*.txt
PROBLEMS := $(subst /,$S,$(wildcard $(subst $S,/,$(PROBLEMSPAT))))
#$(info PROBLEMSPAT=<$(PROBLEMSPAT)>)
#$(info $$(subst $$S,/,$$(PROBLEMSPAT))=<$(subst $S,/,$(PROBLEMSPAT))>)
#$(info $$(wildcard $$(subst $$S,/,$$(PROBLEMSPAT)))=<$(wildcard $(subst $S,/,$(PROBLEMSPAT)))>)
#$(info PROBLEMS=<$(PROBLEMS)>)

#HACK: This is so ugly.
ifneq ($(FIXCOLOURS),)
PROBLEMS := $(PROBLEMS:%__COMBINED__.txt=%__COMBINED__fixed_colours.txt)
endif

UNREDUCEDPROBLEMS := $(filter-out %.reduced.txt,$(PROBLEMS))

ifneq ($(REDUCE),)
# Be careful to avoid creating abc.reduced.reduced.reduced.txt etc...
REDUCEDPROBLEMS := $(UNREDUCEDPROBLEMS:%.txt=%.reduced.txt)
PROBLEMS := $(UNREDUCEDPROBLEMS) $(REDUCEDPROBLEMS)
endif

SORTEDREDUCEDPROBLEMS := $(addsuffix .sorted,$(filter %.reduced.txt,$(PROBLEMS)))
#$(info SORTEDREDUCEDPROBLEMS=<$(SORTEDREDUCEDPROBLEMS)>)

ILPS := $(PROBLEMS:.txt=.ilp)
GRBSOLUTIONS := $(ILPS:.ilp=.grb.sol)
LPSSOLUTIONS := $(ILPS:.ilp=.lps.sol)
SOLUTIONS := $(GRBSOLUTIONS) $(LPSSOLUTIONS)

GRBCANONS := $(GRBSOLUTIONS:=.canon)
LPSCANONS := $(LPSSOLUTIONS:=.canon)

SUMMARIES := $(UNREDUCEDPROBLEMS:.txt=.summary.md5)
#ALLRESULTFILES := $(REDUCEDPROBLEMS) $(ILPS) $(SOLUTIONS) $(GRBCANONS) $(LPSCANONS) $(SUMMARIES)
ALLRESULTFILES := $(REDUCEDPROBLEMS) $(REDUCEDPROBLEMS:.txt=.txt.ftr.log) $(SORTEDREDUCEDPROBLEMS) $(ILPS) $(SOLUTIONS) $(SOLUTIONS:.sol=.log) $(GRBCANONS) $(LPSCANONS) $(SUMMARIES)
#HACK: The $(wildcard ...) call selects just those output files that *already exist*, so we don't go building anything just to move it...
#MOVEDRESULTFILES := $(addprefix $(MOVEDIR)$S,$(ALLRESULTFILES))
MOVEDRESULTFILES := $(addprefix $(MOVEDIR)$S,$(subst /,$S,$(wildcard $(subst $S,/,$(ALLRESULTFILES)))))
#$(info ALLRESULTFILES=<$(ALLRESULTFILES)>)
$(info MOVEDRESULTFILES=<$(MOVEDRESULTFILES)>)

# Default rule
all: all_gurobi all_lp_solve

all_summaries: $(SUMMARIES)

all_gurobi: $(GRBCANONS)

all_lp_solve: $(LPSCANONS)

sol_gurobi: $(GRBSOLUTIONS)

sol_lp_solve: $(LPSSOLUTIONS)

all_ilps: $(ILPS)

all_reductions: $(REDUCEDPROBLEMS)

all_colour_fixes: $(PROBLEMS)

all_reductions_sorted: $(SORTEDREDUCEDPROBLEMS)

#move_results:
#	-mkdir $(addprefix $(MOVEDIR),$(sort $(dir $(ALLRESULTFILES))))
move_results: $(MOVEDRESULTFILES)

#HACK: Ugh.
#$(MOVEDIR)%: %
#MOVED$S%: %
#MOVED\%: %
#MOVED\\%: %
#$(MOVEDIR)$S%: %
$(MOVEDIR)$(ES)%: %
	-$(SAFEMKDIR) $(dir $@)
	$(MOVECMD) $< $@
#	$(MOVECMD) $(basename $<)*.log $(dir $@)

# All remaining rules are implicit rules, not static pattern rules, so that (a) they can be applied to any given input file, and
# (b) you can apply them to .ilp files without needing to produce that .ilp file from a precursor graph (.txt) file.

# We can't actually use the $$ shell variable sensibly because make starts a new shell for each recipe line!
# We could dig out make's own $$, but that would break if we tried to use make -j.  Using $@ in a recursively-expanded variable is simple and 100% safe :)
#GRBWORKDIR := grb_work_dir.`hostname`.$$$$
GRBWORKDIR = $@.grb_work_dir
GRBTIMELIMIT := 7200
# In GB.  The only way you can actually tell Gurobi to limit its memory usage is with the NodefileStart parameter, apparently.
GRBWORKMEM := 11

# $(RECORDCMD) will run the command line in $(CMD) and save info about it to $@.run.
# Gurobi will still output a (generally suboptimal) solution if it runs out of time, so detect that case by grepping the logfile and (for debugging purposes) moving any .sol file to .sol.SUBOPT.
# I considered writing a one-line message to the .sol file to prevent a second run of make wasting time, but this is more likely to lead to problems than save time.
# Also don't forget to delete any existing logfile first, since Gurobi always APPENDS!
%.grb.sol: CMD=$(TIME) gurobi_cl ResultFile=$@ LogFile=$(@:.sol=.log) Threads=$(GRBTHREADS) TimeLimit=$(GRBTIMELIMIT) NodefileStart=$(GRBWORKMEM) NodefileDir=$(GRBWORKDIR) $< > $(@:.sol=.stdout) 2> $(@:.sol=.stderr)
%.grb.sol: %.ilp
#	$(TIME) gurobi_cl ResultFile=$@ LogFile=$(@:.sol=.log) $<
#	$(TIME) gurobi_cl ResultFile=$@ LogFile=$(@:.sol=.log) Threads=$(GRBTHREADS) $< > $(@:.sol=.stdout) 2> $(@:.sol=.stderr)
	$(SAFEMKDIR) $(GRBWORKDIR)
	$(DELCMD) $(@:.sol=.log)
#	$(TIME) gurobi_cl ResultFile=$@ LogFile=$(@:.sol=.log) Threads=$(GRBTHREADS) TimeLimit=$(GRBTIMELIMIT) NodefileStart=$(GRBWORKMEM) NodefileDir=$(GRBWORKDIR) $< > $(@:.sol=.stdout) 2> $(@:.sol=.stderr)
	$(RECORDCMD)
	grep -q "^Optimal solution found" $(@:.sol=.log) || ($(MOVECMD) $@ $@.SUBOPT && exit 1)
	-rmdir $(GRBWORKDIR)

%.lps.sol: %.ilp
	$(TIME) lp_solve -rxli xli_CPLEX $< > $@ 2> $(@:.sol=.log)

# lp_solve's -presolve actually made things slightly slower, but what the hell...
%.lps.-presolve.sol: %.ilp
	$(TIME) lp_solve -presolve -rxli xli_CPLEX $< > $@ 2> $(@:.sol=.log)

# Kai's latest ilp.jar (from his email at 1:20am on 14/8/2013) actually produces output that is directly usable by both Gurobi and lp_solve,
# so I don't need to do any further processing.
# 13/10/2013: Kai's NONcombined instances now have a special "x + y" line to hold their score, which gives valuable info but
# causes ilp.jar to choke.  Strip it out en route.
# WTJW 21/3/2014: Java by default grabs more than the 4GB of RAM that SGE allows by default!
#%.bad_ilp: %.txt
%.ilp: %.txt
#	java -jar ilp.jar < $< > $@
#	perl -lpe "s/\s*\+.*//" < $< | java -jar ilp.jar > $@
	perl -lpe "s/\s*\+.*//" < $< | java -Xmx2G -jar ilp.jar > $@

# Rules for reducing problems
%.reduced.txt: %.txt
#	$(TIME) $(FTREDUCE) --reduce-mode vub --vertex-ubs --vertex-seb-ubs < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "* ( tim-vertex-ubs seb-vertex-ubs reduce-vub ) reduce-unreach unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "* ( tim-vertex-ubs seb-vertex-ubs reduce-vub ) reduce-unreach * ( seb-vertex-ubs tim-vertex-ubs calc-rec-slide-lbs reduce-slide reduce-domedge ) unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "* ( clear-vertex-ubs seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach seb-vertex-ubs tim-vertex-ubs clear-slide-lbs calc-rec-slide-lbs reduce-slide clear-slide-lbs calc-rec-slide-lbs reduce-domedge ) unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "* ( clear-vertex-ubs seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach seb-vertex-ubs tim-vertex-ubs clear-slide-lbs calc-rec-slide-lbs reduce-dompath ) unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "* ( clear-vertex-ubs seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach seb-vertex-ubs tim-vertex-ubs clear-slide-lbs calc-rec-slide-lbs calc-implied-anc reduce-dompath ) unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "enable-seb-vub-strength * ( clear-vertex-ubs seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach seb-vertex-ubs tim-vertex-ubs clear-slide-lbs calc-rec-slide-lbs calc-implied-anc reduce-dompath ) reduce-combcol unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "enable-seb-vub-strength * ( * ( clear-vertex-ubs seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) clear-slide-lbs calc-rec-slide-lbs calc-implied-anc reduce-dompath ) reduce-combcol unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "enable-seb-vub-strength * ( * ( clear-vertex-ubs seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) clear-slide-lbs calc-rec-slide-lbs calc-implied-anc reduce-dompath calc-implied-edges ) reduce-combcol unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "enable-seb-vub-strength * ( * ( clear-vertex-ubs seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) clear-slide-lbs calc-rec-slide-lbs calc-implied-anc reduce-slide reduce-domedge calc-implied-edges ) reduce-combcol unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) --time "read renumber-verts enable-seb-vub-strength * ( * ( clear-vertex-ubs seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) clear-slide-lbs calc-rec-slide-lbs calc-implied-anc reduce-slide reduce-domedge calc-implied-edges ) reduce-combcol unrenumber-verts write" < $< > $@ 2> $@.ftr.log
	$(TIME) $(FTREDUCE) --compare-mode ascasc --time "read renumber-verts enable-seb-vub-strength * ( * ( seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) calc-rec-slide-lbs DEBUG-calc-anchor-lbs reduce-slide-strong reduce-unreach reduce-dompath2 calc-implied-edges ) reduce-combcol unrenumber-verts write" < $< > $@ 2> $@.ftr.log
#	$(TIME) $(FTREDUCE) read renumber-verts "enable-seb-vub-strength * ( * ( clear-vertex-ubs seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) clear-slide-lbs calc-rec-slide-lbs calc-implied-anc reduce-slide reduce-domedge calc-implied-edges ) reduce-combcol unrenumber-verts write" < $< > $@

# Specific reductions for the paper
%.reduced1.txt: CMD=$(TIME) $(FTREDUCE) --compare-mode ascasc --time "read renumber-verts * ( tim-vertex-ubs reduce-vub reduce-unreach ) unrenumber-verts write" < $< > $@ 2> $@.ftr.log
%.reduced1.txt: %.txt
	$(RECORDCMD)

%.reduced2.txt: CMD=$(TIME) $(FTREDUCE) --compare-mode ascasc --time "read renumber-verts * ( seb-vertex-ubs reduce-vub reduce-unreach ) unrenumber-verts write" < $< > $@ 2> $@.ftr.log
%.reduced2.txt: %.txt
	$(RECORDCMD)

%.reduced3.txt: CMD=$(TIME) $(FTREDUCE) --compare-mode ascasc --time "read renumber-verts * ( seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) unrenumber-verts write" < $< > $@ 2> $@.ftr.log
%.reduced3.txt: %.txt
	$(RECORDCMD)

%.reduced4.txt: CMD=$(TIME) $(FTREDUCE) --compare-mode ascasc --time "read renumber-verts enable-seb-vub-strength * ( seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) unrenumber-verts write" < $< > $@ 2> $@.ftr.log
%.reduced4.txt: %.txt
	$(RECORDCMD)

%.reduced5.txt: CMD=$(TIME) $(FTREDUCE) --compare-mode ascasc --time "read renumber-verts enable-seb-vub-strength * ( * ( seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) calc-anc-col-lbs calc-rec-slide-lbs DEBUG-calc-anchor-lbs reduce-slide-strong reduce-unreach reduce-dompath2 calc-implied-edges ) reduce-combcol unrenumber-verts write" < $< > $@ 2> $@.ftr.log
%.reduced5.txt: %.txt
	$(RECORDCMD)

%.reduced6.txt: CMD=$(TIME) $(FTREDUCE) --compare-mode ascasc --time "read renumber-verts enable-seb-vub-strength * ( * ( seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) DEBUG-calc-anchor-lbs reduce-colsubtree-adv reduce-unreach calc-anc-col-lbs calc-rec-slide-lbs DEBUG-calc-anchor-lbs reduce-slide-strong reduce-unreach reduce-dompath2 calc-implied-edges ) reduce-combcol unrenumber-verts write" < $< > $@ 2> $@.ftr.log
%.reduced6.txt: %.txt
	$(RECORDCMD)

# Rules for grouping together various combinations of reductions and solvers
# I don't understand why, but these implicit rules will not even be seen unless I
# add an empty recipe line (a single blank) -- even though I didn't have to do that for "all", etc.
%.reduced046.orig_cplex.sol.canon: %.orig_cplex.sol.canon %.reduced4.orig_cplex.sol.canon %.reduced6.orig_cplex.sol.canon
	

%.reduced046.cut_cplex.sol.canon: %.cut_cplex.sol.canon %.reduced4.cut_cplex.sol.canon %.reduced6.cut_cplex.sol.canon
	

%.both_cplex.sol.canon: %.orig_cplex.sol.canon %.cut_cplex.sol.canon
	
#	touch $@

# And some for Gurobi...
%.reduced046.grb.sol.canon: %.grb.sol.canon %.reduced4.grb.sol.canon %.reduced6.grb.sol.canon
	

%.cplex_and_grb.sol.canon: %.both_cplex.sol.canon %.grb.sol.canon
	

# Test new $(RECORDCMD) trick
%.RECORDED_reduced4.txt: CMD=$(TIME) $(FTREDUCE) --compare-mode ascasc --time "read renumber-verts enable-seb-vub-strength * ( seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) unrenumber-verts write" < $< > $@ 2> $@.ftr.log
%.RECORDED_reduced4.txt: %.txt
	$(RECORDCMD)

# Put solutions in a canonical form that can easily be compared with e.g. diff
# Ridiculous $Q and $(PERCENT) needed because Windows and Unix treat ", ' and % characters differently...
#%.sol.canon: %.sol
#	grep -E "^#|[^0-9]1$$" < $< | sort > $@
%.grb.sol.canon: %.grb.sol
#	grep -E "^#|[^0-9]1$$" < $< | sort > $@
#	perl -lne "print if s/^# Objective value = (.+)/sprintf q{# Objective value (to 6 d.p.) = %%.6f}, $$1/e or s/\s+1$$/ 1/" < $< | sort > $@
	perl -lne $Qprint if s/^# Objective value = (.+)/sprintf q{# Objective value (to 6 d.p.) = $(PERCENT).6f}, $$1/e or s/\s+1$$/ 1/$Q < $< | sort > $@

%.lps.sol.canon: %.lps.sol
	perl -lne $Qprint if s/^Value of objective function: (.+)/sprintf q{# Objective value (to 6 d.p.) = $(PERCENT).6f}, $$1/e or s/\s+1$$/ 1/$Q < $< | sort > $@

# NOTE: We don't list any prerequisites because we want to use all pre-existing
# This means that you should run "make all_summaries" only *after* a previous make has finished doing all the actual work you want to do.
%.summary.md5:
	md5sum $*.*.sol.canon > $@

# Kai's combined instances give the super-root colour 0, the same colour as all the sub-roots...  This causes problems for my reduction.
# This also applies to the instances Stephan emailed me on 19/11/2013 -- these contain the right colour for the root, but the colour count is too small by 1!
#HACK: In order for this rule to work in both the current directory and subdirs, it works for ANY prefix...  Could be a feature?  :-P
#%__COMBINED__fixed_colours.txt: %__COMBINED__.txt
%fixed_colours.txt: %.txt
	perl recolour_root_vertex.pl < $< > $@

# Currently used for canonicalising reduced problems.
%.sorted: %
	sort < $< > $@

# Force make to delete a target if its recipe fails.
.DELETE_ON_ERROR:

# Force make to keep all files (e.g. .ilp files) built during chains of implicit rules.
.SECONDARY:
