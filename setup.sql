create schema frag_tree_reduction;

-- You need to close and reopen the DB connection after the following (or alternatively use SET search_path = frag_tree_reduction, "$user", public).
ALTER DATABASE wtwhite set search_path = frag_tree_reduction, "$user", public;

-- Should probably make other tables have a FK to this, but it's not crucial.
create table reduction_types (
	reduction varchar primary key,
	description varchar
);

insert into reduction_types (reduction, description)
values
('', 'Original problem instance'),
('.reduced3', 'read renumber-verts * ( seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) unrenumber-verts write'),
('.reduced4', 'read renumber-verts enable-seb-vub-strength * ( seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) unrenumber-verts write'),
('.reduced6', 'read renumber-verts enable-seb-vub-strength * ( * ( seb-vertex-ubs tim-vertex-ubs reduce-vub reduce-unreach ) DEBUG-calc-anchor-lbs reduce-colsubtree-adv reduce-unreach calc-anc-col-lbs calc-rec-slide-lbs DEBUG-calc-anchor-lbs reduce-slide-strong reduce-unreach reduce-dompath2 calc-implied-edges ) reduce-combcol unrenumber-verts write');

-- Now we actually encode the parent_path info by joining to the reductions table.  This is more consistent with the solutions table, which also has input_path and output_path fields
-- (although the latter's output_path field does not link to anything in the instances table).
create table instances (
	dataset varchar not null,	-- Comes from the path and basename.  path/to/xyz.reduced4.txt becomes path/to/xyz.
	reduction varchar not null,	-- '' for original instance, '.reduced4' or '.reduced6'.  Must be enough to make this unique in combination with the dataset field.
	data_group varchar,			-- Breaks DB normalisation technically, since it's a function of dataset.  But we need it as a separate field for speed, especially for indexing.
	is_combined boolean,		-- Breaks DB normalisation again, since it's a function of dataset.  Anyway...
	path varchar /*not null*/,	-- Really "just in case".  So all instances must have distinct filenames, but this is no longer the PK.  Include directory paths if necessary.
--	parent_path varchar references instances (path),	-- Non-NULL only for reduced instances
	nverts integer /*not null*/,
	nedges integer /*not null*/,
	ncols integer /*not null*/,
	opt_sol_val numeric,		-- Contains the known optimal solution value (or the first number if 2 are separated by a '+') if this is specified as line 4 in the file
	opt_sol_extra_val numeric,	-- Contains the number after the '+' (the weight of the edge to this subtree) if line 4 specifies a known optimal solution value in 2 parts
	create_time timestamp not null default current_timestamp,
	primary key (dataset, reduction),
	check (data_group = substring(dataset from '([^/]+)/(?:b\.[0-9a-f][0-9a-f]/)?[^/]+$')),		-- We need to allow the bucket to be absent for the fmm1 pseudo-dataset.
	check (is_combined = (dataset ~ '/__COMBINED__fixed_colours$'))
);

-- Only instance rows that have all no "problems" (NULL values in columns I don't want to have NULLs in).  Don't need to check PK fields.
create or replace view ok_instances as
select * from instances
where true
and path is not null
and nverts is not null
and nedges is not null
and ncols is not null;

create or replace view bad_instances as
select * from instances
where false
or path is null
or nverts is null
or nedges is null
or ncols is null;

-- Could use a separate table and even add a FK constraint, but it's unnecessary.
create or replace view datasets_vw as
select distinct dataset as dataset from instances;

-- Slight loss of expressiveness: we don't allow a reduction to operate on an already-reduced dataset.  Big deal.
create table reductions (
	dataset varchar not null,	-- Comes from the path and basename.  path/to/xyz.reduced4.txt becomes path/to/xyz.
	reduction varchar not null,	-- '.reduced4' or '.reduced6'.  Must be enough to make this unique in combination with the dataset field.  No rows for ''.
	input_path varchar /*not null*/,		-- Refers to the path of the instance that was input to the reduction.
	output_path varchar /*not null*/,	-- Identifies the reduced instance.
	job_id varchar /*not null*/,		-- Generally integer, but doesn't have to be.  Could use my own bespoke string IDs if I want.  Mainly used to allow multiple reps.
	command_line varchar /*not null*/,		-- Should contain the reduction "program"
	hostname varchar /*not null*/,
	completed_time timestamp /*not null*/,
	elapsed_secs numeric /*not null*/,
	elapsed_nonio_secs numeric /*not null*/,	-- Includes all time except that used by "read" and "write" actions
	user_secs numeric,
	system_secs numeric,
	percent_cpu numeric,
	max_mem_kb numeric /*not null*/,
	create_time timestamp not null default current_timestamp,
	primary key (dataset, reduction)
);

-- Only reduction rows that have all no "problems" (NULL values in columns I don't want to have NULLs in).  Don't need to check PK fields.
create or replace view ok_reductions as
select * from reductions
where true
and input_path is not null
and output_path is not null
and job_id is not null
and command_line is not null
and hostname is not null
and completed_time is not null
and elapsed_secs is not null
and elapsed_nonio_secs is not null
and max_mem_kb is not null;

create or replace view bad_reductions as
select * from reductions
where false
or input_path is null
or output_path is null
or job_id is null
or command_line is null
or hostname is null
or completed_time is null
or elapsed_secs is null
or elapsed_nonio_secs is null
or max_mem_kb is null;

-- Just contains the info common to all solution methods.
create table solutions (
	dataset varchar not null,	-- Comes from the path and basename.  path/to/xyz.reduced4.txt becomes path/to/xyz.
	reduction varchar not null,	-- '', '.reduced4' or '.reduced6'.  Must be enough to make this unique in combination with the dataset field.  No rows for ''.
	solver varchar not null,	-- 'orig_cplex', 'cut_cplex' or 'grb'.
	input_path varchar /*not null*/,		-- Refers to the path of the instance that was solved.
	output_path varchar /*not null*/,	-- Name of the .sol (not .sol.canon) file produced.
	job_id varchar /*not null*/,		-- Generally integer, but doesn't have to be.  Could use my own bespoke string IDs if I want.  Mainly used to allow multiple reps.
	command_line varchar /*not null*/,		-- May optionally have a " #" followed by e.g. CPLEX/Gurobi version string or other guff to distinguish e.g. runs with different parameter settings.
	hostname varchar /*not null*/,
	completed_time timestamp /*not null*/,
	nthreads integer /*not null*/,
--	found_opt boolean /*not null*/,
	result_type varchar /*not null*/,		-- 'OPT' if optimal, otherwise an error message in any format
	orig_result_type varchar /*not null*/,	--HACK: The original value of result_type, before I converted '[CPLEX solutionStatusString: integer optimal, tolerance]OPT' to 'OPT'...
	opt_sol_val numeric,				-- Can only be NULL if result_type <> 'OPT'
	opt_bound numeric,					-- Gurobi can report the optimum LB.  Should be equal to opt_sol_val to within numerical error if result_type = 'OPT'.
	elapsed_secs numeric /*not null*/,
	elapsed_nonio_secs numeric /*not null*/,
	ncuts_generated integer /*not null*/,	-- Always 0 for unstrengthened ILPs
	user_secs numeric,
	system_secs numeric,
	percent_cpu numeric,
	max_mem_kb numeric /*not null*/,
	time_limit_secs numeric,
	mem_limit_kb numeric,
	canon_first_line varchar,			-- The first line from the .canon file
	canon_tail_md5 varchar,				-- The hex MD5 hash of the .canon file, except for the first line (which sometimes has extra guff)
	create_time timestamp not null default current_timestamp,
	primary key (dataset, reduction, solver),
	foreign key (dataset, reduction) references instances (dataset, reduction),
	check (result_type <> 'OPT' or opt_sol_val is not null)
);

-- Only solution rows that have all no "problems" (NULL values).  Don't need to check PK fields.  Also note that the result_type must be 'OPT'!
create or replace view ok_solutions as
select * from solutions
where true
and result_type = 'OPT'
and input_path is not null
and output_path is not null
and job_id is not null
and command_line is not null
and hostname is not null
and completed_time is not null
and nthreads is not null
and result_type is not null
and orig_result_type is not null
and elapsed_secs is not null
and elapsed_nonio_secs is not null
and ncuts_generated is not null
and max_mem_kb is not null
and canon_first_line is not null
and canon_tail_md5 is not null;

create or replace view bad_solutions as
select * from solutions
where false
or result_type <> 'OPT'
or input_path is null
or output_path is null
or job_id is null
or command_line is null
or hostname is null
or completed_time is null
or nthreads is null
or result_type is null
or orig_result_type is null
or elapsed_secs is null
or elapsed_nonio_secs is null
or ncuts_generated is null
or max_mem_kb is null
or canon_first_line is null
or canon_tail_md5 is null;

---- Currently not actually used, since we do all solving via Stephan's solve_cplex or solve_gurobi programs.
--create table solutions_gurobi (
--	output_path varchar not null,	-- Name of the .sol (not .sol.canon) file produced.
--	job_id varchar not null,		-- Generally integer, but doesn't have to be.  Could use my own bespoke string IDs if I want.  Mainly used to allow multiple reps.
--	solver varchar not null,		-- e.g. Gurobi version string.  Can add extra guff to distinguish e.g. runs with different parameter settings.
--	nrows_before_presolve integer,		-- All remaining items are only specified for solvers that report this info (e.g. Gurobi)
--	ncols_before_presolve integer,
--	nnonzeros_before_presolve integer,
--	nrows_after_presolve integer,
--	ncols_after_presolve integer,
--	nnonzeros_after_presolve integer,
--	nnodes_explored integer,
--	tolerance numeric,
--	opt_bound numeric,
--	primary key (output_path, job_id),
--	foreign key (output_path, job_id) references solutions (output_path, job_id)
--);


-- This is getting paper-specific, now
create or replace view nice_instances as
select i.*
from ok_instances i
join nice_data_groups g using (data_group)
where dataset like 'fmm1/%';


-- Tidying up after scrape.pl
update solutions
set result_type = 'OPT'
where result_type = '[CPLEX solutionStatusString: integer optimal, tolerance]OPT';


-- Teardown
begin;
drop view nice_instances;
drop view ok_solutions, ok_reductions, ok_instances;
drop view bad_solutions, bad_reductions, bad_instances;
drop view datasets_vw;
drop table solutions;
drop table reductions;
drop table instances;
drop table tmp_solutions;
drop table tmp_reductions;
drop table tmp_instances;
