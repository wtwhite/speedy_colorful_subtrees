-- Show the fraction of edges remaining in each reduced problem
select
	o.path,
	cast(r.nedges as numeric) / o.nedges as fraction,
	o.opt_sol_val - r.opt_sol_val		-- Actually this is misleading, because it's just what it says in the input file, not what e.g. Gurobi actually found!
from
	instances r,
	instances o
where o.path = r.parent_path;

-- Show the fraction of edges remaining in each reduced problem, and elapsed time
select
	o.path,
	o.nedges - r.nedges as nedges_deleted,
	cast(o.nedges - r.nedges as numeric) / o.nedges as fraction_deleted,
	red.elapsed_time
from
	instances r,
	instances o,
	reductions red
where o.path = r.parent_path
and red.path = r.path;

-- Compare direct solution times with solution times via reduction.
-- Still gives results for reductions, even if the original hasn't been solved.
-- Also gives results for just the reductions, even if the reduced instance hasn't been solved either.
select
	oi.path,
	so.elapsed_time as total_orig_sol_time,
	oi.nedges - ri.nedges as nedges_deleted,
	cast(oi.nedges - ri.nedges as numeric) / oi.nedges as fraction_deleted,
	red.elapsed_time as reduction_time,
	sr.elapsed_time as red_sol_time,
	red.elapsed_time + sr.elapsed_time as total_red_sol_time
from
	reductions red,
	instances ri
	cross join instances oi				--HACK: This cross join is restricted by WHERE conditions later...
	left join solutions sr on sr.path = ri.path
	left join solutions so on
		so.path = oi.path
		and so.solver = sr.solver		-- Restricts the full join to a more useful subset in the case where we're using multiple solvers
where ri.parent_path = oi.path
and red.path = ri.path;
