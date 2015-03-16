with data as (
	select
		solver,
		reduction,
		sum(coalesce(r.elapsed_secs, 0)) as red_secs,
		sum(s.elapsed_secs) as ilp_secs,
		sum(coalesce(r.elapsed_secs, 0) + s.elapsed_secs) as tot_secs,		/* This is actually the total time, but GLE needs this title for the key */
		min(coalesce(r.elapsed_secs, 0) + s.elapsed_secs) as min_secs,
		max(coalesce(r.elapsed_secs, 0) + s.elapsed_secs) as max_secs,
		count(*)
	from ok_solutions s
	left join ok_reductions r using (dataset, reduction)
	where dataset like 'graphs100/b.__/mpos%'
	group by 1, 2
	/*order by case when solver = 'orig_cplex' then '1' when solver = 'cut_cplex' then '2' else solver end, 2*/
	/*order by regexp_replace(regexp_replace(substring(solver, 1, 1) || case when reduction = '' then '' else right(reduction, 1) end, '^g', 'zzz'), '^o', 'aaa');*/
)
select
	case
		when solver = 'orig_cplex' then 'CPLEX'
		when solver = 'cut_cplex' then 'CPLEX + Cuts'
		when solver = 'grb' then 'Gurobi'
	end as solver,
	r0.ilp_secs as r0_ilp_secs,
	r0.tot_secs as r0_tot_secs,
	r3.ilp_secs as r3_ilp_secs,
	r3.tot_secs as r3_tot_secs,
	r4.ilp_secs as r4_ilp_secs,
	r4.tot_secs as r4_tot_secs,
	r6.ilp_secs as r6_ilp_secs,
	r6.tot_secs as r6_tot_secs
from (select * from data where reduction = '') r0
join (select * from data where reduction = '.reduced3') r3 using (solver)
join (select * from data where reduction = '.reduced4') r4 using (solver)
join (select * from data where reduction = '.reduced6') r6 using (solver)
order by case when solver = 'orig_cplex' then '1' when solver = 'cut_cplex' then '2' else solver end, 2;
