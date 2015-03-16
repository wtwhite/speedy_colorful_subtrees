select
	regexp_replace(i.dataset, '^.*/mpos', 'm') as "ID",
	i.nverts as "Vertices",
	i.ncols as "Colors",
	i.nedges as "Edges",
	r3.nedges as "R1 Edges",
	r4.nedges as "R2 Edges",
	r6.nedges as "R3 Edges",
	slow.elapsed_nonio_secs as "o Time (s)",
	fastred.elapsed_nonio_secs + fast.elapsed_nonio_secs as "cR2 Time (s)"
from instances i
join instances r3 using (dataset)
join instances r4 using (dataset)
join instances r6 using (dataset)
join solutions slow using (dataset)
join solutions fast using (dataset)
join reductions fastred using (dataset)
where i.reduction = ''
and i.dataset like 'graphs100/b.__/mpos%'
and r3.reduction = '.reduced3'
and r4.reduction = '.reduced4'
and r6.reduction = '.reduced6'
and slow.solver = 'orig_cplex'
and slow.reduction = ''
and fast.solver = 'cut_cplex'
and fast.reduction = '.reduced3'
and fastred.reduction = '.reduced3'
order by 1;
