select
	row_number() over (order by r0.nedges) as x,		/*HACK: This ORDER BY clause is pretty disgusting, but "ORDER BY 2" has no effect! */
	r0.nedges as "\Delta\ edges R1 \rightarrow\ unreduced",
	r3.nedges as "\Delta\ edges R2 \rightarrow\ R1",
	r4.nedges as "\Delta\ edges R3 \rightarrow\ R2",
	r6.nedges as "# edges after R3",
	r0.ncols
from ok_instances r0
join ok_instances r3 using (dataset)
join ok_instances r4 using (dataset)
join ok_instances r6 using (dataset)
where dataset like 'graphs100/b.__/mpos%'
and r0.reduction = ''
and r3.reduction = '.reduced3'
and r4.reduction = '.reduced4'
and r6.reduction = '.reduced6'
order by r0.nedges;
