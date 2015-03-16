select
	1 - max(cast(r3.nedges as numeric) / r0.nedges) as worst_r3_red,
	1 - max(cast(r4.nedges as numeric) / r0.nedges) as worst_r4_red,
	1 - max(cast(r6.nedges as numeric) / r0.nedges) as worst_r6_red
from ok_instances r0
join ok_instances r3 using (dataset)
join ok_instances r4 using (dataset)
join ok_instances r6 using (dataset)
where dataset like 'graphs100/b.__/mpos%'
and r0.reduction = ''
and r3.reduction = '.reduced3'
and r4.reduction = '.reduced4'
and r6.reduction = '.reduced6';
