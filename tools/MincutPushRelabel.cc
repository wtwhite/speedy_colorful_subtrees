
/* Maximal flow - Push-Relabel algorithm */
/* Queue with gap and global updates */
/* Boris Cherkassky - cher@theory.stanford.edu, on.cher@zib-berlin.de */
/* Andrew V. Goldberg - goldberg@cs.stanford.edu */
#include "MincutPushRelabel.h"

#include <limits.h>

#define BIGGEST_FLOW (numeric_limits<long>::max())
#define GLOB_UPDT_FREQ 1
#define WHITE 0
#define GREY 1
#define BLACK 2

#define FLOW(a) (cap[a - arcs] - a->r_cap)

/*
 * important: the graph has to be directed, i.e. each pair of adjacent
 * nodes u,v has two edges (u,v) and (v,u)
 *
 */
MincutPushRelabel::MincutPushRelabel(	const ogdf::Graph& g,
										const ogdf::NodeArray<int>& nodeToID,
										const ogdf::EdgeArray<double>& capacities,
										ogdf::node s,
										ogdf::node t)
{
	long n, /* internal number of nodes */
	     node_min, /* minimal no of node  */
	     node_max, /* maximal no of nodes */
	     *arc_first, /* internal array for holding
			    - node degree
			    - position of the first outgoing arc */
	     *arc_tail, /* internal array: tails of the arcs */
	     source, /* no of the source */
	     sink, /* no of the sink */
	     /* temporary variables carrying no of nodes */
	     head, tail, i;

	long m, /* internal number of arcs */
	     /* temporary variables carrying no of arcs */
	     last, arc_num, arc_new_num;

	mf_node *nodes, /* pointer to the node structure */
		*head_p,
		*ndp;

	arc *arcs, /* pointer to the arc structure */
	    *arc_current,
	    *arc_new,
	    *arc_tmp;

	double *acap, /* array of capasities */
	       cap; /* capacity of the current arc */

	long pos_current=0; /* 2*no_alines */

	n = g.numberOfNodes();
	m = g.numberOfEdges();

	//	std::cout << "DheaMaxFlow::DheaMaxFlow(): n==" << n << ", m==" << m << endl << flush;

	/* allocating memory for 'nodes', 'arcs'  and internal arrays */
	nodes = (mf_node*) calloc (n+2, sizeof(mf_node));
	this->nodearray = nodes;
	arcs = (arc*) calloc (2*m+1, sizeof(arc));
	arc_tail = (long*) calloc (2*m, sizeof(long));
	arc_first= (long*) calloc (n+2, sizeof(long));
	acap = (double*) calloc (2*m, sizeof(double));
	/* arc_first [ 0 .. n+1 ] = 0 - initialized by calloc */
	queue = (mf_node**) calloc (n, sizeof (mf_node*));

	if (!nodes || !arcs || !arc_tail || !arc_first || !acap || !queue) {
		std::cerr << "Couldn't allocate memory!" << std::endl;
	}

	/* setting pointer to the first arc */
	arc_current = arcs;

	source = nodeToID[s];
	sink = nodeToID[t];

	node_max = 0;
	node_min = n;
	int my_arc_num = 0;

	ogdf::edge e;
	forall_edges(e, g) {
		//			std::cout << "DheaMaxFlow::DheaMaxFlow(): edge " << e << ", twin == " << twin[e] << flush;

		cap = capacities[e];

		//			std::cout << ", cap==" << cap << ", sistercap==" << sistercap << endl << flush;

		head = nodeToID[e->target()];
		tail = nodeToID[e->source()];

		//			std::cout << "DheaMaxFlow::DheaMaxFlow(): head==" << head << ", tail==" << tail << endl << flush;

		/* no of arcs incident to node i is stored in arc_first[i+1] */
		arc_first[tail + 1] ++;
		arc_first[head + 1] ++;

		/* storing information about the arc */
		arc_tail[pos_current] = tail;
		arc_tail[pos_current+1] = head;
		arc_current->head = nodes + head;
		arc_current->r_cap = cap;
		arc_current->sister = arc_current + 1;
		arc_current->ogdfEdge = e;
		arc_current->ogdfSrc = 0;
		arc_current->ogdfTgt = 0;
		arc_current->rev_edge = false;
		arc_current->head->id = e->target()->index();
		arc_current->id = my_arc_num;

		(arc_current + 1)->head = nodes + tail;
		(arc_current + 1)->r_cap = 0;
		(arc_current + 1)->sister = arc_current;
		(arc_current + 1)->ogdfEdge = 0;
		(arc_current + 1)->ogdfSrc = 0;
		(arc_current + 1)->ogdfTgt = 0;
		(arc_current + 1)->rev_edge = true;
		(arc_current + 1)->head->id = e->source()->index();
		(arc_current + 1)->id = my_arc_num + 1;

		//			std::cout << " edge " << my_arc_num << ":(" << arc_current->head->id << "," << (arc_current + 1)->head->id << ")" << endl << flush;
		//			std::cout << " rev-edge " << my_arc_num+1 << ":(" << arc_current->head->id << "," << (arc_current + 1)->head->id << ")" << endl << flush;

		/* searching minimumu and maximum node */
		if (head < node_min) node_min = head;
		if (tail < node_min) node_min = tail;
		if (head > node_max) node_max = head;
		if (tail > node_max) node_max = tail;

		arc_current += 2;
		pos_current += 2;
		my_arc_num += 2;
	}

	/********** ordering arcs - linear time algorithm ***********/

	/* first arc from the first node */
	(nodes + node_min)->first = arcs;

	/* before below loop arc_first[i+1] is the number of arcs outgoing from i;
	   after this loop arc_first[i] is the position of the first
	   outgoing from node i arcs after they would be ordered;
	   this value is transformed to pointer and written to node.first[i]
	 */

	for (i = node_min + 1; i <= node_max + 1; i ++) {
		arc_first[i] += arc_first[i-1];
		(nodes + i)->first = arcs + arc_first[i];
	}

	for (i = node_min; i < node_max; i ++) { /* scanning all the nodes exept the last*/
		last = ((nodes + i + 1)->first) - arcs;
		/* arcs outgoing from i must be cited
		   from position arc_first[i] to the position
		   equal to initial value of arc_first[i+1]-1  */

		for (arc_num = arc_first[i]; arc_num < last; arc_num ++) {
			tail = arc_tail[arc_num];

			while (tail != i)
				/* the arc no  arc_num  is not in place because arc cited here
				   must go out from i;
				   we'll put it to its place and continue this process
				   until an arc in this position would go out from i */

			{	arc_new_num = arc_first[tail];
				arc_current = arcs + arc_num;
				arc_new = arcs + arc_new_num;

				/* arc_current must be cited in the position arc_new
				   swapping these arcs:                                 */

				head_p = arc_new->head;
				arc_new->head = arc_current->head;
				arc_current->head = head_p;

				cap = arc_new->r_cap;
				arc_new->r_cap = arc_current->r_cap;
				arc_current->r_cap = cap;

				ogdf::edge e = arc_new->ogdfEdge;
				arc_new->ogdfEdge = arc_current->ogdfEdge;
				arc_current->ogdfEdge = e;

				ogdf::node currNode = arc_new->ogdfSrc;
				arc_new->ogdfSrc = arc_current->ogdfSrc;
				arc_current->ogdfSrc = currNode;

				currNode = arc_new->ogdfTgt;
				arc_new->ogdfTgt = arc_current->ogdfTgt;
				arc_current->ogdfTgt = currNode;

				int rev_edge = arc_new->rev_edge;
				arc_new->rev_edge = arc_current->rev_edge;
				arc_current->rev_edge = rev_edge;

				if (arc_new != arc_current->sister) {
					arc_tmp = arc_new->sister;
					arc_new->sister = arc_current->sister;
					arc_current->sister = arc_tmp;

					(arc_current->sister)->sister = arc_current;
					(arc_new->sister)->sister = arc_new;
				}

				arc_tail[arc_num] = arc_tail[arc_new_num];
				arc_tail[arc_new_num] = tail;

				/* we increase arc_first[tail]  */
				arc_first[tail] ++;

				tail = arc_tail[arc_num];
			}
		}
		/* all arcs outgoing from  i  are in place */
	}

	/* -----------------------  arcs are ordered  ------------------------- */

	/*----------- constructing lists ---------------*/

	for (ndp = nodes + node_min; ndp <= nodes + node_max; ndp ++)
		ndp->first = (arc*) NULL;

	for (arc_current = arcs + (2*m-1); arc_current >= arcs; arc_current --) {
		arc_num = arc_current - arcs;
		tail = arc_tail [arc_num];
		ndp = nodes + tail;
		arc_current->next = ndp->first;
		ndp->first = arc_current;
	}

	/* ----------- assigning output values ------------*/
	this->narcs = m;
	this->n = node_max - node_min + 1;
	this->nsource = nodes + source;
	this->nsink = nodes + sink;
	this->nodes = nodes + node_min;
	this->arcs = arcs;
	this->cap = acap;
	this->first_node = node_min;

	for (arc_current = arcs, arc_num = 0; arc_num < 2*m; arc_current ++, arc_num ++) {
		acap [ arc_num ] = arc_current->r_cap;
		this->cap [ arc_num ] = arc_current->r_cap;
	}

	/* free internal memory */
	free (arc_first);
	free (arc_tail);

}

MincutPushRelabel::~MincutPushRelabel()
{
	if (nodearray) free(nodearray);
	if (arcs) free(arcs);
	if (cap) free(cap);
	if (queue) free(queue);
}


void MincutPushRelabel::update(const ogdf::EdgeArray<double>& capacities, int s, int t)
{
	nsource = nodes + s;
	nsink = nodes + t;

	arc* arc_current;
	int arc_num;

	for (arc_current = arcs, arc_num = 0;
			arc_num < 2*this->narcs;
			arc_current ++, arc_num ++)
	{
		if (arc_current->ogdfEdge == 0)
			arc_current->r_cap = 0;
		else
			arc_current->r_cap = capacities[arc_current->ogdfEdge];
		this->cap [ arc_num ] = arc_current->r_cap;
	}
}


double
MincutPushRelabel::max_flow(double r_cap[], double eps)
{
	int i=0;
	double f = 0;
	int e=0;

	if ((e = prflow(&f, eps))) {
		std::cout << "Error in max_flow:" << e << std::endl;
		return 0;
	}

	for (i=0; i<narcs; i++)
		r_cap[i] = 0;//arcs[i].r_cap;

	return f;
}

double
MincutPushRelabel::min_cut(double border, int cut[], double eps)
{
	double f = 0;
	int e=0;

	if ((e = prflow(&f, eps))) {
		std::cout << "Error in max_flow:" << e << std::endl;
		return 0;
	}

	for (int j=0; j<n; j++)
		cut[j] = 0;

	if (f < border) {
		ogdf::List<int> Q;

		Q.pushBack(nsource-nodes);
		cut[nsource - nodes + first_node] = 1;

		while (!Q.empty()) {
			mf_node akt = nodes[Q.popFrontRet()];

			arc * a = akt.first;

			while (a) {
				//  cout <<  cap[a - arcs] << "/" << a->r_cap << " ";
				//if (cap[a - arcs] > f && !cut[a->head - nodes + first_node])
				if (a->r_cap > 0 && cut[a->head - nodes + first_node] == 0) {
					cut[a->head - nodes + first_node] = 1;
					Q.pushBack(a->head - nodes);
				}
				a = a->next;
			}
		}

		Q.pushBack(nsink-nodes);
		cut[nsink - nodes + first_node] = 2;

		while (!Q.empty()) {
			mf_node akt = nodes[Q.popFrontRet()];

			arc * a = akt.first;

			while (a) {
				//  cout <<  cap[a - arcs] << "/" << a->r_cap << " ";
				//if (cap[a - arcs] > f && !cut[a->head - nodes + first_node])
				if (a->sister->r_cap > 0 && cut[a->head - nodes + first_node] == 0) {
					cut[a->head - nodes + first_node] = 2;
					Q.pushBack(a->head - nodes);
				}
				a = a->next;
			}
		}
		//    cout << endl;
	}

	return f;
}


/*--- initialization */

int MincutPushRelabel::pr_init ()
{
	mf_node  *i;        /* current node */
	/*
	   arc * a;

	   for (i = nodes; i < nodes + n; i++)
	   {
	   cout << "node " << (i-nodes) << endl;
	   for (a=i->first; a; a=a->next)
	   cout << a->r_cap << endl;
	   }
	 */

	for (i = nodes; i < nodes + n; i++)
		i->excess = 0;

	nsource->excess = BIGGEST_FLOW;

	lmax = n-1;

	return 0;

} /* end of initialization */


/*--- global rank update - breadth first search */

void MincutPushRelabel::def_ranks ()
{
	mf_node  *i, *j; //, *jn;  /* current nodes */
	arc   *a;           /* current arc   */
	long  j_rank;       /* rank of node j */

	n_up ++; /* statistics */

	/* initialization */

	for (i = nodes; i < nodes + n; i ++)
		i   ->rank = n;

	nsink->rank = 0;

	*queue = nsink;

	qp_first = qp_last = NULL;

	lmax = 0;

	/* breadth first search */

	for (q_read = queue, q_write = queue + 1;
			q_read != q_write;
			q_read ++
	   )
	{ /* scanning arcs incident to node i */

		i = *q_read;
		j_rank = (i->rank) + 1;

		for (a = i->first; a != NULL; a = a->next) {
			j = a->head;

			if (j->rank == n)
				/* j is not labelled */

				if (((a->sister)->r_cap) > 0)
				{ /* arc (j, i) is not saturated */

					j->rank    = j_rank;
					j->current = j->first;

					if (j_rank > lmax) lmax = j_rank;

					if ((j->excess) > 0) {
						j->q_next     = qp_first;
						if (qp_first == NULL) qp_last = j;
						qp_first        = j;
					}

					*q_write = j; q_write ++; /* put j  to scanning queue */
				}
		} /* node "i" is scanned */
	} /* end of scanning queue */
} /* end of global update */

/*--- pushing flow from node  i  */

int MincutPushRelabel::push (mf_node * i)
{
	mf_node  *j;                /* sucsessor of i */
	long  j_rank;            /* rank of the next layer */
	arc   *a;                /* current arc (i,j) */
	double  fl;                /* flow to push through the arc */

	j_rank = (i->rank) - 1;

	/* scanning arcs outgoing from  i  */

	for (a = i->current; a != NULL; a = a->next) {
		if (a->r_cap > 0) { /* "a" is not saturated */
			j = a->head;

			if (j->rank == j_rank)
			{ /* j belongs to the next layer */

				fl = min(i->excess, a->r_cap);

				a            ->r_cap -= fl;
				(a->sister)->r_cap += fl;
				n_push ++; /* statistics */

				if (j_rank > 0) {
					if (j->excess == 0)
					{ /* before current push  j  had zero excess */

						/* put  j  to the push-list */

						if (qp_first != NULL)
							qp_last->q_next = j;
						else
							qp_first = j;


						qp_last = j;
						j->q_next = NULL;

					} /* j->excess == 0 */

				} /* j->rank > 0 */

				j->excess += fl;
				i->excess -= fl;

				if (i->excess == 0) break;

			} /* j belongs to the next layer */
		} /* a  is not saturated */
	} /* end of scanning arcs from  i */

	i->current = a;

	return (a == NULL) ? 1 : 0;
} /* end of push */

/*--- relabelling node i */

long MincutPushRelabel::relabel (mf_node * i)
{
	mf_node  *j;        /* sucsessor of i */
	long  j_rank;    /* minimal rank of a node available from j */
	arc   *a;        /* current arc */
	arc   *a_j=NULL;      /* an arc which leads to the node with minimal rank */

	n_rel ++; /* statistics */

	i->rank = j_rank = n;

	/* looking for a node with minimal rank available from i */

	for (a = i->first; a != NULL; a = a->next) {
		if (a->r_cap > 0) {
			j = a->head;

			if (j->rank < j_rank) {
				j_rank = j->rank;
				a_j    = a;
			}
		}
	}

	j_rank++;
	if (j_rank < n) {
		/* siting  i  into the manual */

		i->rank    = j_rank;
		i->current = a_j;

		if (j_rank > lmax) lmax = j_rank;

	} /* end of j_rank < n */

	return j_rank;
} /* end of relabel */


/*--- organizer */

int MincutPushRelabel::prflow (double * fl, double eps)
{
	mf_node   *i;          /* current node */
	long   i_rank;      /* rank of  i */
	long   n_r;         /* the number of relabels */
	int    cc;          /* condition code */

	//int ii=0;

	cc = pr_init ();

	if (cc) return cc;

	def_ranks ();

	n_r = 0;

	/* queue method */

	while (qp_first != NULL) { /* main loop */
		/* checking the necessity of global update */
		if (n_r > GLOB_UPDT_FREQ * (float) n)
		{ /* it is time for global update */
			def_ranks ();
			n_r = 0;
			if (qp_first == NULL) break;
		}

		i = qp_first;
		qp_first = qp_first->q_next;
		if (qp_first == NULL) qp_last = NULL;

		i_rank = i->rank;

		while (i_rank < n) {
			/* untill i will be free from excess or beyond the gap */

			cc = push (i);
			if (cc == 0) break;

			/* i must be relabeled */

			i_rank = relabel (i);
			n_r ++;

		} /* end of scanning i */

	} /* end of the main loop */

	*fl += nsink->excess;

	prefl_to_flow (eps);

	return 0;
} /* end of constructing flow */


/*--- removing excessive flow - second phase of PR-algorithm */
void MincutPushRelabel::prefl_to_flow (double eps)
/*
   do dsf in the reverse flow graph from nodes with excess
   cancel cycles if found
   return excess flow in topological order
*/

/*
   rank is used for dfs labels
   nl_prev is used for DSF tree
   q_next is used for topological order list
*/
{
	mf_node *i, *j, *tos=NULL, *bos, *restart, *r;
	arc *a;
	double delta;

	/* deal with self-loops */
	for (i = nodes; i < nodes + n; i++) {
		for (a = i->first; a != NULL; a = a->next)
			if (a->head == i) {
				a->r_cap = cap[a - arcs];
			}
	}

	/* initialize */
	bos = NULL;
	for (i = nodes; i < nodes + n; i++) {
		i->rank = WHITE;
		i->nl_prev = NULL;
		i->current = i->first;
	}

	for (i = nodes; i < nodes + n; i++)
		if ((i->rank == WHITE) && (i->excess > 0) &&
				(i != nsource) && (i != nsink)) {
			r = i;
			r->rank = GREY;
			do {
				for (; i->current != NULL; i->current = i->current->next) {
					a = i->current;
					if ((cap[a - arcs] == 0) && (a->r_cap > 0) &&
							(a->head != nsource) && (a->head != nsink)) {
						j = a->head;
						if (j->rank == WHITE) {
							/* start scanning j */
							j->rank = GREY;
							j->nl_prev = i;
							i = j;
							break;
						}
						else
							if (j->rank == GREY) {
								/* find minimum flow on the cycle */
								delta = a->r_cap;
								while (1) {
									delta = min(delta, j->current->r_cap);
									if (j == i)
										break;
									else
										j = j->current->head;
								}

								/* remove delta flow units */
								j = i;
								while (1) {
									a = j->current;
									a->r_cap -= delta;
									a->sister->r_cap += delta;
									j = a->head;
									if (j == i)
										break;
								}

								/* back DFS to the first zerod arc */
								restart = i;
								for (j = i->current->head; j != i; j = a->head) {
									a = j->current;
									if ((j->rank == WHITE) || (a->r_cap == 0)) {
										j->current->head->rank = WHITE;
										if (j->rank != WHITE)
											restart = j;
									}
								}

								if (restart != i) {
									i = restart;
									i->current = i->current->next;
									break;
								}
							}
					}
				}

				if (i->current == NULL) {
					/* scan of i complete */
					i->rank = BLACK;
					if (i != nsource) {
						if (bos == NULL) {
							bos = i;
							tos = i;
						}
						else {
							i->q_next = tos;
							tos = i;
						}
					}

					if (i != r) {
						i = i->nl_prev;
						i->current = i->current->next;
					}
					else
						break;
				}
			} while (1);
		}


	/* return excesses */
	/* note that sink is not on the stack */
	if (bos != NULL) {
		i = tos;
		do {
			a = i->first;

			if (a) {
				while (i->excess > eps) {
					if (!a) {
						// epsilon was too small, just fix it
						i->excess = 0;
						break;
					}

					if ((cap[a - arcs] <= eps) && (a->r_cap > eps)) {
						delta = min(i->excess, a->r_cap);
						a->r_cap -= delta;
						a->sister->r_cap += delta;
						i->excess -= delta;
						a->head->excess += delta;
					}
					a = a->next;
				}
			}

			if (i == bos)
				break;
			else
				i = i->q_next;
		} while (1);
	}
}
