/*
 * provides class for separation of directed cuts
 *
 * based on the mincut-code by Cherkassy and Goldberg
 * see MincutPushRelabel.cpp
 */
#ifndef MINCUT_PUSH_RELABEL_H
#define MINCUT_PUSH_RELABEL_H

#include "ogdf/basic/Graph_d.h"
#include "ogdf/basic/EdgeArray.h"
#include "ogdf/basic/NodeArray.h"

typedef  /* arc */
	struct arc_st
{
	double		r_cap;           /* residual capacity */
	ogdf::edge	ogdfEdge;        /* corresponding edge in ogdf graph,
					     0 for backedges */
	bool 		rev_edge;	 /* we don't need this varibale for the
					    STP separation, but for the
					    GSEC-separation */
	ogdf::node	ogdfSrc;	 /* we don't need this varibale for the
					    STP separation, but for the
					    GSEC-separation */
	ogdf::node	ogdfTgt;	 /* we don't need this varibale for the
					    STP separation, but for the
					    GSEC-separation */
	struct node_st   *head;           /* head */
	struct arc_st    *sister;         /* opposite arc */
	struct arc_st    *next;           /* next arc with the same tail */
	int		id;		 /* for debugging, testing, etc.
					  * i.e., not really necessary */
}
  arc;

typedef  /* node */
	struct node_st
{
	arc              *first;           /* first outgoing arc */
	arc              *current;         /* current incident arc */
	double           excess;           /* excess of the node */
	long             rank;             /* distance from the sink */
	struct node_st   *q_next;          /* next node in queue */
	struct node_st   *nl_prev;         /* used by prefl_to_flow */
	int 		id;		   /* for debugging, testing, etc.
					    * i.e., not really necessary */
} mf_node;


typedef /* layer */
	struct layer_st
{
	mf_node             *push_first;      /* 1st mf_node with pisitive excess */
	mf_node             *trans_first;     /* 1st mf_node with zero excess */
} layer;

class MincutPushRelabel
{
public:
	MincutPushRelabel(const ogdf::Graph& g,
				const ogdf::NodeArray<int>& nodeToID,
				const ogdf::EdgeArray<double>& capacities,
				ogdf::node s,
				ogdf::node t);

	~MincutPushRelabel();

	void update(const ogdf::EdgeArray<double>& capacities, int s, int t);

	double max_flow(double r_cap[], double eps);
	double min_cut(double border, int cut[], double eps);

protected:
	int pr_init();
	void def_ranks();
	int push(mf_node *i);
	long relabel(mf_node *i);
	int prflow(double *fl, double eps);
	void prefl_to_flow(double eps);

	long      n;                    /* number of nodes */
	long      narcs;                /* number of arcs */
	long      first_node;
	mf_node   *nodes;               /* array of nodes */
	mf_node   *nodearray;           /* pointer to allocated memory */
	arc       *arcs;                /* array of arcs */
	double    *cap;                 /* array of capasities */
	mf_node   *nsource;             /* origin */
	mf_node   *nsink;               /* destination */
	mf_node   **queue;              /* queue for storing nodes */
	mf_node   **q_read, **q_write;  /* queue pointers */
	long      lmax;                 /* maximal layer */
	mf_node   *qp_first, *qp_last;  /* beg and end of push-queue */

	/* Statistics */
	long n_up;
	long n_push;
	long n_rel;
};
#endif
