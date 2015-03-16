#include <ogdf/basic/Graph.h>
#include <ogdf/basic/HashArray.h>

#include <sstream>
#include <map>
#include <set>

#define SAME_COLOR_EDGES_POSSIBLE
//#define REDUCE_GRAPH_WITH_SOURCE_INFORMATION
// FIXME: the code with REDUCE_GRAPH_WITH_SOURCE_INFORMATION is currently not valid,
// try: speedy/100/mpos2978C18H47N11O6PS2_44.txt speedy/100/mpos2978C19H48N9O7P2S_41.txt speedy/100/mpos2978C21H50N6O8P2S_27.txt

using namespace ogdf;
using std::map;
using std::set;
using std::pair;

#include "WeightedNodeColoredGraph.h"

static void
removeUnreachableNodes(WeightedNodeColoredGraph &G)
{
	const node root = G.firstNode();
	NodeArray<bool> seen(G.theGraph(), false);
	List<node> todo;
	todo.pushBack(root);
	seen[root] = true;
	while (!todo.empty()) { // BFS
		node v = todo.popFrontRet();
		adjEntry adj;
		forall_adj(adj, v) {
			node u = adj->twinNode();
			if (adj->theEdge()->source() == v
			 && !seen[u]) {
				todo.pushBack(u);
				seen[u] = true;
			}
		}
	}

	// remove unseen nodes
	int count = 0;
	node next;
	for (node v = G.firstNode(); v; v = next) {
		next = v->succ();
		if (!seen[v]) {
			G.delNode(v);
			++count;
		}
	}
	cout << "  -- removed " << count << " unreachable nodes.\n";
}

// reduction count enum
enum {
	cntIsolatedNode = 0,
	cntNPosSinkEdge,
	cntOneNonPosSourceEdge,
	cntMergeNodes,
	cntParallelEdges,
#ifdef SAME_COLOR_EDGES_POSSIBLE
	cntSameColorEdge,
#endif
#ifdef REDUCE_GRAPH_WITH_SOURCE_INFORMATION
	cntSmallWeightEdge,
#endif
	// cntLimit must be the last one
	cntLimit
};

#ifdef REDUCE_GRAPH_WITH_SOURCE_INFORMATION
static void
reduceByEdge(WeightedNodeColoredGraph &wcG, edge e, int *count)
{
	// NOTE: this reduction rule may only be applied if we know that u is not a root!
	const node u = e->source();
	const node v = e->target();
	double min = numeric_limits<double>::max();
	bool minUsed = false;
	adjEntry uAdj, vAdj;
	forall_adj(uAdj, u) {
		const node w = uAdj->theEdge()->source();
		if (w != u) {
			forall_adj(vAdj, v) {
				if (vAdj->theEdge()->source() == w) { // w -> u && w -> v
					// find min weight(wv) - weight(wu)
					double alpha = wcG.weight(vAdj->theEdge());
					const double beta = wcG.weight(uAdj->theEdge());
					if (beta > 0) {
						alpha -= beta;
					} // alpha := min(alpha, alpha - beta)
					if (min > alpha) {
						min = alpha;
						minUsed = true;
					}
					break;
				}
			}
		}
	}
	if (minUsed && min > wcG.weight(e)) {
		++count[cntSmallWeightEdge];
		wcG.delEdge(e);
	}
}
#endif

static void
reduceParallelEdges(WeightedNodeColoredGraph &wcG, node v, int &count)
{
	// merge parallel edges that are leaving v
	if (v->outdeg() >= 2) {
		NodeArray<edge> hit(wcG.theGraph(), NULL);
		adjEntry nextAdj;
		for (adjEntry adj = v->firstAdj(); adj; adj = nextAdj) {
			nextAdj = adj->succ();
			const edge e = adj->theEdge();
			const node w = e->target();
			if (w != v) { // outgoing edges
				if (!hit[w]) {
					hit[w] = e;
				} else {
					if (wcG.weight(hit[w]) < wcG.weight(e)) {
						wcG.delEdge(hit[w]);
						hit[w] = e;
					} else {
						wcG.delEdge(e);
					}
					++count;
				}
			}
		}
	}
}

static void
reduceByNode(WeightedNodeColoredGraph &wcG, node v, int *count)
{
	// delete non-positive edges to sinks
	if (v->outdeg() == 0) { // v is a sink
		adjEntry next_adj;
		for (adjEntry adj = v->firstAdj(); adj; adj = next_adj) {
			next_adj = adj->succ();
			edge e = adj->theEdge();
			if (wcG.weight(e) <= 0) {
				wcG.delEdge(e);
				++count[cntNPosSinkEdge];
			}
		}
	}

	/* if:
	 *   - the color of v is unique, and
	 *   - all the parents of v that are connected over non-negative edges
	 *     have the same color, and
	 *   - there may be negative edges going into v,
	 * then:
	 *   - we merge v into each parent node that is connected over a
	 *     non-negative edge
	 *     but *keep* v itself, just removing its incoming non-negative
	 *     edges (because their weights are already added to the incoming
	 *     weights of the parent node),
	 *     and
	 *   - we change the color of v to the (unique) color of the parents
	 *     we have merged into (otherwise the color constraint may not be
	 *     satisfied if v is in the tree over a negative edge).
	 *
	 * This is a reduction that may generate parallel edges!
	 */
	if (wcG.color(v) >= 0 && wcG.numberOfEquallyColoredNodes(v) == 0) { // the color is unique
		adjEntry adj;
		List<edge> nonnegEdges;
		bool otherIncomingEdges = false;
		forall_adj(adj, v) {
			const edge e = adj->theEdge();
			if (e->source() == wcG.firstNode()) {
				continue;
			}
			if (e->target() == v) { // incoming edge
				if (wcG.weight(e) >= 0) {
					if (nonnegEdges.empty()
					 || wcG.color(e->source()) == wcG.color(nonnegEdges.front()->source())) {
						nonnegEdges.pushBack(e);
					} else {
						otherIncomingEdges = true;
						break;
					}
				} else { // negative edge
					otherIncomingEdges = true;
					break;
				}
			}
		}
		if (!nonnegEdges.empty() && !otherIncomingEdges) {
			const int color = wcG.color(nonnegEdges.front()->source());
			List<edge>::iterator it2;
			for (List<edge>::iterator it = nonnegEdges.begin(); it.valid(); it = it2) {
				it2 = it.succ();
				edge e = *it;
				// add weight of edge to all parent node's incoming edges
				forall_adj(adj, e->source()) {
					if (adj->theEdge()->target() == adj->theNode()) {
						wcG.weight(adj->theEdge()) += wcG.weight(e);
					}
				}

				// contract (copy outgoing edges of v to e->source())
				bool checkParallel = false;
				forall_adj(adj, v) {
					if (adj->theEdge()->target() != v) {
						wcG.newEdge(e->source(), adj->twinNode(), wcG.weight(adj->theEdge()));
						checkParallel = true;
					}
				}

				// now check if parallel edges occured
				if (checkParallel) {
					reduceParallelEdges(wcG, e->source(), count[cntParallelEdges]);
				}

				// finally remove the old edge
				wcG.delEdge(e);
			}
			wcG.updateColor(v, color);
			++count[cntMergeNodes];
			if (!otherIncomingEdges) {
				wcG.delNode(v);
				return;
			}
		}
	}

	// reduceParallel(wcG, v, count[cntParallelEdges]);

#if 0
	// delete source edges with only one non-positive out-edge
	if (v->indeg() == 0
	 && v->outdeg() == 1
	 && wcG.weight(v) + wcG.weight(v->firstAdj()->theEdge()) <= 0) {
		wcG.delEdge(v->firstAdj()->theEdge());
		++count[cntOneNonPosSourceEdge];
	}
#endif

	// delete isolated nodes
	if (v->degree() == 0) {
		wcG.delNode(v);
		++count[cntIsolatedNode];
	}
}

static void
reduceGraph(WeightedNodeColoredGraph &wcG)
{
	int count[cntLimit];

	cout << "  -- reducing input graph:\n";

	do {
		int maxc = sizeof(count)/sizeof(int);
		memset(count, 0, sizeof(count));

#ifdef SAME_COLOR_EDGES_POSSIBLE
		edge nextEdge;
		for (edge e = wcG.firstEdge(); e; e = nextEdge) {
			nextEdge = e->succ();
			if (wcG.color(e->source()) == wcG.color(e->target())) {
				++count[cntSameColorEdge];
				wcG.delEdge(e);
			}
		}
#endif // SAME_COLOR_EDGES_POSSIBLE

#ifdef REDUCE_GRAPH_WITH_SOURCE_INFORMATION
#ifndef SAME_COLOR_EDGES_POSSIBLE
		edge nextEdge;
#endif
		for (edge e = wcG.firstEdge(); e; e = nextEdge) {
			nextEdge = e->succ();
			if (!wcG.theGraph().searchEdge(wcG.firstNode(), e->source())) {
				reduceByEdge(wcG, e, count);
			}
		}
#endif

		node nextNode;
		for (node v = wcG.firstNode()->succ(); v; v = nextNode) { // omit root (firstNode)
			nextNode = v->succ();
			reduceByNode(wcG, v, count);
		}

#ifdef SAME_COLOR_EDGES_POSSIBLE
		if (count[cntSameColorEdge]) {
			cout
			  << "     -- "
			  << count[cntSameColorEdge]
			  << " edges connecting nodes of the same color removed.\n";
		}
#endif // SAME_COLOR_EDGES_POSSIBLE
#ifdef REDUCE_GRAPH_WITH_SOURCE_INFORMATION
		if (count[cntSmallWeightEdge]) {
			cout
			  << "     -- "
			  << count[cntSmallWeightEdge]
			  << " edges with real small weights removed.\n";
		}
#endif // REDUCE_GRAPH_WITH_SOURCE_INFORMATION
		if (count[cntIsolatedNode]) {
			cout
			  << "     -- "
			  << count[cntIsolatedNode]
			  << " isolated nodes removed.\n";
		}
		if (count[cntNPosSinkEdge]) {
			cout
			  << "     -- "
			  << count[cntNPosSinkEdge]
			  << " negative edges to sinks removed.\n";
		}
		if (count[cntOneNonPosSourceEdge]) {
			cout
			  << "     -- "
			  << count[cntOneNonPosSourceEdge]
			  << " single non-positive source edges removed.\n";
		}
		if (count[cntMergeNodes]) {
			cout
			  << "     -- "
			  << count[cntMergeNodes]
			  << " nodes merged.\n";
		}
		if (count[cntParallelEdges]) {
			cout
			  << "     -- "
			  << count[cntParallelEdges]
			  << " parallel edges removed.\n";
		}

		while (--maxc) {
			count[maxc-1] += count[maxc];
		}
	} while (count[cntIsolatedNode] > 0);

	cout
	  << "  -- reduced graph has "
	  << wcG.numberOfNodes() << " nodes, "
	  << wcG.numberOfEdges() << " edges, "
	  << wcG.numberOfColors() << " colors.\n";
}

static void
usage(char *prog)
{
	cout
	  << "Usage: "
	  << prog << " [-k] <filename> [<node> ...]\n"
	     " Deletes the given nodes and reduces the graph. Write output to output.dat\n"
	     "  -k    keep the given nodes and delete all others\n";
}

int
main(int argc, char *argv[])
{
	int argi = 1;
	bool keepMode = false;
	while (argv[argi] && argv[argi][0] == '-') {
		switch (argv[argi][1]) {
		case 'k':
			keepMode = true;
			break;
		default:
			usage(argv[0]);
			return -1;
		}
		++argi;
	}

	if (argc - argi < 1) {
		usage(argv[0]);
		return -1;
	}

	WeightedNodeColoredGraph G;
	cout << argv[argi] << ":\n";
	map<int,node> indexNode;
	if (!G.readMCAfile(argv[argi], indexNode)) {
		cerr << "Error reading file (perhaps check number of nodes and number of edges).\n";
		return 1;
	}
	cout
	 << "  -- "
	 << G.numberOfNodes() << " nodes, "
	 << G.numberOfEdges() << " edges, "
	 << G.numberOfColors() << " colors\n";

	++argi;
	if (keepMode) {
		NodeArray<bool> keep(G.theGraph(), false);
		for (; argi < argc; ++argi) {
			int index = atoi(argv[argi]);
			cout << "  -- keeping node " << index << "\n";
			map<int,node>::iterator f = indexNode.find(index);
			if (f == indexNode.end()) {
				cout << "   ** warning: node does not exist!\n";
			} else {
				keep[f->second] = true;
			}
		}
		node w;
		for (node v = G.firstNode()->succ(); v; v = w) { // omit root
			w = v->succ();
			if (!keep[v]) {
				G.delNode(v);
			}
		}
	} else {
		for (; argi < argc; ++argi) {
			int index = atoi(argv[argi]);
			cout << "  -- removing node " << index << "\n";
			map<int,node>::iterator f = indexNode.find(index);
			if (f == indexNode.end()) {
				cout << "   ** warning: node does not exist!\n";
			} else {
				G.delNode(f->second);
				indexNode.erase(f);
			}
		}
	}

	// reduce the input graph
	removeUnreachableNodes(G);
	reduceGraph(G);

	G.writeMCAfile("output.dat");
	return 0;
}
