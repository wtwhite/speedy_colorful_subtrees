#include <ogdf/basic/Graph.h>
#include <ogdf/basic/NodeSet.h>
#include <ogdf/basic/Stopwatch.h>
#include <ogdf/basic/HashArray.h>

#include <sstream>
#include <map>
#include <set>

using namespace ogdf;
using std::map;
using std::set;
using std::pair;

#include "WeightedNodeColoredGraph.h"

static void
findReachableNodes(const WeightedNodeColoredGraph &G, node source, NodeSet &reachable)
{
	NodeArray<bool> seen(G.theGraph(), false);
	SListPure<node> todo;
	todo.pushBack(source);
	seen[source] = true;
	while (!todo.empty()) { // BFS
		node v = todo.popFrontRet();
		reachable.insert(v);
		adjEntry adj;
		forall_adj(adj, v) {
			node w = adj->twinNode();
			if (adj->theEdge()->source() == v
			 && !seen[w]) {
				todo.pushBack(w);
				seen[w] = true;
			}
		}
	}
}

static void
writeInducedSubgraph(WeightedNodeColoredGraph &G, NodeSet &reachable, const string &filename)
{
	reachable.insert(G.firstNode()); // insert root
	for (edge f, e = G.firstEdge(); e; e = f) {
		f = e->succ();
		if (!reachable.isMember(e->source())
		 || !reachable.isMember(e->target())) {
			G.theGraph().hideEdge(e);
		}
	}
	G.writeMCAfile(filename);
	G.theGraph().restoreAllEdges();
}

static int
computeUnionCardinality(const NodeSet &set1, const NodeSet &set2)
{
	OGDF_ASSERT(set1.size() >= set2.size());
	int cnt = set1.size();
	forall_listiterators(node, it, set2.nodes()) {
		if (!set1.isMember(*it)) {
			++cnt;
		}
	}
	return cnt;
}

static void
run(WeightedNodeColoredGraph &G, int limit = 0)
{
	node root = G.firstNode();

	StopwatchWallClock timer;

	cout << "Obtaining sets of reachable nodes for " << root->outdeg() << " sets..." << flush; // O(n^2)
	timer.start();
	Array<NodeSet *> reachable(root->outdeg());
	Array< SListPure<node> > sources(root->outdeg());
	ListPure< Prioritized<int, int> > sortReachable;
	int i = 0;
	for (adjEntry adj = root->firstAdj(); adj; adj = adj->succ()) {
		const node v = adj->twinNode();
		OGDF_ASSERT(root == adj->theEdge()->source());
		OGDF_ASSERT(v->indeg() == 1);
		reachable[i] = new NodeSet(G.theGraph());
		sources[i].pushBack(v);
		findReachableNodes(G, v, *reachable[i]);
		sortReachable.pushBack(Prioritized<int, int>(i, -reachable[i]->size()));
		++i;
	}
	timer.stop();
	cout << " -> " << (timer.milliSeconds() * 1e-3) << " seconds\n";

	cout << "Sorting sets..." << flush; // O(n log n)
	timer.start(true);
	sortReachable.quicksort();
	timer.stop();
	cout << " -> " << (timer.milliSeconds() * 1e-3) << " seconds\n";

	if (!limit) {
		limit = max(750, 50 - sortReachable.front().priority());
	}
	cout << "Partitioning with instance node bound " << limit << "..." << flush;
	timer.start(true);
	for (ListIterator< Prioritized<int, int> > it1 = sortReachable.begin(); it1 != sortReachable.rbegin();) {
		int minUC = G.numberOfNodes() + 1;
		ListIterator< Prioritized<int, int> > minIt = NULL;
		i = (*it1).item();
		NodeSet &set1 = *reachable[i];
		for (ListIterator< Prioritized<int, int> > it2 = it1.succ(); it2.valid(); ++it2) {
			int j = (*it2).item();
			int sd = computeUnionCardinality(set1, *reachable[j]);
			if (sd < minUC) {
				minUC = sd;
				minIt = it2;
			}
		}
		OGDF_ASSERT(minIt.valid());
		if (minUC <= limit) {
			NodeSet *set2 = reachable[(*minIt).item()];
			forall_listiterators(node, it, set2->nodes()) {
				set1.insert(*it);
			}
			sources[i].conc(sources[(*minIt).item()]);
			sortReachable.del(minIt);
		} else {
			//cout << "\t" << "filled instance with " << set1.size() << " nodes (gap: " << limit - set1.size() << ")\n";
			++it1;
		}
	}
	timer.stop();
	cout << " -> " << (timer.milliSeconds() * 1e-3) << " seconds\n";

	cout << "Infering new " << sortReachable.size() << " instances and writing each output file..." << flush;
	timer.start(true);
	i = 0;
	for (ListIterator< Prioritized<int, int> > it = sortReachable.begin(); it != sortReachable.rbegin(); ++it) {
		int j = (*it).item();
		writeInducedSubgraph(G, *reachable[j], "output-" + to_string(++i) + ".dat");
	}
	timer.stop();
	cout << " -> " << (timer.milliSeconds() * 1e-3) << " seconds\n";

	cout << "Cleaning up..." << flush;
	timer.start(true);
	for (i = 0; i < reachable.size(); ++i) {
		delete reachable[i];
	}
	timer.stop();
	cout << " -> " << (timer.milliSeconds() * 1e-3) << " seconds\n";
}

static void
usage(char *prog)
{
	cout
	  << "Usage: "
	  << prog << " <filename> <number of edges>\n"
	     " Partition global file into separate problems.\n"
	     " Write the separate problems to output-{1,2,3,...}.dat\n" // TODO
	     " in the `global' file format, i.e. the uncolored super-root is included\n";
}

int
main(int argc, char *argv[])
{
	int argi = 1;
	while (argv[argi] && argv[argi][0] == '-') {
		switch (argv[argi][1]) {
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
	if (!G.readMCAfile(argv[argi])) {
		cerr << "Error reading file (perhaps check number of nodes and number of edges).\n";
		return 1;
	}
	cout
	 << "  -- "
	 << G.numberOfNodes() << " nodes, "
	 << G.numberOfEdges() << " edges, "
	 << G.numberOfColors() << " colors\n";

	++argi;

	run(G);

	return 0;
}
