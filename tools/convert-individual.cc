#include <ogdf/basic/Graph.h>
#include <ogdf/basic/HashArray.h>

#include <sstream>
#include <map>
#include <set>

using namespace ogdf;
using std::map;
using std::set;
using std::pair;

#include "WeightedNodeColoredGraph.h"

static bool
readMCAfile(WeightedNodeColoredGraph &G, const string &filename)
{
	map<int,node> indexNode;
	ifstream in(filename.c_str());
	if (!in.is_open()) {
		cerr << "Error: could not open file `" << filename << "'.\n";
		return false;
	}

	int numNodes, numEdges, numColors;
	in >> numNodes >> numEdges >> numColors;

	in >> G.solution;
	//cout << "  ** ILP solution is " << G.solution << "\n";
	double offset(0);
	std::istream::sentry s(in);
	if (in.peek() == '+') {
		string plus;
		in >> plus >> offset;
		//cout << "  ** ILP offset is " << offset << "\n";
	}

	for (int i = 0; i < numNodes; ++i) {
		int nodeId, colorId;
		if (!(in >> nodeId >> colorId)) {
			return false;
		}
		indexNode[nodeId] = G.newNode(colorId);
		if (i == 0) { // first node is source
			// connect super-root to source
			OGDF_ASSERT(nodeId == 0);
			G.newEdge(G.firstNode(), indexNode[nodeId], offset);
		}
	}

	for (int i = 0; i < numEdges; ++i) {
		int node1, node2;
		double weight;
		if (!(in >> node1 >> node2 >> weight)) {
			return false;
		}
		G.newEdge(indexNode[node1], indexNode[node2], weight);
	}
	in.close();

	return true;
}

static void
usage(char *prog)
{
	cout
	  << "Usage: "
	  << prog << " <infile> <outfile>\n"
	     " Converts <infile> from any old individual format to global format\n"
	     " and write result to <outfile>\n";
}

int
main(int argc, char *argv[])
{
	if (argc != 3) {
		usage(argv[0]);
		return -1;
	}

	WeightedNodeColoredGraph G;
	G.newNode(-1); // super-root (by G.firstNode())
	if (!readMCAfile(G, argv[1])) {
		cerr << "Error reading file (perhaps check number of nodes and number of edges).\n";
		return 1;
	}

	G.writeMCAfile(argv[2]);
	return 0;
}
