// Note: include after using namespace ogdf;
#ifndef _WEIGHTED_NODE_COLORED_GRAPH_H
#define _WEIGHTED_NODE_COLORED_GRAPH_H

// WTJW: All this needed just to handle line 4 of input files properly...
#include <sstream>
#include <iterator>
#include <algorithm>
#include <vector>

class GraphCompositor
{
protected:
	Graph G;

public:
	const Graph &theGraph() const
	{
		return G;
	}

	Graph &theGraph() // for special things like hideEdge/restoreEdge
	{
		return G;
	}

	void delEdge(edge e)
	{
		G.delEdge(e);
	}

	int numberOfNodes() const
	{
		return G.numberOfNodes();
	}

	int numberOfEdges() const
	{
		return G.numberOfEdges();
	}

	node firstNode() const
	{
		return G.firstNode();
	}

	node lastNode() const
	{
		return G.lastNode();
	}

	edge firstEdge() const
	{
		return G.firstEdge();
	}

	edge lastEdge() const
	{
		return G.lastEdge();
	}
};

typedef HashArray< int, set<node> > ColorSetMap;
typedef pair<edge, double> EdgeValPair;

#define SUPERTARGETCOLOR(x)	(- x - 3)

class WeightedNodeColoredGraph : public GraphCompositor
{
protected:
	EdgeArray<double> edgeWeight;
	NodeArray<int> nodeColor;
	ColorSetMap coloredNodes;
	List<node> m_superTargets;
	int numColors;

	//! Set color of a node, initially
	void setColor(node v, int color)
	{
		nodeColor[v] = color;
		if (!coloredNodes.isDefined(color)) {
			// color not yet used for a node
			++numColors;
		}
		coloredNodes[color].insert(v);
	}

public:
	double solution;

	WeightedNodeColoredGraph()
	 : GraphCompositor()
	 , edgeWeight(G)
	 , nodeColor(G)
	 , coloredNodes()
	 , numColors(0)
	 , solution(-0.0)
	{}

	//! Returns the number of (actually used!) colors
	int numberOfColors() const
	{
		return numColors;
	}

	//! Returns the number of nodes with the same color
	int numberOfEquallyColoredNodes(node v) // const
	{
		return coloredNodes[nodeColor[v]].size() - 1;
	}

	//! Update color of a node
	void updateColor(node v, int color)
	{
		// remove old color
		set<node> &oldColorSet = coloredNodes[nodeColor[v]];
		set<node>::iterator itOld = oldColorSet.find(v);
		OGDF_ASSERT(itOld != oldColorSet.end());
		oldColorSet.erase(itOld);
		if (oldColorSet.empty()) {
			--numColors;
		}
		setColor(v, color);
	}

	//! Creates a new node with specified color and returns it.
	node newNode(int color)
	{
		node v = G.newNode();
		setColor(v, color);
		return v;
	}

#if 0
	//! Creates a new node with unique color and returns it.
	node newNode()
	{
		node v = G.newNode();
		// find new color id
		int color = 0;
		for (ColorSetMap::const_iterator icpair = coloredNodes.begin(); icpair != coloredNodes.end(); icpair++) {
			if ((*icpair).first > color) {
				color = (*icpair).first;
			}
		}
		setColor(v, ++color);
		return v;
	}
#endif

	//! Creates a new edge (\a v, \a w) with specified weight and returns it.
	edge newEdge(node v, node w, double weight)
	{
		edge e = G.newEdge(v, w);
		edgeWeight[e] = weight;
		return e;
	}

	//! Removes node \a v and all incident edges from the graph.
	void delNode(node v)
	{
		set<node> &colorSet = coloredNodes[nodeColor[v]];
		colorSet.erase(colorSet.find(v));
		if (colorSet.empty()) {
			--numColors;
		}
		G.delNode(v);
	}

	int color(node v) const
	{
		return nodeColor[v];
	}

	double weight(edge e) const
	{
		return edgeWeight[e];
	}

	double &weight(edge e)
	{
		return edgeWeight[e];
	}

	const ColorSetMap &colorSetMap() const
	{
		return coloredNodes;
	}

	const set<node> &colorSet(int c) const
	{
		OGDF_ASSERT(coloredNodes.isDefined(c));
		return coloredNodes[c];
	}

	void addSuperTargets()
	{
		map<int,node> colorTarget;
		for (node v = G.firstNode(); v; v = v->succ()) {
			const int c = nodeColor[v];
			if (c >= 0) {
				map<int,node>::iterator st = colorTarget.find(c);
				node w;
				if (st == colorTarget.end()) {
					w = newNode(SUPERTARGETCOLOR(c));
					colorTarget[c] = w;
					m_superTargets.pushBack(w);
				} else {
					w = (*st).second;
				}
				newEdge(v, w, 0);
			}
		}
	}

	const List<node> &superTargets() const
	{
		return m_superTargets;
	}

	bool readMCAfile(const string &filename, map<int,node> &indexNode)
	{
		ifstream in(filename.c_str());
		if (!in.is_open()) {
			cerr << "Error: could not open file `" << filename << "'.\n";
			return false;
		}

		// WTJW: Line 4 can take various formats: A single FP number giving the score of the optimal subtree;
		// two FP numbers separated by " + ", giving the score of the optimal subtree and (I think...) the score of the edge from the superroot
		// in the combined dataset to this candidate dataset's root node; or even omitted completely (i.e. we begin listing vertex colours right
		// away), in which case it will consist of two integers separated by whitespace.
		// All this means we can't just read stuff in "nicely", we have to process it line by line.
		int nNodes, nEdges, nColors;
//		in >> nNodes >> nEdges >> nColors;
//
//		in >> solution;
		string line;
		bool alreadyHaveFirstColor = false;
		int firstNodeId, firstColorId;		// Only get values if line 4 is "missing", i.e. alreadyHaveFirstColor == true.
		{
			getline(in, line);
			std::istringstream iss(line);
			iss >> nNodes;
		}
		{
			getline(in, line);
			std::istringstream iss(line);
			iss >> nEdges;
		}
		{
			getline(in, line);
			std::istringstream iss(line);
			iss >> nColors;
		}
		{
			// Now the tricky line 4...
			getline(in, line);
			std::istringstream iss(line);
			std::vector<string> tokens;
			std::copy(std::istream_iterator<string>(iss), std::istream_iterator<string>(), std::back_inserter(tokens));		// Read all tokens as strings into tokens[]
			if (tokens.size() == 2) {
				std::istringstream issTok1(tokens[0]), issTok2(tokens[1]);
				issTok1 >> firstNodeId;
				issTok2 >> firstColorId;
				alreadyHaveFirstColor = true;
			} else {
				std::istringstream issTok(tokens[0]);
				issTok >> solution;
				if (tokens.size() == 3 && tokens[1] == "+") {
					std::istringstream issTok2(tokens[2]);
					double superRootEdgeWeight;
					issTok2 >> superRootEdgeWeight;
#ifdef DEBUG
					cout << "  ** Ignoring specified weight of edge from superroot of " << superRootEdgeWeight << "\n";
#endif
				} else if (tokens.size() != 1) {
					cerr << "Error: Expected line 4 to have 1, 2 or 3 tokens with the middle one being '+'!\n";
					return false;
				}
			}
		}

#ifdef DEBUG
		cout << "  ** ILP solution should be " << solution << "\n";
#endif

		for (int i = 0; i < nNodes; ++i) {
			int nodeId, colorId;
			if (alreadyHaveFirstColor && i == 0) {
				// There was no objective score on line 4, so we already read this data in.
				nodeId = firstNodeId;
				colorId = firstColorId;
			} else {
				if (!(in >> nodeId >> colorId)) {
					return false;
				}
			}
			indexNode[nodeId] = newNode(colorId);
		}
		// assign the right color to the super-root
		updateColor(G.firstNode(), -1);
		// XXX: this is not needed if we use the sane global format

		for (int i = 0; i < nEdges; ++i) {
			int node1, node2;
			double weight;
			if (!(in >> node1 >> node2 >> weight)) {
				return false;
			}
			newEdge(indexNode[node1], indexNode[node2], weight);
		}
		in.close();

		return true;
	}

	bool readMCAfile(const string &filename)
	{
		map<int,node> indexNode;
		return readMCAfile(filename, indexNode);
	}

	bool writeMCAfile(const string &filename) const
	{
		ofstream of(filename.c_str());
		if (!of.is_open()) {
			cerr << "Error: could not open file `" << filename << "'.\n";
			return false;
		}

		of.precision(16);
		of << G.numberOfNodes() << "\n"
		   << G.numberOfEdges() << "\n"
		   << numColors << "\n"
		   << solution << "\n";

		node v;
		int i = 0, c = 0;
		NodeArray<int> nodeIndex(theGraph());
		map<int, int> colorIndex;
		forall_nodes(v, G) {
			nodeIndex[v] = i++;
			map<int, int>::iterator colorIt = colorIndex.find(color(v));
			int col;
			if (colorIt == colorIndex.end()) {
				col = c;
				colorIndex[color(v)] = c++;
			} else {
				col = (*colorIt).second;
			}
			of << nodeIndex[v] << " " << col << "\n";
		}
		edge e;
		forall_edges(e, G) {
			of << nodeIndex[e->source()] << " " << nodeIndex[e->target()] << " " << edgeWeight[e] << "\n";
		}
		of.close();
		return true;
	}

}; // class WeightedNodeColoredGraph
#endif // _WEIGHTED_NODE_COLORED_GRAPH_H
