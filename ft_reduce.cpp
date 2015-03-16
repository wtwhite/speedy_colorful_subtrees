#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
#include <cfloat>
#include <map>
#include <set>
#include <numeric>
#include <fstream>
#include <sstream>
#include <functional>
#include <cmath>		// Needed for FP abs()!
#include <climits>		// Needed for INT_MAX
#include <ctime>
#include <iterator>

// This turns out to not be especially helpful, and surprisingly it makes things quite a bit slower.
//#define USE_OC 1		// Add an "opportunity cost" field to each edge

using namespace std;

struct edge {
	int u, v;
	double w;
#ifdef USE_OC
	double oc;			// A LB on the opportunity cost of this edge -- i.e. if this edge is present, then we miss out on another set of edges weighing at least oc which we could otherwise add.
#endif	// USE_OC
};

// Used for unique().
struct is_equal_on_both_endpoints {
	bool operator()(edge a, edge b) const {
		return a.u == b.u && a.v == b.v;
	}
};

enum { DONTCARE, ASCBYSTART, DESCBYSTART, ASCBYWEIGHT, DESCBYWEIGHT, ASCASC, ASCDESC, DESCASC, DESCDESC } compareMode = DONTCARE;		//DEBUG

char *inputFName = NULL;		// NULL means read from stdin.
char *outputFName = NULL;		// NULL means write to stdout.

bool shouldCheckPreconds = false;	// If true, various actions will perform tests of various preconditions (e.g. to check that vertices are numbered in topological order)
bool shouldSaveTables = false;		// If true, a variety of functions will write out their results to particular files.
vector<edge> watchedEdges;			// A list of edges that we will check for existence after every action.  Used for debugging, for finding the reduction that is deleting an edge when it shouldn't be.

struct action_info {
	int nIterations;
	clock_t elapsed;
	int nEdgesDeleted;
};

bool timeActions = false;
//map<void (*)(), pair<int, clock_t> > actionTimes;		// Key is pointer to action function (HACK: yes, this is totally weird, it should be the name but that duplicates info in the "AST"); first elem of pair is execution count, second is total elapsed time.
map<void (*)(), action_info> actionTimes;		// Key is pointer to action function (HACK: yes, this is totally weird, it should be the name but that duplicates info in the "AST"); first elem of pair is execution count, second is total elapsed time.

//HACK: We "should" use numeric_limits<double>::min() instead of DBL_MIN etc., but it's not a constant, so it forces us to put it in
// a variable or face a function call each time...  yech.
// Another reason why we don't define these to use the built-in FP "infinity" values is because equality comparisons don't necessarily work
// with them.  It's useful to be able to use one of these values as a "seen yet?" indicator and to test for this with equality.
//#define INFINITY 1e12
#define INF DBL_MAX
// Yup, INT_MIN means the smallest (i.e. most negative integer), but of course DBL_MIN means the smallest *positive* double!
#define NEGINF (-DBL_MAX)

map<string, int> dumpTableLastSuffixFor;		// Used by dumpTable()

// The following guff is to enable functions to be conveniently loaded into a global map and then called by name.
map<string, void (*)()> actions;

clock_t lastClock = clock();		// For timing actions inside run()

struct action_inserter {
	action_inserter(void (*func)(), string name) {
		actions[name] = func;
	}
};

//HACK: we use __LINE__ to get "unique" names for these objects, because otherwise you can't have more than one name for an action.
// This also means that you should not put more than one ADD_ACTION() on a given source line!
#define ADD_ACTION(f, name) \
	void f(); /* Forward declaration may be necessary */ \
	namespace ACTION_NAMESPACE { \
		action_inserter f##_##__LINE__(f, name); \
	}

// Dumps a vector of any type to a uniquely-named file having some given prefix.  Each element will be preceded by its 0-based index and a tab, and terminated with a newline.
template <typename T>
void dumpTable(string fnPrefix, vector<T> const&v) {
	if (shouldSaveTables) {
		ostringstream oss;
		oss << fnPrefix << ++dumpTableLastSuffixFor[fnPrefix];
		string fName(oss.str());
		cerr << "Writing file '" << fName << "'...\n";
		
		ofstream f(fName.c_str());
		for (int i = 0; i < v.size(); ++i) {
			f << i << '\t' << v[i] << '\n';
		}
	}
}

ostream& operator<<(ostream& os, edge e) {
	return os << '(' << e.u << ", " << e.v << ": w=" << e.w << ')';
}

// Used mainly for dumping 2D arrays with dumpTable().
template <typename T>
ostream& operator<<(ostream& os, vector<T> const& v) {
	for (int i = 0; i < v.size(); ++i) {
		if (i > 0) {
			os << ' ';
		}
		
		os << v[i];
	}
	
	return os;
}

struct compare_by_endpoint_incr {
	bool operator()(edge a, edge b) const {
		//DEBUG: The line below is all I believe we need, but to do some basic testing let's change the order around a bit...
		// DAMMIT: On graphs\pos44247C16H41N11O11P.txt, sorting ascending by start point leaves 7536 edges remaining, as does sorting descending by start point... but the sets of edges
		// are different -- even after sorting the output files!  E.g. the edge "0 289 -12.8022" appears only in the desc file, and no other edge with that weight appears there :-(
		// Also: some edges (like "300 371 -1.1674") appear twice in the desc file!
		//return a.v < b.v;
		//return a.v < b.v || (a.v == b.v && a.u < b.u);
		//return a.v < b.v || (a.v == b.v && a.u > b.u);
		switch (compareMode) {
		case DONTCARE:
			return a.v < b.v;
			
		case ASCBYSTART:
			return a.v < b.v || (a.v == b.v && a.u < b.u);
		
		case DESCBYSTART:
			return a.v < b.v || (a.v == b.v && a.u > b.u);
		
		case ASCBYWEIGHT:
			return a.v < b.v || (a.v == b.v && a.w < b.w);
		
		case DESCBYWEIGHT:
			return a.v < b.v || (a.v == b.v && a.w > b.w);
		
		case ASCASC:
			//if (a.v < b.v) return true;
			//if (a.u < b.u) return true;
			//return a.w < b.w;
			if (a.v != b.v) return a.v < b.v;
			if (a.u != b.u) return a.u < b.u;
			return a.w < b.w;
		
		case ASCDESC:
			//if (a.v < b.v) return true;
			//if (a.u < b.u) return true;
			//return a.w > b.w;
			if (a.v != b.v) return a.v < b.v;
			if (a.u != b.u) return a.u < b.u;
			return a.w > b.w;
		
		case DESCASC:
		//	if (a.v < b.v) return true;
		//	if (a.u > b.u) return true;
		//	return a.w < b.w;
			if (a.v != b.v) return a.v < b.v;
			if (a.u != b.u) return a.u > b.u;
			return a.w < b.w;
		
		case DESCDESC:
			//if (a.v < b.v) return true;
			//if (a.u > b.u) return true;
			//return a.w > b.w;
			if (a.v != b.v) return a.v < b.v;
			if (a.u != b.u) return a.u > b.u;
			return a.w > b.w;
		}
	}
};

int nVerts, nEdges, nCol;
vector<edge> edges;			// Keeps edges, usually in sorted order by endpoint, so that all edges to a particular vertex appear contiguously.
vector<int> startOfEdgesTo;	// startOfEdgesTo[i] is the position in edges[] where all edges ending at vertex i start.
vector<int> c;		// c[i] is the colour of vertex i.  0 <= c[i] < nCol for all i.

// Used for renaming vertices during topological sort
int nextVertexId = 0;
vector<int> newIdFor;		// newIdFor[i] is the new name for vertex i after topological sorting.  -1 if not yet visited.

// Use this for renaming vertices after topological sorting
//// Elements of this struct will be used in an array indexed by the endpoint of
//// the directed edge; v will contain the edge's start point, and w its weight.
//struct halfEdge {
//	int v;
//	double w;
//};
//
//vector<halfEdge> edgesEndingAt;		// edgesEndingAt[i] is a list of directed edges ending at vertex i.

vector<double> ub;			// ub[i] is an upper bound on the weight of the maximum subtree rooted at vertex i.  Always >= 0.
vector<edge> edgesBySource;	// Analogous to edges[], but contains edges in order of their start points instead of their end points.
vector<int> startOfEdgesFrom;		// startOfEdgesFrom[i] is the position within edgesBySource[] where vertex i's outgoing edges begin.

//// lb[u][v] is a lower bound on the score change obtainable by attaching v to the solution tree, given that
//// the "anchor" u is already in the solution tree.  We might attach v to u itself, or to some ancestor of u.
//// Only defined for u <= v.
// lb[v][u] is a lower bound on the score change obtainable by attaching v to the solution tree, given that
// the "anchor" u is already in the solution tree.  We might attach v to u itself, or to some ancestor of u.
// Only defined for u <= v.
vector<vector<double> > lb;

// Similar to lb[][], but with subscript order reversed.  See description under calc-anc-col-lbs action.
vector<vector<double> > lbCol;

edge reverse(edge uv) {
	edge vu;
	vu.u = uv.v;
	vu.v = uv.u;
	vu.w = uv.w;
	return vu;
}

struct compare_by_startpoint_incr {
	//bool operator()(edge a, edge b) const {
	//	return a.u < b.u;
	//}
	//HACK: This way, we don't have to reimplement all that ghastly "switch (compareMode)" stuff...
	bool operator()(edge a, edge b) const {
		return compare_by_endpoint_incr()(reverse(a), reverse(b));
	}
};

// The variables below are used for calculating my vertex UB.
vector<bool> seenUbChild;		// seenUbChild[i] is true iff vertex i has been visited yet during computation of my vertex upper bound.
int nTotalInternalVerticesWithZeroChildUpperBound;
int nTotalEdgesDeletedDueToZeroChildUpperBound;
double highestChildUpperBoundEver;

// The variables below are used for calculating Sebastian's ColourFul Forest vertex UB.
vector<bool> seenUbColForest;		// seenUbColForest[i] is true iff vertex i has been visited yet during computation of the Colourful Forest vertex upper bound.
//vector<vector<int> > reachable;		// reachable[i] is a sorted list of all vertices reachable from vertex i.
int nTotalInternalVerticesWithZeroColForestUpperBound;
int nTotalEdgesDeletedDueToZeroColForestUpperBound;
double highestColForestUpperBoundEver;
int nOrphanColours;		// Number of colours that have no in-edges but at least one out-edge

bool shouldStrengthenColForestUpperBound = false;		// Can be enabled with "enable-col-forest-vub-strength".
double totalStrengtheningColForestUpperBound;
double totalSneakyStrengtheningColForestUpperBound;

vector<vector<int> > impliedAnc;		// impliedAnc[i] is a sorted list of ancestors of i that must be present if i is present.  Includes i itself, and also always includes the root.

//vector<vector<int> > edgesImpliedByVertex;		// edgesImpliedByVertex[i] is an unsorted list of indices into edgesBySource[] of the edges implied by vertex i.
vector<vector<edge> > edgesImpliedByVertex;		// edgesImpliedByVertex[i] is an unsorted list of edges implied by vertex i (through it's max-weight in-edge).

// Here we rely on the fact that all vertices of the same colour form a contiguous block.
vector<int> firstVertOfSameColourAs;
vector<int> firstVertOfNextColour;

// Sort edge list in edges[] by destination vertex, so that all in-edges for a given vertex appear contiguously.
// Then update the array of indices that allow each destination vertex's in-edge list to be quickly found.
void sortEdgesAndCalcVertIndices(bool sortNeeded) {
	if (sortNeeded) {
		// Important!  nEdges is what we trust.  This was different than edges.size() for a while there...
		cerr << "Sorting "<< nEdges << " edges by destination vertex...\n";
		assert(nEdges == edges.size());
		sort(edges.begin(), edges.begin() + nEdges, compare_by_endpoint_incr());
	}
	
	startOfEdgesTo[0] = 0;
	int j = 0;
	for (int i = 1; i <= nVerts; ++i) {		// Note we do fill out startOfEdgesTo[nVerts], so this array must be of size nVerts + 1.
		while (j < nEdges && edges[j].v < i) {
			++j;
		}
		
		startOfEdgesTo[i] = j;
	}
}

//TODO: Move this.
// We have to map colours to colours -- I can't think of an O(n) way to map vertices to colours!
//vector<int> newColFor;			// newColFor[i] is the colour that old *vertex* (not *colour*) i will be mapped to.
vector<int> newColFor;			// newColFor[i] is the colour that old *colour* (not *vertex*) i will be mapped to.

void checkVerticesAreTopSorted(string errMsg) {
	for (int i = 0; i < nEdges; ++i) {
		if (edges[i].u >= edges[i].v) {
			cerr << "Edge (" << edges[i].u << ", " << edges[i].v << ") exists from a higher-numbered vertex to a lower-numbered vertex!\n";
			cerr << errMsg << "\n";
			exit(1);
		}
		
		// Also check colours.
		if (c[edges[i].u] >= c[edges[i].v]) {
			cerr << "Edge (" << edges[i].u << ", " << edges[i].v << ") exists from a higher-numbered colour (" << c[edges[i].u] << ") to a lower-numbered colour (" << c[edges[i].v] << ")!\n";
			cerr << errMsg << "\n";
			exit(1);
		}
	}
}

// The new topSort() respects colours as well: every vertex of a parent's colour will be recursively renamed first.
// Having thought about it a lot, this is necessary!

//// Find new IDs for v and all its parents, so that they will obey (x, y) in graph => x < y.
//void topSort(int v) {
//	if (newIdFor[v] == -1) {
//		// Haven't visited this vertex yet.
//		// Recurse to rename all parents first.
//		for (int i = startOfEdgesTo[v]; i < startOfEdgesTo[v + 1]; ++i) {
//			topSort(edges[i].u);
//		}
//		
//		newIdFor[v] = nextVertexId++;
//		
//		// We figure out new colours by recording the lowest-numbered new vertex ID for each colour.
//		// We'll update this (by squeezing it down to the range 0 .. nCol-1) in a second pass.
//		//assert(c[v] != 3);		//DEBUG
//		//DEBUG
//		//if (c[v] == 90) {
//		if (c[v] == 50) {
//			cerr << "c[v]=" << c[v] << ": Previous newColFor=" << newColFor[c[v]] << ", new contender = " << newIdFor[v] << "\n";
//		}
//		newColFor[c[v]] = min(newColFor[c[v]], newIdFor[v]);
//	}
//}
//BUG: Actually this approach is busted because it can still assign -- see i_hate_colour_renaming.txt for a mathematical counterexample!
// Find new IDs for v and all its parents, so that they will obey (x, y) in graph => x < y.
//vector<vector<int> > topSortVertsOfColour;		//HACK: This is used only by topSort(), though it could in fact be used elsewhere and this might save time.
//void topSort(int v) {
//	if (newIdFor[v] == -1) {
//		// Haven't visited this vertex yet.
//		// Find the colours of all its parents, and recurse to rename all vertices of any of these colours first.
//		vector<bool> coloursToProcessFirst(nCol, false);
//		for (int i = startOfEdgesTo[v]; i < startOfEdgesTo[v + 1]; ++i) {
//			coloursToProcessFirst[c[edges[i].u]] = true;
//		}
//		
//		for (int i = 0; i < nCol; ++i) {
//			if (coloursToProcessFirst[i]) {
//				for (int j = 0; j < topSortVertsOfColour[i].size(); ++j) {
//					topSort(topSortVertsOfColour[i][j]);
//				}
//			}
//		}
//		newIdFor[v] = nextVertexId++;
//		
//		// We figure out new colours by recording the lowest-numbered new vertex ID for each colour.
//		// We'll update this (by squeezing it down to the range 0 .. nCol-1) in a second pass.
//		//assert(c[v] != 3);		//DEBUG
//		//DEBUG
//		//if (c[v] == 90) {
//		if (c[v] == 50) {
//			cerr << "c[v]=" << c[v] << ": Previous newColFor=" << newColFor[c[v]] << ", new contender = " << newIdFor[v] << "\n";
//		}
//		newColFor[c[v]] = min(newColFor[c[v]], newIdFor[v]);
//	}
//}

// The (only?) right way is to actually forget about topologically sorting vertices, and topologically sort colours instead.
// To store these colour-edges, we *could* use a flat array plus indexes into it a la startOfEdgesTo[], but it's easier and still
// fast enough to just use an array of set<>s.
// Much simpler than the previous approach!  :)
vector<set<int> > topSortColourEdgesToColour;		// topSortColourEdgesToColour[i][j] is the jth colour (in no particular order) that has an edge to colour i.
vector<vector<int> > topSortVertsOfColour;		//HACK: This is used only by topSort(), though it could in fact be used elsewhere and this might save time.
int nextColourId;									// Used for renumbering colours

void topSort(int col) {
	//cerr << "topSort(col=" << col << ") called.\n";		//DEBUG
	if (newColFor[col] == -1) {
		// Haven't visited this colour yet.
		// Find the colours of all its parents (a colour x is a parent of a colour y if there is an edge from any x-coloured vertex to any y-coloured vertex)
		// and recurse to process them first.
		for (set<int>::const_iterator i(topSortColourEdgesToColour[col].begin()); i != topSortColourEdgesToColour[col].end(); ++i) {
			topSort(*i);
		}
		newColFor[col] = nextColourId++;
		//cerr << "Mapping old colour " << col << " to new colour " << newColFor[col] << ".\n";		//DEBUG
		
		// Now renumber all vertices of this colour.  It's safe to just assign consecutive integers,
		// since we know that no two vertices of the same colour can have any edges between each other.
		for (int i = 0; i < topSortVertsOfColour[col].size(); ++i) {
			newIdFor[topSortVertsOfColour[col][i]] = nextVertexId++;
		}
		
		//DEBUG
		//cerr << "Renumbered " << topSortVertsOfColour[col].size() << " vertices of colour " << col << " and renumbered their colour to " << newColFor[col] << ".\n";
	}
}


// This action has the rare distinction of being runnable *before* vertex renaming.
// We figure out which colours are reachable from each given colour, and output a graph containing this, with colours in the input
// graph encoded as vertices, with an edge (u, v) if colour v is reachable from colour u in the original graph.
// The graph output has dummy colours.
// IMPORTANT: I mean reachable in the sense that you can get from a vertex of colour u to a vertex of colour v by following "colour-edges",
// i.e. there exists a path
// Actually let's not try to calculate reachability here.  Let's just write out a graph in which we show which colours you can get
// to *directly* from each colour.  This is much simpler, and the cycle-testing can be done by the second invocation anyway.
ADD_ACTION(writeColourGraphDebug, "debug-write-colour-graph")
void writeColourGraphDebug() {
	char const *fn = "colourgraph.txt";
	cerr << "Writing colour graph to file " << fn << "...\n";
	ofstream f(fn);
	
	vector<vector<bool> > seen(nCol, vector<bool>(nCol, false));		// Just do it the dumbest possible way
	vector<edge> colEdges;
	for (int i = 0; i < nEdges; ++i) {
		int cU = c[edges[i].u];
		int cV = c[edges[i].v];
		if (!seen[cU][cV]) {
			if (seen[cV][cU]) {		//DEBUG: we can report this special case nice and early
				cerr << "There's an edge directly from colour " << cU << " to colour " << cV << " and also one in the opposite direction!\n";
			}
			
			//edge ce(cU, cV, 1.0);		// Dummy edge weight
			edge ce = { cU, cV, 1.0 };		// Dummy edge weight
			colEdges.push_back(ce);
			seen[cU][cV] = true;
		}
	}
	
	// Write out header and dummy colours
	f << nCol << "\n" << colEdges.size() << "\n" << nCol << "\n";
	for (int i = 0; i < nCol; ++i) {
		f << i << ' ' << i << "\n";
	}
	
	// Write out unique edges between colours
	for (int i = 0; i < colEdges.size(); ++i) {
		f << colEdges[i].u << ' ' << colEdges[i].v << " 1.0\n";
	}
	
	cerr << "File written.\n";
}

double knownOptSolVal = -999;		//HACK: Only meaningfull if hasKnownOptSolVal == true, which is set inside readInput().
bool hasKnownOptSolVal;

// Needed by slideLb() and reduce-dompath2.
void calcFirstVertOfSameColour() {
	firstVertOfSameColourAs.resize(nVerts);
	firstVertOfNextColour.resize(nVerts);
	
	int firstVert = 0;
	//int nextVert = 0;
	for (int i = 1; i < nVerts; ++i) {
		if (c[i - 1] != c[i]) {
			assert(c[i] == c[i - 1] + 1);		//DEBUG: We depend on this I think...
			for (int j = firstVert; j < i; ++j) {
				firstVertOfNextColour[j] = i;
			}
			firstVert = i;
		}
		
		firstVertOfSameColourAs[i] = firstVert;
	}
	
	for (int j = firstVert; j < nVerts; ++j) {
		firstVertOfNextColour[j] = nVerts;
	}
}

void readInput(istream& is) {
	// Read in number of vertices, edges and colours
	is >> nVerts >> nEdges >> nCol;
	
	c.resize(nVerts);
	
	// A known optimal solution value may or may not appear next.  We can tell
	// if it is there because it will appear by itself on the next line; if
	// it is absent, the next line will contain a pair of numbers instead (the
	// first line of data for assignment of colours to vertices).
	string line;
	is >> ws;		// Need to skip end of previous line!
	getline(is, line);
	cerr << "line=<" << line << ">.\n";		//DEBUG
	istringstream iss(line);
	iss >> knownOptSolVal;
	int dummy;
	iss >> dummy;
	if (!iss) {
		hasKnownOptSolVal = true;
		cerr << "Input file records known optimal solution value of " << knownOptSolVal << ".\n";
	} else {
		hasKnownOptSolVal = false;
		cerr << "Input file does not record any known optimal solution value.\n";
		cerr << "knownOptSolVal=" << knownOptSolVal << ", dummy=" << dummy << ".\n";		//DEBUG
		c[static_cast<int>(knownOptSolVal)] = dummy;		// Reinterpret the line as the first colour assignment
	}
	
	// Read in colours
	//for (int i = 0; i < nVerts; ++i) {
	for (int i = !hasKnownOptSolVal; i < nVerts; ++i) {		// Read one fewer item if we already handled the first colour assignment
		int j;
		is >> j;
		assert(j == i);		// Sanity check.  Maybe we should allow vertex colours to be listed in any order, but there's no reason to give them in any other order, and this makes it easy to check that we have them all.
		is >> c[j];
	}
	
	// Read in weighted edges
	edges.resize(nEdges);
	for (int i = 0; i < nEdges; ++i) {
		edge e;
		is >> e.u >> e.v >> e.w;
		edges[i] = e;
		//int u;			// Edge endpoint
		//halfEdge he;	// Edge start point and weight
		//is >> he.v >> u >> he.w;
		
	}
	
	startOfEdgesTo.resize(nVerts + 1);		// Need an extra item at the end to get the number of edges to the last vertex
	sortEdgesAndCalcVertIndices(true);
}

//HACK: This function's name is now a misnomer, since we only read from stdin if -i was not specified...
ADD_ACTION(readInputFromStdin, "read")
void readInputFromStdin() {
	if (inputFName) {
		cerr << "Reading input from file \"" << inputFName << "\"...\n";
		ifstream in(inputFName);
		readInput(in);
	} else {
		cerr << "Reading input from stdin...\n";
		readInput(cin);
	}
	
	cerr << "Read in graph with " << nVerts << " vertices, " << nEdges << " edges and " << nCol << " colours.\n";
}

// Now both renames vertices using newIdFor[] AND (if calcInvert is true) turns newIdFor[] into its inverse (for a later call to revertVertices()).
void renameVertices(bool calcInvert) {
	// Actually rename vertices.
	for (int i = 0; i < nEdges; ++i) {
		edges[i].u = newIdFor[edges[i].u];
		edges[i].v = newIdFor[edges[i].v];
	}
	
	// Now fix up the edge list again.
	sortEdgesAndCalcVertIndices(true);
	
	// Similarly rename any watched edges.
	for (int i = 0; i < watchedEdges.size(); ++i) {
		watchedEdges[i].u = newIdFor[watchedEdges[i].u];
		watchedEdges[i].v = newIdFor[watchedEdges[i].v];
	}
	
	//// Fix up vertex colours.  This could be done in-place instead, but it's hard work (cycle chasing), will be slower and the extra space we use isn't important anyway.
	//vector<int> newCol(nVerts);
	//for (int i = 0; i < nVerts; ++i) {
	//	newCol[newIdFor[i]] = c[i];
	//}
	//c = newCol;
	// Rename vertex colours.
	vector<int> tmpCol(nVerts, -1);
	for (int i = 0; i < nVerts; ++i) {
		//c[i] = newColFor[c[i]];
		//c[newIdFor[i]] = newColFor[c[i]];
		tmpCol[newIdFor[i]] = newColFor[c[i]];
	}
	c = tmpCol;
	
	dumpTable("renameVertices", newIdFor);
	dumpTable("renameVertices_colours", newColFor);
	dumpTable("renameVertices_colours_of_vertices", c);
	
	if (calcInvert) {
		if (shouldCheckPreconds) {
			checkVerticesAreTopSorted("Somehow the vertices are not top-sorted, even though renameVertices(calcInvert=true) was called!");		//HACK: We only check this if we are forward-converting, which usually corresponds to calcInvert == true.
		}
		
		// Might as well compute the back-renaming function now instead of at the call to revertVertices() --
		// this will be useful for debugging.
		vector<int> tmp(newIdFor);
		for (int i = 0; i < nVerts; ++i) {
			newIdFor[tmp[i]] = i;
		}
		
		vector<int> tmpCol(newColFor);
		for (int i = 0; i < nCol; ++i) {
			newColFor[tmpCol[i]] = i;
		}
		
		calcFirstVertOfSameColour();			// Of general utility, so do it now
	}
}

//HACK: For convenience with debugging etc.
edge unrenumber(edge uv) {
	uv.u = newIdFor[uv.u];
	uv.v = newIdFor[uv.v];
	return uv;
}

//HACK: Currently this is safe to call only once, due to dependence on global variables etc.
ADD_ACTION(calcNewVertexIds, "renumber-verts")
void calcNewVertexIds() {
	cerr << "Renumbering vertices in topological order, so that an edge (u, v) implies v > u...\n";
	
	//newColFor.assign(nCol, nVerts);		// A too-large value that will be reduced by later min() calls
	newColFor.assign(nCol, -1);		// -1 indicates "not processed yet"
	
	// First find all the "colour-edges"
	topSortColourEdgesToColour.assign(nCol, set<int>());
	for (int i = 0; i < nEdges; ++i) {
		// Sanity checking.  A fuller test would be for absence of cycles, but this was enough to find a colour assignment problem with Kai's combined instances.
		if (c[edges[i].u] == c[edges[i].v]) {
			cerr << "Both endpoints of edge " << edges[i] << " have colour " << c[edges[i].u] << "!\n";
			exit(1);
		}
		topSortColourEdgesToColour[c[edges[i].v]].insert(c[edges[i].u]);
	}
	
	nextColourId = 0;
	nextVertexId = 0;
	//newIdFor.resize(nVerts, -1);
	newIdFor.assign(nVerts, -1);		// Fixes a potential bug: resize() won't reset elements to -1 if calcNewVertexIds() is called more than once.
	topSortVertsOfColour.assign(nCol, vector<int>());
	for (int i = 0; i < nVerts; ++i) {
		topSortVertsOfColour[c[i]].push_back(i);
	}
	
	//// Perform a topological sort to figure out new IDs for vertices so that for every edge (u, v), u < v.
	//for (int i = 0; i < nVerts; ++i) {
	//	cerr << "About to call topSort(" << i << ")...\n";		//DEBUG
	//	topSort(i);
	//}
	// Perform a topological sort to figure out new IDs for colours and vertices so that for every edge (u, v), u < v, and the same for every "colour-edge".
	for (int i = 0; i < nCol; ++i) {
		cerr << "About to call topSort(" << i << ")...\n";		//DEBUG
		topSort(i);
	}
	
	assert(nextColourId == nCol);
	assert(nextVertexId == nVerts);
	
	//// Now "compress" the new colours back down to the range 0 .. nCol-1.
	//
	//dumpTable("calcNewVertexIds_colours_before_compression", newColFor);		//DEBUG
	//dumpTable("calcNewVertexIds_colours_for_vertices_before_compression", newColFor);		//DEBUG
	//
	//// First build the compression mapping
	//vector<int> sortedColours(newColFor);
	//sort(sortedColours.begin(), sortedColours.end());
	//vector<int> colourRank(nVerts, -1);
	//for (int i = 0; i < nCol; ++i) {
	//	colourRank[sortedColours[i]] = i;
	//}
	//
	//// Now map the colours
	//for (int i = 0; i < nCol; ++i) {
	//	newColFor[i] = colourRank[newColFor[i]];
	//	assert(newColFor[i] != -1 && newColFor[i] != nVerts);		// Sanity
	//}
	
	renameVertices(true);
}

void writeOutput(ostream& os) {
	os << nVerts << '\n' << nEdges << '\n' << nCol << '\n';
	
	if (hasKnownOptSolVal) {
		int oldPrec = os.precision();		// Ugh, I hate iostreams
		os.precision(30);		//HACK: guesstimated from looking at Kai's files...
		
		os << knownOptSolVal << '\n';		// Only write this out if it was in the input
		os.precision(oldPrec);
	}
	
	// Write out colours
	for (int i = 0; i < nVerts; ++i) {
		os << i << ' ' << c[i] << '\n';
	}
	
	int oldPrec = os.precision();		// Ugh, I hate iostreams
	os.precision(30);		//HACK: guesstimated from looking at Kai's files...
	
	// Write out weighted edges
	for (int i = 0; i < nEdges; ++i) {
		os << edges[i].u << ' ' << edges[i].v << ' ' << edges[i].w << '\n';
	}
	
	os.precision(oldPrec);
}

//HACK: This function's name is now a misnomer, since we only write to stdout if -o was not specified...
ADD_ACTION(writeOutputToStdout, "write")
void writeOutputToStdout() {
	if (outputFName) {
		cerr << "Writing output to file \"" << outputFName << "\"...\n";
		ofstream out(outputFName);
		writeOutput(out);
	} else {
		cerr << "Writing output to stdout...\n";
		writeOutput(cout);
	}
	
	cerr << "Wrote out graph with " << nVerts << " vertices, " << nEdges << " edges and " << nCol << " colours.\n";
}

vector<vector<double> > h;		// h[i][j] stores the maximum, over all paths from i to j EXCEPT the single-edge path (i, j), of the minimum-weight subpath of length at least 1.
vector<vector<double> > e;		// e[i][j] stores the maximum, over all paths from i to j EXCEPT the single-edge path (i, j), of the minimum-weight subpath ending at j of length 0 or more.
vector<vector<double> > b;		// b[i][j] stores the maximum, over all paths from i to j, of the minimum-weight subpath of length at least 1.

double bestEverH = NEGINF;	//DEBUG

// The edge reduction made possible by these DP calculations is actually invalid for the Maximum **Colourful** Subtree problem,
// because it doesn't consider the possibility that the optimal solution contains a vertex of the same colour as the vertex x
// that we (may) introduce.  It's still valid for the uncoloured variant of the problem, where we are always allowed to introduce
// a vertex that is not already in the solution if added edge will not create an (unordered) cycle.
void calcDpTablesForUncolouredReduction() {
	if (h.empty()) {
		cerr << "Calculating DP tables for uncoloured reduction...\n";
		
		if (shouldCheckPreconds) {
			checkVerticesAreTopSorted("Vertices must be numbered in topological order to calculate the DP tables for uncoloured reduction.");
		}
		
		//HACK: The memory allocations below are horribly inefficient, at least in C++98.  Might be better just to new up a plain old 2D array.
		h.resize(nVerts, vector<double>(nVerts));
		e.resize(nVerts, vector<double>(nVerts));
		b.resize(nVerts, vector<double>(nVerts));
		
		// Putting j in the outside loop should improve memory locality of accesses to edges[].
		// Actually it might be necessary for correctness too!  (Since when calculating h[i][j],
		// we might need h[a][b] for *any* a, but only for any b < j -- meaning we must have
		// already calculated this h[a][b] value on an earlier pass of the outer loop.)
		for (int j = 0; j < nVerts; ++j) {
			for (int i = 0; i < j; ++i) {
				double bestH = NEGINF;
				double bestE = 0;
				double wIJ = NEGINF;
				for (int ki = startOfEdgesTo[j]; ki < startOfEdgesTo[j + 1]; ++ki) {
					int k = edges[ki].u;
					double w = edges[ki].w;
					
					if (k == i) {
						wIJ = w;
					} else if (k > i) {
						double hx = min(b[i][k], e[i][k] + w);
						bestH = max(bestH, hx);
						//DEBUG
						if (bestH > bestEverH) {
							bestEverH = bestH;
							cerr << "Found new best-ever h-value for h[" << i << "][" << j << "] = " << bestH << "!\n";
						}
						
						double ex = min(0.0, e[i][k] + w);
						bestE = max(bestE, ex);
					}
				}
				
				h[i][j] = bestH;
				e[i][j] = bestE;
				b[i][j] = max(bestH, wIJ);
			}
		}
	}
}

ADD_ACTION(clearDpTablesForUncolouredReduction, "clear-uncoloured-dp-tables")
void clearDpTablesForUncolouredReduction() {
	if (!h.empty()) {
		cerr << "Clearing DP tables for uncoloured reduction...\n";
		h.clear();
		e.clear();
		b.clear();
	}
}

// See the comments above calcDpTablesForUncolouredReduction() for why this is unfortunately invalid on the coloured
// variant of the problem.
ADD_ACTION(uncolouredReduceEdges, "reduce-uncoloured")
void uncolouredReduceEdges() {
	cerr << "Reducing edges using the uncoloured reduction...\n";
	
	calcDpTablesForUncolouredReduction();
	
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		edges[j] = edges[i];
		
		if (h[edges[i].u][edges[i].v] <= edges[i].w) {		// Could actually use '<' without producing a worse solution, but '<=' should guarantee that the same solution (not just a solution of the same score) can be returned.
			++j;		// This edge lives...
		}
	}
	
	nEdges = j;
	edges.erase(edges.begin() + nEdges, edges.end());		// Keep the container size in sync with nEdges
	
	// Now fix up the edge list again.  But it's already in sorted order, so don't bother resorting it.
	sortEdgesAndCalcVertIndices(false);
	edgesBySource.clear();		// So that anyone who needs these edges in future knows to recalculate them
}

//DEBUG
// Just testing edge deletion using something that will be reliable across vertex renamings: is the last digit before the decimal point of the weight a 7?  If so, delete!
// Even THIS produces different results for different asc/desc sorts of pos44247C16H41N11O11P!  E.g. "0 614 -10.6208" appears in desc but not asc, while "0 616 -10.6208" appears in desc but not asc!
// There also appear to be other edge *weights* that appear in just one or the other...
// In the original file: "0 614 -13.009059286941962", "0 616 -8.335963514068549"
//
// OK: They are the same after sorting when using --reduce-mode silly --skip-rename --skip-dp.  Phew -- a starting point.
// YES! They are DIFFERENT when using just --reduce-mode silly --skip-dp.  So it must be the renaming step that is causing problems.  Also, importantly, the files are identical to the files
// produced when --skip-dp is turned off.  :-)
//
// OK, have now fixed a BUG in sillyReduceEdges() :-P  We actually do get rid of edges whose last integer digit is a 7 now.  This leaves 50151 of the 56941 edges in both the following cases:
// --compare-mode ascbystart --reduce-mode silly --skip-rename --skip-dp
// --compare-mode descbystart --reduce-mode silly --skip-rename --skip-dp
// Also a test with grep "7\." finds no lines on either one.  Also the output files are identical after sorting :-)
//
// When we drop "--skip-rename", suddenly (a) the files differ and (b) 453 lines containing "7." appear in asc, 454 in desc!
//
// The puzzling part is that getting rid of all reductions causes the outputs to be the same again...  the problem is with sorting THE ENTIRE CONTAINER, not just the first nEdges!!!
ADD_ACTION(sillyReduceEdges, "reduce-silly")
void sillyReduceEdges() {
	cerr << "Reducing edges sillily... ;-)\n";
	
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		edges[j] = edges[i];
		
		//if (static_cast<int>(abs(edges[i].w)) % 10 != 7) {
		if (static_cast<int>(static_cast<double>(abs(edges[i].w))) % 10 != 7) {
			++j;		// This edge lives...
		}
	}
	
	nEdges = j;
	
	// Now fix up the edge list again.  But it's already in sorted order, so don't bother resorting it.
	sortEdgesAndCalcVertIndices(false);
	edgesBySource.clear();		// So that anyone who needs these edges in future knows to recalculate them
}

ADD_ACTION(revertVertexIds, "unrenumber-verts")
void revertVertexIds() {
	cerr << "Renumbering vertices back to their original IDs...\n";
	
	// All the work is now done by the previous call to renameVertices(true).
	renameVertices(false);
}

void calcChildUpperBoundFor(int v) {
	//cerr << "calcChildUpperBoundFor(v=" << v << ") called.  Already processed: " << (seenUbChild[v] ? "yes" : "no") << ".\n";		//DEBUG
	if (!seenUbChild[v]) {
		// Haven't already calculated this.
		
		// First make sure we have calculated UBs for each of our children.
		// Also accumulate children into sets based on their colours.
		//HACK: Might be slightly faster to just use a map<int, double> bestEdgeToColour instead and accumulate the max in this same loop...
		map<int, vector<edge> > verticesOfColour;		// Technically we only need the endpoint and the weight, but who cares...
		//cerr << "   " << (startOfEdgesFrom[v + 1] - startOfEdgesFrom[v]) << " children to process for " << v << ".\n";		//DEBUG
		for (int i = startOfEdgesFrom[v]; i < startOfEdgesFrom[v + 1]; ++i) {
			calcChildUpperBoundFor(edgesBySource[i].v);
			verticesOfColour[c[edgesBySource[i].v]].push_back(edgesBySource[i]);
		}
		
		// We can sum UBs across colours.
		double x = 0.0;
		for (map<int, vector<edge> >::const_iterator i(verticesOfColour.begin()); i != verticesOfColour.end(); ++i) {
			// For each colour, we can choose the best UB of any child of that colour, or not to include any child of this colour at all if doing so would add a negative score.
			double best = 0.0;
			for (int j = 0; j < (*i).second.size(); ++j) {
				best = max(best, (*i).second[j].w + ub[(*i).second[j].v]);
			}
			
			x += best;
		}
		
		// Is the existing (e.g. Colourful Forest) UB for this vertex better?  If so, grab it!
		x = min(x, ub[v]);
		
		// Compute some interesting stats
		if (x == 0.0 && startOfEdgesFrom[v] != startOfEdgesFrom[v + 1]) {
			++nTotalInternalVerticesWithZeroChildUpperBound;
			nTotalEdgesDeletedDueToZeroChildUpperBound += startOfEdgesFrom[v + 1] - startOfEdgesFrom[v];
		}
		highestChildUpperBoundEver = max(highestChildUpperBoundEver, x);
		
		ub[v] = x;
		seenUbChild[v] = true;
	}
	//cerr << "calcUpperBoundFor(v=" << v << ") finished.\n";		//DEBUG
}

struct compare_by_startpoint_incr_then_endpoint_incr {
	bool operator()(edge a, edge b) const {
		return a.u < b.u || (a.u == b.u && a.v < b.v);
	}
};

//struct compare_by_endpoint_incr_then_startpoint_incr {
//	bool operator()(edge a, edge b) const {
//		return a.v < b.v || (a.v == b.v && a.u < b.u);
//	}
//};

//HACK: This function and the following should be more "symmetrical" to sortEdgesAndCalcVertIndices()...
void calcVertIndicesForEdgesBySource() {
	startOfEdgesFrom.resize(nVerts + 1);		// Need room for one more entry at the end
	startOfEdgesFrom[0] = 0;
	int j = 0;
	for (int i = 1; i <= nVerts; ++i) {		// Note we do fill out startOfEdgesTo[nVerts], so this array must be of size nVerts + 1.
		while (j < nEdges && edgesBySource[j].u < i) {
			++j;
		}
		
		startOfEdgesFrom[i] = j;
	}
	assert(j == nEdges);
}

// If not already calculated, will populate edgesBySource[] and startOfEdgesFrom[] arrays from edges[].
void calcEdgesFrom() {
	if (edgesBySource.empty()) {
		cerr << "Sorting "<< nEdges << " edges by source vertex...\n";
		
		//HACK: Pretty disgusting code near-duplication here...
		// Sort edge list in edgesBySource[] by source vertex, so that all out-edges for a given vertex appear contiguously.
		// Then update the array of indices that allow each destination vertex's in-edge list to be quickly found.
		assert(nEdges == edges.size());
		edgesBySource = edges;
		sort(edgesBySource.begin(), edgesBySource.begin() + nEdges, compare_by_startpoint_incr());
		//sort(edgesBySource.begin(), edgesBySource.begin() + nEdges, compare_by_startpoint_incr_then_endpoint_incr());
		
		calcVertIndicesForEdgesBySource();
	}
}

// Calculates my vertex UBs, *assuming the ub[] already holds other valid UBs if it is non-empty*.
ADD_ACTION(calcChildVertexUpperBounds, "calc-child-vertex-ubs")
void calcChildVertexUpperBounds() {
	cerr << "Calculating child vertex upper bounds (using any vertex upper bounds already computed)...\n";
	
	calcEdgesFrom();
	
	ub.resize(nVerts, INF);		// Crucially this will not disrupt ub[] if it already contains useful UBs.
	seenUbChild.assign(nVerts, false);
	nTotalInternalVerticesWithZeroChildUpperBound = 0;
	nTotalEdgesDeletedDueToZeroChildUpperBound = 0;
	highestChildUpperBoundEver = NEGINF;
	for (int i = 0; i < nVerts; ++i) {
		//cerr << "********** Starting calcVertexUpperBounds() loop iter " << i << "**********\n";		//DEBUG
		calcChildUpperBoundFor(i);
	}
	
	dumpTable("calcChildVertexUpperBounds", ub);
	
	cerr << nTotalInternalVerticesWithZeroChildUpperBound << " internal vertices have an upper bound of 0, resulting in " << nTotalEdgesDeletedDueToZeroChildUpperBound << " potential edge deletions.\n";
	cerr << "Worst upper bound for any vertex: " << highestChildUpperBoundEver << "\n";
}

// What it used to be called, left in for backcompat.
ADD_ACTION(calcTimVertexUpperBounds, "tim-vertex-ubs")
void calcTimVertexUpperBounds() {
	calcChildVertexUpperBounds();
}

// You can run this action to force the vertex UBs to be recalculated from scratch by a following call to calcVertex...UpperBounds().
// (E.g. after some edges have been deleted, these UBs are still valid so the vector containing them will not be
// automatically cleared -- but they may now be loose.)
ADD_ACTION(clearVertexUpperBounds, "clear-vertex-ubs")
void clearVertexUpperBounds() {
	if (!ub.empty()) {
		cerr << "Clearing vertex upper bounds.\n";
		ub.clear();
	}
}

//HACK: I know, I know, I should use Boost.  But I don't want to bring all that in now...
// This is necessary because arrays aren't copyable or assignable in C++, so you can't use them inside a vector.
template <typename T, size_t N>
struct array {
	T operator[](size_t i) const {
		return _array[i];
	}
	
	T& operator[](size_t i) {
		return _array[i];
	}
	
	size_t size() const {
		return N;
	}
	
private:
	T _array[N];
};

struct inEdgeToColourInfo {
	edge edges[2];		// edges[0] is the best edge, edges[1] is the second-best
	double maxInEdge;	// The weight of the greatest reachable edge ending at the start of edges[0]
};

inline void colForestUbMergeEdge(inEdgeToColourInfo& bestInEdgeToColour, edge newEdge, double newEdgeMaxInEdge) {
	if (newEdge.w >= bestInEdgeToColour.edges[0].w) {		// Important to make this >=, so 2 edges of equal, max weight edges will become 1st and 2nd best
		// IMPORTANT: We need to check that the edge being added is not already here!  This can happen when the same edge is recorded as the best in-edge
		// for a colour via 2 different children, because the ">=" (as opposed to ">") will allow it!  This fixes a longstanding bug.
		if (newEdge.u == bestInEdgeToColour.edges[0].u && newEdge.v == bestInEdgeToColour.edges[0].v) {
			// IMPORTANT: Even though the edge is the same, it might itself now be accessible from another, better edge.
			//HACK: This seems fragile as we are relying on the same edge to have its weight compare <= itself in order to perform this necessary-for-correctness update...  Dunno if FP is actually sane enough to guarantee this...
			bestInEdgeToColour.maxInEdge = max(bestInEdgeToColour.maxInEdge, newEdgeMaxInEdge);
			return;
		}
		
		bestInEdgeToColour.edges[1] = bestInEdgeToColour.edges[0];
		bestInEdgeToColour.edges[0] = newEdge;
		bestInEdgeToColour.maxInEdge = newEdgeMaxInEdge;
	} else if (newEdge.w > bestInEdgeToColour.edges[1].w) {
		bestInEdgeToColour.edges[1] = newEdge;		// We don't record any maxInEdge for the 2nd-best edge
	}
	
	// Check we don't have 2 copies of the exact same edge!
	assert(bestInEdgeToColour.edges[0].u == -1 || bestInEdgeToColour.edges[0].u != bestInEdgeToColour.edges[1].u || bestInEdgeToColour.edges[0].v != bestInEdgeToColour.edges[1].v);
}

// Used for running accumulate() on edges to sum their weights.
struct accum_weights_from_first {
	double operator()(double total, inEdgeToColourInfo e2) const {
		return total + e2.edges[0].w;
	}
};

ADD_ACTION(enableColForestUpperBoundStrengthening, "enable-col-forest-vub-strength")
void enableColForestUpperBoundStrengthening() {
	cerr << "Enabled strengthening for Colourful Forest vertex UB.\n";
	shouldStrengthenColForestUpperBound = true;
}

// What it used to be called, left in for backcompat.
ADD_ACTION(enableSebUpperBoundStrengthening, "enable-seb-vub-strength")
void enableSebUpperBoundStrengthening() {
	enableColForestUpperBoundStrengthening();
}

// My strengthening of Sebastian's original bound.
// The idea is that if the graph we have is not a subtree, then we can find a bound
// on how many edges need to be changed, and from that, find a bound on how much we can safely reduce the score by.
//double calcStrengthenedColForestUpperBoundFor(int v, vector<vector<inEdgeToColourInfo> > const& colForestUbBestInEdges) {
double calcStrengthenedColForestUpperBoundFor(int v, vector<inEdgeToColourInfo> const& bestInEdges) {
	// Calculate a tentative UB (same as unstrengthened).
	double x = accumulate(bestInEdges.begin(), bestInEdges.end(), 0.0, accum_weights_from_first());
	double DEBUGorigSum = x;		//DEBUG
	
	// Now look for violations of the subtree property.
	vector<vector<int> > coloursFromColour(nCol);
	for (int i = 0; i < nCol; ++i) {
		if (bestInEdges[i].edges[0].u != -1) {
			int fromCol = c[bestInEdges[i].edges[0].u];
			coloursFromColour[fromCol].push_back(i);
		}
	}
	
	vector<bool> safeToMessWithInEdgeToColour(nCol, false);
	
	// We can consider all out-edges leaving a particular colour as an independent problem to fix.
	for (int i = c[v]; i < nCol; ++i) {		// Only need to consider colours >= v's colour.
		//HACK: May be a bit expensive to allocate an nVerts-size array in this loop...
		//vector<double> repairCostsFromVertex(nVerts, 0.0);		// repairCostsFromVertex[i] is the cost to repair all edges from vertex i to some other vertex that are in the current solution.
		map<int, double> repairCostsFromVertex;		// Should be faster than a complete array, since there will only be at most nCol entries.  Has the advantage that a newly-created entry gets value 0.0.
		double totalRepairCost = 0.0;
		for (int j = 0; j < coloursFromColour[i].size(); ++j) {
			// Note that although this loop seems to add an nCol^2 factor to the runtime, it doesn't because the total number
			// of iterations across all values of i is just the number of colours.
			
			// We know that bestInEdges[coloursFromColour[i][j]].edges[1].w cannot be NEGINF because we know there is at least 1 edge to coloursFromColour[i][j], so it will be min(0.0, that).
			double cost = bestInEdges[coloursFromColour[i][j]].edges[0].w - bestInEdges[coloursFromColour[i][j]].edges[1].w;
			repairCostsFromVertex[bestInEdges[coloursFromColour[i][j]].edges[0].u] += cost;
			totalRepairCost += cost;
		}
		
		double bestRepairCost = totalRepairCost;		// Corresponds to repairing *all* outgoing edges.
		//cerr << "calcColForestUpperBoundFor(v=" << v << "): Trying to delete all " << coloursFromColour[i].size() << " outgoing edges from col " << i << ", with cost " << totalRepairCost << ".\n";		//DEBUG
		
		// Try making each i-coloured vertex u the parent of all these edges.
		// We only bother trying the i-coloured vertices that are actually the starting point of
		// an out-edge, since inEdgeRepairCost must be >= 0 and it will cost totalRepairCost + inEdgeRepairCost
		// to move all these out-edges to some different i-coloured vertex, and we already consider
		// repairing all out-edges for a total cost of totalRepairCost.
		for (int j = 0; j < coloursFromColour[i].size(); ++j) {
			int u = bestInEdges[coloursFromColour[i][j]].edges[0].u;
			
			double inEdgeRepairCost;
			double maxInEdge = bestInEdges[coloursFromColour[i][j]].maxInEdge;		//HACK: maybe get rid of this var altogether.
			
			// We know !coloursFromColour[i].empty(), because otherwise this loop would not have executed!
			if (i > c[v] && bestInEdges[i].edges[0].u == -1) {
				// We have no in-edges to this colour so far but at least one out-edge, so we can add in the best edge to this vertex (which will be negative, since otherwise it would already be present)
				//HACK: We no longer compute this exactly because it's too expensive -- instead we just take
				// the second-best edge to this colour, which is a UB.
				//inEdgeRepairCost = -maxInEdgeToVertex[u];
				//inEdgeRepairCost = -bestInEdges[c[u]][1].w;
				inEdgeRepairCost = -maxInEdge;
			} else {
				if (i > c[v] && bestInEdges[i].edges[0].u != -1 && bestInEdges[i].edges[0].v != u && safeToMessWithInEdgeToColour[i]) {
					//inEdgeRepairCost = bestInEdge[i].w - maxInEdgeToVertex[u];
					//inEdgeRepairCost = bestInEdges[i][0].w - bestInEdges[c[u]][1].w;
					inEdgeRepairCost = bestInEdges[i].edges[0].w - maxInEdge;
					totalSneakyStrengtheningColForestUpperBound += inEdgeRepairCost;			// To answer the question: Does this do anything at all?
				} else {
					inEdgeRepairCost = 0.0;
				}
			}
			
			double tryUsingU = totalRepairCost - repairCostsFromVertex[u] + inEdgeRepairCost;		// This will work out to be 0 in the valid case, i.e. when all out-edges leave from the same vertex, and this vertex has its colour's in-edge.
			assert(tryUsingU >= 0.0);		//DEBUG: Sanity check.  We should never be able to improve things!
			//cerr << "calcColForestUpperBoundFor(v=" << v << "): Trying to change vertex for col " << i << " to " << u << " (start of best in-edge for colour " << coloursFromColour[i][j] << "; has end " << bestInEdges[coloursFromColour[i][j]].edges[0].v << "), with cost " << tryUsingU << ".\n";		//DEBUG
			bestRepairCost = min(bestRepairCost, tryUsingU);
		}
		
		if (bestRepairCost == 0.0) {		// I know, I know -- FP equality comparisons are The Devil.  But an accidental nonzero here doesn't compromise safety
			// Everything about this colour looks good.  So mark all colours we get to from its out-edges
			// as safe to mess with.  (We'll be processing these later.)
			// Also note this only takes O(nCol) time across all iterations.
			for (int j = 0; j < coloursFromColour[i].size(); ++j) {
				safeToMessWithInEdgeToColour[coloursFromColour[i][j]] = true;
			}
		}
		
		x -= bestRepairCost;
		totalStrengtheningColForestUpperBound += bestRepairCost;
	}
	
	assert(x <= DEBUGorigSum);		//DEBUG: Sanity check -- we can't do worse than just summing everything!
	return x;
}

// Used for initialising colForestUbBestInEdges[][] in calcColForestVertexUpperBounds, and an equivalent data structure in reduceColSubtreeAdvantage().
inEdgeToColourInfo initInEdgeToColourInfo = {
	{
		{
			-1,		// edges[0].u
			-1,		// edges[0].v
			0.0		// edges[0].w
#ifdef USE_OC
			, 0.0		// edges[0].oc
#endif	// USE_OC
		},
		{
			-1,		// edges[1].u
			-1,		// edges[1].v
			NEGINF	// edges[1].w
#ifdef USE_OC
			, 0.0		// edges[1].oc
#endif	// USE_OC
		}
	},
	NEGINF		// maxInEdge
};

// Calculates Colourful Forest (Sebastian's) vertex UBs.  If ub[] contains other UBs already, they'll be overwritten.
// We now do this using a "reaching" DP that uses iteration instead of recursion, relying on the
// fact that vertices have been topologically sorted, so that the last-numbered vertex has no children.
ADD_ACTION(calcColourfulForestVertexUpperBounds, "calc-col-forest-vertex-ubs")
void calcColourfulForestVertexUpperBounds() {
	cerr << "Calculating Colourful Forest vertex upper bounds using reaching DP...\n";
	
	if (shouldCheckPreconds) {
		checkVerticesAreTopSorted("Vertices must be numbered in topological order to calculate Colourful Forest vertex UB!  Use renumber-verts (and later unrenumber-verts) to do this, or modify the code to filter edges starting from the first instead of startOfEdges[v].");
	}
	
	calcEdgesFrom();
	
	ub.resize(nVerts, INF);		// Crucially this will not disrupt ub[] if it already contains useful UBs.
	nTotalInternalVerticesWithZeroColForestUpperBound = 0;
	nTotalEdgesDeletedDueToZeroColForestUpperBound = 0;
	highestColForestUpperBoundEver = NEGINF;
	totalStrengtheningColForestUpperBound = 0.0;
	totalSneakyStrengtheningColForestUpperBound = 0.0;
	nOrphanColours = 0;
	
	// Simpler but burns a lot of memory to allocate all at once; let's allocate as we go instead.
	//vector<vector<inEdgeToColourInfo> > colForestUbBestInEdges(nVerts, vector<inEdgeToColourInfo>(nCol, initInEdgeToColourInfo));	// colForestUbBestInEdges[v][c].edges[0] is the best in-edge to colour c for the subtree rooted at v, and colForestUbBestInEdges[v][c].edges[1] is the second-best in-edge.  colForestUbBestInEdges[v][c].maxInEdge is the highest-weight in-edge to the start of ~.edges[0] that is reachable from v.
	vector<vector<inEdgeToColourInfo> > colForestUbBestInEdges(nVerts);	// colForestUbBestInEdges[v][c].edges[0] is the best in-edge to colour c for the subtree rooted at v, and colForestUbBestInEdges[v][c].edges[1] is the second-best in-edge.  colForestUbBestInEdges[v][c].maxInEdge is the highest-weight in-edge to the start of ~.edges[0] that is reachable from v.
	
	// For each child vertex v
	for (int v = nVerts - 1; v >= 0; --v) {
		// We know that colForestUbBestInEdges[v] is up-to-date, because all children of v have already been processed and pushed their info up to v.
		//cerr << "Processing in-edges for vertex " << v << "...\n";		//DEBUG
		
		// Calculate the new UB for v from this info.
		double x;
		
		if (colForestUbBestInEdges[v].empty()) {
			// This can only mean that v is a leaf, in which case its UB is trivially 0.  Leave it empty, and work around it later.
			x = 0.0;
		} else if (shouldStrengthenColForestUpperBound) {
			x = calcStrengthenedColForestUpperBoundFor(v, colForestUbBestInEdges[v]);
		} else {
			x = accumulate(colForestUbBestInEdges[v].begin(), colForestUbBestInEdges[v].end(), 0.0, accum_weights_from_first());
		}
		ub[v] = min(ub[v], x);
		
		// Compute some interesting stats
		if (x == 0.0 && startOfEdgesFrom[v] != startOfEdgesFrom[v + 1]) {
			++nTotalInternalVerticesWithZeroColForestUpperBound;
			nTotalEdgesDeletedDueToZeroColForestUpperBound += startOfEdgesFrom[v + 1] - startOfEdgesFrom[v];
		}
		highestColForestUpperBoundEver = max(highestColForestUpperBoundEver, x);
		
		// For each parent u of v
		for (int iu = startOfEdgesTo[v]; iu < startOfEdgesTo[v + 1]; ++iu) {
			int u = edges[iu].u;
			//cerr << "Processing in-edge " << iu << ": " << edges[iu] << ".\n";		//DEBUG
			
			// Lazily initialise info for u
			if (colForestUbBestInEdges[u].empty()) {
				colForestUbBestInEdges[u].resize(nCol, initInEdgeToColourInfo);
			}
			
			// Only bother merging in results for vertices that are not leaves.  This lets us avoid allocating
			if (!colForestUbBestInEdges[v].empty()) {
				// Merge the info about the best in-edges to each colour for v into its parent u
				for (int k = c[v] + 1; k < nCol; ++k) {			// There can't be any entries for v's own colour or lower
					colForestUbMergeEdge(colForestUbBestInEdges[u][k], colForestUbBestInEdges[v][k].edges[0], colForestUbBestInEdges[v][k].maxInEdge);
					colForestUbMergeEdge(colForestUbBestInEdges[u][k], colForestUbBestInEdges[v][k].edges[1], 0.0);		//HACK: maxInEdge parameter is meaningless here, because it can never be used
					
					if (colForestUbBestInEdges[u][k].edges[0].u == v) {
						// The best in-edge for colour k is actually the one from u's child v -- so update the maxInEdge for this edge.
						// This is (the only place) where we calculate new values for maxInEdge.
						colForestUbBestInEdges[u][k].maxInEdge = max(colForestUbBestInEdges[u][k].maxInEdge, edges[iu].w);
					}
				}
			}
			
			// Now merge in the edge (u, v) to u's info.
			colForestUbMergeEdge(colForestUbBestInEdges[u][c[v]], edges[iu], NEGINF);
		}
		
		// We have applied everything known about v to all its parents, so we no longer need v's info.
		colForestUbBestInEdges.resize(v);			// Shrinks the vector by 1, effectively deleting v's info
	}
	
	dumpTable("calcColForestVertexUpperBounds", ub);
	
	cerr << "Worst col-forest-vertex-UB for any vertex: " << highestColForestUpperBoundEver << "\n";
	cerr << nTotalInternalVerticesWithZeroColForestUpperBound << " internal vertices with ColForest UB of 0; this would result in " << nTotalEdgesDeletedDueToZeroColForestUpperBound << " edge deletions.\n";
	if (shouldStrengthenColForestUpperBound) {
		cerr << "Total strengthening achieved: " << totalStrengtheningColForestUpperBound << "\n";
		cerr << "Total sneakiness (not all or perhaps any of this was actually used for strengthening): " << totalSneakyStrengtheningColForestUpperBound << "\n";
		cerr << "Number of colours with out-edges but no in-edges: " << nOrphanColours << "\n";
	}
}

// What it used to be called, left in for backcompat.
ADD_ACTION(calcSebVertexUpperBounds, "seb-vertex-ubs")
void calcSebVertexUpperBounds() {
	calcColourfulForestVertexUpperBounds();
}

ADD_ACTION(vertexUbReduceEdges, "reduce-vub")
void vertexUbReduceEdges() {
	cerr << "Reducing edges using vertex upper bounds...\n";
	
	if (ub.empty()) {
		cerr << "Must compute vertex upper bounds (using calc-child-vertex-ubs and/or calc-col-forest-vertex-ubs) before you can apply reduce-vub!\n";
		exit(1);
	}
	
	int j = 0;
	double scaredEdge = INF;
	for (int i = 0; i < nEdges; ++i) {
		edges[j] = edges[i];
		
		if (edges[i].w + ub[edges[i].v] > 0) {
		// The following version doesn't actually help, since any edge that it would delete would have already been deleted by the
		// stronger reduce-slide-strong or reduce-dompath2 reduction that produced the oc value in the first place.
		//if (edges[i].w - edges[i].oc + ub[edges[i].v] > 0) {
			++j;		// This edge lives...
			scaredEdge = min(scaredEdge, edges[i].w + ub[edges[i].v]);
		}
	}
	
	cerr << "Deleted " << (nEdges - j) << " edges.\n";
	cerr << "The most scared edge survives by " << scaredEdge << ".\n";		// Gives us a clue how close we are to pruning the next edge
	
	nEdges = j;
	edges.erase(edges.begin() + nEdges, edges.end());		// Keep container size in sync with nEdges
	
	sortEdgesAndCalcVertIndices(false);
	edgesBySource.clear();		// So that anyone who needs these edges in future knows to recalculate them
}

ADD_ACTION(reduceNegativePendantEdges, "reduce-negpend")
void reduceNegativePendantEdges() {
	cerr << "Reducing negative pendant edges...\n";
	
	int nEdgesDeleted = 0;
	
	calcEdgesFrom();
	
	for (int iPass = 1; iPass == 1 || nEdgesDeleted > 0; ++iPass) {
		cerr << "reduceNegativePendantEdges() pass " << iPass << "...\n";
		
		int j = 0;
		for (int i = 0; i < nEdges; ++i) {
			edgesBySource[j] = edgesBySource[i];
			
			// Note: it doesn't matter that we destroy the information in edgesBySource[] being "pointed to" by startOfEdgesFrom[], since all we care about is the number of children.
			if (edgesBySource[i].w > 0 || startOfEdgesFrom[edgesBySource[i].v] != startOfEdgesFrom[edgesBySource[i].v + 1]) {
				++j;		// This edge lives..
			}
		}
		
		nEdgesDeleted = nEdges - j;
		cerr << "Deleted " << nEdgesDeleted << " edges.\n";
		
		nEdges = j;
		edgesBySource.erase(edgesBySource.begin() + nEdges, edgesBySource.end());		// Keep container size in sync with nEdges
		
		calcVertIndicesForEdgesBySource();
	}
	
	// Update "main" edge list
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);
}

ADD_ACTION(reduceUnreachableEdges, "reduce-unreach")
void reduceUnreachableEdges() {
	cerr << "Reducing unreachable edges...\n";
	
	if (shouldCheckPreconds) {
		checkVerticesAreTopSorted("Vertices must be numbered in topological order to remove unreachable edges in a single pass!  Use renumber-verts (and later unrenumber-verts) to do this.");
	}
	
	calcEdgesFrom();
	
	// Only need 1 pass if vertices have been renamed so that (u, v) => u < v and we scan edges in source-vertex order :-)
	vector<int> nInEdges(nVerts, 0);
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		edgesBySource[j] = edgesBySource[i];
		
		if (edgesBySource[i].u == 0 || nInEdges[edgesBySource[i].u]) {		// NOTE: the first vertex is allowed to have no in-edges!
			// At least 1 in-edges to this edge: keep its out-edges
			++nInEdges[edgesBySource[i].v];
			++j;		// This edge lives...
		}
	}
	
	// The number of unreachable vertices is an interesting stat.
	int nUnreachableVertices = 0;
	for (int i = 0; i < nVerts; ++i) {
		if (!nInEdges[i]) {
			++nUnreachableVertices;
		}
	}
	
	int nEdgesDeleted = nEdges - j;
	cerr << "reduceUnreachableEdges(): Deleted " << nEdgesDeleted << " edges.  " << nUnreachableVertices << " unreachable vertices (including the root) remain; these will not be deleted, to preserve vertex numbering.\n";
	
	nEdges = j;
	edgesBySource.erase(edgesBySource.begin() + nEdges, edgesBySource.end());		// Keep container size in sync with nEdges
	
	calcVertIndicesForEdgesBySource();
	
	// Update "main" edge list
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);
}

ADD_ACTION(nop, "nop")
void nop() {
	// What it says on the tin.
}

// m_[i][j - firstVertOfSameColourAs[j]], for i and j of the same colour, is a lower bound on the score change that could result by deleting the subtree rooted at i and adding a subtree rooted at j.  Always <= 0.
// Access it through slideLb(i, j), which does the subtraction.
//vector<vector<double> > m;				// m[i][j], for i and j of the same colour, is a lower bound on the score change that could result by deleting the subtree rooted at i and adding a subtree rooted at j.  Always <= 0.
vector<vector<double> > m_;
vector<double> worstSlideLb;		// worstSlideLb[u] is the worst value of slideLb(u, v) for any v of the same colour as u -- i.e. a LB if we don't know which vertex we'll slide u to.

inline double slideLb(int u, int v) {
	return m_[u][v - firstVertOfSameColourAs[v]];
}

// Populate m[][] using the simplest, crappiest way: whenever u != v, just take the -ve of the UB for u, and for u == v, take 0.
ADD_ACTION(calcCrappySlideLowerBound, "calc-crappy-slide-lbs")
void calcCrappySlideLowerBound() {
	if (ub.empty()) {
		cerr << "Must compute vertex upper bounds (using calc-child-vertex-ubs and/or calc-col-forest-vertex-ubs) before you can apply calc-crappy-slide-lbs!\n";
		exit(1);
	}
	
	cerr << "Calculating crappy slide lower bounds...\n";
	vector<vector<int> > vertsOfColour(nCol);
	for (int i = 0; i < nVerts; ++i) {
		vertsOfColour[c[i]].push_back(i);
	}
	
	//m.resize(nVerts, vector<double>(nVerts, NEGINF));
	//m_.resize(nVerts);
	m_.assign(nVerts, vector<double>());
	for (int u = 0; u < nVerts; ++u) {
		m_[u].resize(vertsOfColour[c[u]].size(), NEGINF);
	}
	
	for (int col = 0; col < nCol; ++col) {
		//cerr << "Calculating slide bound for the " << vertsOfColour[col].size() << " vertices of colour " << col << "...\n";
		for (int i = 0; i < vertsOfColour[col].size(); ++i) {
			for (int j = 0; j < vertsOfColour[col].size(); ++j) {
				if (i == j) {
					//m[vertsOfColour[col][i]][vertsOfColour[col][j]] = 0.0;
					m_[vertsOfColour[col][i]][vertsOfColour[col][j] - firstVertOfSameColourAs[vertsOfColour[col][j]]] = 0.0;
				} else {
					//m[vertsOfColour[col][i]][vertsOfColour[col][j]] = -ub[vertsOfColour[col][i]];
					m_[vertsOfColour[col][i]][vertsOfColour[col][j] - firstVertOfSameColourAs[vertsOfColour[col][j]]] = -ub[vertsOfColour[col][i]];
				}
			}
		}
	}
}

struct compare_by_endpoint_colour_asc_then_endpoint_asc {
	bool operator()(edge a, edge b) const {
		if (c[a.v] != c[b.v]) {
			return c[a.v] < c[b.v];
		}
		
		return a.v < b.v;
	}
};

ADD_ACTION(calcSlideLowerBound, "calc-slide-lbs")
void calcSlideLowerBound() {
	if (ub.empty()) {
		cerr << "Must compute vertex upper bounds (using calc-child-vertex-ubs and/or calc-col-forest-vertex-ubs) before you can apply calc-slide-lbs!\n";
		exit(1);
	}
	
	cerr << "Calculating slide lower bounds...\n";
	vector<vector<int> > vertsOfColour(nCol);
	for (int i = 0; i < nVerts; ++i) {
		vertsOfColour[c[i]].push_back(i);
	}
	
	//m.resize(nVerts, vector<double>(nVerts, NEGINF));
	//m.assign(nVerts, vector<double>(nVerts, NEGINF));
	m_.assign(nVerts, vector<double>());
	for (int u = 0; u < nVerts; ++u) {
		m_[u].resize(vertsOfColour[c[u]].size(), NEGINF);
	}
	
	for (int col = 0; col < nCol; ++col) {
		//cerr << "Calculating slide bound for the " << vertsOfColour[col].size() << " vertices of colour " << col << "...\n";
		for (int i = 0; i < vertsOfColour[col].size(); ++i) {
			int u = vertsOfColour[col][i];
			vector<edge> uChildren(edgesBySource.begin() + startOfEdgesFrom[u], edgesBySource.begin() + startOfEdgesFrom[u + 1]);
			sort(uChildren.begin(), uChildren.end(), compare_by_endpoint_colour_asc_then_endpoint_asc());		// Sort by colour of destination vertex, then by destination vertex
			
			for (int j = 0; j < vertsOfColour[col].size(); ++j) {
				int v = vertsOfColour[col][j];
				
				double best = -ub[u];		// This will sometimes turn out to be better than the more complicated cases below, because they can sum edges many times via different children.
				if (u == v) {
					best = 0.0;
				} else {
					vector<edge> vChildren(edgesBySource.begin() + startOfEdgesFrom[v], edgesBySource.begin() + startOfEdgesFrom[v + 1]);
					sort(vChildren.begin(), vChildren.end(), compare_by_endpoint_colour_asc_then_endpoint_asc());		// Sort by colour of destination vertex, then by destination vertex
					
					// Build a table giving the highest-weight edge to each colour from v.  If no positive-weight edge to some colour exists, use 0.
					vector<double> bestKVforColour(nCol, 0.0);		// Will hold
					for (int k = 0; k < vChildren.size(); ++k) {
						bestKVforColour[c[vChildren[k].v]] = max(bestKVforColour[c[vChildren[k].v]], vChildren[k].w);
					}
					
					// Now merge the two lists of children in linear time, deciding the best option for each child of u, and the worst child for each colour.
					// Note that we allow the case where kV == vChildren.size(), so we have to check for this before any access to vChildren[vK].
					int kV = 0;
					double worstForColour = INF;
					double total = 0.0;
					for (int kU = 0; kU < uChildren.size(); ++kU) {
						double bestForKid;
						
						// Skip any irrelevant-coloured children of v
						while (kV < vChildren.size() && c[uChildren[kU].v] > c[vChildren[kV].v]) {
							++kV;
						}
						
						if (kV == vChildren.size() || c[uChildren[kU].v] < c[vChildren[kV].v]) {
							// There are no more children of v to consider for the current colour.
							bestForKid = min(0.0, -uChildren[kU].w - ub[uChildren[kU].v]);
						} else {
							// At this point, we must have a u-child and a v-child of the same colour.
							
							// Skip any irrelevant children of v
							while (kV < vChildren.size() && uChildren[kU].v > vChildren[kV].v) {
								++kV;
							}
							
							// At this point, we either have uChildren[kU].v == vChildren[kV].v (u and v share this child), or
							// (kV == vChildren.size() || uChildren[kU].v < vChildren[kV].v) (v does not have this child of u).
							
							// Option 1: Delete the edge from u to this child and any subtree hanging under it, then add in the best
							// possible edge from v to any vertex of this child's colour (or 0 if there is no such positive-weight edge from v).
							// This can always be done.  Note that it's also possible that the edge (u, x) does not exist in the subtree
							// rooted at u, meaning that we are not allowed to use a vertex of its colour after all (since this vertex, or another vertex
							// of its colour, may be used by a path somewhere else in the tree), meaning we can never go above 0.
							bestForKid = min(0.0, bestKVforColour[c[uChildren[kU].v]] - uChildren[kU].w - ub[uChildren[kU].v]);
							if (kV < vChildren.size() && uChildren[kU].v == vChildren[kV].v) {
								// v has this child too, so there's a second option: delete the edge from u to this child, and add the edge from v to it.
								// As before, we have to consider the case that (u, x) is not in the tree, so we cap at 0.
								bestForKid = max(bestForKid, min(0.0, vChildren[kV].w - uChildren[kU].w));
							}
						}
						
						// Is this the worst child for this colour so far?
						worstForColour = min(worstForColour, bestForKid);
						
						// Have we finished processing this colour?
						if (kU == uChildren.size() - 1 || c[uChildren[kU + 1].v] > c[uChildren[kU].v]) {
							// Yes we have.
							total += worstForColour;
							worstForColour = INF;
						}
					}
					
					best = max(best, total);
				}
				
				//m[u][v] = max(m[u][v], best);
				m_[u][v - firstVertOfSameColourAs[v]] = max(m_[u][v - firstVertOfSameColourAs[v]], best);
			}
		}
	}
}

struct compare_by_endpoint_colour_asc {
	bool operator()(edge a, edge b) const {
		return c[a.v] < c[b.v];
	}
};

ADD_ACTION(calcRecursiveSlideLowerBound, "calc-rec-slide-lbs")
void calcRecursiveSlideLowerBound() {
	if (ub.empty()) {
		cerr << "Must compute vertex upper bounds (using calc-child-vertex-ubs and/or calc-col-forest-vertex-ubs) before you can apply calc-rec-slide-lbs!\n";
		exit(1);
	}
	
	if (lbCol.empty()) {
		cerr << "Must compute anchor-to-colour lower bounds (using calc-anc-col-lbs) before you can apply calc-rec-slide-lbs!\n";
		exit(1);
	}
	
	if ( firstVertOfSameColourAs.empty() ) {
		calcFirstVertOfSameColour();
	}
	
	cerr << "Calculating recursive slide lower bounds...\n";
	vector<vector<int> > vertsOfColour(nCol);
	for (int i = 0; i < nVerts; ++i) {
		vertsOfColour[c[i]].push_back(i);
	}
	
	//TODO: Actually the comments below are out of date -- nowadays colours DO obey this condition!
	// So we're doing a small amount of unnecessary bookkeeping
	
	// We need to process colours in order so that whenever a vertex of colour c
	// is reachable from a vertex of colour d, c is processed before d.  The colours
	// that we currently have do not necessarily obey this condition.
	// We can do this by looking up the colours of vertices in reverse vertex order.
	// Alternatively, we could have the vertex-renaming code also rename colours to obey this condition (this is probably the nicer way actually).
	vector<int> colourOrder;
	vector<bool> seenColour(nCol, false);
	for (int i = nVerts - 1; i >= 0; --i) {
		if (!seenColour[c[i]]) {
			colourOrder.push_back(c[i]);
			seenColour[c[i]] = true;
		}
	}
	
	//cerr << "There appear to be vertices of " << colourOrder.size() << " distinct colours.\n";		//DEBUG
	assert(colourOrder.size() <= nCol);		// There may be fewer colours than declared if there have been vertex deletions, though this program doesn't (currently) perform any.
	
	//m.resize(nVerts, vector<double>(nVerts, NEGINF));
	//m.assign(nVerts, vector<double>(nVerts, NEGINF));
	m_.assign(nVerts, vector<double>());
	for (int u = 0; u < nVerts; ++u) {
		m_[u].resize(vertsOfColour[c[u]].size(), NEGINF);
	}
	
	worstSlideLb.assign(nVerts, INF);		// Will be overwritten with something lower for each vertex
	
	int nInitImprovementsUsingLbCol = 0, nOverallImprovementsUsingLbCol = 0;		//DEBUG
	
	for (int iCol = 0; iCol < colourOrder.size(); ++iCol) {
		int col = colourOrder[iCol];		// Note the difference from calcSlideLowerBound
		
		// First, for each vertex of this colour, sort all its out-edges by colour.  We'll undo this at the end.
		for (int i = 0; i < vertsOfColour[col].size(); ++i) {
			int u = vertsOfColour[col][i];
			sort(edgesBySource.begin() + startOfEdgesFrom[u], edgesBySource.begin() + startOfEdgesFrom[u + 1], compare_by_endpoint_colour_asc());
		}
		
		//cerr << "Calculating recursive slide bound for the " << vertsOfColour[col].size() << " vertices of colour " << col << "...\n";
		for (int i = 0; i < vertsOfColour[col].size(); ++i) {
			int u = vertsOfColour[col][i];
			for (int j = 0; j < vertsOfColour[col].size(); ++j) {
				int v = vertsOfColour[col][j];
				
				double best = -ub[u];		// This will sometimes turn out to be better than the more complicated cases below, because they can sum edges many times via different children.
				if (u == v) {
					best = 0.0;
				} else {
					int startKV = startOfEdgesFrom[v];
					double worstForColour = 0.0;		// Corresponds to case in which no child of this colour appears under u.
					double total = 0.0;
					for (int kU = startOfEdgesFrom[u]; kU < startOfEdgesFrom[u + 1]; ++kU) {
						// Skip any v-children of irrelevant colours
						while (startKV < startOfEdgesFrom[v + 1] && c[edgesBySource[kU].v] > c[edgesBySource[startKV].v]) {
							++startKV;
						}
						
						// Always try just deleting this edge and any subtree under it.
						double bestForKid = -ub[edgesBySource[kU].v];
						if (startKV == startOfEdgesFrom[v + 1] || c[edgesBySource[kU].v] < c[edgesBySource[startKV].v]) {
							// There are no more children of v to consider for the current colour.
						} else {
							// v has children of the same colour as this child of u.  Try all of them.
							bestForKid += max(0.0, lbCol[v][c[edgesBySource[startKV].v]]);		//HACK: Can't add this at the definition of bestForKid because there may be no kids!
							bestForKid = max(bestForKid, lbCol[v][c[edgesBySource[startKV].v]] + worstSlideLb[edgesBySource[kU].v]);		// Alternatively, we can try adding a new edge to some unknown c[x]-coloured vertex, and shifting any tree underneath x across to it.		//TODO: Add this to the -strong version too.
							if (lbCol[v][c[edgesBySource[startKV].v]] > 0.0) {		//DEBUG
								++nInitImprovementsUsingLbCol;
							}
							for (int k = startKV; k < startOfEdgesFrom[v + 1] && c[edgesBySource[k].v] == c[edgesBySource[kU].v]; ++k) {
								bestForKid = max(bestForKid, edgesBySource[k].w + slideLb(edgesBySource[kU].v, edgesBySource[k].v));
							}
							
							//DEBUG
							if (-ub[edgesBySource[kU].v] < bestForKid && -ub[edgesBySource[kU].v] + lbCol[v][c[edgesBySource[startKV].v]] == bestForKid) {
								++nOverallImprovementsUsingLbCol;
							}
						}
						
						// Is this the worst child for this colour so far?
						worstForColour = min(worstForColour, bestForKid - edgesBySource[kU].w);
						
						// Have we finished processing this colour?
						if (kU == startOfEdgesFrom[u + 1] - 1 || c[edgesBySource[kU + 1].v] > c[edgesBySource[kU].v]) {
							// Yes we have.
							total += worstForColour;
							worstForColour = 0.0;
						}
					}
					
					best = max(best, total);
				}
				
				m_[u][v - firstVertOfSameColourAs[v]] = max(m_[u][v - firstVertOfSameColourAs[v]], best);
				worstSlideLb[u] = min(worstSlideLb[u], m_[u][v - firstVertOfSameColourAs[v]]);
			}
		}
	}
	
	// Now undo the sorting that we did before.
	//HACK: Because we don't use a stable sort, it's possible that the new order will be different than the original!
	for (int u = 0; u < nVerts; ++u) {
		sort(edgesBySource.begin() + startOfEdgesFrom[u], edgesBySource.begin() + startOfEdgesFrom[u + 1], compare_by_startpoint_incr());
	}
	
	//DEBUG
	cerr << "Improvements to initial LB due to lbCol = " << nInitImprovementsUsingLbCol << "; improvements to final LB = " << nOverallImprovementsUsingLbCol << ".\n";
	dumpTable("calcRecursiveSlideLowerBound", m_);
}

// Like calc-rec-slide-lbs, except that instead of only considering sliding children of u to vertices that are children
// of v, we slide each child to all same-coloured vertices, and use lb[][] to try to connect it to (or above) v.
ADD_ACTION(calcRecursiveSlideLowerBoundStrong, "calc-rec-slide-lbs-strong")
void calcRecursiveSlideLowerBoundStrong() {
	if (ub.empty()) {
		cerr << "Must compute vertex upper bounds (using calc-child-vertex-ubs and/or calc-col-forest-vertex-ubs) before you can apply calc-rec-slide-lbs-strong!\n";
		exit(1);
	}
	
	if (lb.empty()) {
		cerr << "You must calculate anchor lower bounds using calc-anchor-lbs before running calc-rec-slide-lbs-strong.\n";
		exit(1);
	}
	
	cerr << "Calculating strong recursive slide lower bounds...\n";
	vector<vector<int> > vertsOfColour(nCol);
	for (int i = 0; i < nVerts; ++i) {
		vertsOfColour[c[i]].push_back(i);
	}
	
	//TODO: Actually the comments below are out of date -- nowadays colours DO obey this condition!
	// So we're doing a small amount of unnecessary bookkeeping
	
	// We need to process colours in order so that whenever a vertex of colour c
	// is reachable from a vertex of colour d, c is processed before d.  The colours
	// that we currently have do not necessarily obey this condition.
	// We can do this by looking up the colours of vertices in reverse vertex order.
	// Alternatively, we could have the vertex-renaming code also rename colours to obey this condition (this is probably the nicer way actually).
	vector<int> colourOrder;
	vector<bool> seenColour(nCol, false);
	for (int i = nVerts - 1; i >= 0; --i) {
		if (!seenColour[c[i]]) {
			colourOrder.push_back(c[i]);
			seenColour[c[i]] = true;
		}
	}
	
	//cerr << "There appear to be vertices of " << colourOrder.size() << " distinct colours.\n";		//DEBUG
	assert(colourOrder.size() <= nCol);		// There may be fewer colours than declared if there have been vertex deletions, though this program doesn't (currently) perform any.
	
	//m.resize(nVerts, vector<double>(nVerts, NEGINF));
	//m.assign(nVerts, vector<double>(nVerts, NEGINF));
	m_.assign(nVerts, vector<double>());
	for (int u = 0; u < nVerts; ++u) {
		m_[u].resize(vertsOfColour[c[u]].size(), NEGINF);
	}
	
	for (int iCol = 0; iCol < colourOrder.size(); ++iCol) {
		int col = colourOrder[iCol];		// Note the difference from calcSlideLowerBound
		
		// First, for each vertex of this colour, sort all its out-edges by colour.  We'll undo this at the end.
		for (int i = 0; i < vertsOfColour[col].size(); ++i) {
			int u = vertsOfColour[col][i];
			sort(edgesBySource.begin() + startOfEdgesFrom[u], edgesBySource.begin() + startOfEdgesFrom[u + 1], compare_by_endpoint_colour_asc());
		}
		
		//cerr << "Calculating recursive slide bound for the " << vertsOfColour[col].size() << " vertices of colour " << col << "...\n";
		for (int i = 0; i < vertsOfColour[col].size(); ++i) {
			int u = vertsOfColour[col][i];
			for (int j = 0; j < vertsOfColour[col].size(); ++j) {
				int v = vertsOfColour[col][j];
				
				double best = -ub[u];		// This will sometimes turn out to be better than the more complicated cases below, because they can sum edges many times via different children.
				if (u == v) {
					best = 0.0;
				} else {
					int startKV = startOfEdgesFrom[v];
					double worstForColour = 0.0;		// Corresponds to case in which no child of this colour appears under u.
					double total = 0.0;
					for (int kU = startOfEdgesFrom[u]; kU < startOfEdgesFrom[u + 1]; ++kU) {
						// Skip any v-children of irrelevant colours
						while (startKV < startOfEdgesFrom[v + 1] && c[edgesBySource[kU].v] > c[edgesBySource[startKV].v]) {
							++startKV;
						}
						
						// Always try just deleting this edge and any subtree under it.
						double bestForKid = -ub[edgesBySource[kU].v];
						for (int ix = 0; ix < vertsOfColour[c[edgesBySource[kU].v]].size(); ++ix) {
							int x = vertsOfColour[c[edgesBySource[kU].v]][ix];
							bestForKid = max(bestForKid, lb[x][v] + slideLb(edgesBySource[kU].v, x));
						}
						
						// Is this the worst child for this colour so far?
						worstForColour = min(worstForColour, bestForKid - edgesBySource[kU].w);
						
						// Have we finished processing this colour?
						if (kU == startOfEdgesFrom[u + 1] - 1 || c[edgesBySource[kU + 1].v] > c[edgesBySource[kU].v]) {
							// Yes we have.
							total += worstForColour;
							worstForColour = 0.0;
						}
					}
					
					best = max(best, total);
				}
				
				//m[u][v] = max(m[u][v], best);
				m_[u][v - firstVertOfSameColourAs[v]] = max(m_[u][v - firstVertOfSameColourAs[v]], best);
			}
		}
	}
	
	// Now undo the sorting that we did before.
	//HACK: Because we don't use a stable sort, it's possible that the new order will be different than the original!
	for (int u = 0; u < nVerts; ++u) {
		sort(edgesBySource.begin() + startOfEdgesFrom[u], edgesBySource.begin() + startOfEdgesFrom[u + 1], compare_by_startpoint_incr());
	}
}

ADD_ACTION(clearSlideLowerBounds, "clear-slide-lbs")
void clearSlideLowerBounds() {
	//if (!m.empty()) {
	if (!m_.empty()) {
		cerr << "Clearing slide lower bounds...\n";
		//m.clear();
		m_.clear();
	}
}

// Given a vertex u that is in the solution and a colour i that has no vertices in the solution,
// what is the largest weight we can add to the graph by adding an edge from a vertex on the (unknown) path
// from the root to u (including u itself as a possibility), to some vertex of colour i?  That will be stored in lbCol[u][i].
// In theory we could even allow chains of edges, but for now, we allow only a single edge to be added,
// making this very similar to lb[][].
//
// Dealing with root-connectivity is tricky.  We could try saying "lbCol[u][iCol] must be NEGINF whenever u is disconnected", but
// that's not enough: without extra info, there's no way for u to tell whether a particular parent p is reporting
// lbCol[p][iCol] == NEGINF for some iCol because it's connected but there's no edge from p to any iCol-coloured vertex, or because p is disconnected.
// In the former case, we should include p in the minimum across all parents (which will therefore go to NEGINF too) but allow
// out-edges from u to iCol-coloured vertices to increase the overall maximum, while in the latter case
// we should not include p in the minimum across all parents (because that edge will never be part of any solution anyway), but
// if no other connected parents exist we should force the overall maximum to NEGINF even if out-edges to iCol-coloured vertices
// do exist.
//			
// Therefore, although there can of course be unreachable vertices (since we never delete a vertex), we require that there are
// no unreachable *edges*.  If such an edge (u, v) exists then the value of lbCol[u][c[v]] may be greater than NEGINF.
//
// This is computed from scratch each time, so there's no need to write a clear-anc-col-lbs action.
ADD_ACTION(calcAnchorToColourLowerBounds, "calc-anc-col-lbs")
void calcAnchorToColourLowerBounds() {
	cerr << "Calculating anchor-to-colour lower bounds...\n";
	
	calcEdgesFrom();
	
	lbCol.assign(nVerts, vector<double>(nCol, NEGINF));
	for (int u = 0; u < nVerts; ++u) {
		// The only unreachable vertices we handle correctly are the "directly" unreachable ones -- those having no in-edges.
		// So ensure that reduce-unreach has been called beforehand!
		if (startOfEdgesTo[u] != startOfEdgesTo[u + 1]) {
			lbCol[u][c[u]] = 0.0;		// Not sure if this is ever useful but it might be.
			
			// Update scores for each colour using all outgoing edges
			for (int iv = startOfEdgesFrom[u]; iv < startOfEdgesFrom[u + 1]; ++iv) {
				int v = edgesBySource[iv].v;
				lbCol[u][c[v]] = max(lbCol[u][c[v]], edgesBySource[iv].w);
			}
			
			// For each colour, see whether we can do better by connecting to it from some vertex on every path from the root to u.
			for (int iCol = 0; iCol < nCol; ++iCol) {
				double worst = INF;		// We know that worst will drop below this since we already checked that we have at least one parent.
				for (int ip = startOfEdgesTo[u]; ip < startOfEdgesTo[u + 1]; ++ip) {
					worst = min(worst, lbCol[edges[ip].u][iCol]);
				}
				
				lbCol[u][iCol] = max(lbCol[u][iCol], worst);
			}
		}
	}
}

#ifdef USE_OC
ADD_ACTION(clearEdgeOpCosts, "clear-edge-op-costs")
void clearEdgeOpCosts() {
	// Note that this might be applied when there are only edges in one of the two lists.
	for (int i = 0; i < edges.size(); ++i) {
		edges[i].oc = 0.0;
	}
	for (int i = 0; i < edgesBySource.size(); ++i) {
		edgesBySource[i].oc = 0.0;
	}
	
	cerr << "Cleared edge opportunity costs.\n";
}
#endif	// USE_OC

void DEBUGtestEdge(edge e) {
	//// This is a sorted list of edges missing from pos44247C21H42N10O6S2 after applying the following reduction program:
	//// ft_reduce read renumber-verts * ( tim-vertex-ubs seb-vertex-ubs reduce-vub ) reduce-unreach * ( seb-vertex-ubs tim-vertex-ubs calc-rec-slide-lbs reduce-slide reduce-domedge ) unrenumber-verts write
	//static edge missingEdges[] = {
	//	{ 271, 367, -999.424242 },
	//	{ 337, 420, -999.424242 },
	//	{ 337, 474, -999.424242 },
	//	{ 367, 380, -999.424242 },
	//	{ 367, 433, -999.424242 },
	//	{ 367, 664, -999.424242 },
	//	{ 380, 442, -999.424242 },
	//	{ 420, 511, -999.424242 },
	//	{ 420, 553, -999.424242 }
	//};
	//
	//edge renamedBack(e);
	//renamedBack.u = newIdFor[renamedBack.u];
	//renamedBack.v = newIdFor[renamedBack.v];
	//
	//edge* ep = lower_bound(missingEdges, missingEdges + sizeof missingEdges / sizeof missingEdges[0], renamedBack, compare_by_startpoint_incr_then_endpoint_incr());
	//if (ep < missingEdges + sizeof missingEdges / sizeof missingEdges[0] && ep->u == renamedBack.u && ep->v == renamedBack.v) {
	//	cerr << "FOUND IT!  About to delete edge " << renamedBack << "!\n";
	//	exit(2);
	//}
}

ADD_ACTION(reduceWithSlide, "reduce-slide")
void reduceWithSlide() {
	cerr << "Reducing edges with slide reduction...\n";
	
	if (m_.empty()) {
		cerr << "Must compute slide lower bounds (e.g. with calc-crappy-slide-lbs) before using reduce-slide!\n";
		exit(1);
	}
	
	int nContenders = 0;
	int nDecentContenders = 0;
	double biggestEdgeWeightGap = NEGINF;
	double bestEver = NEGINF;
	vector<edge> newEdges;			// We make a copy instead of the usual in-place deletion of edges because we need to preserve earlier edges for testing later ones.
	for (int i = 0; i < nEdges; ++i) {
		bool keep = true;
		for (int k = startOfEdgesFrom[edgesBySource[i].u]; k < startOfEdgesFrom[edgesBySource[i].u + 1]; ++k) {
			assert(edgesBySource[k].u == edgesBySource[i].u);		// Sanity
			if (c[edgesBySource[k].v] == c[edgesBySource[i].v] && edgesBySource[k].v != edgesBySource[i].v) {
				// Edge k is to a vertex of the same colour as edge i's destination vertex.
				++nContenders;
				if (edgesBySource[k].w > edgesBySource[i].w) {
					++nDecentContenders;
				}
				biggestEdgeWeightGap = max(biggestEdgeWeightGap, edgesBySource[k].w - edgesBySource[i].w);
				//bestEver = max(bestEver, edgesBySource[k].w - edgesBySource[i].w + m[edgesBySource[i].v][edgesBySource[k].v]);
				bestEver = max(bestEver, edgesBySource[k].w - edgesBySource[i].w + slideLb(edgesBySource[i].v, edgesBySource[k].v));
				
				//if (edgesBySource[k].w - edgesBySource[i].w + m[edgesBySource[i].v][edgesBySource[k].v] > 0) {
				if (edgesBySource[k].w - edgesBySource[i].w + slideLb(edgesBySource[i].v, edgesBySource[k].v) > 0) {
					keep = false;			// This edge sux yo!
					break;
				}
			}
		}
		
		if (keep) {
			newEdges.push_back(edgesBySource[i]);
		//} else {		//DEBUG
		//	DEBUGtestEdge(edgesBySource[i]);
		}
	}
	
	cerr << "Deleted " << (nEdges - newEdges.size()) << " edges.  (There were " << nContenders << " contenders, of which " << nDecentContenders << " were decent, with the best ever score difference being " << bestEver << " and the biggest gap in edge weights being " << biggestEdgeWeightGap << ".)\n";
	
	nEdges = newEdges.size();
	edgesBySource = newEdges;
	calcVertIndicesForEdgesBySource();
	
	// Update "main" edge list
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);
}

// Like reduce-slide, except we require lb[][] to be populated, and we allow edges to *any* vertex of the same colour as the
// endpoint, and it can begin at an ancestor of the start point.
ADD_ACTION(reduceWithSlideStrong, "reduce-slide-strong")
void reduceWithSlideStrong() {
	cerr << "Reducing edges with strong slide reduction...\n";
	
	if (m_.empty()) {
		cerr << "Must compute slide lower bounds (e.g. with calc-crappy-slide-lbs) before using reduce-slide-strong!\n";
		exit(1);
	}
	
	if (lb.empty()) {
		cerr << "You must calculate anchor lower bounds using calc-anchor-lbs before running reduce-slide-strong.\n";
		exit(1);
	}
	
	int startZ = 0;
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		int u = edges[i].u;
		int v = edges[i].v;			// v will strictly increase
		
		//cerr << "Considering edge " << edges[i] << "...\n";		//DEBUG
		
		// Find the start of v's colour (remember all vertices of the same colour now form a contiguous block)
		for (; startZ < nVerts && c[startZ] < c[v]; ++startZ) {
		}
		
		bool keep = true;
		for (int z = startZ; z < nVerts && c[z] == c[v]; ++z) {
#ifdef USE_OC
			if (z != v) {
				edges[i].oc = max(edges[i].oc, lb[z][u]);		// Update the opportunity cost LB as a nice side-effect; shame we can't get the next-best ancestor when z == v though.
			}
#endif	// USE_OC
			if (lb[z][u] + slideLb(v, z) > edges[i].w) {
				// We can delete the edge (u, v) because it's better to connect z to one of u's ancestors.
				keep = false;
				break;
			}
		}
		
		edges[j] = edges[i];
		j += keep;
	}
	
	cerr << "Deleted " << (nEdges - j) << " edges.\n";
	edges.resize(j);
	nEdges = j;
	sortEdgesAndCalcVertIndices(false);
	edgesBySource.clear();		// So that anyone who needs these edges in future knows to recalculate them
}

//vector<vector<double> > w;		// w[i][j] is the weight of the edge from vertex i to vertex j, or NEGINF if no such edge exists.
//
//void calcDpTablesForColouredReduction() {
//	assert(0, "This is massively incomplete...  I started working on the simpler bounds first.");		//TODO
//	
//	// Initialise lookup table w[][] for edge weights
//	w.resize(nVerts, vector<double>(nVerts, NEGINF));
//	for (int i = 0; i < nEdges; ++i) {
//		w[edges[i].u][edges[i].v] = edges[i].w;
//	}
//}

ADD_ACTION(reduceDominatingEdge, "reduce-domedge")
void reduceDominatingEdge() {
	cerr << "Reducing edges with the dominating edge rule...\n";
	
	if (m_.empty()) {
		cerr << "Must compute a slide lower bound (e.g. with calc-rec-slide-lbs) before running reduce-domedge!\n";
		exit(1);
	}
	
	// Vertex 0 must be the root: sanity-check that this at least could be true!
	assert(startOfEdgesTo[0] == startOfEdgesTo[1] && startOfEdgesFrom[0] != startOfEdgesFrom[1]);
	
	//TODO: We should just sort edges[] by end point then by start point at the outset -- that would save having to do this now.
	// Note: reordering the edges to a vertex will not alter any startOfEdges[] entries.
	cerr << "Sorting edges to each vertex by start vertex...\n";
	for (int i = 0; i < nVerts; ++i) {
		sort(edges.begin() + startOfEdgesTo[i], edges.begin() + startOfEdgesTo[i + 1], compare_by_startpoint_incr());
	}
	
	vector<vector<int> > vertsOfColour(nCol);
	for (int i = 0; i < nVerts; ++i) {
		vertsOfColour[c[i]].push_back(i);
	}
	
	vector<bool> keep(nEdges, true);
	int nEdgeImprovements = 0;
	
	// For every vertex x
	cerr << "Processing edges...\n";
	//for (int x = 0; x < nVerts; ++x) {
	for (int x = 1; x < nVerts; ++x) {			// The root cannot be an intermediate vertex
		//cerr << "Calculating s(x) and r(x) for dominating edges starting from " << x << "...\n";
		
		// Compute s(x) and r(x).  We don't need to store these in a table since the outermost loop is on x!
		double r = INF;
		double s = INF;
		
		// For all vertices y != x of the same colour as x
		for (int i = 0; i < vertsOfColour[c[x]].size(); ++i) {
			int y = vertsOfColour[c[x]][i];
			
			if (y != x) {
				// For every parent p of y, update s(x)
				for (int j = startOfEdgesTo[y]; j < startOfEdgesTo[y + 1]; ++j) {
					//int p = edges[j].u;		// Not actually used!
					//s = min(s, m[y][x] - edges[j].w);		//HACK: m[y][x] doesn't vary in this loop, so this could be better optimised.
					s = min(s, slideLb(y, x) - edges[j].w);		//HACK: m[y][x] doesn't vary in this loop, so this could be better optimised.
				}
				
				//TODO: FIX THIS BUG!
				// We're currently calculating r wrong: we're assuming that p must be restricted to parents of c(x)-coloured vertices
				// that also have an edge to x, but of course p is not so restricted -- it could well be that the parent p of the
				// c(x)-coloured vertex y in the current solution (tree) has no edge to x!  In that case, we don't have the
				// option of attaching x to this parent, as the calculation of condition C currently assumes.
				// To actually get a valid non-NEGINF value for r, it must be that *every* parent of *every* c(x)-coloured vertex has
				// an edge to x: in this case we can safely choose the minimum of all these.  But I doubt this will ever be greater
				// than s + w(ux).  (Well, maybe it will, e.g. if there are only 1 or 2 other vertices of colour c(x).)
				//assert(!"This is a BUG!  We are calculating r wrong: see the comment lines above here.");
				// For every parent p of y that is also a parent of x, update r(x)
				// Use a linear-time merge to find parents shared by x and y (remember edge start points for each end point have been sorted by start point)
				int j = startOfEdgesTo[x];
				int k = startOfEdgesTo[y];
				while (j < startOfEdgesTo[x + 1] && k < startOfEdgesTo[y + 1]) {
					if (edges[j].u < edges[k].u) {
						++j;
					} else if (edges[j].u > edges[k].u) {
						r = NEGINF;		// We only arrive here if there is some parent p of y that does not have an edge to x, meaning we cannot attach x to p!
						++k;
					} else {
						// This vertex is a parent of both x and y.
						//r = min(r, edges[j].w - edges[k].w + m[y][x]);
						r = min(r, edges[j].w - edges[k].w + slideLb(y, x));
						++j;
						++k;
					}
				}
			}
		}
		
		// This is all now irrelevant, since we set r to NEGINF if there is *any* other vertex of colour c(x) that has a parent p
		// without an edge to x -- which will happen a lot.
		//// Tricky: if r == INF, then it means that x does not share a parent with any of its same-coloured vertices.
		//// In this case, we want to "disable" the option of creating the edge px (since no such edge exists).
		//// OTOH, if s == INF, then it means that none of its same-coloured vertices has any in-edges (probably
		//// because x is the only vertex of its colour), meaning it's always safe to add x to a solution that doesn't
		//// already contain it, so it's OK to satisfy C by letting it remain at INF.
		//if (r == INF) {
		//	r = NEGINF;
		//}
		
		// For all pairs (x, v) (not necessarily an edge in the graph), compute the best
		// vertex z such that c(z) == c(v), (x, z) is an edge and maximises w(x, z) + m(v, z).
		// Actually we use "reaching" to calculate this more easily: for every edge (x, z), compare with
		// all bestZ[x][v] for v the same colour as z and update where necessary
		vector<double> bestZ(nVerts, NEGINF);
		for (int i = startOfEdgesFrom[x]; i < startOfEdgesFrom[x + 1]; ++i) {
			int z = edgesBySource[i].v;
			
			for (int k = 0; k < vertsOfColour[c[z]].size(); ++k) {
				int v = vertsOfColour[c[z]][k];
				//bestZ[v] = max(bestZ[v], edgesBySource[i].w + m[v][z]);
				bestZ[v] = max(bestZ[v], edgesBySource[i].w + slideLb(v, z));
			}
		}
		
		// Calculate the cheapest way to insert x using a single edge, given that u is in the tree.
		// I.e. take the minimum over all paths from the root to u of the best edge on this path to x.
		// Because of our processing order, we can get away with a 1D array :)
		//HACK: May actually be faster to calculate using memoisation (instead of tabulation) if in-degrees are low...
		vector<double> edgeToXFrom(x, NEGINF);				// A 1D slice of the adjacency matrix, giving the weight of an edge from any given vertex < x to x, or NEGINF if there is no such edge.
		for (int i = startOfEdgesTo[x]; i < startOfEdgesTo[x + 1]; ++i) {
			edgeToXFrom[edges[i].u] = edges[i].w;
		}
		
		vector<double> bestEdgeToXGiven(x, NEGINF);		// Don't need the full nVerts indices since any parent of x must be numbered lower than x
		bestEdgeToXGiven[0] = edgeToXFrom[0];
		for (int u = 0; u < x; ++u) {
			if (u > 0) {
				double rest = INF;
				for (int j = startOfEdgesTo[u]; j < startOfEdgesTo[u + 1]; ++j) {
					rest = min(rest, bestEdgeToXGiven[edges[j].u]);
				}
				
				bestEdgeToXGiven[u] = max(edgeToXFrom[u], rest);
			}
			
			if (bestEdgeToXGiven[u] > edgeToXFrom[u]) {
				++nEdgeImprovements;		// Just for my interest
			}
			
			// We used to end this loop, and have a separate loop start below that looped u through all of x's parents,
			// but now we let u be any possible vertex above x.
			if (bestEdgeToXGiven[u] == NEGINF) {
				continue;
			}
			
			// For every child v of u potentially reachable from x
			for (int j = startOfEdgesFrom[u]; j < startOfEdgesFrom[u + 1]; ++j) {
				int v = edgesBySource[j].v;
				
				if (v < x || c[v] == c[x]) {
					continue;		// v cannot possibly be reachable from x if it is lower-numbered or the same colour
				}
				
				if (keep[j]) {
					// We haven't deleted (u, v) yet.
					
					// Does it satisfy criteria A, B and C?
					//if (bestZ[x][v] - edgesBySource[j].w > 0 &&								// A
					//	bestZ[x][v] - edgesBySource[j].w + bestEdgeToXGiven[u] > 0 &&					// B
					//	bestZ[x][v] - edgesBySource[j].w + max(bestEdgeToXGiven[u] + s, r) > 0)		// C
					if (bestZ[v] - edgesBySource[j].w > 0 &&								// A
						bestZ[v] - edgesBySource[j].w + bestEdgeToXGiven[u] > 0 &&					// B
						bestZ[v] - edgesBySource[j].w + max(bestEdgeToXGiven[u] + s, r) > 0)		// C
					{
						keep[j] = false;		// We can safely delete (u, v).
					}
				}
			}
		}
	}
	
	// Actually delete the edges
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		edgesBySource[j] = edgesBySource[i];
		j += keep[i];		// Only advance if we want to keep this edge
	}
	
	cerr << "Deleted " << (nEdges - j) << " edges.\n";
	cerr << nEdgeImprovements << " edges were either improved or considered when they were not parents of x.\n";
	nEdges = j;
	edgesBySource.erase(edgesBySource.begin() + nEdges, edgesBySource.end());
	calcVertIndicesForEdgesBySource();
	
	// Update "main" edge list
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);
}

ADD_ACTION(calcImpliedAncestors, "calc-implied-anc")
void calcImpliedAncestors() {
	cerr << "Calculating implied ancestor vertices...\n";
	impliedAnc.assign(nVerts, vector<int>());
	impliedAnc[0] = vector<int>(1, 0);			// The root is its own ancestor, and the only one.
	
	for (int i = 1; i < nVerts; ++i) {
		if (startOfEdgesTo[i] != startOfEdgesTo[i + 1]) {
			impliedAnc[i] = impliedAnc[edges[startOfEdgesTo[i]].u];
		}
		
		for (int j = startOfEdgesTo[i] + 1; j < startOfEdgesTo[i + 1]; ++j) {
			int u = edges[j].u;
			
			// Intersect with u's list of implied ancestor vertices
			int from = 0;
			int to = 0;
			for (int ji = 0; ji < impliedAnc[u].size() && from < impliedAnc[i].size(); ++ji) {
				while (from < impliedAnc[i].size() && impliedAnc[i][from] < impliedAnc[u][ji]) {
					++from;
				}
				
				if (from == impliedAnc[i].size()) {
					break;
				}
				
				if (impliedAnc[i][from] == impliedAnc[u][ji]) {
					// Keep this one
					impliedAnc[i][to++] = impliedAnc[i][from++];
				}
			}
			
			impliedAnc[i].resize(to);
		}
		
		impliedAnc[i].push_back(i);		// A vertex is considered to be an ancestor of itself
		
		if (impliedAnc[i].size() > 2) {
			cerr << "Vertex " << i << " has " << impliedAnc[i].size() << " implied ancestors!\n";
		}
	}
}

ADD_ACTION(reduceDominatingPath, "reduce-dompath")
void reduceDominatingPath() {
	cerr << "Reducing edges with the dominating path rule...\n";
	
	if (m_.empty()) {
		cerr << "Must compute a slide lower bound (e.g. with calc-rec-slide-lbs) before running reduce-dompath!\n";
		exit(1);
	}
	
	if (impliedAnc.empty()) {
		cerr << "Must compute implied ancestor vertices with calc-implied-anc before running reduce-dompath!\n";
		exit(1);
	}
	
	vector<vector<int> > vertsOfColour(nCol);
	for (int i = 0; i < nVerts; ++i) {
		vertsOfColour[c[i]].push_back(i);
	}
	
	//vector<vector<double> > ra;			// DON'T ACTUALLY NEED TO STORE THIS.  ra[i][j] will contain a LB on the score change needed to connect j to the solution, assuming i is already in it.
	vector<vector<double> > rb;			// rb[i][j] will contain a LB on the score change needed to connect i to j via some path, assuming no c(j)-coloured vertex is in the solution already.
	//vector<vector<double> > rc;			// rc[i][j] will contain a LB on the score change needed to connect j to the solution, assuming i is already in it, and that some c(j)-coloured vertex (possibly j itself) is in the solution already.
	vector<vector<double> > rcp;			// rcp[i][j] will contain a LB on the score change needed to connect j to the solution, assuming i is already in it, and that some *other* c(j)-coloured vertex is in the solution already.
	//vector<vector<double> > re;			// DON'T ACTUALLY NEED TO STORE THIS.  re[i][j] will contain a LB on the score change needed to connect i to j via some path, assuming that some c(j)-coloured vertex (possibly j itself) is in the solution already.
	//vector<vector<double> > rf;			// DON'T ACTUALLY NEED TO STORE THIS.  rf[i][j] will contain a LB on the score change needed to connect i to the solution via p, assuming that some c(j)-coloured vertex (possibly j itself) is in the solution already, and has a parent p.
	//vector<double> rg;						// rg[i] will contain a LB on the score change needed to make it possible to insert i into the solution, assuming that some c(i)-coloured vertex (possibly i itself) is already in the solution and therefore the edge from it to its parent must be deleted, and its subtree either deleted or "slid" across to i.
	vector<double> rgp;						// rg[i] will contain a LB on the score change needed to make it possible to insert i into the solution, assuming that some *other* c(i)-coloured vertex is already in the solution and therefore the edge from it to its parent must be deleted, and its subtree either deleted or "slid" across to i.
	vector<double> rh;						// rh[i] will contain the heaviest in-edge for i, or NEGINF if none exists.
	
	cerr << "Calculating rb[][]...\n";
	rb.assign(nVerts, vector<double>(nVerts, NEGINF));
	for (int i = 0; i < nVerts; ++i) {
		rb[i][i] = 0.0;
		for (int j = i + 1; j < nVerts; ++j) {
			if (c[j] != c[i]) {		// We can ignore path endpoints of the same colour as the start point, since they cannot exist
				for (int ki = startOfEdgesTo[j]; ki < startOfEdgesTo[j + 1]; ++ki) {
					int k = edges[ki].u;
					rb[i][j] = max(rb[i][j], min(rb[i][k] + edges[ki].w, edges[ki].w));
				}
			}
		}
	}
	
	cerr << "Calculating rh[]...\n";
	rh.assign(nVerts, NEGINF);
	for (int i = 0; i < nVerts; ++i) {
		for (int ji = startOfEdgesTo[i]; ji < startOfEdgesTo[i + 1]; ++ji) {
			rh[i] = max(rh[i], edges[ji].w);
		}
	}
	
	cerr << "Calculating rgp[]...\n";
	rgp.assign(nVerts, INF);
	for (int i = 0; i < nVerts; ++i) {
		for (int ji = 0; ji < vertsOfColour[c[i]].size(); ++ji) {
			int y = vertsOfColour[c[i]][ji];
			if (y != i) {
				//rgp[i] = min(rgp[i], m[y][i] - rh[y]);
				rgp[i] = min(rgp[i], slideLb(y, i) - rh[y]);
			}
		}
	}
	
	cerr << "Filtering edges...\n";
	//HACK: The following can probably be sped up somehow...
	// NB: We modify the edgesBySource[] array in-place, but leave edges[], which we need to read from meanwhile, intact
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		edgesBySource[j] = edgesBySource[i];
		
		int u = edgesBySource[i].u;
		int v = edgesBySource[i].v;
		double best = NEGINF;
		for (int ki = 0; ki < vertsOfColour[c[v]].size(); ++ki) {
			int k = vertsOfColour[c[v]][ki];
			for (int pi = startOfEdgesTo[k]; pi < startOfEdgesTo[k + 1]; ++pi) {
				int p = edges[pi].u;
				
				assert(!"This implementation of the dominating path reduction is flawed in at least 2 ways -- see comments on rev 27159.");
				if (p != i) {			//BUG: This HAS to be wrong!  p is a vertex number, while i is an index into the list of edges by source!
					//best = max(best, m[v][k] + edges[pi].w + /* ra[i][p] == */ min(rb[i][p], /* rc[i][p] == */ max(/* re[i][p] == */ min(0.0, rb[i][p]), /* rf[p] == */ NEGINF) + rg[p]));
					//best = max(best, m[v][k] + edges[pi].w + /* ra[i][p] == */ min(rb[u][p], /* rc[i][p] == */ max(/* re[i][p] == */ min(0.0, rb[u][p]), /* rf[p] == */ NEGINF) + rg[p]));
					//best = max(best, m[v][k] + edges[pi].w + /* ra[i][p] == */ min(rb[0][p], /* rc[i][p] == */ max(/* re[i][p] == */ min(0.0, rb[0][p]), /* rf[p] == */ NEGINF) + rg[p]));	// Try root too!
					//best = max(best, m[v][k] + edges[pi].w + /* rap[i][p] == */ min(0.0, min(rb[u][p], /* rcp[i][p] == */ max(rb[u][p], /* rfp[p] == */ NEGINF) + rgp[p])));		// Assumes p is reachable from the root
					//best = max(best, m[v][k] + edges[pi].w + /* rap[i][p] == */ min(0.0, min(rb[0][p], /* rcp[i][p] == */ max(rb[0][p], /* rfp[p] == */ NEGINF) + rgp[p])));		// Assumes p is reachable from the root.  Try root too!
					for (int a = 0; a < impliedAnc[u].size(); ++a) {
						//best = max(best, m[v][k] + edges[pi].w + /* rap[i][p] == */ min(0.0, min(rb[impliedAnc[u][a]][p], /* rcp[i][p] == */ max(rb[impliedAnc[u][a]][p], /* rfp[p] == */ NEGINF) + rgp[p])));		// Assumes p is reachable from the root
						//double x = m[v][k] + edges[pi].w + /* rap[i][p] == */ min(0.0, min(rb[impliedAnc[u][a]][p], /* rcp[i][p] == */ max(rb[impliedAnc[u][a]][p], /* rfp[p] == */ NEGINF) + rgp[p]));		// Assumes p is reachable from the root
						double x = slideLb(v, k) + edges[pi].w + /* rap[i][p] == */ min(0.0, min(rb[impliedAnc[u][a]][p], /* rcp[i][p] == */ max(rb[impliedAnc[u][a]][p], /* rfp[p] == */ NEGINF) + rgp[p]));		// Assumes p is reachable from the root
						best = max(best, x);		// Assumes p is reachable from the root
						if (edges[pi].u == 1499 && edges[pi].v == 1515) {		//DEBUG
							cerr << "Trying ancestor " << impliedAnc[u][a] << " of vertex " << u << " gives path of score " << x << "; new best is now " << best << ".\n";
						}
					}
				}
			}
		}
		
		edge DEBUGedge(edgesBySource[i]);		//DEBUG: This exists purely because of limitations in VS++ 2008's expression evaluator...
		if (best <= edgesBySource[i].w) {
			++j;		// This edge lives...
		}
	}
	
	cerr << "Deleted " << (nEdges - j) << " edges.\n";
	nEdges = j;
	edgesBySource.erase(edgesBySource.begin() + nEdges, edgesBySource.end());
	calcVertIndicesForEdgesBySource();
	
	// Update "main" edge list
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);
}

ADD_ACTION(reduceCombineColours, "reduce-combcol")
void reduceCombineColours() {
	cerr << "Combining colours...\n";
	
	vector<vector<int> > vertsOfColour(nCol);
	for (int i = 0; i < nVerts; ++i) {
		vertsOfColour[c[i]].push_back(i);
	}
	
	vector<double> maxInEdgeForColour(nCol, NEGINF);
	vector<int> soleOutColourForColour(nCol, -1);
	for (int i = 0; i < nEdges; ++i) {
		int vCol = c[edges[i].v];
		maxInEdgeForColour[vCol] = max(maxInEdgeForColour[vCol], edges[i].w);
		
		int uCol = c[edges[i].u];
		if (soleOutColourForColour[uCol] == -1) {
			soleOutColourForColour[uCol] = vCol;
		} else if (soleOutColourForColour[uCol] != vCol) {
			soleOutColourForColour[uCol] = -2;		// -2 means "more than one distinct colour"
		}
	}
	
	vector<vector<int> > sourceColours(nCol);
	int nColoursWithSoleOutColour = 0;
	int nColoursWithOnlyNegInEdges = 0;
	for (int i = 0; i < nCol; ++i) {
		if (maxInEdgeForColour[i] < 0 && soleOutColourForColour[i] >= 0) {
			cerr << "Colour " << i << ", with " << vertsOfColour[i].size() << " vertices, has only negative in-edges and all out-edges go to colour " << soleOutColourForColour[i] << ".\n";
			sourceColours[soleOutColourForColour[i]].push_back(i);
		}
		
		// For info only
		if (maxInEdgeForColour[i] < 0) {
			++nColoursWithOnlyNegInEdges;
		}
		
		if (soleOutColourForColour[i] >= 0) {
			++nColoursWithSoleOutColour;
		}
	}
	
	cerr << nColoursWithSoleOutColour << " colours have out-edges to just a single colour; " << nColoursWithOnlyNegInEdges << " colours have only negative-weight in-edges.\n";
	
	int nVertsRecoloured = 0;
	int nColoursChanged = 0;
	for (int i = 0; i < nCol; ++i) {
		for (int j = 1; j < sourceColours[i].size(); ++j) {
			cerr << "Changing all vertices of colour " << sourceColours[i][j] << " to colour " << sourceColours[i][0] << ".\n";
			for (int k = 0; k < vertsOfColour[sourceColours[i][j]].size(); ++k) {
				c[vertsOfColour[sourceColours[i][j]][k]] = sourceColours[i][0];
			}
			
			nVertsRecoloured += vertsOfColour[sourceColours[i][j]].size();
			++nColoursChanged;
		}
	}
	
	cerr << nVertsRecoloured << " vertices across " << nColoursChanged << " colours recoloured.\n";
}

// Frustrating that the C++ Standard Library doesn't already have this!
template <typename It, typename F>
//It in_place_intersect(It begin1, It end1, It begin2, It end2, F comp = typename less<It>::value_type()) {
It in_place_intersect(It begin1, It end1, It begin2, It end2, F comp) {
	It out(begin1);
	
	while (begin1 != end1 && begin2 != end2) {
		if (comp(*begin1, *begin2)) {
			++begin1;
		} else if (comp(*begin2, *begin1)) {
			++begin2;
		} else {
			// Must be equal
			*out = *begin1;
			++begin1;
			++begin2;
			++out;
		}
	}
	
	return out;
}

template <typename It>
It in_place_intersect(It begin1, It end1, It begin2, It end2) {
	return in_place_intersect(begin1, end1, begin2, end2, less<typename It::value_type>());
}

// Likewise, this should be in the C++ Standard Library.
template <typename It, typename F>
bool has_nonempty_intersection(It begin1, It end1, It begin2, It end2, F comp) {
	while (begin1 != end1 && begin2 != end2) {
		if (comp(*begin1, *begin2)) {
			++begin1;
		} else if (comp(*begin2, *begin1)) {
			++begin2;
		} else {
			// Must be equal
			return true;
		}
	}
	
	return false;
}

template <typename It>
bool has_nonempty_intersection(It begin1, It end1, It begin2, It end2) {
	return has_nonempty_intersection(begin1, end1, begin2, end2, less<typename It::value_type>());
}

// Compares indices into edgesBySource[].
struct compare_by_edge_weight_plus_ub_desc {
	bool operator()(int a, int b) const {
		return edgesBySource[b].w + ub[edgesBySource[b].v] < edgesBySource[a].w + ub[edgesBySource[a].v];
	}
};

//TODO: These get out of date as soon as any reduction after calcImpliedEdges() deletes an edge!
vector<vector<edge> > impliedEdges;			// impliedEdges[i][j] is the jth descendant edge implied by edge edgesBySource[i].  Edges are not stored in any particular order.
vector<vector<int> > coloursImpliedByEdge;	// coloursImpliedByEdge[i] is a list of colours (in ascending order) that are implied by edge edgesBySource[i].

ADD_ACTION(calcImpliedEdges, "calc-implied-edges")
void calcImpliedEdges() {
	if (ub.empty()) {
		cerr << "Must compute vertex upper bounds (using calc-child-vertex-ubs and/or calc-col-forest-vertex-ubs) before you can apply calc-implied-edges!\n";
		exit(1);
	}
	
	cerr << "Calculating implied descendant edges and colours...\n";
	cerr << "//TODO: Lots of debugging cruft in here...\n";		//HACK
	
	impliedEdges.assign(nEdges, vector<edge>());		//HACK: Yeah, this could get pretty big...
	coloursImpliedByEdge.assign(nEdges, vector<int>());
	
	// Need to do this to get the single colour "implied" by each leaf vertex!
	for (int i = 0; i < nEdges; ++i) {
		coloursImpliedByEdge[i].push_back(c[edgesBySource[i].v]);
	}
	
	int nOnlyChild = 0;
	int nEdgesDeletedDueToKnownForcedEdges = 0;
	int nEdgesDeletedDueToUnknownForcedEdges = 0;
	int nEdgesDeletedDueToBetterUBs = 0;
	int nEdgesDeletedDueToPigeonholePrinciple = 0;
	int nSiblingsIgnored = 0;
	vector<bool> keep(nEdges, true);		// We will delete edges on a second pass.
	
	// Process edges in modified "postorder" -- i.e. when an edge (u, v) is processed, all edges (v, w) will have already been processed.
	for (int i = nEdges - 1; i >= 0; --i) {
#ifdef USE_OC
		double debt = edgesBySource[i].w - edgesBySource[i].oc;			// If we cannot find a subtree having at least this cost then edgesBySource[i] cannot be part of any optimal solution.
#else	// not USE_OC
		double debt = edgesBySource[i].w;			// If we cannot find a subtree having at least this cost then edgesBySource[i] cannot be part of any optimal solution.
#endif	// not USE_OC
		//if (edgesBySource[i].oc > 0) {
		//	cerr << "calc-implied-edges: " << edgesBySource[i] << " has debt " << debt << ".\n";		//DEBUG
		//}
		if (debt < 0) {
			// This edge is negative.  It may therefore imply some child edges.
			int v = edgesBySource[i].v;
			double allEdgesUB = 0.0;
			for (int j = startOfEdgesFrom[v]; j < startOfEdgesFrom[v + 1]; ++j) {
				assert(edgesBySource[j].u == v);		// Sanity check
				if (keep[j]) {			// Only process child edges we haven't already marked for deletion
					double x = edgesBySource[j].w + ub[edgesBySource[j].v];
					// You need to run reduce-vub first, to ensure there are no such edges.
					// IMPORTANT: It *is* OK to run some other reduction in between, since that will only (possibly) remove edges, not change this property of any surviving edges.
					assert(x >= 0.0);
					allEdgesUB += x;
				}
			}
			
			vector<int> coloursImpliedBySpecificEdge;
			coloursImpliedBySpecificEdge.push_back(c[v]);	// Note that colours are NOT ordered!
			vector<int> DEBUGcoloursImpliedBySpecificEdgeWithoutOpCost;		//DEBUG
			DEBUGcoloursImpliedBySpecificEdgeWithoutOpCost.push_back(c[v]);	//DEBUG: Note that colours are NOT ordered!
			vector<int> remainingEdgeIndices;		// Will eventually contain indices in edgesBySource[] of all child edges that are not forced in by this edge.
			double forcedInEdgesUB = 0.0;		// Will eventually be the maximum that the edges that must be present could contribute
			for (int j = startOfEdgesFrom[v]; j < startOfEdgesFrom[v + 1]; ++j) {
				if (keep[j]) {			// Only process child edges we haven't already marked for deletion
					double uBwithoutThisEdge = allEdgesUB - (edgesBySource[j].w + ub[edgesBySource[j].v]);
					assert(uBwithoutThisEdge >= 0.0);
					
					if (debt + uBwithoutThisEdge < 0.0) {
						// Even adding all other child edges with their UBs doesn't recoup the cost paid by the edge to v, so this (jth) child edge is forced in.
						// NOTE: If edge i winds up being kept, then it must be that impliedEdges[i] does not contain any duplicate
						// edges -- since if it did, this would imply duplicate colours, implying that edge i would have been deleted.
						// (And if edge i isn't kept, we don't care what goes in impliedEdges[i] because we'll delete it anyway.)
						// Note that impliedEdges[i] is not sorted.
						impliedEdges[i].push_back(edgesBySource[j]);
						impliedEdges[i].insert(impliedEdges[i].end(), impliedEdges[j].begin(), impliedEdges[j].end());		//HACK: Burns more space than any implicit representation relying on recursion, but easier to work with
						coloursImpliedBySpecificEdge.insert(coloursImpliedBySpecificEdge.end(), coloursImpliedByEdge[j].begin(), coloursImpliedByEdge[j].end());
						if (startOfEdgesFrom[v] + 1 == startOfEdgesFrom[v + 1]) {
							++nOnlyChild;		//HACK: Works because the loop running around this will only run once in this case :)
						}
						forcedInEdgesUB += edgesBySource[j].w + ub[edgesBySource[j].v];
					} else {
						remainingEdgeIndices.push_back(j);		// These will be separately considered later.
					}
					
					//DEBUG
					if (edgesBySource[i].w + uBwithoutThisEdge < 0.0) {
						DEBUGcoloursImpliedBySpecificEdgeWithoutOpCost.insert(DEBUGcoloursImpliedBySpecificEdgeWithoutOpCost.end(), coloursImpliedByEdge[j].begin(), coloursImpliedByEdge[j].end());
					}
				}
			}
			
			sort(coloursImpliedBySpecificEdge.begin(), coloursImpliedBySpecificEdge.end());
			vector<int>::iterator newEnd = unique(coloursImpliedBySpecificEdge.begin(), coloursImpliedBySpecificEdge.end());
			if (newEnd != coloursImpliedBySpecificEdge.end()) {
				// More than one forced-in child edge implies the same colour: this edge cannot be present.
				if (!watchedEdges.empty() && watchedEdges[0].u == edgesBySource[i].u && watchedEdges[0].v == edgesBySource[i].v) {
					//DEBUG
					cerr << "Edge #" << i << ", " << edgesBySource[i] << " deleted due to a colour being implied by 2 or more known forced-in child edges!\n";
					exit(1);
				}
				
				//DEBUG
				sort(DEBUGcoloursImpliedBySpecificEdgeWithoutOpCost.begin(), DEBUGcoloursImpliedBySpecificEdgeWithoutOpCost.end());
				vector<int>::iterator DEBUGnewEnd = unique(DEBUGcoloursImpliedBySpecificEdgeWithoutOpCost.begin(), DEBUGcoloursImpliedBySpecificEdgeWithoutOpCost.end());
				if (DEBUGnewEnd == DEBUGcoloursImpliedBySpecificEdgeWithoutOpCost.end()) {
					cerr << "Edge #" << i << ", " << edgesBySource[i] << ": deleted due to a colour being implied by 2 or more known forced-in child edges only due to opcost.\n";
				}
				
				keep[i] = false;
				++nEdgesDeletedDueToKnownForcedEdges;
			} else {
				// We can still consider colours that are implied by *combinations* of the remaining children.
				
				// First, we can ignore any child edge that implies any colour already implied by a forced-in edge.
				int k = 0;
				for (int ji = 0; ji < remainingEdgeIndices.size(); ++ji) {
					remainingEdgeIndices[k] = remainingEdgeIndices[ji];
					
					int j = remainingEdgeIndices[ji];
					k += !has_nonempty_intersection(coloursImpliedBySpecificEdge.begin(), coloursImpliedBySpecificEdge.end(), coloursImpliedByEdge[j].begin(), coloursImpliedByEdge[j].end());		// I am so cool
				}
				
				nSiblingsIgnored += remainingEdgeIndices.size() - k;
				remainingEdgeIndices.resize(k);
				
				vector<int> colours;
				if (debt + forcedInEdgesUB < 0.0) {
					// At least 1 of the remaining edges must be forced in, though we don't know which.
					if (remainingEdgeIndices.empty()) {
						// This can happen when vertex UBs are loose (i.e. when one or more reductions have deleted edges since the last time reduce-XXX-vub was run till convergence).
						// E.g. ub[v] could be high due to child edges that used to be present, even though the only remaining child edges aren't enough to keep (u, v) in the graph.
						// It can also happen if child edges were ignored due to implying colours already implied by forced-in edges.
						// More than one forced-in child edge implies the same colour: this edge cannot be present.
						if (!watchedEdges.empty() && watchedEdges[0].u == edgesBySource[i].u && watchedEdges[0].v == edgesBySource[i].v) {
							//DEBUG
							cerr << "Edge #" << i << ", " << edgesBySource[i] << " deleted due to better bounds!  edgesBySource[i].w=" << edgesBySource[i].w << ", debt=" << debt << ", forcedInEdgesUB=" << forcedInEdgesUB << ", debt+forcedInEdgesUB=" << (debt + forcedInEdgesUB) << " < 0 and there are no remaining child edges to add!\n";
							exit(1);
						}
						
						//DEBUG
						if (edgesBySource[i].w + forcedInEdgesUB >= 0.0) {
							cerr << "Edge #" << i << ", " << edgesBySource[i] << ": deleted due to better bounds only due to opcost.\n";
						}
						
						keep[i] = false;
						++nEdgesDeletedDueToBetterUBs;
					} else {
						// How many of the remaining edges must be included?  We can figure this out by trying to add in the best k edges,
						// for successive values of k.  This won't tell us exactly which edges must be present, but gives us a LB
						// on the *number* that must be.
						// This code correctly handles the case when none of the remaining edges are forced in.
						sort(remainingEdgeIndices.begin(), remainingEdgeIndices.end(), compare_by_edge_weight_plus_ub_desc());
						
						double total = 0.0;
						int nRemainingEdgesForcedIn = 0;
						int DEBUGnRemainingEdgesForcedInWithoutOpCost = 0;
						for (; nRemainingEdgesForcedIn < remainingEdgeIndices.size() && debt + forcedInEdgesUB + total < 0.0; ++nRemainingEdgesForcedIn) {
							//DEBUG
							if (edgesBySource[i].w + forcedInEdgesUB + total < 0.0) {
								++DEBUGnRemainingEdgesForcedInWithoutOpCost;
							}
							total += edgesBySource[remainingEdgeIndices[nRemainingEdgesForcedIn]].w + ub[edgesBySource[remainingEdgeIndices[nRemainingEdgesForcedIn]].v];
						}
						
						assert(nRemainingEdgesForcedIn > 0);
						
						//TODO: Apply the "full" pigeonhole principle, i.e. where we count minimum numbers of total edges and the union of all reachable colours.
						//if (keep[i]) {
							// Finally try using the pigeonhole principle to detect when there are more implied child edges than colours spanning their endpoints.
							// If the previous test was restricted to just u's children (and not colours they imply recursively) then this test will trigger whenever the it does,
							// and sometimes when it doesn't (i.e. it is strictly stronger), but it's still useful to perform the previous test because
							// it enables us to detect *which* colours have a single edge forced in, which can strengthen the tests for ancestors of u.
							vector<bool> coloursSeen(nCol, false);
							int nDistinctColours = 0;
							for (int ji = 0; ji < remainingEdgeIndices.size(); ++ji) {
								int j = remainingEdgeIndices[ji];
								nDistinctColours += !coloursSeen[c[edgesBySource[j].v]];
								coloursSeen[c[edgesBySource[j].v]] = true;
							}
							
							//cerr << "For " << edgesBySource[i] << ": nRemainingEdges=" << remainingEdgeIndices.size() << ", nDistinctColours=" << nDistinctColours << ", nRemainingEdgesForcedIn=" << nRemainingEdgesForcedIn << ".\n";	//DEBUG
							//DEBUG
							if (DEBUGnRemainingEdgesForcedInWithoutOpCost < nRemainingEdgesForcedIn) {
								cerr << "Edge #" << i << ", " << edgesBySource[i] << ": nRemainingEdges=" << remainingEdgeIndices.size() << ", nDistinctColours=" << nDistinctColours << ", nRemainingEdgesForcedIn=" << nRemainingEdgesForcedIn << ", DEBUGnRemainingEdgesForcedInWithoutOpCost=" << DEBUGnRemainingEdgesForcedInWithoutOpCost << "\n";
							}
							
							if (nDistinctColours < nRemainingEdgesForcedIn) {
								keep[i] = false;
								++nEdgesDeletedDueToPigeonholePrinciple;
							}
//						}
						
						// Which colours are missing in at most (nRemainingEdgesForcedIn-1) child edges?  All such colours must be implied by any solution having
						// at least nRemainingEdgesForcedIn child edges.
						//map<int, int> colourFreqs;						// colourFreqs[i] is the number of child edges that imply colour i.
						vector<int> colourFreqs(nCol, 0);
						// Adjust colour counts so that we correctly detect when there is a colour that is implied by unknown child edges and also by a known edge
						for (int k = 0; k < coloursImpliedBySpecificEdge.size(); ++k) {
							colourFreqs[coloursImpliedBySpecificEdge[k]] = 1;
						}
						
						//for (int ji = 0; ji < remainingEdgeIndices.size(); ++ji) {
						for (int ji = 0; ji < remainingEdgeIndices.size() && keep[i]; ++ji) {
							int j = remainingEdgeIndices[ji];
							
							////HACK: Probably not the fastest way to do it, but easy and not too bad.
							//for (int k = 0; k < coloursImpliedByEdge[j].size(); ++k) {
							//	++colourFreqs[coloursImpliedByEdge[j][k]];
							//}
							for (int k = 0; k < coloursImpliedByEdge[j].size(); ++k) {
								//if (++colourFreqs[coloursImpliedByEdge[j][k]] == remainingEdgeIndices.size() - nRemainingEdgesForcedIn + 1) {
								++colourFreqs[coloursImpliedByEdge[j][k]];
								if (colourFreqs[coloursImpliedByEdge[j][k]] == remainingEdgeIndices.size() - nRemainingEdgesForcedIn + 1) {
									// We have a colour that must appear (at least) once.
									colours.push_back(coloursImpliedByEdge[j][k]);
								} else if (colourFreqs[coloursImpliedByEdge[j][k]] == remainingEdgeIndices.size() - nRemainingEdgesForcedIn + 2) {
									// Bingo -- we have a colour that appears often enough to force a conflict.
									if (!watchedEdges.empty() && watchedEdges[0].u == edgesBySource[i].u && watchedEdges[0].v == edgesBySource[i].v) {
										//DEBUG
										cerr << "Edge #" << i << ", " << edgesBySource[i] << " deleted due to a colour being implied by 2 or more unknown forced-in child edges!\n";
										exit(1);
									}
									keep[i] = false;			// Also breaks out of outer loop
									++nEdgesDeletedDueToUnknownForcedEdges;
									break;
								}
							}
						}
						
						//// Now gather up any colours that occur the necessary number of times.
						//for (map<int, int>::const_iterator j(colourFreqs.begin()); j != colourFreqs.end(); ++j) {
						//	if ((*j).second > remainingEdgeIndices.size() - nRemainingEdgesForcedIn) {
						//		colours.push_back((*j).first);
						//	}
						//}
						//
						//if (!colours.empty()) {
						//	keep[i] = false;
						//	++nEdgesDeletedDueToUnknownForcedEdges;
						//}
					}
				}
				
				if (keep[i]) {
					// We only bother keeping implied colour lists for edges we are not about to delete.
					coloursImpliedByEdge[i].clear();		// Since we now populate this at the start with the endpoint's colour.
					// The colours implied by known forced-in edges and unknown forced-in edges must be disjoint!
					//set_union(colours.begin(), colours.end(), coloursImpliedBySpecificEdge.begin(), coloursImpliedBySpecificEdge.end(), back_inserter(coloursImpliedByEdge[i]));
					//assert(colours.empty());		//DEBUG
					sort(colours.begin(), colours.end());
					merge(colours.begin(), colours.end(), coloursImpliedBySpecificEdge.begin(), coloursImpliedBySpecificEdge.end(), back_inserter(coloursImpliedByEdge[i]));
				}
			}
		}
	}
	
	// Show some interesting stats
	vector<int> nEdgesImplied(nEdges, 0);
	for (int i = 0; i < nEdges; ++i) {
		nEdgesImplied[i] = impliedEdges[i].size();
	}
	
	sort(nEdgesImplied.begin(), nEdgesImplied.end(), greater<int>());
	
	int lastChange = -1;
	int nEdgesWithAtLeastOneForcedChild = 0;
	for (int i = 0; i < nEdges; ++i) {
		if (i == nEdges - 1 || nEdgesImplied[i] != nEdgesImplied[i + 1]) {
			cerr << (i - lastChange) << " edges have " << nEdgesImplied[i] << " implied edges.\n";
			
			if (nEdgesImplied[i] == 0) {
				nEdgesWithAtLeastOneForcedChild = nEdges - (i - lastChange);		//TODO: Is this correct??
			}
			
			lastChange = i;
		}
	}
	
	cerr << nOnlyChild << " of the " << nEdgesWithAtLeastOneForcedChild << " edges with at least one implied child edge had only one child edge.\n";
	cerr << "We can delete " << (nEdgesDeletedDueToKnownForcedEdges + nEdgesDeletedDueToUnknownForcedEdges + nEdgesDeletedDueToBetterUBs + nEdgesDeletedDueToPigeonholePrinciple) << " edges due to implied colour conflicts (" << nEdgesDeletedDueToKnownForcedEdges << " from known forced edges, " << nEdgesDeletedDueToUnknownForcedEdges << " from unknown, " << nEdgesDeletedDueToBetterUBs << " from (implied) better UBs, " << nEdgesDeletedDueToPigeonholePrinciple << " from the pigeonhole principle).\n";
	cerr << "We were able to ignore " << nSiblingsIgnored << " edges due to forced-in siblings that implied the same colour.\n";
	
	vector<int> nColoursImplied(nEdges, 0);
	for (int i = 0; i < nEdges; ++i) {
		nColoursImplied[i] = coloursImpliedByEdge[i].size();
	}
	
	sort(nColoursImplied.begin(), nColoursImplied.end(), greater<int>());
	
	lastChange = -1;
	for (int i = 0; i < nEdges; ++i) {
		if (i == nEdges - 1 || nColoursImplied[i] != nColoursImplied[i + 1]) {
			cerr << (i - lastChange) << " edges have " << nColoursImplied[i] << " implied colours.\n";
			
			lastChange = i;
		}
	}
	
	// Actually delete the edges!
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		edgesBySource[j] = edgesBySource[i];
		impliedEdges[j] = impliedEdges[i];
		coloursImpliedByEdge[j] = coloursImpliedByEdge[i];
		j += keep[i];
	}
	
	cerr << "Deleted " << (nEdges - j) << " edges.\n";
	nEdges = j;
	edgesBySource.resize(nEdges);
	calcVertIndicesForEdgesBySource();
	
	// Update "main" edge list
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);
}

ADD_ACTION(reduceImpliedEdges, "reduce-implied-edges")
void reduceImpliedEdges() {
	assert(!"Don't call this method -- the reduction is actually done inside calc-implied-edges now!");		//HACK
	
	
	
	if (impliedEdges.empty()) {
		cerr << "Must compute implied edges using calc-implied-edges before you can apply reduce-implied-edges!\n";
		exit(1);
	}
	
	cerr << "Reducing edges using implied edges and colour constraints...\n";
	
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		edgesBySource[j] = edgesBySource[i];
		
		vector<bool> seenColour(nCol, false);
		
		bool keep = true;
		for (int k = 0; k < impliedEdges[i].size(); ++k) {
			if (seenColour[c[impliedEdges[i][k].v]]) {
				keep = false;
				break;
			}
			
			seenColour[c[impliedEdges[i][k].v]] = true;
		}
		
		j += keep;
	}
	
	cerr << "Deleted " << (nEdges - j) << " edges.\n";
	nEdges = j;
	edgesBySource.resize(nEdges);
	calcVertIndicesForEdgesBySource();
	
	// Update "main" edge list
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);
}

ADD_ACTION(reduceCombineColoursStrong, "reduce-combcol-strong")
void reduceCombineColoursStrong() {
	cerr << "Combining colours strongly...\n";
	
	vector<vector<int> > coloursImpliedByColour(nCol);
	vector<bool> seenVertex(nVerts, false);
	for (int i = 0; i < nEdges; ++i) {
		int v = edgesBySource[i].v;
		
		if (!seenVertex[v]) {
			//coloursImpliedByColour[c[v]] = coloursImpliedByEdge[i];
			remove_copy(coloursImpliedByEdge[i].begin(), coloursImpliedByEdge[i].end(), back_inserter(coloursImpliedByColour[c[v]]), c[v]);		// Exclude its own colour
			seenVertex[v] = true;
		} else {
			vector<int> temp;
			set_intersection(coloursImpliedByColour[c[v]].begin(), coloursImpliedByColour[c[v]].end(), coloursImpliedByEdge[i].begin(), coloursImpliedByEdge[i].end(), back_inserter(temp));
			coloursImpliedByColour[c[v]] = temp;
		}
	}
	
	for (int i = 0; i < nCol; ++i) {
		if (!coloursImpliedByColour[i].empty()) {
			cerr << "Colour " << i << " implies " << coloursImpliedByColour[i].size() << " other colours:";
			for (int j = 0; j < coloursImpliedByColour[i].size(); ++j) {
				cerr << ' ' << coloursImpliedByColour[i][j];
			}
			cerr << "\n";
		}
	}
	
	//TODO: The whole, you know, reduction.  In particular we still need to check that there are no edges from one colour to another colour that we're trying to combine it with.
}

// Populate lb[][].
//HACK: We currently even calculate lb[u][v] where v is above u, although I don't think we currently use this anywhere.
// Not in the new dom-path reduction anyway.
ADD_ACTION(calcAnchorLowerBounds, "calc-anchor-lbs")
void calcAnchorLowerBounds() {
	assert(!"This function is BROKEN in at least 2 ways: it hasn't swapped the order of subscripts to lb[][], and it seems not to work anyway!");		//TODO
	
	
	
	cerr << "Calculating anchor lower bounds for vertex insertion...\n";
	
	calcEdgesFrom();
	
	lb.assign(nVerts, vector<double>(nVerts, NEGINF));
	for (int u = 0; u < nVerts; ++u) {
		lb[u][0] = 0.0;
	}
	
	int nImprovements = 0;			// Just for diagnostics
	int nEdgeImprovements = 0;		// Likewise
	//for (int v = 1; v < nVerts; ++v) {
	//	// Leave all lb[.][v] at NEGINF if v is unreachable
	//	if (startOfEdgesTo[v] != startOfEdgesTo[v + 1]) {
	//		// Now see whether any of these LBs can be improved via the edges leaving u.
	//		for (int iv = startOfEdgesTo[v]; iv < startOfEdgesTo[v + 1]; ++iv) {
	//			lb[edges[iv].u][v] = max(lb[edges[iv].u][v], edges[iv].w);
	//			//cerr << "D: (u=" << u << ", v=" << edgesBySource[iu].v << "): edge=" << edgesBySource[iu].w << ", new lb=" << lb[u][edgesBySource[iu].v] << "\n";		//DEBUG
	//		}
	//		
	//		for (int u = 0; u < nVerts; ++u) {
	//			if (startOfEdgesTo[u] != startOfEdgesTo[u + 1]) {			// We specifically exclude the root as well -- the only non-NEGINF scores you can get with this anchor are direct edges.
	for (int u = 0; u < nVerts; ++u) {
		// Leave lb[u][v] at NEGINF if u is unreachable
		//if (u == 0 || startOfEdgesTo[u] != startOfEdgesTo[u + 1]) {
		if (startOfEdgesTo[u] != startOfEdgesTo[u + 1]) {			// We specifically exclude the root as well -- the only non-NEGINF scores you can get with this anchor are direct edges.
			for (int v = 1; v < nVerts; ++v) {
				// Leave all lb[.][v] at NEGINF if v is unreachable
				if (startOfEdgesTo[v] != startOfEdgesTo[v + 1]) {
					if (u == v) {
						//cerr << "(u=" << u << ", v=" << v << "): 0.0 (same)\n";		//DEBUG
						lb[u][v] = 0.0;
					} else {
						double worstParent = INF;
						//double directEdge = NEGINF;
						for (int ip = startOfEdgesTo[u]; ip < startOfEdgesTo[u + 1]; ++ip) {
							worstParent = min(worstParent, lb[edges[ip].u][v]);
							//if (edges[ip].u == u) {
							//	directEdge = edges[ip].w;
							//}
						}
						
						//cerr << "(u=" << u << ", v=" << v << "): worstParent=" << worstParent << ", directEdge=" << directEdge << "\n";		//DEBUG
						//if (directEdge < worstParent) {
						//	++nImprovements;
						//	if (directEdge != NEGINF) {
						//		++nEdgeImprovements;
						//	}
						//}
						
						//lb[u][v] = max(directEdge, worstParent);
						lb[u][v] = max(lb[u][v], worstParent);
						//cerr << "P: (u=" << u << ", v=" << v << "): worstParent=" << worstParent << ", new lb=" << lb[u][v] << "\n";		//DEBUG
					}
				}
			}
			
			// Now see whether any of these LBs can be improved via the edges leaving u.
			for (int iu = startOfEdgesFrom[u]; iu < startOfEdgesFrom[u + 1]; ++iu) {
				lb[u][edgesBySource[iu].v] = max(lb[u][edgesBySource[iu].v], edgesBySource[iu].w);
				//cerr << "D: (u=" << u << ", v=" << edgesBySource[iu].v << "): edge=" << edgesBySource[iu].w << ", new lb=" << lb[u][edgesBySource[iu].v] << "\n";		//DEBUG
			}
		}
	}
	
	// Purely for diagnostics.
	for (int i = 0; i < nEdges; ++i) {
		if (edges[i].w < lb[edges[i].u][edges[i].v]) {
			++nEdgeImprovements;
		}
	}
	
	//cerr << nImprovements << " vertex pairs (including " << nEdgeImprovements << " edges) had their bounds improved.\n";
	cerr << nEdgeImprovements << " edges had their bounds improved.\n";
}

// memo[] and seen[] are indexed by u only; v stays fixed throughout the recursion.
//double DEBUGanchorLowerBoundFor(int u, int v, vector<double>& memo, vector<bool>& seen) {
//	cerr << "DEBUGanchorLowerBoundFor(u=" << u << ", v=" << v << ") called.\n";		//DEBUG
//	if (!seen[u]) {
//		double x = INF;
//		
//		if (v == 0) {
//			x = 0.0;
//		} else {
//			for (int iup = startOfEdgesTo[u]; iup < startOfEdgesTo[u + 1]; ++iup) {
//				int p = edges[iup].u;
//				x = min(x, DEBUGanchorLowerBoundFor(p, v, memo, seen));
//			}
//			
//			cerr << "DEBUGanchorLowerBoundFor(u=" << u << ", v=" << v << "): worst parent = " << x << "\n";		//DEBUG
//			
//			// And is there actually an edge directly from u to v?  We could also choose that.
//			// We do this the hard way...
//			for (int ivp = startOfEdgesTo[v]; ivp < startOfEdgesTo[v + 1]; ++ivp) {
//				if (edges[ivp].u == u) {
//					x = max(x, edges[ivp].w);
//					cerr << "DEBUGanchorLowerBoundFor(u=" << u << ", v=" << v << "): found direct in-edge " << edges[ivp] << " with weight " << edges[ivp].w << "\n";		//DEBUG
//					break;
//				}
//			}
//		}
//		
//		// At this point, x will be INF iff u is disconnected from the root.
//		
//		memo[u] = x;
//		seen[u] = true;
//	}
//	
//	cerr << "DEBUGanchorLowerBoundFor(u=" << u << ", v=" << v << ") returning " << memo[u] << ".\n";		//DEBUG
//	return memo[u];
//}
double DEBUGanchorLowerBoundFor(int u, int v, vector<double>& memo, vector<bool>& seen) {
	//cerr << "DEBUGanchorLowerBoundFor(u=" << u << ", v=" << v << ") called.\n";		//DEBUG
	if (u > v) {
		return NEGINF;			// Only handle u <= v
	}
	
	if (!seen[u]) {
		// We require that there be no unreachable edges in the graph.
		//// This is assured if every vertex is either (a) the root, (b) has at least 1 in-edge, or
		//// (c) has no in-edges and no out-edges.  The following assert checks this.
		//assert(u == 0 || startOfEdgesTo[u] != startOfEdgesTo[u + 1] || (startOfEdgesTo[u] == startOfEdgesTo[u + 1] && startOfEdgesFrom[u] == startOfEdgesFrom[u + 1]));
		//double x = NEGINF;
		// This is assured if every vertex with 1 or more out-edges is either (a) the root, or (b) has at least 1 in-edge.
		// The following assert checks this.
		assert(startOfEdgesFrom[u] == startOfEdgesFrom[u + 1] || u == 0 || startOfEdgesTo[u] != startOfEdgesTo[u + 1]);
		
		double x = NEGINF;
		
		if (v == 0) {
			x = 0.0;
		} else {
			double north;
			if (startOfEdgesTo[u] == startOfEdgesTo[u + 1]) {
				// This vertex has no in-edges: Either it's the root, or it has no out-edges either.
				// Either way, there's no possibility to connect v north of it.
				north = NEGINF;
			} else {
				// This vertex is reachable from the root and has at least 1 in-edge.
				north = INF;			// We know there is at least 1 term in the min(), so this will come down.
				for (int iup = startOfEdgesTo[u]; iup < startOfEdgesTo[u + 1]; ++iup) {
					int p = edges[iup].u;
					north = min(north, DEBUGanchorLowerBoundFor(p, v, memo, seen));
				}
			}
			
			//cerr << "DEBUGanchorLowerBoundFor(u=" << u << ", v=" << v << "): north = " << north << "\n";		//DEBUG
			
			double direct = NEGINF;
			
			// And is there actually an edge directly from u to v?  We could also choose that.
			// We do this the hard way...
			for (int ivp = startOfEdgesTo[v]; ivp < startOfEdgesTo[v + 1]; ++ivp) {
				if (edges[ivp].u == u) {
					direct = edges[ivp].w;
					//cerr << "DEBUGanchorLowerBoundFor(u=" << u << ", v=" << v << "): found direct in-edge " << edges[ivp] << " with weight " << direct << "\n";		//DEBUG
					break;
				}
			}
			
			// If u is root-reachable and v == u, then we trivially can force v in for 0.
			if (u == v && startOfEdgesTo[v] != startOfEdgesTo[v + 1]) {
				direct = 0.0;
			}
			
			x = max(north, direct);
		}
		
		memo[u] = x;
		seen[u] = true;
	}
	
	//cerr << "DEBUGanchorLowerBoundFor(u=" << u << ", v=" << v << ") returning " << memo[u] << ".\n";		//DEBUG
	return memo[u];
}

// Populate lb[][] in the easiest possible way.  For comparison with calc-anchor-lbs, which I'm not 100% sure of.
//HACK: We currently even calculate lb[u][v] where v is above u, although I don't think we currently use this anywhere.
// Not in the new dom-path reduction anyway.
ADD_ACTION(DEBUGcalcAnchorLowerBounds, "DEBUG-calc-anchor-lbs")
void DEBUGcalcAnchorLowerBounds() {
	cerr << "Calculating anchor lower bounds for vertex insertion (DEBUG version)...\n";
	
	calcEdgesFrom();
	
	//lb.assign(nVerts, vector<double>(nVerts, NEGINF));
	//lb.assign(nVerts, vector<double>(nVerts, 4200.0069));		//DEBUG: Hopefully this number is sufficiently weird-looking
	lb.resize(nVerts);
	for (int v = 0; v < nVerts; ++v) {
		lb[v].assign(v + 1, NEGINF);			// lb[v][u] is defined only for u <= v
	}
	//for (int u = 0; u < nVerts; ++u) {
	//	lb[u][0] = 0.0;
	//}
	
	int nInfsSeen = 0;			// Purely for diagnostics
	int nNegInfsSeen = 0;
	int nNonInfsSeen = 0;
	for (int v = 0; v < nVerts; ++v) {
		vector<bool> seen(nVerts, false);
		vector<double> memo(nVerts);		// Initial value doesn't matter
		//for (int u = 0; u < nVerts; ++u) {
		for (int u = 0; u <= v; ++u) {
			//lb[u][v] = DEBUGanchorLowerBoundFor(u, v, memo, seen);
			lb[v][u] = DEBUGanchorLowerBoundFor(u, v, memo, seen);
			
			////DEBUG
			//if (lb[v][u] == INF) {
			//	++nInfsSeen;
			//	//lb[u][v] = NEGINF;
			//} else if (lb[v][u] == NEGINF) {
			//	++nNegInfsSeen;
			//} else {
			//	++nNonInfsSeen;
			//	//cerr << "Saw non-INF LB: lb[" << u << "][" << v << "]=" << lb[u][v] << "!\n";		//DEBUG
			//}
		}
	}
	
	//cerr << "Anchor lower bound calculation (DEBUG version) finished.  " << nInfsSeen << " INFs seen; " << nNegInfsSeen << " NEGINFs seen; " << nNonInfsSeen << " non-[NEG]INFs seen.\n";
	cerr << "Anchor lower bound calculation (DEBUG version) finished.\n";
	dumpTable("DEBUGcalcAnchorLowerBounds", lb);
}

// The appendix of the paper reads better without "DEBUG"...  I know.  I know.
ADD_ACTION(calcAnchorLowerBoundsNew, "calc-anchor-lbs-new")
void calcAnchorLowerBoundsNew() {
	DEBUGcalcAnchorLowerBounds();
}

int DEBUGreduceDomPath2StopAfterN = INT_MAX;		// Stop after deleting this many edges in a row (also used for reduce-colsubtree-adv)

// Dominating Path II: The Redomination!  (I.e. this time it actually works.)
ADD_ACTION(reduceDominatingPath2, "reduce-dompath2")
void reduceDominatingPath2() {
	cerr << "Reducing edges using dominating path 2 reduction...\n";
	
	if (ub.empty()) {
		cerr << "Must compute vertex upper bounds (using calc-child-vertex-ubs and/or calc-col-forest-vertex-ubs) before you can apply reduce-dompath2!\n";
		exit(1);
	}
	
	if (m_.empty()) {
		cerr << "Must compute a slide lower bound (e.g. with calc-rec-slide-lbs) before running reduce-dompath2!\n";
		exit(1);
	}
	
	if (lb.empty()) {
		cerr << "You must calculate anchor lower bounds using calc-anchor-lbs before running reduce-dompath2.\n";
		exit(1);
	}
	
	calcEdgesFrom();
	
	cerr << "Calculating max in-edges to each vertex...\n";
	vector<double> maxInEdgeTo(nVerts, NEGINF);
	for (int i = 0; i < nEdges; ++i) {
		maxInEdgeTo[edges[i].v] = max(maxInEdgeTo[edges[i].v], edges[i].w);		// edges[] better for locality than edgesBySource[]
	}
	
	// severAndSlideTo[u][v - firstVertOfSameColourAs[v]] is the maximum it could cost to delete an in-edge to any other vertex u of that colour and slide the subtree below u to being below v.
	// The indexing system drops us from O(V^2) space to O(VD), where D is the maximum number of vertices of any colour.
	// (Note: we require u and v to be of the same colour, so firstVertOfSameColourAs[v] == firstVertOfSameColourAs[u].)
	vector<double> severAndSlideTo(nVerts, INF);
	
	// Now actually populate the array.
	cerr << "Calculating max sever-and-slide costs for each vertex...\n";		// "max cost" = "min score change", and we record the latter.
	for (int v = 0; v < nVerts; ++v) {
		for (int u = firstVertOfSameColourAs[v]; u < firstVertOfNextColour[v]; ++u) {
			severAndSlideTo[v] = min(severAndSlideTo[v], slideLb(u, v) - maxInEdgeTo[u]);		// Can be > 0 if all in-edges to c[v] are negative.
			if ( v == 663 ) {
				cerr << "next vert is: " << firstVertOfNextColour[v] << "\n";
			}
		}
	}
	
	dumpTable("severAndSlide", severAndSlideTo);
	
	cerr << "Reducing edges...\n";
	vector<bool> keep(nEdges, true);			// keep[i] == true iff edgesBySource[i] should be retained.
	int nEdgesDeleted = 0;
	for (int y = 0; y < nVerts; ++y) {
		// Lazily calculate dplb[a] for all a, given that y is in the solution tree.
		vector<double> dplb(nVerts);		// dplb[i] is the minimum cost to force vertex i into the solution, given that y is already in it
		int nextVertToProcess = firstVertOfNextColour[y];			// All lower-numbered vertices have already been calculated.
		
		for (int ix = startOfEdgesFrom[y]; ix < startOfEdgesFrom[y + 1]; ++ix) {
			int x = edgesBySource[ix].v;
			
			// Lazily calculate dplb[] up to what we need to handle x.
			for (; nextVertToProcess < firstVertOfSameColourAs[x]; ++nextVertToProcess) {
				// Calculate bestConn(a, y)
				// The next line gets executed O(V^2) times.
				double bestConn = lb[nextVertToProcess][y];
				for (int ip = startOfEdgesTo[nextVertToProcess]; ip < startOfEdgesTo[nextVertToProcess + 1]; ++ip) {
					// We're only allowed to try forcing in vertices that could not force y out of the graph (e.g. by removing all paths from the root to it).
					// Vertices of colours greater than c[y] is a safe set of vertices to play with.
					// (Actually that's not true -- we could try forcing in vertices of colours <= c[y] if we calculated lb[v][u] for all pairs of vertices instead of
					// only for pairs with u < v.  In fact we could even try forcing in a different vertex of colour c[y], and this would be *cheaper*
					// in general because we know which c[y]-coloured vertex is already in the graph (it's y!) -- but in this case we can still hope that all in-edges
					// to y to get killed later on anyway, which will lead to (y, x) being deleted by reduce-unreach.)
					if (c[edges[ip].u] > c[y]) {
						bestConn = max(bestConn, dplb[edges[ip].u] + edges[ip].w);
					}
				}
				
				dplb[nextVertToProcess] = min(0.0, min(bestConn, bestConn + severAndSlideTo[nextVertToProcess]));
			}
			
			// For each vertex z of the same colour as x:
			for (int z = firstVertOfSameColourAs[x]; z < firstVertOfNextColour[x] && keep[ix]; ++z) {
				//if (z != x) {
					for (int iv = startOfEdgesTo[z]; iv < startOfEdgesTo[z + 1]; ++iv) {
						int v = edges[iv].u;
						
						if (c[v] > c[y]) {
							double oc = edges[iv].w + dplb[v];
#ifdef USE_OC
							edgesBySource[ix].oc = max(edgesBySource[ix].oc, oc);		// Update the opportunity cost LB as a nice side-effect :)
#endif	// USE_OC
							if (oc + slideLb(x, z) > edgesBySource[ix].w) {
								// It's possible to get a better solution by deleting (y, x), inserting (v, z), and connecting v to an ancestor of y somehow.
								//DEBUG
								if (!watchedEdges.empty() && watchedEdges[0].u == y && watchedEdges[0].v == x) {
									//DEBUG
									cerr << "Edge #" << ix << ", " << edgesBySource[ix] << " is about to be deleted!  Writing out the current instance...\n";
									cerr << "(v, z) = " << edges[iv] << "\n";
									cerr << "w(v, z) (edges[iv].w) = " << edges[iv].w << "\n";
									cerr << "w(y, x) (edgesBySource[ix].w) = " << edgesBySource[ix].w << "\n";
									cerr << "dplb(y, v) (dplb[v]) = " << dplb[v] << "\n";
									//cerr << "lb[y][v] = " << lb[y][v] << "\n";
									cerr << "lb[v][y] = " << lb[v][y] << "\n";
									cerr << "edges[iv].w - edgesBySource[ix].w + dplb[v] = " << (edges[iv].w - edgesBySource[ix].w + dplb[v]) << "\n";
									cerr << ">>>>> BUT: <<<<<\n";
									//cerr << "m[x][z] = " << m[x][z] << "\n";
									cerr << "slideLb(x, z) = " << slideLb(x, z) << "\n";
									//cerr << "edges[iv].w - edgesBySource[ix].w + dplb[v] + m[x][z] = " << (edges[iv].w - edgesBySource[ix].w + dplb[v] + m[x][z]) << "\n";
									cerr << "edges[iv].w - edgesBySource[ix].w + dplb[v] + slideLb(x, z) = " << (edges[iv].w - edgesBySource[ix].w + dplb[v] + slideLb(x, z)) << "\n";
									cerr << ">>>>> HMMM <<<<<\n";
									cerr << "severAndSlideTo[v] = " << severAndSlideTo[v] << "\n";
									cerr << "maxInEdgeTo[v] = " << maxInEdgeTo[v] << "\n";
									//cerr << "Therefore
									
									writeOutputToStdout();
									cerr << "And now we resume, but presumably --watch-edges will kick in in a sec...\n";
									//exit(1);
								}
								
								keep[ix] = false;
								if (++nEdgesDeleted == DEBUGreduceDomPath2StopAfterN) {
									cerr << "Stopping early due to hitting " << DEBUGreduceDomPath2StopAfterN << " edge deletions.\n";		//DEBUG
									goto FINISHED;			//HACK: For debugging
								}
								
								break;		// Test on keep[ix] in outer loop ensures we break out of it too
							}
						}
					}
				//}
			}
		}
	}
	
FINISHED:		//HACK: For debugging with DEBUGreduceDomPath2StopAfterN
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		edgesBySource[j] = edgesBySource[i];
		j += keep[i];
	}
	
	cerr << (nEdges - j) << " edges deleted; " << j << " edges remain.\n";
	
	nEdges = j;
	edgesBySource.resize(nEdges);
	calcVertIndicesForEdgesBySource();
	
	// Update "main" edge list
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);
}

ADD_ACTION(calcVertexImpliedEdges, "calc-vertex-implied-edges")
void calcVertexImpliedEdges() {
	if (impliedEdges.empty()) {
		cerr << "You must calculate edges implied by edges using calc-implied-edges before running calc-vertex-implied-edges.\n";
		exit(1);
	}
	
	cerr << "Calculating edges implied by vertices...\n";
	
	// First find the index of the maximum-weight edge into each vertex.
	// All other in-edges will imply at least the edges that this edge implies.
	vector<int> maxInEdgeIndex(nVerts, -1);
	for (int i = 0; i < nEdges; ++i) {
		int v = edgesBySource[i].v;
		if (maxInEdgeIndex[v] == -1 || edgesBySource[i].w > edgesBySource[maxInEdgeIndex[v]].w) {
			maxInEdgeIndex[v] = i;
		}
	}
	
	// Now assign that edge's list of implied edges to the vertex.
	edgesImpliedByVertex.assign(nVerts, vector<edge>());
	for (int i = 0; i < nVerts; ++i) {
		if (maxInEdgeIndex[i] != -1) {
			edgesImpliedByVertex[i] = impliedEdges[maxInEdgeIndex[i]];
		}
	}
	
	// Show some interesting stats
	//vector<int> nEdgesImpliedFreq(nEdges + 1, 0);		// Can't imply more than nEdges edges
	vector<int> nEdgesImpliedFreq;
	for (int i = 0; i < nVerts; ++i) {
		if (nEdgesImpliedFreq.size() <= edgesImpliedByVertex[i].size()) {
			nEdgesImpliedFreq.resize(edgesImpliedByVertex[i].size() + 1);		// Just grow it as necessary
		}
		
		++nEdgesImpliedFreq[edgesImpliedByVertex[i].size()];
	}
	
	for (int i = nEdgesImpliedFreq.size() - 1; i >= 0; --i) {
		if (nEdgesImpliedFreq[i]) {
			cerr << nEdgesImpliedFreq[i] << " vertices have " << i << " implied edges.\n";
		}
	}
}

//// Insert info for all edges reachable from x.  We make use of the assumption that some edge (u, v) is known to be in the tree already
//// (we don't need to know v).
//// bestEdges[i] = info about the best edge into colour i from any vertex reachable from v.
//// seen[u] = true iff all edges leaving u have already been processed.
////TODO: Currently testing by ignoring the weight of uv -- should give identical results to applying seb-vertex-ubs then reduce-vub.
//void colSubtreeAdvantageGatherEdges(int x, int u, vector<inEdgeToColourInfo>& bestEdges, vector<bool>& seen, double inEdgeWeight) {
//	// Even if x has already been seen, we still need to update maxInEdge for all of its outgoing edges.
//	for (int iy = startOfEdgesFrom[x]; iy < startOfEdgesFrom[x + 1]; ++iy) {
//		int y = edgesBySource[iy].v;
//		
//		//// Lazily initialise for this vertex
//		//if (bestEdges[y].empty()) {
//		//	bestEdges[y].resize(nCol, initInEdgeToColourInfo);
//		//}
//		
//		if (!seen[x]) {
//			//colForestUbMergeEdge(bestEdges[c[y]], edgesBySource[iy], inEdgeWeight);
//			edge temp(edgesBySource[iy]);
//			temp.w -= max(0.0, lbCol[u][c[y]]);
//			colForestUbMergeEdge(bestEdges[c[y]], temp, inEdgeWeight);
//			
//			// Recurse to gather rest of reachable edges
//			//colSubtreeAdvantageGatherEdges(y, u, bestEdges, seen, edgesBySource[iy].w);
//			colSubtreeAdvantageGatherEdges(y, u, bestEdges, seen, temp.w);
//		}
//		
//		bestEdges[c[y]].maxInEdge = max(bestEdges[c[y]].maxInEdge, inEdgeWeight);
//	}
//	
//	seen[x] = true;
//}

int reduceColSubtreeAdvantageNEdgesIgnored;
int reduceColSubtreeAdvantageNBetterEdgesNorth;
int reduceColSubtreeAdvantageNPositiveEdgesNorth;
int reduceColSubtreeAdvantageNEdgesDeleted;		//DEBUG: Used for stopping early

// Insert info for all edges reachable from x.  We make use of the assumption that some edge (u, v) is known to be in the tree already
// (we don't need to know v).
// bestEdges[i] = info about the best edge into colour i from any vertex reachable from v.
// seen[u] = true iff all edges leaving u have already been processed.
// If lb[y][u] > w(x, y) then we know that (x, y) cannot be part of any optimal tree anyway, so we can ignore it.
//TODO: Currently testing by ignoring the weight of uv -- should give identical results to applying seb-vertex-ubs then reduce-vub.
// bestConn will normally be lb[x][u], but when x == v it will be set to the minimum of lb[x][p] over all u's parents p.  That's
// to avoid the possibility of trying to add (u, v) itself back in as a way of connecting x!  Bit of a hack.
void colSubtreeAdvantageGatherEdges(int x, int u, vector<inEdgeToColourInfo>& bestEdges, vector<bool>& seen, double inEdgeWeight, double bestConn) {
	// Even if x has already been seen, we still need to update maxInEdge for all of its outgoing edges.
	for (int iy = startOfEdgesFrom[x]; iy < startOfEdgesFrom[x + 1]; ++iy) {
		int y = edgesBySource[iy].v;
		
		if (!seen[x]) {
			if (lb[y][u] < edgesBySource[iy].w) {
				// It's possible that (x, y) is part of some optimal (u, v)-containing solution.
				//colForestUbMergeEdge(bestEdges[c[y]], edgesBySource[iy], inEdgeWeight);
				edge temp(edgesBySource[iy]);
				//temp.w -= max(0.0, lbCol[u][c[y]]);
				
				// We can either delete the edge, for a cost of -w(x, y), or we can try attaching
				// x to an ancestor of u, for a cost of lb[x][u] -- but note that in the latter case,
				// if lb[x][u] is positive, then we need to divide by a lower bound on the number of
				// paths that might force this edge in to ensure that we never wind up
				// counting it more than once in total.  If we allowed this edge to be forced in by deeper
				// edges, then a suitable lower bound would be the number of remaining colours;
				// since we currently only allow x's immediate children to force this edge in,
				// we can use either that or the number of children of x, whichever is smaller.
				// Note that the sense of cost is reversed here.
				// Also note that x == v is a special case, which we handle by calculating bestConnUb in the caller.
				if (bestConn > 0) {
					bestConn /= min(nCol - c[x], startOfEdgesFrom[x + 1] - startOfEdgesFrom[x]);
					//bestConn = 0;		//DEBUG: Very conservative!
					++reduceColSubtreeAdvantageNPositiveEdgesNorth;
				}
				
				if (-bestConn < temp.w) {
					++reduceColSubtreeAdvantageNBetterEdgesNorth;
				}
				temp.w = min(temp.w, -bestConn);
				
				colForestUbMergeEdge(bestEdges[c[y]], temp, inEdgeWeight);
				
				// Recurse to gather rest of reachable edges
				colSubtreeAdvantageGatherEdges(y, u, bestEdges, seen, temp.w, lb[y][u]);		// It's always safe to connect a non-child descendant of u directly to u
			} else {
				++reduceColSubtreeAdvantageNEdgesIgnored;
			}
		}
		
		if (bestEdges[c[y]].edges[0].u == x) {
			bestEdges[c[y]].maxInEdge = max(bestEdges[c[y]].maxInEdge, inEdgeWeight);
		}
	}
	
	seen[x] = true;
}

void doColSubtreeAdvantageReductionFor(int u, vector<bool>& keep) {
	for (int iv = startOfEdgesFrom[u]; iv < startOfEdgesFrom[u + 1]; ++iv) {
		int v = edgesBySource[iv].v;
		
		vector<inEdgeToColourInfo> bestEdges(nCol, initInEdgeToColourInfo);		// Since it's only a 1D array, might as well initialise it right away.
		vector<bool> seen(nVerts, false);
		
		// First, because we don't store a table containing L_{anc'}(u, v), we need to calculate that
		// on-the-fly for (u, v).  This is a special case.
		double bestConn = INF;
		for (int ip = startOfEdgesTo[u]; ip < startOfEdgesTo[u + 1]; ++ip) {
			bestConn = min(bestConn, lb[v][edges[ip].u]);
		}
		
		if (u == 0) {
			bestConn = NEGINF;
		}
		
		// Now gather up info about the best, second-best and maxInEdge-to-the-best for all edges reachable from v,
		// given that we know u is in the solution.  This is what makes the complexity O(E^2*C)...
		colSubtreeAdvantageGatherEdges(v, u, bestEdges, seen, NEGINF, bestConn);
		
		// Can we safely delete edge (u, v)?
		double total;
		if (shouldStrengthenColForestUpperBound) {
			total = calcStrengthenedColForestUpperBoundFor(v, bestEdges);
		} else {
			total = accumulate(bestEdges.begin(), bestEdges.end(), 0.0, accum_weights_from_first());
		}
		
		if (total + edgesBySource[iv].w < 0.0) {
			keep[iv] = false;			// This edge can be deleted.
			
			//DEBUG
			if (u == 0 && v == 67) {
				cerr << "Discarding (u=" << u << ", v=" << v << ") by a margin of total=" << total << " + edgesBySource[iv].w=" << edgesBySource[iv].w << " = " << (total + edgesBySource[iv].w) << "!\n";
			}
			if (++reduceColSubtreeAdvantageNEdgesDeleted == DEBUGreduceDomPath2StopAfterN) {
				return;
			}
		}
	}
}

// Adapted from seb-vertex-ubs.  When considering a particular edge uv, instead of using "raw" edge costs for an edge xy reachable from v,
// we subtract the weight of the best edge that we know we could add to a vertex of colour c[y].
// (This means we have to do a full recursion through every reachable edge, like the original (very inefficient) version
// of seb-vertex-ubs.)
ADD_ACTION(reduceColSubtreeAdvantage, "reduce-colsubtree-adv")
void reduceColSubtreeAdvantage() {
	cerr << "Performing colourful subtree advantage reduction...\n";
	
	//if (lbCol.empty()) {
	//	cerr << "Must compute anchor-to-colour lower bounds (using calc-anc-col-lbs) before you can apply reduce-colsubtree-adv!\n";
	//	exit(1);
	//}
	if (lb.empty()) {
		cerr << "Must compute anchor lower bounds (using calc-anchor-lbs) before you can apply reduce-colsubtree-adv!\n";
		exit(1);
	}
	
	calcEdgesFrom();
	
	reduceColSubtreeAdvantageNEdgesIgnored = 0;
	reduceColSubtreeAdvantageNPositiveEdgesNorth = 0;
	reduceColSubtreeAdvantageNEdgesDeleted = 0;		//DEBUG
	vector<bool> keep(nEdges, true);
	for (int u = 0; u < nVerts; ++u) {
		doColSubtreeAdvantageReductionFor(u, keep);
		
		//DEBUG
		if (reduceColSubtreeAdvantageNEdgesDeleted == DEBUGreduceDomPath2StopAfterN) {
			break;
		}
	}
	
	int j = 0;
	for (int i = 0; i < nEdges; ++i) {
		edgesBySource[j] = edgesBySource[i];
		j += keep[i];
	}
	
	cerr << (nEdges - j) << " edges deleted. " << reduceColSubtreeAdvantageNEdgesIgnored << " edges were ignored due to better edges from (an anc of) u. " << reduceColSubtreeAdvantageNBetterEdgesNorth << " edges from (an anc of) u were better than deleting the edge; of these, " << reduceColSubtreeAdvantageNPositiveEdgesNorth << " were positive.\n";
	nEdges = j;
	edgesBySource.resize(nEdges);
	calcVertIndicesForEdgesBySource();
	
	// Update "main" edge list
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);
}

// Have any of the edges list with --watch-edge x y been deleted?  If so, sound the alarm.
void checkWatchedEdges() {
	for (int i = 0; i < watchedEdges.size(); ++i) {
		// Use O(n) linear search to find the right in-edge
		//HACK: If we maintained edges in edges[] or edgesBySource[] in a completely sorted list, then we could just use binary_search() directly,
		// which would be faster, but this should be fast enough because we only do a brute-force search through the in-edges of the correct vertex.
		int j = startOfEdgesTo[watchedEdges[i].v];
		for (; j < startOfEdgesTo[watchedEdges[i].v + 1]; ++j) {
			if (edges[j].u == watchedEdges[i].u) {
				break;		// Found it
			}
		}
		
		if (j == startOfEdgesTo[watchedEdges[i].v + 1]) {
			cerr << "Edge " << watchedEdges[i] << " has been deleted!\n";
			exit(1);
		}
	}
}

// Parses a very simple language consisting of a sequence of terms, in which each term must be either an action name or a parenthesised sequence of action names,
// optionally preceded by either an integer repetition count or "*" (meaning "repeat until no more edges have been deleted").
// The nice thing about this language is that the parser doesn't ever have to backtrack, so we can read from an arbitrary istream
// without unholy complications...
// Note that all terms, including parentheses, must be separated from each other with whitespace.
// Using a tagged union instead of subclasses is not very C++ish, but who cares :-P
class Command {
	enum command_type { ACTION, SEQUENCE };
	
	command_type type;
	union {
		struct {
			void (*func)();
		} act;
		struct {
			Command* first;
			Command* second;
			int repeatCount;		// Applies to the first command.  -1 means "repeat until no more edges have been deleted".
		} seq;
	} data;
	
public:
	static Command* parse(istream& is) {
		string token;
		is >> token;
		
		if (!is || token == ")") {
			return NULL;
		}
		
		Command* cmd = new Command;
		cmd->type = SEQUENCE;
		cmd->data.seq.repeatCount = 1;
		if (token == "*") {
			// Next command should be repeated until no more edges are deleted.
			cmd->data.seq.repeatCount = -1;
			is >> token;
			if (!is || token == ")") {
				cerr << "Unexpected ')' or end of stream after '*'!\n";
				exit(1);
			}
		} else if (token.find_first_not_of("0123456789") == string::npos) {
			cmd->data.seq.repeatCount = atoi(token.c_str());
			is >> token;
			if (!is || token == ")") {
				cerr << "Unexpected ')' or end of stream after repeat count of " << cmd->data.seq.repeatCount << "!\n";
				exit(1);
			}
		}
		
		// At this point, we have just read the command into token.
		Command* firstCmd;
		if (token == "(") {
			//firstCmd = new Command(SEQUENCE, repeatCount, parse(is));		// Recurse to process parenthesised expression
			firstCmd = parse(is);		// Recurse to process parenthesised expression, which will gobble the trailing ")"
		} else {
			// It must be an action name.
			map<string, void (*)()>::const_iterator iter(actions.find(token));
			if (iter == actions.end()) {
				cerr << "Unrecognised command '" << token << "'.\n";
				exit(1);
			} else {
				firstCmd = new Command;
				firstCmd->type = ACTION;
				firstCmd->data.act.func = (*iter).second;
			}
		}
		
		cmd->data.seq.first = firstCmd;
		cmd->data.seq.second = parse(is);		// Recurse to process rest of input
		return cmd;
	}
	
	static void destroy(Command* cmd) {
		if (!cmd) {
			return;
		}
		
		if (cmd->type == SEQUENCE) {
			destroy(cmd->data.seq.first);
			destroy(cmd->data.seq.second);
		}
		
		delete cmd;
	}
	
	void run() {
		if (type == ACTION) {
			// Run the action!
			//pair<int, clock_t>* pt;
			action_info* pt;
			if (timeActions) {
				pt = &actionTimes[data.act.func];		// Creates a new entry if none exists
				//++pt->first;
				//pt->second -= clock();
				++pt->nIterations;
				pt->elapsed -= lastClock;				// Only call clock() once, to avoid "cracks" in time
				pt->nEdgesDeleted += nEdges;
			}
			
			(*data.act.func)();
			
			if (timeActions) {
				//pt->second += clock();
				lastClock = clock();
				pt->elapsed += lastClock;
				pt->nEdgesDeleted -= nEdges;
			}
			
			// If necessary, perform any desired checks.
			// It would be more elegant to build this into the action directly, but C++ isn't really made for dealing with first-order functions.
			checkWatchedEdges();
		} else {
			// SEQUENCE type
			assert(type == SEQUENCE);
			
			if (data.seq.repeatCount == -1) {
				int nEdgesBefore = INT_MAX;		// INT_MAX ensures we get at least 1 iteration.
				while (nEdges < nEdgesBefore) {
					nEdgesBefore = nEdges;
					data.seq.first->run();
				}
			} else {
				for (int i = 0; i < data.seq.repeatCount; ++i) {
					data.seq.first->run();
				}
			}
			
			if (data.seq.second) {
				data.seq.second->run();
			}
		}
	}
};

vector<int> verticesToKeepBelow;		// Contains a list of vertex numbers; if keep-subgraph-below-verts is run, then only edges below any of these vertices will be retained.

// Called by keepSubgraphBelow().  Accumulates edges into edges[], which should start out empty.
void addEdgesBelow(int v, vector<bool>& seen) {
	//cerr << "addEdgesBelow(v=" << v << ") called.  seen=" << seen[v] << ".\n";		//DEBUG
	if (!seen[v]) {
		seen[v] = true;
		
		for (int i = startOfEdgesFrom[v]; i < startOfEdgesFrom[v + 1]; ++i) {
			edges.push_back(edgesBySource[i]);
			addEdgesBelow(edgesBySource[i].v, seen);
		}
	}
}

// Discard all edges in the graph, except those reachable from vertices in verticesToKeepBelow.
//HACK: The vertices could be specified much more elegantly if our parser allowed actions to take parameters...
ADD_ACTION(keepSubgraphBelow, "keep-subgraph-below-verts")
void keepSubgraphBelow() {
	//cerr << "Deleting all edges except those reachable from the " << verticesToKeepBelow.size() << " vertices specified with --keep-edges-below...\n";
	calcEdgesFrom();						// edgesBySource[] probably hasn't even been initialised yet!
	edges.clear();		// We will use this to hold the new list of edges, and recurse through edgesBySource.
	vector<bool> seen(nVerts, false);
	for (int i = 0; i < verticesToKeepBelow.size(); ++i) {
		addEdgesBelow(verticesToKeepBelow[i], seen);
	}
	
	cerr << "Deleted " << (nEdges - edges.size()) << " edges: " << edges.size() << " edges remain.\n";
	nEdges = edges.size();
	
	sortEdgesAndCalcVertIndices(true);		// Put edges[] in a proper state
	edgesBySource.clear();		// calcEdgesFrom() needs edgesBySource[] to be empty first; this is inefficient, but who cares.
	calcEdgesFrom();						// Put edgesBySource[] in a proper state
}

vector<int> verticesToKeepAbove;		// Contains a list of vertex numbers; if keep-subgraph-above-verts is run, then only edges on a path to any of these vertices will be retained.

// Called by keepSubgraphAbove().  Accumulates edges into edgesBySource[], which should start out empty.
void addEdgesAbove(int v, vector<bool>& seen) {
	//cerr << "addEdgesAbove(v=" << v << ") called.  seen=" << seen[v] << ".\n";		//DEBUG
	if (!seen[v]) {
		seen[v] = true;
		
		for (int i = startOfEdgesTo[v]; i < startOfEdgesTo[v + 1]; ++i) {
			edgesBySource.push_back(edges[i]);
			addEdgesAbove(edges[i].u, seen);
		}
	}
}

ADD_ACTION(keepSubgraphAbove, "keep-subgraph-above-verts")
void keepSubgraphAbove() {
	cerr << "Deleting all edges except those on some path from the root to " << verticesToKeepAbove.size() << " vertices specified with --keep-edges-above...\n";
	edgesBySource.clear();		// We will use this to hold the new list of edges, and recurse through edges.
	vector<bool> seen(nVerts, false);
	for (int i = 0; i < verticesToKeepAbove.size(); ++i) {
		addEdgesAbove(verticesToKeepAbove[i], seen);
	}
	
	cerr << "Deleted " << (nEdges - edgesBySource.size()) << " edges: " << edgesBySource.size() << " edges remain.\n";
	nEdges = edgesBySource.size();
	
	edges = edgesBySource;
	sortEdgesAndCalcVertIndices(true);		// Put edges[] in a proper state
	edgesBySource.clear();		// So that anyone who needs these edges in future knows to recalculate them
}

int dottyPrecision = -1;		// Number of digits to use for edge labels when write-dotty is called; -1 (the default) means no edge labels.

void writeDottyOutput(ostream& os) {
	calcEdgesFrom();						// edgesBySource[] probably hasn't even been initialised yet!
	
	os << "digraph {\n";
	//os << "clustermode = local\n";
	os << "newrank = true\t//HACK: This undocumented property, found on http://www.graphviz.org/content/rank-between-cluster, enables us to both put boxes around clusters AND have each colour on a distinct row; without it, you can only get one or the other!\n";
	
	int oldPrec = os.precision();		// Ugh, I hate iostreams
	if (dottyPrecision >= 0) {
		os.precision(dottyPrecision);
	}
	
	for (int i = 0; i < nEdges; ++i) {
		os << edges[i].u << " -> " << edges[i].v;
		if (dottyPrecision >= 0) {
			os << " [label = " << edges[i].w << "]";
		}
		
		os << "\n";
	}
	
	// Add "subgraphs" to group all nodes of the same colour into the same row.
	// We only mention a node if it has at least one edge arriving at it or leaving it, and
	// we only mention a colour if it has at least one such node.
	int aVertexFromTheLastColourThatActuallyAppears = -1;
	vector<int> verticesThatActuallyAppear;
	for (int v = 0; v < nVerts; ++v) {
		if (startOfEdgesTo[v] != startOfEdgesTo[v + 1] || startOfEdgesFrom[v] != startOfEdgesFrom[v + 1]) {
			verticesThatActuallyAppear.push_back(v);
		}
		
		if (v == nVerts - 1 || c[v + 1] != c[v]) {
			if (!verticesThatActuallyAppear.empty()) {
				os << "subgraph cluster_col" << c[v] << " {\nrank = same\nlabel = \"Colour " << c[v] << "\"\n";
				for (int i = 0; i < verticesThatActuallyAppear.size(); ++i) {
					os << verticesThatActuallyAppear[i] << "\n";
				}
				
				os << "}\n";
				
				// Also add an "invisible edge" from the previous actually-appearing colour to this one,
				// to make sure that vertices of different colours appear on different rows (otherwise this won't happen!)
				if (aVertexFromTheLastColourThatActuallyAppears != -1) {
					os << aVertexFromTheLastColourThatActuallyAppears << " -> " << verticesThatActuallyAppear[0] << " [style = invis]\n";
				}
				
				aVertexFromTheLastColourThatActuallyAppears = verticesThatActuallyAppear[0];
				verticesThatActuallyAppear.clear();
			}
		}
	}
	
	os << "}\n";
	
	os.precision(oldPrec);
}

// Write out a dotty-format file for drawing a directed graph.  All nodes of a given colour will appear on the same row.
ADD_ACTION(writeDottyAction, "write-dotty")
void writeDottyAction() {
	if (outputFName) {
		cerr << "Writing dotty output to file \"" << outputFName << "\"...\n";
		ofstream out(outputFName);
		writeDottyOutput(out);
	} else {
		cerr << "Writing dotty output to stdout...\n";
		writeDottyOutput(cout);
	}
	
	cerr << "Wrote out dotty graph with " << nVerts << " vertices, " << nEdges << " edges and " << nCol << " colours.\n";
}

int main(int argc, char **argv) {
	string programString;
	
	//HACK: Needed to find bugs!
	cerr.precision(30);
	cout.precision(30);
	
	// Process command-line arguments
	for (int i = 1; i < argc; ++i) {
		if (string(argv[i]) == "--compare-mode") {
			++i;
			if (string(argv[i]) == "dontcare") {
				compareMode = DONTCARE;
			} else if (string(argv[i]) == "ascbystart") {
				compareMode = ASCBYSTART;
			} else if (string(argv[i]) == "descbystart") {
				compareMode = DESCBYSTART;
			} else if (string(argv[i]) == "ascbyweight") {
				compareMode = ASCBYWEIGHT;
			} else if (string(argv[i]) == "descbyweight") {
				compareMode = DESCBYWEIGHT;
			} else if (string(argv[i]) == "ascasc") {
				compareMode = ASCASC;
			} else if (string(argv[i]) == "ascdesc") {
				compareMode = ASCDESC;
			} else if (string(argv[i]) == "descasc") {
				compareMode = DESCASC;
			} else if (string(argv[i]) == "descdesc") {
				compareMode = DESCDESC;
			} else {
				cerr << "Unrecognised compare mode '" << argv[i] << "'.\n";
				return 1;
			}
		} else if (string(argv[i]) == "--save-tables") {
			shouldSaveTables = true;
		} else if (string(argv[i]) == "--check-preconds") {
			shouldCheckPreconds = true;
		} else if (string(argv[i]) == "--watch-edge") {
			edge e;
			e.u = atoi(argv[++i]);
			e.v = atoi(argv[++i]);
			e.w = -1;		// Ignored
			watchedEdges.push_back(e);
			cerr << "Watching edge " << e << ".\n";		// For debugging.
		} else if (string(argv[i]) == "--keep-edges-below") {	// Mostly for debugging.
			int v = atoi(argv[++i]);
			verticesToKeepBelow.push_back(v);
			cerr << "Will keep subgraph below vertex " << v << " when keep-subgraph-below-verts is run.\n";		//HACK: Yes, this is an ugly hack forced by the fact that we can't easily supply arguments to actions...
		} else if (string(argv[i]) == "--keep-edges-above") {	// Mostly for debugging.
			int v = atoi(argv[++i]);
			verticesToKeepAbove.push_back(v);
			cerr << "Will keep subgraph above vertex " << v << " when keep-subgraph-above-verts is run.\n";		//HACK: Yes, this is an ugly hack forced by the fact that we can't easily supply arguments to actions...
		} else if (string(argv[i]) == "--dotty-precision") {
			dottyPrecision = atoi(argv[++i]);	// -1 means "no edge labels"
		} else if (string(argv[i]) == "--time") {
			timeActions = true;
		} else if (string(argv[i]) == "--DEBUG-dompath2-stop-after") {
			DEBUGreduceDomPath2StopAfterN = atoi(argv[++i]);
		} else if (string(argv[i]) == "-i") {
			inputFName = argv[++i];			//HACK: Needed for debugging under VS++ 2008, since it doesn't support redirection with debugging even though http://msdn.microsoft.com/en-us/library/kcw4dzyf%28v=vs.90%29.aspx says it does!
		} else if (string(argv[i]) == "-o") {
			outputFName = argv[++i];			//HACK: Needed for debugging under VS++ 2008, since it doesn't support redirection with debugging even though http://msdn.microsoft.com/en-us/library/kcw4dzyf%28v=vs.90%29.aspx says it does!
		} else if (argv[i][0] != '-') {
			// A "command" in the "program" describing what to do.  Add it to the program strig to be parsed later.
			if (!programString.empty()) {
				programString += " ";
			}
			programString += argv[i];
		} else {
			cerr << "Unrecognised argument '" << argv[i] << "'.\n";
			return 1;
		}
	}
	
	cerr << "About to parse the following program: <" << programString << ">.\n";
	istringstream iss(programString);
	Command* prog = Command::parse(iss);
	if (prog) {
		cerr << "About to run the program...\n";
		prog->run();
	} else {
		cerr << "No program!\n";
	}
	
	Command::destroy(prog);
	
	//int nIterations;
	//clock_t elapsed;
	//int nEdgesDeleted;
	//if (timeActions) {
	//	//HACK: First we need to be able to convert back from function pointers to names...  Yes, this is thoroughly weird...
	//	map<void (*)(), string> nameFor;
	//	for (map<string, void (*)()>::const_iterator i(actions.begin()); i != actions.end(); ++i) {
	//		nameFor[(*i).second] = (*i).first;			// If it happens that more than one action name maps to the same function, that function will get the last name.
	//	}
	//	
	//	int oldPrec = cerr.precision();		// Ugh, I hate iostreams
	//	cerr.precision(4);
	//	
	//	double total = 0;
	//	cerr << "Iters\tTotSecs\tAvgSecs\tAction\n";		// Put action last to preserve columns as much as possible
	//	for (map<void (*)(), pair<int, clock_t> >::const_iterator i(actionTimes.begin()); i != actionTimes.end(); ++i) {
	//		double t = static_cast<double>((*i).second.second) / CLOCKS_PER_SEC;
	//		total += t;
	//		cerr << (*i).second.first << '\t' << t << '\t' << (t / (*i).second.first) << '\t' << nameFor[(*i).first] << '\n';
	//	}
	//	
	//	cerr << '\t' << total << "\t\t===== Total =====\n";
	//	cerr.precision(oldPrec);
	//}
	if (timeActions) {
		//HACK: First we need to be able to convert back from function pointers to names...  Yes, this is thoroughly weird...
		map<void (*)(), string> nameFor;
		for (map<string, void (*)()>::const_iterator i(actions.begin()); i != actions.end(); ++i) {
			nameFor[(*i).second] = (*i).first;			// If it happens that more than one action name maps to the same function, that function will get the last name.
		}
		
		int oldPrec = cerr.precision();		// Ugh, I hate iostreams
		cerr.precision(4);
		
		double total = 0;
		int totalEdgesDeleted = 0;
		cerr << "Iters\tTotSecs\tAvgSecs\tEdgeDel\tAction\n";		// Put action last to preserve columns as much as possible
		//for (map<void (*)(), pair<int, clock_t> >::const_iterator i(actionTimes.begin()); i != actionTimes.end(); ++i) {
		for (map<void (*)(), action_info>::const_iterator i(actionTimes.begin()); i != actionTimes.end(); ++i) {
			double t = static_cast<double>((*i).second.elapsed) / CLOCKS_PER_SEC;
			total += t;
			totalEdgesDeleted += (*i).second.nEdgesDeleted;
			cerr << (*i).second.nIterations << '\t' << t << '\t' << (t / (*i).second.nIterations) << '\t' << (*i).second.nEdgesDeleted << '\t' << nameFor[(*i).first] << '\n';
		}
		
		cerr << '\t' << total << "\t\t" << totalEdgesDeleted << "\t===== Total =====\n";
		cerr.precision(oldPrec);
	}
	
	return 0;
}
