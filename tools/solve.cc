#define DEBUG

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/Stopwatch.h>
#include <ogdf/basic/HashArray.h>

// configuration for the separation oracle
#define ORIGINAL_GOLDBERG
#ifdef ORIGINAL_GOLDBERG
#  include "./MincutPushRelabel.h"
#  undef USE_BACK_CUTS
#else
#  undef USE_EDMONDS_KARP
#  ifdef USE_EDMONDS_KARP
#    include "../MaxFlowEdmondsKarp.h"
#  else
#    include "../MaxFlowGoldbergTarjan.h"
#  endif
#endif
#undef USE_CARDINALITY_HEURISTIC
#undef OMIT_SOME_FLOWS // also try with this being defined

// set LP_SOLVER_GUROBI to use GUROBI, otherwise CPLEX is used
#ifdef LP_SOLVER_GUROBI
#  include <gurobi_c++.h>
#else
#  include <ilcplex/ilocplex.h>
#endif

#define SILENT_LP_SOLVER
#undef SOLVE_RELAXATION

#include <sstream>
#include <cstdlib> // rand
#include <map>
#include <set>
#include <iterator>		// ostream_iterator
#include <algorithm>	// copy()

using namespace ogdf;
using std::map;
using std::set;
using std::pair;
using std::flush;
using std::copy;
using std::ostream_iterator;

#include "WeightedNodeColoredGraph.h"

class MinCut
{
	const Graph *m_G;
#ifdef ORIGINAL_GOLDBERG
	NodeArray<int> m_nodeIDs;
	MincutPushRelabel *m_mincut;
	int *m_cutset;
	List<edge> m_frontEdges;
	int m_last;
#else
	MaxFlowModule<double> *m_mincut;
#endif
	int m_count;

public:
	MinCut()
	{
	}

	MinCut(const Graph &G)
	{
		init(G);
	}

	void
	init(const Graph &G)
	{
		m_G = &G;
		m_mincut = NULL;
#ifdef ORIGINAL_GOLDBERG
		m_nodeIDs.init(G);
		int i = 0;
		for (node v = G.firstNode(); v; v = v->succ()) {
			m_nodeIDs[v] = i++;
		}
		m_cutset = new int[m_G->numberOfNodes()];
		m_last = -1;
#endif
		m_count = 0;
	}

	~MinCut()
	{
		if (m_mincut) {
			delete m_mincut;
		}
#ifdef ORIGINAL_GOLDBERG
		delete[] m_cutset;
#endif
	}

	double
	call(const EdgeArray<double> &cap, node source, node target, double epsilon)
	{
		++m_count;
#ifdef ORIGINAL_GOLDBERG
		if (!m_mincut) {
			m_mincut = new MincutPushRelabel(*m_G, m_nodeIDs, cap, source, target);
		} else {
			m_mincut->update(cap, m_nodeIDs[source], m_nodeIDs[target]);
		}
		return m_mincut->min_cut(1, m_cutset, epsilon);
#else // ORIGINAL_GOLDBERG
		if (m_mincut) {
			delete m_mincut;
		}
# ifdef USE_EDMONDS_KARP
		m_mincut = new MaxFlowEdmondsKarp<double>(*m_G, cap, source, target, epsilon);
# else
		m_mincut = new MaxFlowGoldbergTarjan<double>(*m_G, cap, source, target, epsilon);
# endif
		return m_mincut->flow();
#endif // ORIGINAL_GOLDBERG
	}

#ifdef ORIGINAL_GOLDBERG
	int
	cutset(node v) const
	{
		return m_cutset[m_nodeIDs[v]];
	}
#endif

	const List<edge> &
	frontEdges()
	{
#ifdef ORIGINAL_GOLDBERG
		if (m_last < m_count) {
			m_frontEdges.clear();
			m_last = m_count;
			for (edge e = m_G->firstEdge(); e; e = e->succ()) {
				if (cutset(e->source()) != 2
				 && cutset(e->target()) == 2) {
					m_frontEdges.pushBack(e);
				}
			}
		}
		return m_frontEdges;
#else
		return m_mincut->cutSourceEdges();
#endif
	}

#ifdef USE_BACK_CUTS
	// only available for ORIGINAL_GOLDBERG
	double
	backEdgeVal(List<edge> &list, const EdgeArray<double> &cap)
	{
		double cutVal = 0;
		for (edge e = m_G->firstEdge(); e; e = e->succ()) {
			if (cutset(e->source()) == 1
			 && cutset(e->target()) != 1) {
				cutVal += cap[e];
				list.pushBack(e);
			}
		}
		if (list.empty()) {
			return 1;
		}
		return cutVal;
	}
#endif

	int
	count() const
	{
		return m_count;
	}
}; // class MinCut

#ifdef SOLVE_RELAXATION
class CallbackSplitAndJoin
#else
class CallbackSplitAndJoin
#  ifndef LP_SOLVER_GUROBI
   : public IloCplex::UserCutCallbackI
#  else
   : public GRBCallback
#  endif
#endif
{
#ifndef LP_SOLVER_GUROBI
	IloEnv env;
#  ifdef SOLVE_RELAXATION
	IloCplex *cplex;
#  endif
#else
	GRBModel &model;
#endif
	int separationCount;
	const WeightedNodeColoredGraph &G;
	const node root;
#ifndef LP_SOLVER_GUROBI
	IloBoolVarArray &xVars;
	const EdgeArray<IloInt> &index;
#else
	const EdgeArray<GRBVar> &xVars;
#endif
	MinCut minCut;
	const double eps;

public:
#ifdef SOLVE_RELAXATION
#  ifndef LP_SOLVER_GUROBI
	CallbackSplitAndJoin(IloCplex *_cplex, const WeightedNodeColoredGraph &_G, const node _root, IloBoolVarArray &_xVars, const EdgeArray<IloInt> &_index, double _eps)
	  : env(_cplex->getEnv())
	  , cplex(_cplex)
#  else
	CallbackSplitAndJoin(GRBModel &_model, const WeightedNodeColoredGraph &_G, const node _root, EdgeArray<GRBVar> &_xVars, double _eps)
	  : model(_model)
#  endif
	  ,
#else
#  ifndef LP_SOLVER_GUROBI
	IloCplex::CallbackI* duplicateCallback() const
	{
		return (new (env) CallbackSplitAndJoin(*this));
	}

	CallbackSplitAndJoin(IloEnv _env, const WeightedNodeColoredGraph &_G, const node _root, IloBoolVarArray &_xVars, const EdgeArray<IloInt> &_index, double _eps)
	  : IloCplex::UserCutCallbackI(_env)
	  , env(_env)
	  ,
#  else
	CallbackSplitAndJoin(const WeightedNodeColoredGraph &_G, const node _root, const EdgeArray<GRBVar> &_xVars, double _eps)
	  :
#  endif
#endif
	    separationCount(0)
	  , G(_G)
	  , root(_root)
	  , xVars(_xVars)
#ifndef LP_SOLVER_GUROBI
	  , index(_index)
#endif
	  , minCut(G.theGraph())
	  , eps(_eps)
	{
	}

#ifdef DEBUG
	~CallbackSplitAndJoin()
	{
		cout << "  ** number of max flow computations: " << minCut.count() << "\n";
	}
#endif

	// split-and-join constraint (separated!)
	// for all colors c and all all S \subseteq V with V_c \subseteq S: sum_(in-edge e of V_c) x_e <= sum_(in-edge e of S){x_e}
	void addSplitJoinConstraint(const List<edge> &inEdges, int c);

	bool separateMinCut(const node v, const EdgeArray<double> &capacity, double flowEpsilon, double cardEpsilon);

#ifndef LP_SOLVER_GUROBI
	bool doCallback(const IloNumArray &sol);
#else
	bool doCallback(const EdgeArray<double> &sol);
#endif

#ifdef SOLVE_RELAXATION
	bool call()
	{
#  ifndef LP_SOLVER_GUROBI
		cout << "  -* Separation " << separationCount++ << " result: " << cplex->getObjValue() << "\n";
		IloNumArray sol(env);
		cplex->getValues(sol, xVars);
#  else
		cout
		  << "  -* Separation " << separationCount++
		  << " result: " << -model.get(GRB_DoubleAttr_ObjVal) << "\n";
		EdgeArray<double> sol(G.theGraph());
		for (edge e = G.firstEdge(); e; e = e->succ()) {
			sol[e] = xVars[e].get(GRB_DoubleAttr_X);
		}
#  endif
		return doCallback(sol);
	}
#else
#  ifndef LP_SOLVER_GUROBI
	void main()
	{
		if (isAfterCutLoop()) {
			// nothing to separate after cut loop
			return;
		}
		cout << "  -* Separation " << separationCount++ << " result: " << getObjValue() << "\n";
		IloNumArray sol(env);
		getValues(sol, xVars);
		doCallback(sol);
	}
#  else
	void callback()
	{
		if (where == GRB_CB_MIPNODE
		 && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL) {
			cout
			  << "  -* Separation " << separationCount++
			  << " result: " << getDoubleInfo(GRB_CB_MIPNODE_OBJBND)
			  << " [best: " << getDoubleInfo(GRB_CB_MIPNODE_OBJBST)
			  << ", feasible: " << getIntInfo(GRB_CB_MIPNODE_SOLCNT)
			  << ", #nodes: " << getDoubleInfo(GRB_CB_MIPNODE_NODCNT)
			  << "]\n";
			EdgeArray<double> sol(G.theGraph());
			for (edge e = G.firstEdge(); e; e = e->succ()) {
				sol[e] = getNodeRel(xVars[e]);
			}
			doCallback(sol);
		}
	}
#  endif
#endif
}; // class CallbackSplitAndJoin

#if (!defined LP_SOLVER_GUROBI) && (!defined SOLVE_RELAXATION)
IloCplex::Callback callbackSplitAndJoin(IloEnv env, const WeightedNodeColoredGraph &G, const node root, IloBoolVarArray &xVars, const EdgeArray<IloInt> &index, double eps)
{
	return (IloCplex::Callback(new (env) CallbackSplitAndJoin(env, G, root, xVars, index, eps)));
}
#endif

void
CallbackSplitAndJoin::addSplitJoinConstraint(const List<edge> &inEdges, int color)
{
	OGDF_ASSERT(!inEdges.empty());
#ifndef LP_SOLVER_GUROBI
	IloExpr expr(env);
	forall_listiterators(edge, it, inEdges) {
		expr += xVars[index[*it]];
	}
#else
	GRBLinExpr expr = 0.0;
	forall_listiterators(edge, it, inEdges) {
		expr += xVars[*it];
	}
#endif
	for (set<node>::const_iterator it = G.colorSet(color).begin(); it != G.colorSet(color).end(); ++it) {
		const node v = *it;
		for (adjEntry adj = v->firstAdj(); adj; adj = adj->succ()) {
			const edge e = adj->theEdge();
			if (e->target() == v) { // incoming edges
#ifndef LP_SOLVER_GUROBI
				expr -= xVars[index[e]];
#else
				expr -= xVars[e];
#endif
			}
		}
	}
#ifdef SOLVE_RELAXATION
#  ifndef LP_SOLVER_GUROBI
	cplex->getModel().add(expr >= 0);
#  else
	model.addConstr(expr >= 0);
#  endif
#else
#  ifndef LP_SOLVER_GUROBI
	add(expr >= 0).end();
#  else
	addCut(expr >= 0);
#  endif
#endif
}

bool
CallbackSplitAndJoin::separateMinCut(const node target, const EdgeArray<double> &capacity, double flowEpsilon, double cardEpsilon)
{
	double maximum(0);
	double threshold(0);
	int color = SUPERTARGETCOLOR(G.color(target));
	for (adjEntry adjSuper = target->firstAdj(); adjSuper; adjSuper = adjSuper->succ()) {
		OGDF_ASSERT(adjSuper->theEdge()->target() == target); // all edges are incoming
		const node w = adjSuper->twinNode();
		OGDF_ASSERT(G.color(w) == color);
		for (adjEntry adj = w->firstAdj(); adj; adj = adj->succ()) {
			const edge e = adj->theEdge();
			if (e->target() == w) { // incoming edge
#ifdef USE_CARDINALITY_HEURISTIC
				const double cap = capacity[e] - cardEpsilon;
				threshold += cap;
#  ifdef OMIT_SOME_FLOWS
				maximum = max(maximum, cap);
#  endif // OMIT_SOME_FLOWS
# else
				threshold += capacity[e];
#  ifdef OMIT_SOME_FLOWS
				maximum = max(maximum, capacity[e]);
#  endif // OMIT_SOME_FLOWS
#endif
			}
		}
	}
	if (threshold <= maximum + eps) {
		return false;
	}

	OGDF_ASSERT(threshold < 1 + eps);

	double cutVal = minCut.call(capacity, root, target, flowEpsilon);
#ifdef USE_CARDINALITY_HEURISTIC
	cutVal -= minCut.frontEdges().size() * cardEpsilon;
#endif
	if (cutVal < threshold - eps
	 && !minCut.frontEdges().empty()) {
		addSplitJoinConstraint(minCut.frontEdges(), color);
#ifdef USE_BACK_CUTS
		List<edge> edgeList;
		cutVal = minCut.backEdgeVal(edgeList, capacity);
#  ifdef USE_CARDINALITY_HEURISTIC
		cutVal -= edgeList.size() * cardEpsilon;
#  endif // USE_CARDINALITY_HEURISTIC
		if (cutVal < threshold - eps) {
			addSplitJoinConstraint(edgeList, color);
		}
#endif // USE_BACK_CUTS
		return true;
	}

	return false;
}

bool
#ifndef LP_SOLVER_GUROBI
CallbackSplitAndJoin::doCallback(const IloNumArray &sol)
#else
CallbackSplitAndJoin::doCallback(const EdgeArray<double> &sol)
#endif
{
	int cutsFound = 0;
	const double inputEpsilon = 1e-16;

	EdgeArray<double> capacity(G.theGraph());
#ifdef USE_CARDINALITY_HEURISTIC
	double minCap = numeric_limits<double>::max();
#endif
	for (edge e = G.firstEdge(); e; e = e->succ()) {
#ifndef LP_SOLVER_GUROBI
		const double cap = sol[index[e]];
#else
		const double cap = sol[e];
#endif
#ifdef USE_CARDINALITY_HEURISTIC
		if (cap >= inputEpsilon
		 && cap < minCap) {
			minCap = cap;
		}
#endif
		capacity[e] = cap;
	}
#ifdef USE_CARDINALITY_HEURISTIC
	const double cardEpsilon = minCap / G.numberOfEdges();
	for (edge e = G.firstEdge(); e; e = e->succ()) {
		capacity[e] += cardEpsilon;
	}
	const double flowEpsilon = cardEpsilon / 100;
#else
	const double flowEpsilon = inputEpsilon, cardEpsilon = inputEpsilon;
#endif

	forall_listiterators(node, it, G.superTargets()) {
		const node v = *it;
		for (adjEntry adjSuper = v->firstAdj(); adjSuper; adjSuper = adjSuper->succ()) {
			OGDF_ASSERT(adjSuper->theEdge()->target() == v); // all edges are incoming
#ifdef USE_CARDINALITY_HEURISTIC
			capacity[adjSuper->theEdge()] = 1 + cardEpsilon;
#else
			capacity[adjSuper->theEdge()] = 1;
#endif
		}
		if (separateMinCut(v, capacity, flowEpsilon, cardEpsilon)) {
			++cutsFound;
		}
	}

#ifdef DEBUG
	cout << "  ** Cuts found (by MinCut): " << cutsFound << "\n";
#endif
	return (cutsFound > 0);
}

class LinearProgram
{
public:
	typedef enum {
		LP_ORIGINAL = 0,
		LP_CUT = 4,
	} lp_type;

protected:
	const WeightedNodeColoredGraph &G;
	const node root;
	const lp_type lp;
#ifndef LP_SOLVER_GUROBI
	IloEnv env;
	IloModel model;
	IloCplex cplex;
	IloBoolVarArray xVars;
	EdgeArray<IloInt> index;
	IloCplex::UserCutCallbackI *callbackImplSplitAndJoin;
#else
	GRBEnv env;
	GRBModel model;
	EdgeArray<GRBVar> xVars;
	GRBCallback *callbackImplSplitAndJoin;
#endif

	// add objective (without root node weight)
	void addObjective();

	// color constraints:
	// for all colors c: sum_(in-edges e into nodes with color c){x_e} <= 1
	void addColorConstraint(const set<node> &nodeset);
	void addIndegree1Constraint(node v) {
		set<node> nodeset;
		nodeset.insert(v);
		addColorConstraint(nodeset);
	}

	// connectivity constraint (by out-edge e)
	// for all out-edges e of v != root: x_e <= sum_(in-edge f of v){x_f}
	void addConnectivityConstraint(edge e);

	void addConstraints();

public:
	LinearProgram(const WeightedNodeColoredGraph &wcG, const node _root, lp_type _lp, int timeLimitSecs, int memLimitMb, char* workDir)
	 : G(wcG)
	 , root(_root)
	 , lp(_lp)
	 , env()
	 , model(env)
#ifndef LP_SOLVER_GUROBI
	 , cplex(model)
	 , xVars(env, G.numberOfEdges())
	 , index(G.theGraph())
#else
	 , xVars(G.theGraph())
#endif
	 , callbackImplSplitAndJoin(NULL)
	{
#ifndef LP_SOLVER_GUROBI
#  ifdef SILENT_LP_SOLVER
		cplex.setOut(env.getNullStream());
#  else
		cplex.setOut(cerr);
#  endif
		cplex.setParam(IloCplex::IntParam::Threads, 1);
		if (timeLimitSecs != -1) {
			cplex.setParam(IloCplex::TiLim, timeLimitSecs);
		}
		if (memLimitMb != -1) {
			cplex.setParam(IloCplex::WorkMem, memLimitMb);
		}
		if (workDir) {
			cplex.setParam(IloCplex::WorkDir, workDir);
		}
#  if 0
		cplex.setParam(IloCplex::MIPEmphasis, CPX_MIPEMPHASIS_OPTIMALITY);
		cplex.setParam(IloCplex::RootAlg, CPX_ALG_DUAL);
		cplex.setParam(IloCplex::NodeAlg, CPX_ALG_DUAL);
#  endif
#  ifndef SOLVE_RELAXATION
		if (lp >= LP_CUT) {
			cplex.setParam(IloCplex::MIPSearch, CPX_MIPSEARCH_TRADITIONAL);
			callbackImplSplitAndJoin = new (env) CallbackSplitAndJoin(env, G, root, xVars, index, cplex.getParam(IloCplex::EpRHS));
			cplex.use(IloCplex::Callback(callbackImplSplitAndJoin));
		}
#  endif
#else
		GRBEnv menv = model.getEnv();
#  ifdef SILENT_LP_SOLVER
		menv.set(GRB_IntParam_OutputFlag, 0);
#  else
		menv.set(GRB_IntParam_OutputFlag, 1);
#  endif
		menv.set(GRB_IntParam_Threads, 1);
		menv.set(GRB_IntParam_Method, 1); // use dual simplex
		if (timeLimitSecs != -1) {
			assert(!"Don't support setting a time limit for Gurobi yet!");
		}
		if (memLimitMb != -1) {
			assert(!"Don't support setting a memory limit for Gurobi yet!");
		}
		if (workDir) {
			assert(!"Don't support setting a working directory for Gurobi yet!");
		}
		// TODO: try GRB_IntParam_MIPFocus \in {1,2,3}
#  ifndef SOLVE_RELAXATION
		if (lp >= LP_CUT) {
			menv.set(GRB_IntParam_PreCrush, 1);
			callbackImplSplitAndJoin = new CallbackSplitAndJoin(G, root, xVars, 1e-6);
			model.setCallback(callbackImplSplitAndJoin);
		}
#  endif
#endif // !LP_SOLVER_GUROBI

		addObjective();
		addConstraints();
	}

	~LinearProgram()
	{
		if (callbackImplSplitAndJoin) {
			delete callbackImplSplitAndJoin;
		}
#ifndef LP_SOLVER_GUROBI
		env.end();
#endif
	}

	// solve LP
	double solve();

	void writeSolution(char* fname);
};

void
LinearProgram::addObjective()
{
	edge e;
#ifndef LP_SOLVER_GUROBI
	IloExpr expr(env);
	int i = 0;
#endif
	forall_edges(e, G) {
		std::stringstream ss;
		ss << "x_" << e->source() << "_" << e->target();
#ifndef LP_SOLVER_GUROBI
#  ifdef SOLVE_RELAXATION
		model.add(IloConversion(env, xVars[i], ILOFLOAT));
#  endif
		xVars[i].setName(ss.rdbuf()->str().c_str());
		expr += G.weight(e) * xVars[i];
		index[e] = i++;
#else
#  ifdef SOLVE_RELAXATION
		xVars[e] = model.addVar(0.0, 1.0, -G.weight(e), GRB_CONTINUOUS, ss.rdbuf()->str().c_str());
#  else
		xVars[e] = model.addVar(0.0, 1.0, -G.weight(e), GRB_BINARY, ss.rdbuf()->str().c_str());
#  endif
#endif
	}
#ifndef LP_SOLVER_GUROBI
	model.add(IloMaximize(env, expr));
#else
	model.update();
#endif
}

void
LinearProgram::addConstraints()
{
	for (node v = G.firstNode(); v; v = v->succ()) {
		addIndegree1Constraint(v);
		adjEntry adj;
		forall_adj(adj, v) {
			const edge e = adj->theEdge();
			const node w = e->target();
			if (w != v) { // for all outgoing edges
				addConnectivityConstraint(e);
			}
		}
	}

	// add color constraints
	const ColorSetMap colorSetMap = G.colorSetMap();
	for (ColorSetMap::const_iterator _nodeList = colorSetMap.begin(); _nodeList.valid(); ++_nodeList) {
		addColorConstraint(_nodeList.info());
	}
}

void
LinearProgram::addColorConstraint(const set<node> &nodeset)
{
	int count = 0;
#ifndef LP_SOLVER_GUROBI
	IloExpr expr(env);
#else
	GRBLinExpr expr = 0;
#endif
	for (set<node>::const_iterator _v = nodeset.begin(); _v != nodeset.end(); _v++) {
		node v = *_v;
		if (!v) {
			continue;
		}
		adjEntry adj;
		forall_adj(adj, v) {
			edge e = adj->theEdge();
			if (e->target() == v) { // is an incoming edge
#ifndef LP_SOLVER_GUROBI
				expr += xVars[index[e]];
#else
				expr += xVars[e];
#endif
				++count;
			}
		}
	}
	if (count >= 2) { // do not add useless constraints like 0 <= 1 or one variable <= 1
#ifndef LP_SOLVER_GUROBI
		model.add(expr <= 1);
#else
		model.addConstr(expr <= 1);
#endif
	}
}

void
LinearProgram::addConnectivityConstraint(edge e)
{
	const node v = e->source();
	if (v == root) {
		return;
	}

#ifndef LP_SOLVER_GUROBI
	IloExpr expr(env);
	expr += xVars[index[e]];
#else
	GRBLinExpr expr = xVars[e];
#endif
	adjEntry adjIn;
	forall_adj(adjIn, v) {
		edge eIn = adjIn->theEdge();
		if (eIn->target() == v) { // eIn is in-edge
#ifndef LP_SOLVER_GUROBI
			expr -= xVars[index[eIn]];
#else
			expr -= xVars[eIn];
#endif
		}
	}
#ifndef LP_SOLVER_GUROBI
	model.add(expr <= 0);
#else
	model.addConstr(expr <= 0);
#endif
}

double
LinearProgram::solve()
{
#ifndef LP_SOLVER_GUROBI
	try {
#  ifdef SOLVE_RELAXATION
		if (lp >= LP_CUT) {
			CallbackSplitAndJoin cb(&cplex, G, root, xVars, index, cplex.getParam(IloCplex::EpRHS));
			do {
				if (!cplex.solve()) {
					cerr
					  << "Failed to optimize LP! ("
					  << cplex.getCplexStatus() << ", "
					  << cplex.getCplexSubStatus() << ")\n";
					throw(-1);
				}
			} while (cb.call());
		} else
#  endif
		{
			if (!cplex.solve()) {
				cerr
				  << "Failed to optimize LP! ("
				  << cplex.getCplexStatus() << ", "
				  << cplex.getCplexSubStatus() << ")\n";
				throw(-1);
			}
		}
	} catch (IloCplex::Exception &exception)
#else
	try {
#  ifdef SOLVE_RELAXATION
		if (lp >= LP_CUT) {
			CallbackSplitAndJoin cb(model, G, root, xVars, 1e-6);
			do {
				model.optimize();
			} while (cb.call());
		} else
#endif
		{
			model.optimize();
		}
	} catch (GRBException exception)
#endif
	{
		cerr << "Failed to solve model!\n" << exception.getMessage();
		throw(-1);
	}

#ifndef LP_SOLVER_GUROBI
	return cplex.getObjValue();
#else
	return -model.get(GRB_DoubleAttr_ObjVal);
#endif
}

void
LinearProgram::writeSolution(char* fname) {
#ifndef LP_SOLVER_GUROBI
	cplex.writeSolution(fname);
#else
	assert(!"Gurobi does not yet support writing the solution to a file!");
#endif
}

static void usage(char *prog)
{
	cout
	  << "Usage: "
	  << prog << " [options] <filename>\n\n"
#ifdef SOLVE_RELAXATION
	"+++ NOTE: this tool only solves the *relaxation*, not the ILP +++\n\n"
#endif
	"Options are:\n"
	"  -l <mode>   choose the linear program [default: c]\n"
	"  -T <secs>   set the time limit in seconds\n"
	"  -M <mb>     set the memory limit in MB\n"
	"  -D <path>   set the working directory to path (CPLEX only)\n"
	"\n"
	"Available linear programs:\n"
	"   o   Original (colorful forest and connectivity constraints)\n"
	"   c   Cut (Original + Cut constraints)\n"
	;
}

int
main(int argc, char *argv[])
{
	int argi = 1;
	LinearProgram::lp_type lp = LinearProgram::LP_CUT;
	char* solutionFName = NULL;		// Setting this to non-NULL causes the solution to be written to this file.
	int timeLimitSecs = -1;
	int memLimitMb = -1;
	char* workDir = NULL;
	while (argv[argi] && argv[argi][0] == '-') {
		switch (argv[argi][1]) {
		case 'l':
			++argi;
			switch (argv[argi][0]) {
			case 'O': case 'o':
				lp = LinearProgram::LP_ORIGINAL;
				break;
			case 'C': case 'c':
				lp = LinearProgram::LP_CUT;
				break;
			default:
				usage(argv[0]);
				return -1;
			}
			break;
		case 'o':
			solutionFName = argv[++argi];
			break;
		case 'T':
			timeLimitSecs = atoi(argv[++argi]);
			break;
		case 'M':
			memLimitMb = atoi(argv[++argi]);
			break;
		case 'D':
			workDir = argv[++argi];
			break;
		default:
			usage(argv[0]);
			return -1;
		}
		++argi;
	}

	if (argc - argi != 1) {
		usage(argv[0]);
		return -1;
	}

	cout << "solve called as: ";
	copy(argv, argv + argc, ostream_iterator<char*>(cout, " "));		// Will append a trailing space, but who cares
	cout << "\n" << flush;				// flush to ensure that we see at least this output if the solver crashes (e.g. runs out of memory).  Yes, I could write endl but this is clearer!
	
	WeightedNodeColoredGraph G;
	cout << argv[argi] << ":\n" << flush;
	if (!G.readMCAfile(argv[argi])) {
		cerr << "Error reading file (perhaps check number of nodes and number of edges).\n";
		return 1;
	}
	cout
	 << "  -- "
	 << G.numberOfNodes() << " nodes, "
	 << G.numberOfEdges() << " edges, "
	 << G.numberOfColors() << " colors\n";

	cout << "Mode: " << (lp == LinearProgram::LP_ORIGINAL ? "Original\n" : "Cut\n");
	if (solutionFName) {
		cout << "Solution will be written to: " << solutionFName << "\n";
	} else {
		cout << "Solution will not be written to any file.\n";
	}
	cout << "Time limit: ";
	if (timeLimitSecs != -1) {
		cout << timeLimitSecs << "s";
	} else {
		cout << "None (default)";
	}
	cout << "\nMemory limit: ";
	if (memLimitMb != -1) {
		cout << memLimitMb << "MB";
	} else {
		cout << "None (default)";
	}
	cout << "\nWorking directory: ";
	if (workDir) {
		cout << workDir;
	} else {
		cout << "None (default)";
	}
	cout << "\n" << flush;		// flush to ensure that we see at least this output if the solver crashes (e.g. runs out of memory).  Yes, I could write endl but this is clearer!

	StopwatchWallClock timer;
	timer.start();

	if (lp >= LinearProgram::LP_CUT) {
		G.addSuperTargets();
	}

	LinearProgram linprog(G, G.firstNode(), lp, timeLimitSecs, memLimitMb, workDir);

	cout << "  -- initialization time: " << timer.milliSeconds()*1e-3 << " s\n";

	double optSol = linprog.solve();

	timer.stop();

	cout
	  << "\t" << optSol << "\n"
	  << "  -- total time used to solve problem: "
	  << timer.milliSeconds()*1e-3 << " seconds\n";

	if (solutionFName) {
		cout << "Writing out solution to " << solutionFName << "...\n";
		linprog.writeSolution(solutionFName);
	}

	return 0;
}
