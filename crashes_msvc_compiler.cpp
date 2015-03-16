// Running
//    cl crashes_msvc_compiler.cpp
// produces a "Microsoft C/C++ Compiler Driver has encountered a problem and needs to close" box.
// Compiler version: Microsoft (R) 32-bit C/C++ Optimizing Compiler Version 15.00.21022.08 for 80x86
// Running on: my laptop (Core2 Duo, 2Gb RAM, WinXP SP3)

#include <vector>

using namespace std;

vector<vector<double> > m;

void f() {
	// DOESN'T CRASH
	//m.resize(1, vector<int>(1, -42.123));
	//m.resize(1, vector<int>(1, -1.7976931348623158e+307));
	//m.resize(1, vector<int>(1, 1.7976931348623158e+308));
	//m.resize(1, vector<int>(1, 1e+308));
	
	// CRASHES
	//m.resize(1, vector<int>(1, -DBL_MAX));		// Requires #include <cfloat>
	//m.resize(1, vector<int>(1, -1.7976931348623158e+308));
	//m.resize(1, vector<int>(1, -1.0976931348623158e+308));
	m.resize(1, vector<int>(1, -1e+308));
}
