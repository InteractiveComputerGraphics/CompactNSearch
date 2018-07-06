#pragma once

namespace CompactNSearch
{
#ifdef USE_DOUBLE
	using Real = double;
#else
	using Real = float;
#endif
}

#define INITIAL_NUMBER_OF_INDICES   50
#define INITIAL_NUMBER_OF_NEIGHBORS 50

#ifdef _MSC_VER
	#include <ppl.h>
#else
	#include <parallel/algorithm>
#endif
