#pragma once

namespace CompactNSearch
{
	using Real = double;
}

#define INITIAL_NUMBER_OF_INDICES   50
#define INITIAL_NUMBER_OF_NEIGHBORS 50

#ifdef _MSC_VER
	#include <ppl.h>
#else
	#include <parallel/algorithm>
#endif
