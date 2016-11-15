#ifndef __CONFIG_H__
#define __CONFIG_H__

namespace CompactNSearch
{
    using Real = double;
}

#define INITIAL_NUMBER_OF_INDICES 50

#ifdef _MSC_VER
	#include <ppl.h>
#else
	#include <parallel/algorithm>
#endif


#endif
