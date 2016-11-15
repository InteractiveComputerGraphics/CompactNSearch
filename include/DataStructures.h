

#ifndef __DATA_STRUCTURES_H__
#define __DATA_STRUCTURES_H__

#include "Config.h"

#include <atomic>
#include <vector>

namespace CompactNSearch
{
struct PointID
{
	unsigned int point_set_id;
	unsigned int point_id;

	bool operator==(PointID const& other) const
	{
		return point_id == other.point_id && point_set_id == other.point_set_id;
	}
};

struct HashKey
{
	HashKey() = default;
	HashKey(int i, int j, int k)
		: k{i, j, k}
	{
	}

	HashKey& operator=(HashKey const& other)
	{
		k[0] = other.k[0];
		k[1] = other.k[1];
		k[2] = other.k[2];
		return *this;
	}

	bool operator==(HashKey const& other) const
	{
		return 
			k[0] == other.k[0] &&
			k[1] == other.k[1] && 
			k[2] == other.k[2];
	}

	bool operator!=(HashKey const& other) const
	{
		return !(*this == other);
	}

	int k[3];
};

struct HashEntry
{
	HashEntry() 
	{
		indices.reserve(INITIAL_NUMBER_OF_INDICES);
	}

	HashEntry(PointID const& id)
	{
		add(id);
	}

	void add(PointID const& id)
	{
		indices.push_back(id);
	}

	void erase(PointID const& id)
	{
		indices.erase(std::find(indices.begin(), indices.end(), id));
	}

	unsigned int n_indices() const
	{
		return static_cast<unsigned int>(indices.size());
	}

	std::vector<PointID> indices;
};

struct SpatialHasher
{
	std::size_t operator()(HashKey const& k) const
	{
		return 
			73856093 * k.k[0] ^
			19349663 * k.k[1] ^
			83492791 * k.k[2];
	}
};

class Spinlock
{
public:

	void lock()
	{
		while (m_lock.test_and_set(std::memory_order_acquire));
	}

	void unlock()
	{
		m_lock.clear(std::memory_order_release);
	}

	Spinlock() = default;
	Spinlock(Spinlock const& other) {};
	Spinlock& operator=(Spinlock const& other) { return *this; }


private:

	std::atomic_flag m_lock = ATOMIC_FLAG_INIT;
};
}

#endif // __DATA_STRUCTURES_H__
