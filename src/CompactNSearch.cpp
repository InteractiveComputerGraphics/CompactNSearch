
#include <CompactNSearch>
#include "../extern/libmorton/libmorton/include/morton.h"

#include <iostream>
#include <numeric>
#include <array>
#include <cstdint>
#include <limits>
#include <algorithm>


namespace CompactNSearch
{
namespace
{
// Determines Morten value according to z-curve.
inline uint_fast64_t
z_value(HashKey const& key)
{
	return morton3D_64_encode(
			static_cast<uint_fast32_t>(static_cast<int64_t>(key.k[0]) -
									  (std::numeric_limits<int>::lowest() + 1)),
			static_cast<uint_fast32_t>(static_cast<int64_t>(key.k[1]) -
									  (std::numeric_limits<int>::lowest() + 1)),
			static_cast<uint_fast32_t>(static_cast<int64_t>(key.k[2]) -
									  (std::numeric_limits<int>::lowest() + 1))
	);
}
}


NeighborhoodSearch::NeighborhoodSearch(Real r, bool erase_empty_cells)
	: m_r2(r * r), m_inv_cell_size(1.0 / r)
	, m_erase_empty_cells(erase_empty_cells)
	, m_initialized(false)
{
	if (r <=  0.0)
	{
		std::cerr << "WARNING: Neighborhood search may not be initialized with a zero or negative"
			<< " search radius. This may lead to unexpected behavior." << std::endl;
	}
}

// Computes triple index to a world space position x.
HashKey
NeighborhoodSearch::cell_index(Real const* x) const
{
	HashKey ret;
	for (unsigned int i = 0; i < 3; ++i)
	{
		if (x[i] >= 0.0) ret.k[i] = static_cast<int>(m_inv_cell_size * x[i]);
		else ret.k[i] = static_cast<int>(m_inv_cell_size * x[i]) - 1;
	}
	return ret;
}

// Determines permutation table for point array.
void
NeighborhoodSearch::z_sort()
{
	for (PointSet& d : m_point_sets)
	{
		d.m_sort_table.resize(d.n_points());
		std::iota(d.m_sort_table.begin(), d.m_sort_table.end(), 0);

		std::sort(d.m_sort_table.begin(), d.m_sort_table.end(),
			[&](unsigned int a, unsigned int b)
		{
			return z_value(cell_index(d.point(a))) < z_value(cell_index(d.point(b)));
		});
	}
	m_initialized = false;
}


// Build hash table and entry array from scratch.
void
NeighborhoodSearch::init()
{
	m_entries.clear();
	m_map.clear();

	// Determine existing entries.
	std::vector<HashKey> temp_keys;
	for (unsigned int j = 0; j < m_point_sets.size(); ++j)
	{
		PointSet& d = m_point_sets[j];
		d.m_locks.resize(d.n_points());
		for (unsigned int i = 0; i < d.n_points(); i++)
		{
			HashKey const& key = cell_index(d.point(i));
			d.m_keys[i] = d.m_old_keys[i] = key;
			auto it = m_map.find(key);
			if (it == m_map.end())
			{
				m_entries.push_back({{ j, i }});
				temp_keys.push_back(key);
				m_map[key] = static_cast<unsigned int>(m_entries.size() - 1);
			}
			else
			{
				m_entries[it->second].add({j, i});
			}
		}
	}

	m_map.clear();
	for (unsigned int i = 0; i < m_entries.size(); ++i)
	{
		m_map.emplace(temp_keys[i], i);
	}

	m_initialized = true;
}


void
NeighborhoodSearch::increase_point_set_size(unsigned int index, Real const* x, std::size_t size)
{
	PointSet& point_set = m_point_sets[index];
	std::size_t old_size = point_set.n_points();

	// Insert new entries.
	for (unsigned int i = static_cast<unsigned int>(old_size); i < point_set.n_points(); i++)
	{
		HashKey key = cell_index(point_set.point(i));
		point_set.m_keys[i] = point_set.m_old_keys[i] = key;
		auto it = m_map.find(key);
		if (it == m_map.end())
		{
			m_entries.push_back({{ index, i }});
			m_map[key] = static_cast<unsigned int>(m_entries.size() - 1);
		}
		else 
		{
			m_entries[it->second].add({ index, i });
		}
	}

	point_set.resize(x, size);

}

void
NeighborhoodSearch::find_neighbors()
{
	if (!m_initialized)
	{
		init();
		m_initialized = true;
	}

	// Precompute cell indices.
	for (PointSet& d : m_point_sets)
	{
		if (!d.is_dynamic()) continue;
		d.m_keys.swap(d.m_old_keys);
		for (unsigned int i = 0; i < d.n_points(); ++i)
			d.m_keys[i] = cell_index(d.point(i));
	}

	std::vector<unsigned int> to_delete;
	if (m_erase_empty_cells)
	{
		to_delete.reserve(m_entries.size());
	}
	update_hash_table(to_delete);
	if (m_erase_empty_cells)
	{
		erase_empty_entries(to_delete);
	}
	query();
}

void
NeighborhoodSearch::erase_empty_entries(std::vector<unsigned int> const& to_delete)
{
	if (to_delete.empty())
		return;

	// Indicated empty cells.
	m_entries.erase(std::remove_if(m_entries.begin(), m_entries.end(), [](HashEntry const& entry)
	{
		return entry.indices.empty();
	}), m_entries.end());

	{
		auto it = m_map.begin();
		while (it != m_map.end())
		{
			auto& kvp = *it;

			if (kvp.second <= to_delete.front() && kvp.second >= to_delete.back() &&
				std::binary_search(to_delete.rbegin(), to_delete.rend(), kvp.second))
			{
				it = m_map.erase(it);
			}
			else
			{
				++it;
			}
		}
	}

	std::vector<std::pair<HashKey const, unsigned int>*> kvps(m_map.size());
	std::transform(m_map.begin(), m_map.end(), kvps.begin(),
	[](std::pair<HashKey const, unsigned int>& kvp)
	{
		return &kvp;
	});

	// Perform neighborhood search.
#ifdef _MSC_VER
	concurrency::parallel_for_each
#else
	__gnu_parallel::for_each
#endif
	(kvps.begin(), kvps.end(), [&](std::pair<HashKey const, unsigned int>* kvp_)
	{
		auto& kvp = *kvp_;

		for (unsigned int i = 0; i < to_delete.size(); ++i)
		{
			if (kvp.second >= to_delete[i])
			{
				kvp.second -= static_cast<unsigned int>(to_delete.size() - i);
				break;
			}
		}
	});
}

void
NeighborhoodSearch::update_hash_table(std::vector<unsigned int>& to_delete)
{
	// Indicate points changing inheriting cell.
	for (unsigned int j = 0; j < m_point_sets.size(); ++j)
	{
		PointSet& d = m_point_sets[j];
		for (unsigned int i = 0; i < d.n_points(); ++i)
		{
			if (d.m_keys[i] == d.m_old_keys[i]) continue;

			HashKey const& key = d.m_keys[i];
			auto it = m_map.find(key);
			if (it == m_map.end())
			{
				m_entries.push_back({{j, i}});
				m_map.insert({ key, static_cast<unsigned int>(m_entries.size() - 1) });
			}
			else
			{
				HashEntry& entry = m_entries[it->second];
				entry.add({j, i});
			}

			unsigned int entry_index = m_map[d.m_old_keys[i]];
			m_entries[entry_index].erase({j, i});
			if (m_erase_empty_cells)
			{
				if (m_entries[entry_index].n_indices() == 0)
				{
					to_delete.push_back(entry_index);
				}
			}
		}
	}

	to_delete.erase(std::remove_if(to_delete.begin(), to_delete.end(),
	[&](unsigned int index)
	{
		return m_entries[index].n_indices() != 0;
	}), to_delete.end());
	std::sort(to_delete.begin(), to_delete.end(), std::greater<unsigned int>());
}

void
NeighborhoodSearch::query()
{
	for (PointSet& d : m_point_sets)
	{
		if (d.is_neighborsearch_enabled())
		{
			for (auto& n : d.m_neighbors)
			{
				n.clear();
			}
		}
	}

	std::vector<std::pair<HashKey const, unsigned int> const*> kvps(m_map.size());
	std::transform(m_map.begin(), m_map.end(), kvps.begin(),
	[](std::pair<HashKey const, unsigned int> const& kvp)
	{
		return &kvp;
	});

	// Perform neighborhood search.
#ifdef _MSC_VER
	concurrency::parallel_for_each
#else
	__gnu_parallel::for_each
#endif
	(kvps.begin(), kvps.end(), [&](std::pair<HashKey const, unsigned int> const* kvp_)
	{
		auto const& kvp = *kvp_;
		HashEntry const& entry = m_entries[kvp.second];
		HashKey const& key = kvp.first;

		for (unsigned int a = 0; a < entry.n_indices(); ++a)
		{
			PointID const& va = entry.indices[a];
			PointSet& da = m_point_sets[va.point_set_id];
			for (unsigned int b = a + 1; b < entry.n_indices(); ++b)
			{
				PointID const& vb = entry.indices[b];
				PointSet& db = m_point_sets[vb.point_set_id];

				if (!da.is_neighborsearch_enabled() && 
					!db.is_neighborsearch_enabled())
				{
					continue;
				}

				Real const* xa = da.point(va.point_id);
				Real const* xb = db.point(vb.point_id);
				Real tmp = xa[0] - xb[0];
				Real l2 = tmp * tmp;
				tmp = xa[1] - xb[1];
				l2 += tmp * tmp;
				tmp = xa[2] - xb[2];
				l2 += tmp * tmp;

				if (l2 < m_r2)
				{
					if (da.is_neighborsearch_enabled())
					{
						da.m_neighbors[va.point_id].push_back(vb);
					}
					if (db.is_neighborsearch_enabled())
					{
						db.m_neighbors[vb.point_id].push_back(va);
					}
				}
			}
		}
	});


	std::vector<std::array<bool, 27>> visited(m_entries.size(), {false});
	std::vector<Spinlock> entry_locks(m_entries.size());

#ifdef _MSC_VER
	concurrency::parallel_for_each
#else
	__gnu_parallel::for_each
#endif
	(kvps.begin(), kvps.end(), [&](std::pair<HashKey const, unsigned int> const* kvp_)
	{
		auto const& kvp = *kvp_;
		HashEntry const& entry = m_entries[kvp.second];
		HashKey const& key = kvp.first;

		for (int dj = -1; dj <= 1; dj++)
		for (int dk = -1; dk <= 1; dk++)
		for (int dl = -1; dl <= 1; dl++)
		{
			int l_ind = 9 * (dj + 1) + 3 * (dk + 1) + (dl + 1);
			if (l_ind == 13)
			{
				continue;
			}
			entry_locks[kvp.second].lock();
			if (visited[kvp.second][l_ind])
			{
				entry_locks[kvp.second].unlock();
				continue;
			}
			entry_locks[kvp.second].unlock();



			auto it = m_map.find({ key.k[0] + dj, key.k[1] + dk, key.k[2] + dl });
			if (it == m_map.end())
				continue;

			std::array<unsigned int, 2> entry_ids{{kvp.second, it->second}};
			if (entry_ids[0] > entry_ids[1])
				std::swap(entry_ids[0], entry_ids[1]);
			entry_locks[entry_ids[0]].lock();
			entry_locks[entry_ids[1]].lock();

			if (visited[kvp.second][l_ind])
			{
				entry_locks[entry_ids[1]].unlock();
				entry_locks[entry_ids[0]].unlock();
				continue;
			}

			visited[kvp.second][l_ind] = true;
			visited[it->second][26 - l_ind] = true;

			entry_locks[entry_ids[1]].unlock();
			entry_locks[entry_ids[0]].unlock();

			for (unsigned int i = 0; i < entry.n_indices(); ++i)
			{
				PointID const& va = entry.indices[i];
				HashEntry const& entry_ = m_entries[it->second];
				unsigned int n_ind = entry_.n_indices();
				for (unsigned int j = 0; j < n_ind; ++j)
				{
					PointID const& vb = entry_.indices[j];
					PointSet& db = m_point_sets[vb.point_set_id];

					PointSet& da = m_point_sets[va.point_set_id];
					if (!da.is_neighborsearch_enabled() && 
						!db.is_neighborsearch_enabled())
					{
						continue;
					}

					Real const* xa = da.point(va.point_id);
					Real const* xb = db.point(vb.point_id);
					Real tmp = xa[0] - xb[0];
					Real l2 = tmp * tmp;
					tmp = xa[1] - xb[1];
					l2 += tmp * tmp;
					tmp = xa[2] - xb[2];
					l2 += tmp * tmp;
					if (l2 < m_r2)
					{
						if (da.is_neighborsearch_enabled())
						{
							da.m_locks[va.point_id].lock();
							da.m_neighbors[va.point_id].push_back(vb);
							da.m_locks[va.point_id].unlock();
						}
						if (db.is_neighborsearch_enabled())
						{
							db.m_locks[vb.point_id].lock();
							db.m_neighbors[vb.point_id].push_back(va);
							db.m_locks[vb.point_id].unlock();
						}
					}
				}
			}
		
		}
	});
}

}

