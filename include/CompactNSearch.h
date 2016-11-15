
#ifndef COMPACT_N_SEARCH__HPP
#define COMPACT_N_SEARCH__HPP

#include "Config.h"
#include "DataStructures.h"
#include "PointSet.h"

#include <unordered_map>

namespace CompactNSearch
{

/**
* @class NeighborhoodSearch
* Stores point data multiple set of points in which neighborhood information for a fixed
* radius r should be generated.
*/
class NeighborhoodSearch
{

public:

	/**
	* Constructor.
	* Creates a new instance of the neighborhood search class.
	* @param r Search radius. If two points are closer to each other than a distance r they are considered neighbors.
	* @param erase_empty_cells If true. Empty cells in spatial hashing grid are erased if the points move.
	*/
	NeighborhoodSearch(Real r, bool erase_empty_cells = false);

	/**
	* Destructor.
	*/
	virtual ~NeighborhoodSearch() = default;

	/**
	* Get method to access a point set.
	* @param i Index of the point set to retrieve.
	*/
	PointSet const& point_set(unsigned int i) const { return m_point_sets[i];     }

	/**
	* Get method to access a point set.
	* @param i Index of the point set to retrieve.
	*/
	PointSet      & point_set(unsigned int i)       { return m_point_sets[i];     }


	/**
	* Returns the number of point sets contained in the search.
	*/
	std::size_t  n_point_sets()               const { return m_point_sets.size(); }

	/**
	* Get method to access the list of point sets.
	*/
	std::vector<PointSet> const& point_sets() const { return m_point_sets;        }

	/**
	* Get method to access the list of point sets.
	*/
	std::vector<PointSet>      & point_sets()       { return m_point_sets;        }

	/**
	* Increases the size of a point set under the assumption that the existing points remain at
	* the same position.
	* @param i Index of point set that will be resized.
	* @param x Pointer to the point position data. Must point to continguous data of 3 * n
	* real values.
	* @param n Number of points.
	*/
	void increase_point_set_size(unsigned int i, Real const* x, std::size_t n);

	/**
	* Creates and adds a new set of points.
	* @param x Pointer to the point position data. Must point to continguous data of 3 * n
	* real values.
	* @param n Number of points.
	* @param is_dynamic Specifies wether the point positions will change for future queries.
	* @param search_neighbors If true, no neighbors of this point set will be determined during the
	* actual search. However, other point set will still determine neighboring points belonging to
	* this point set.
	* @returns Returns unique identifier in form of an index assigned to the newly created point
	* set.
	*/
	unsigned int add_point_set(Real const* x, std::size_t n, bool is_dynamic = true,
		bool search_neighbors = true) 
	{ 
		m_point_sets.push_back({x, n, is_dynamic, search_neighbors});
		return static_cast<unsigned int>(m_point_sets.size() - 1);
	}

	/**
	* Performs the actual query. This method will assign a list of neighboring points to each point
	* every added point set.
	*/
	void find_neighbors();

	/*
	* Generates a sort table according to a space-filling Z curve. Any array-based per point
	* information can then be reordered using the function sort_field of the PointSet class.
	* Please note that the position data will not be modified by this class, such that the user has
	* to invoke the sort_field function on the position array. Moreover, be aware the the grid has
	* be reinitialized after each sort. Therefore, the points should not be reordered too
	* frequently.
	*/
	void z_sort();

	/*
	* @returns Returns the radius in which point neighbors are searched.
	*/
	double radius() const { return std::sqrt(m_r2); }

	/**
	* Sets the radius in which point point neighbors are searched.
	* @param r Search radius.
	*/
	void set_radius(double r) 
	{ 
		m_r2 = r * r; 
		m_inv_cell_size = 1.0 / r;
		m_initialized = false;
	}

private:

	void init();
	void update_hash_table(std::vector<unsigned int>& to_delete);
	void erase_empty_entries(std::vector<unsigned int> const& to_delete);
	void query();

	HashKey cell_index(Real const* x) const;

private:


	std::vector<PointSet> m_point_sets;

	Real m_inv_cell_size;
	Real m_r2;
	std::unordered_map<HashKey, unsigned int, SpatialHasher> m_map;
	std::vector<HashEntry> m_entries;

	bool m_erase_empty_cells;
	bool m_initialized;
};

}

#endif //COMPACT_N_SEARCH__HPP
