#ifndef __my_rtree_h__
#define __my_rtree_h__

#include "PhysiCell_cell.h"
// Bibliotecas externas para RTree (Boost.Geometry) (NUEVO)
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

using namespace PhysiCell;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Define 3D point type (any differences if using 2D?)
typedef bg::model::point<double, 3, bg::cs::cartesian> PointTy;

// indexable type for R-tree based on Cell position
//typedef std::pair<PointTy, Cell*> ValueTy;
typedef Cell* ValueTy;

// Define the bounding box type for the R-tree
typedef bg::model::box<PointTy> BoxTy;

// Extrae el punto de ValueTy
struct Indexable {
	/*
	using result_type = const PointTy&;
	result_type operator()(const ValueTy& v) const {
		return v.first;
	}
	*/
	const PointTy &operator()(const Cell* agent) const {
			return PointTy(agent->position[0], agent->position[1], agent->position[2]);
	}
};

// Compara si dos celdas (por puntero) son iguales
struct EqualTo {
	bool operator()(const ValueTy& lhs, const ValueTy& rhs) const {
		//return lhs.second == rhs.second;
		return lhs == rhs; // Compara punteros de celdas
	}
};

# endif // __my_rtree_h__