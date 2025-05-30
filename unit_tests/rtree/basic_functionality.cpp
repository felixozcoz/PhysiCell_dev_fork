
/*
Compilar con:


*/


#include "../../core/PhysiCell.h"
#include "../../core/PhysiCell_cell_container.h" // o tu propia clase RTree
#include <cassert>
#include <iostream>

using namespace PhysiCell;

void init_prng_for_tests()
{
    int n_threads = omp_get_max_threads();
    PhysiCell::physicell_random_seeds.resize(n_threads);
    for (int i = 0; i < n_threads; i++)
    {
        PhysiCell::physicell_random_seeds[i] = std::time(nullptr) + i;
    }
}

void test_rtree_basic_operations() {
    std::cout << "[Test] R-tree basic operations...\n";

    // 1. Crear contenedor
    Cell_RTree_Container container;

    // 2. Crear células ficticias
    Cell* cell1 = new Cell();
    cell1->position = {1.0, 1.0, 1.0};

    Cell* cell2 = new Cell();
    cell2->position = {2.0, 2.0, 2.0};

    Cell* cell3 = new Cell();
    cell3->position = {100.0, 200.0, 300.0};

    // 3. Registrar células
    container.register_agent(cell1);
    container.register_agent(cell2);
    container.register_agent(cell3);

    // 4. Consultar punto cercano a cell1
    PointTy query_point(1.0, 1.0, 2.0);

    std::vector<ValueTy> result;
    container.rtree.query(bgi::nearest(query_point, 1), std::back_inserter(result));

    assert(!result.empty());
    assert(result[0] == cell1);
    std::cout << " - Nearest neighbor query PASSED.\n";

    // 5. Consultar bounding box con centro entre cell1 y cell2
    result.clear();
    PointTy min_point(0.0, 0.0, 0.0);
    PointTy max_point(3.0, 3.0, 3.0);

    container.rtree.query(bgi::intersects(BoxTy(min_point, max_point)), std::back_inserter(result));
    assert(result.size() == 2);
    std::cout << "results.size() = " << result.size() << "\n";
    std::cout << " - Bounding box query PASSED.\n";

    // 6. Eliminar cell2
    container.rtree.remove(cell2); //remove_agent(cell2);

    // 6bis. Consultar si cell2 aún está en el contenedor (no debería aparecer)
    result.clear();
    PointTy query_point_cell2(2.0, 2.0, 2.0); // Punto donde estaba cell2

    container.rtree.query(bgi::nearest(query_point_cell2, 1), std::back_inserter(result));

    // Buscar explícitamente cell2 en los resultados
    assert(result[0] != cell2 && "Cell2 should not be found after removal");



    // 7. Consultar de nuevo cerca de cell2
    bg::assign_values(query_point, 15.0, 25.0, 35.0);
    result.clear();
    container.rtree.query(bgi::nearest(query_point, 1), std::back_inserter(result));

    assert(!result.empty());
    assert(result[0] != cell2);
    std::cout << " - Removal PASSED.\n";

    // 7. Liberar memoria
    delete cell1;
    delete cell2;
    delete cell3;

    std::cout << "[Test] All R-tree basic operations PASSED!\n";
}

int main() {
    init_prng_for_tests();
    test_rtree_basic_operations();
    return 0;
}
