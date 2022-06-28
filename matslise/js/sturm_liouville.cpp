#include "module.h"
#include "../liouville.h"

using namespace matslise;
using namespace emscripten;
using namespace std;
using namespace Eigen;

void bind_sturm_lioville() {

    class_<SturmLiouville<double>>("SturmLiouville")
            .constructor(optional_override(
                    [](val p, val q, val r, double min, double max, double tolerance) -> SturmLiouville<double> * {
                        return new SturmLiouville<double>(
                                val2function<double(double)>(p), val2function<double(double)>(q),
                                val2function<double(double)>(r), {min, max}, tolerance);
                    }))
            .function("eigenvaluesByIndex", optional_override(
                    [](const SturmLiouville<double> &m, int Imin, int Imax, const Vector2d &left,
                       const Vector2d &right) -> val {
                        auto es = m.eigenvaluesByIndex(Imin, Imax, Y<>(left, {0, 0}), Y<>(right, {0, 0}));
                        return toValArray(es.begin(), es.end(), [](auto t) { return Eigenvalue(t.first, t.second); });
                    }))
            .function("eigenvalueError", optional_override(
                    [](const SturmLiouville<double> &m, double E, const Vector2d &left, const Vector2d &right,
                       int index = -1) -> double {
                        return m.eigenvalueError(E, Y<>(left, {0, 0}), Y<>(right, {0, 0}), index);
                    }))
            .function("eigenpairsByIndex", optional_override(
                    [](const SturmLiouville<double> &m, int Imin, int Imax, const Vector2d &left,
                       const Vector2d &right) -> val {
                        auto es = m.eigenpairsByIndex(Imin, Imax, Y<>(left, {0, 0}), Y<>(right, {0, 0}));
                        return toValArray(es.begin(), es.end(), [](auto &t) {
                            return Eigenpair(get<0>(t), get<1>(t), wrapEigenfunction(std::move(get<2>(t))));
                        });
                    }));

}
