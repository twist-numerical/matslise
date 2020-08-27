#ifndef MATSLISE_MATSLISE3D_H
#define MATSLISE_MATSLISE3D_H

#include "matslise.h"
#include "util/rectangle.h"

namespace matslise {
    template<typename Scalar>
    class AbstractMatslise3D {
    public:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> ArrayXXs;

        const std::function<Scalar(Scalar, Scalar, Scalar)> potential;
        const Rectangle<Scalar, 3> domain;

        AbstractMatslise3D(
                const std::function<Scalar(Scalar, Scalar, Scalar)> &potential, const Rectangle<Scalar, 3> &domain)
                : potential(potential), domain(domain) {
        }

        virtual std::pair<Scalar, Eigen::Index> eigenvalue(const Scalar &Eguess) const = 0;

        virtual Scalar eigenvalueError(const Scalar &E) const = 0;

        virtual std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvalues(const Scalar &Emin, const Scalar &Emax) const = 0;

        virtual std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvaluesByIndex(Eigen::Index Imin, Eigen::Index Imax) const = 0;

        virtual Eigen::Index estimateIndex(const Scalar &E) const = 0;

        virtual ~AbstractMatslise3D() = default;
    };

    template<typename>
    class Matslise3D;

    template<typename _Scalar>
    class Matslise3DSector : public MatsliseNDSector<_Scalar> {
    public:
        using typename MatsliseNDSector<_Scalar>::Scalar;
        static const bool expensive = true;
        typedef typename AbstractMatslise2D<Scalar>::VectorXs VectorXs;
        typedef typename AbstractMatslise2D<Scalar>::MatrixXs MatrixXs;
        typedef typename AbstractMatslise2D<Scalar>::ArrayXs ArrayXs;
        typedef typename AbstractMatslise2D<Scalar>::ArrayXXs ArrayXXs;

        const Matslise3D<Scalar> *matslise3d;
        std::shared_ptr<matslise::AbstractMatslise2D<Scalar>> matslise2d;
        using MatsliseNDSector<Scalar>::matscs;
        using MatsliseNDSector<Scalar>::min;
        using MatsliseNDSector<Scalar>::max;
        Scalar zbar;
        ArrayXXs vbar;

        std::vector<Scalar> eigenvalues;
        std::vector<Eigenfunction2D < Scalar>> eigenfunctions;
        std::vector<ArrayXXs> eigenfunctions_grid;
        Direction direction = none;

        Matslise3DSector(const Matslise3D<Scalar> *matslise3d) :MatsliseNDSector<Scalar>(), matslise3d(matslise3d) {}

        Matslise3DSector(const Matslise3D<Scalar> *matslise3d, const Scalar &min, const Scalar &max, Direction);

        void setDirection(Direction);

        bool contains(const Scalar &point) const {
            return point <= max && point >= min;
        }

        std::vector<ArrayXXs> basis(const ArrayXs &x, const ArrayXs &y) const;

        static bool compare(const Matslise3DSector<Scalar> &a, const Matslise3DSector<Scalar> &b) {
            return a.matslise2d->estimatePotentialMinimum() < b.matslise2d->estimatePotentialMinimum();
        }

        Matslise3DSector<_Scalar> *
        refine(const Matslise3D<Scalar> *problem, const Scalar &min, const Scalar &max, Direction) const;
    };

    template<typename _Scalar=double>
    class Matslise3D : public AbstractMatslise3D<_Scalar>, private MatsliseND<_Scalar, Matslise3DSector<_Scalar>> {
    public:
        typedef _Scalar Scalar;
        using typename AbstractMatslise3D<Scalar>::VectorXs;
        using typename AbstractMatslise3D<Scalar>::MatrixXs;
        using typename AbstractMatslise3D<Scalar>::ArrayXs;
        using typename AbstractMatslise3D<Scalar>::ArrayXXs;
        static constexpr int order = matslise::Matscs<Scalar>::order;
        static constexpr bool refineSectors = true;

        struct Config {
            Scalar tolerance = 1e-6;
            bool xSymmetric = false;
            bool ySymmetric = false;
            std::optional<SectorBuilder < Matslise < Scalar>, Scalar>> xSectorBuilder;
            std::optional<SectorBuilder < Matslise2D < Scalar>, Scalar>> ySectorBuilder;
            std::optional<SectorBuilder < Matslise3D<Scalar>, Scalar>> zSectorBuilder;
            Eigen::Index zStepsPerSector = 5;
            Eigen::Index yStepsPerSector = 3;
            Eigen::Index xyBasisSize = 12;
            Eigen::Index xBasisSize = 12;
        };

        const Config config;
        ArrayXs grid_x;
        ArrayXs grid_y;


        using Sector = Matslise3DSector<Scalar>;
        using MatsliseND<Scalar, Matslise3DSector<Scalar>>::M;
        using MatsliseND<Scalar, Matslise3DSector<Scalar>>::sectors;
        using MatsliseND<Scalar, Matslise3DSector<Scalar>>::matchIndex;
        using MatsliseND<Scalar, Matslise3DSector<Scalar>>::dirichletBoundary;
        using AbstractMatslise3D<Scalar>::domain;
        using AbstractMatslise3D<Scalar>::potential;
    public:
        Matslise3D(const std::function<Scalar(Scalar, Scalar, Scalar)> &potential,
                   const matslise::Rectangle<Scalar, 3> &domain, const Config &config = Config());

        virtual ~Matslise3D();

        using MatsliseND<Scalar, Sector>::propagate;
        using MatsliseND<Scalar, Sector>::matchingErrorMatrix;
        using MatsliseND<Scalar, Sector>::matchingErrors;
        using MatsliseND<Scalar, Sector>::matchingError;
        using MatsliseND<Scalar, Sector>::eigenvalue;
        using MatsliseND<Scalar, Sector>::eigenvalueError;
        using MatsliseND<Scalar, Sector>::eigenvalues;
        using MatsliseND<Scalar, Sector>::eigenvaluesByIndex;
        using MatsliseND<Scalar, Sector>::estimateIndex;

    public: // overrides
        std::pair<Scalar, Eigen::Index> eigenvalue(const Scalar &Eguess) const override {
            return MatsliseND<Scalar, Sector>::eigenvalue(Eguess);
        }

        Scalar eigenvalueError(const Scalar &E) const override {
            return MatsliseND<Scalar, Sector>::eigenvalueError(E);
        }

        std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvalues(const Scalar &Emin, const Scalar &Emax) const override {
            return MatsliseND<Scalar, Sector>::eigenvalues(Emin, Emax);
        }

        std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvaluesByIndex(Eigen::Index Imin, Eigen::Index Imax) const override {
            return MatsliseND<Scalar, Sector>::eigenvaluesByIndex(Imin, Imax);
        }

        Eigen::Index estimateIndex(const Scalar &E) const override {
            return MatsliseND<Scalar, Sector>::estimateIndex(E);
        }
    };
}

#endif //MATSLISE_MATSLISE3D_H
