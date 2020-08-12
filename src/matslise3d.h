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
        const Rectangle<3, Scalar> domain;

        AbstractMatslise3D(
                const std::function<Scalar(Scalar, Scalar, Scalar)> &potential, const Rectangle<3, Scalar> &domain)
                : potential(potential), domain(domain) {
        }

        virtual Scalar eigenvalue(const Scalar &Eguess) const = 0;

        virtual Scalar eigenvalueError(const Scalar &E) const = 0;

        virtual ~AbstractMatslise3D() = default;
    };

    template<typename _Scalar=double>
    class Matslise3D : public AbstractMatslise3D<_Scalar> {
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
            Eigen::Index xyBasisSize = 12;
            Eigen::Index xBasisSize = 12;
        };

        const Config config;
        ArrayXs grid_x;
        ArrayXs grid_y;

        class Sector;

        using AbstractMatslise3D<Scalar>::domain;
        using AbstractMatslise3D<Scalar>::potential;
        std::vector<MatrixXs> M;
        int matchIndex;

        std::vector<typename Matslise3D<Scalar>::Sector *> sectors;
    public:
        Matslise3D(const std::function<Scalar(Scalar, Scalar, Scalar)> &potential,
                   const matslise::Rectangle<3, Scalar> &domain, const Config &config = Config());

        virtual ~Matslise3D();

        Y <Scalar, Eigen::Dynamic>
        propagate(const Scalar &E, const Y <Scalar, Eigen::Dynamic> &y0, const Scalar &a, const Scalar &b,
                  bool use_h = true) const;

    public:
        std::pair<MatrixXs, MatrixXs> matchingErrorMatrix(const Scalar &E, bool use_h = true) const;

        std::vector<std::pair<Scalar, Scalar>> matchingErrors(const Scalar &E, bool use_h = true) const;

        std::pair<Scalar, Scalar> matchingError(const Scalar &E, bool use_h = true) const;

        Scalar eigenvalue(const Scalar &Eguess, bool use_h) const;

        Scalar eigenvalue(const Scalar &Eguess) const override {
            return eigenvalue(Eguess, true);
        }

        Scalar eigenvalueError(const Scalar &E) const override;

    protected:
        MatrixXs conditionY(Y <Scalar, Eigen::Dynamic> &y) const;

    public:
        class Sector {
        public:
            static const bool expensive = true;

            const Matslise3D<Scalar> *matslise3d;
            std::shared_ptr<matslise::AbstractMatslise2D<Scalar>> matslise2d;
            typename matslise::Matscs<Scalar>::Sector *matscs;
            Scalar min, max;
            Scalar zbar;
            ArrayXXs vbar;

            std::vector<Scalar> eigenvalues;
            std::vector<Eigenfunction2D < Scalar>> eigenfunctions;
            std::vector<ArrayXXs> eigenfunctions_grid;
            Direction direction = none;

            Sector(const Matslise3D<Scalar> *matslise3d) : matslise3d(matslise3d) {}

            Sector(const Matslise3D<Scalar> *matslise3d, const Scalar &min, const Scalar &max, Direction);

            virtual ~Sector();

            void setDirection(Direction);

            template<int r>
            Y <Scalar, Eigen::Dynamic, r> propagate(
                    const Scalar &E, const Y <Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b,
                    bool use_h = true) const;

            bool contains(const Scalar &point) const {
                return point <= max && point >= min;
            }

            std::vector<ArrayXXs> basis(const ArrayXs &x, const ArrayXs &y) const;

            Scalar error() const;

            static bool compare(const Sector &a, const Sector &b) {
                return a.matslise2d->estimatePotentialMinimum() < b.matslise2d->estimatePotentialMinimum();
            }

            Sector *
            refine(const Matslise3D<Scalar> *problem, const Scalar &min, const Scalar &max, Direction) const;
        };
    };
}

#endif //MATSLISE_MATSLISE3D_H
