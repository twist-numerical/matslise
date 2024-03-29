#ifndef SCHRODINGER_MATSLISE_H
#define SCHRODINGER_MATSLISE_H

#include <cmath>
#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <functional>
#include <memory>
#include "util/eigen.h"
#include "util/y.h"
#include "util/sectorbuilder.h"
#include "./formula_constants.h"

namespace matslise {
    enum Direction {
        none, forward, backward
    };

    template<typename Scalar>
    struct Domain {
        Scalar min;
        Scalar max;

        bool contains(const Scalar &x) const {
            return min <= x && x <= max;
        }
    };

    template<typename Scalar>
    class AbstractMatslise {
    public:
        static void checkSymmetry(const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right) {
            if (left != right)
                throw std::runtime_error("Halfrange::checkSymmetry(), left and right sides have to be identical");
        }

    public:
        const std::function<Scalar(Scalar)> potential;
        const Domain<Scalar> domain;

        AbstractMatslise(const std::function<Scalar(Scalar)> &potential, const Domain<Scalar> &domain)
                : potential(potential), domain(domain) {
        }

        class Eigenfunction {
        public:
            virtual Eigen::Array<Scalar, 2, 1> operator()(const Scalar &x) const = 0;

            virtual Eigen::Array<Scalar, Eigen::Dynamic, 2>
            operator()(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) const = 0;

            virtual ~Eigenfunction() = default;
        };

    public:

        virtual Scalar estimatePotentialMinimum() const = 0;

        virtual std::vector<std::pair<int, Scalar>>
        eigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left,
                    const matslise::Y<Scalar> &right) const {
            checkSymmetry(left, right);
            return eigenvalues(Emin, Emax, right);
        }


        virtual std::vector<std::pair<int, Scalar>>
        eigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left) const {
            return eigenvalues(Emin, Emax, left, left);
        }

        virtual std::vector<std::pair<int, Scalar>>
        eigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left,
                           const matslise::Y<Scalar> &right) const {
            checkSymmetry(left, right);
            return eigenvaluesByIndex(Imin, Imax, right);
        }

        virtual std::vector<std::pair<int, Scalar>>
        eigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &side) const {
            return eigenvaluesByIndex(Imin, Imax, side, side);
        };

        virtual std::vector<std::tuple<int, Scalar, std::unique_ptr<Eigenfunction>>>
        eigenpairsByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left,
                          const matslise::Y<Scalar> &right) const;

        virtual std::vector<std::tuple<int, Scalar, std::unique_ptr<Eigenfunction>>>
        eigenpairsByIndex(int Imin, int Imax, const matslise::Y<Scalar> &side) const {
            return eigenpairsByIndex(Imin, Imax, side, side);
        };

        virtual Scalar
        eigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left,
                        const matslise::Y<Scalar> &right, int index = -1) const {
            checkSymmetry(left, right);
            return eigenvalueError(E, right, index);
        }

        virtual Scalar eigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left, int index = -1) const {
            return eigenvalueError(E, left, left, index);
        };

        virtual std::unique_ptr<Eigenfunction> eigenfunction(
                const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                int index = -1) const {
            checkSymmetry(left, right);
            return eigenfunction(E, right, index);
        };

        virtual std::unique_ptr<Eigenfunction> eigenfunction(
                const Scalar &E, const matslise::Y<Scalar> &left, int index = -1) const {
            return eigenfunction(E, left, left, index);
        };

        virtual ~AbstractMatslise() = default;
    };

    template<typename Scalar_=double>
    class Matslise : public AbstractMatslise<Scalar_> {
    public:
        typedef Scalar_ Scalar;
        static constexpr int order = MATSLISE_HMAX_delta - 1;
        static constexpr bool refineSectors = false;

        class Sector;

        using AbstractMatslise<Scalar>::domain;
        using AbstractMatslise<Scalar>::potential;
        typedef typename AbstractMatslise<Scalar>::Eigenfunction Eigenfunction;

        int sectorCount;
        int matchIndex;
        std::vector<std::unique_ptr<Matslise::Sector>> sectors;

        Scalar tolerance;
    public:
        Matslise(std::function<Scalar(const Scalar &)> V, const Scalar &xmin, const Scalar &xmax,
                 const Scalar &tolerance = 1e-8) : Matslise(V, Domain<Scalar>{xmin, xmax}, tolerance) {}

        Matslise(std::function<Scalar(const Scalar &)> V, const Domain<Scalar> &domain,
                 const Scalar &tolerance = 1e-8)
                : Matslise(V, domain, tolerance, sector_builder::AutomaticSectorBuilder<Matslise<Scalar>>(tolerance)) {}

        Matslise(std::function<Scalar(const Scalar &)> V, const Scalar &xmin, const Scalar &xmax,
                 const Scalar &tolerance, const sector_builder::SectorBuilder<Matslise<Scalar>> &sectorBuilder)
                : Matslise(V, Domain<Scalar>{xmin, xmax}, tolerance, sectorBuilder) {}

        Matslise(std::function<Scalar(const Scalar &)> V, const Domain<Scalar> &domain,
                 const Scalar &tolerance, const sector_builder::SectorBuilder<Matslise<Scalar>> &sectorBuilder)
                : AbstractMatslise<Scalar>(V, domain),
                  sectorCount(0), matchIndex(0), sectors(), tolerance(tolerance) {
            auto buildSectors = sectorBuilder(this, domain.min, domain.max);
            sectors = std::move(buildSectors.sectors);
            matchIndex = buildSectors.matchIndex;
            sectorCount = sectors.size();
            Scalar &v0Match = sectors[matchIndex]->vs[0];
            for (auto &s: sectors)
                s->v0Match = v0Match;
        }

        template<int cols = 1>
        std::pair<matslise::Y<Scalar, 1, cols>,
                typename std::conditional<cols == 1, Scalar, Eigen::Array<Scalar, cols, 1>>::type>
        propagate(const Scalar &E, const matslise::Y<Scalar, 1, cols> &y, const Scalar &a, const Scalar &b,
                  bool use_h = true) const;

        std::tuple<Scalar, Scalar, Scalar>
        matchingError(const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                      bool use_h = true) const;

    public: // Override
        Scalar estimatePotentialMinimum() const override;

        using AbstractMatslise<Scalar>::eigenvalues;

        std::vector<std::pair<int, Scalar>>
        eigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left,
                    const matslise::Y<Scalar> &right) const override;

        using AbstractMatslise<Scalar>::eigenvaluesByIndex;

        std::vector<std::pair<int, Scalar>>
        eigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left,
                           const matslise::Y<Scalar> &right) const override;

        using AbstractMatslise<Scalar>::eigenvalueError;

        Scalar eigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left,
                               const matslise::Y<Scalar> &right, int index = -1) const override;

        using AbstractMatslise<Scalar>::eigenfunction;

        std::unique_ptr<Eigenfunction> eigenfunction(const Scalar &E, const matslise::Y<Scalar> &left,
                                                     const matslise::Y<Scalar> &right, int index = -1) const override;

        virtual ~Matslise() = default;

    public:
        class Sector {
        public:
            static const bool expensive = false;

            Eigen::Array<Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign>, MATSLISE_ETA_delta, MATSLISE_HMAX_delta, Eigen::DontAlign> t_coeff;
            Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign> t_coeff_h[MATSLISE_ETA_h];
            Scalar v0Match = 0;
            std::array<Scalar, MATSLISE_N> vs;
            Scalar min, max, h;
            Direction direction = none;

            Sector(const Matslise *problem, const Scalar &min, const Scalar &max, Direction direction);

            void setDirection(Direction);

            void calculateTCoeffs();

            bool contains(const Scalar &point) const {
                return point <= max && point >= min;
            }

            T<Scalar> calculateT(const Scalar &E, bool use_h = true) const;

            T<Scalar> calculateT(const Scalar &E, const Scalar &delta, bool use_h = true) const;

            template<bool withPrufer, int cols = 1>
            typename std::conditional<withPrufer,
                    std::pair<Y<Scalar, 1, cols>, Eigen::Array<Scalar, cols, 1>>,
                    Y<Scalar, 1, cols>
            >::type
            propagate(const Scalar &E, const Y<Scalar, 1, cols> &y0, const Scalar &a, const Scalar &b,
                      bool use_h = true) const;

            Scalar theta0(const Scalar &E, const Y<Scalar> &y0) const;

            template<int col>
            typename std::conditional<col == 1, Scalar, Eigen::Array<Scalar, col, 1>>::type
            prufer(const Scalar &E, const Scalar &delta, const Y<Scalar, 1, col> &y0,
                   const Y<Scalar, 1, col> &y1) const;

            Scalar error() const;

            ~Sector() = default;

            static bool compare(const Sector &a, const Sector &b) {
                return a.vs[0] < b.vs[0];
            }
        };
    };

    template<typename Scalar_ = double>
    class MatsliseHalf : public AbstractMatslise<Scalar_> {
        typedef Scalar_ Scalar;
        typedef typename AbstractMatslise<Scalar>::Eigenfunction Eigenfunction;

    public:
        std::unique_ptr<const Matslise<Scalar>> ms;
    public:
        MatsliseHalf(std::function<Scalar(Scalar)> V, const Scalar &xmax, const Scalar &tolerance)
                : MatsliseHalf(V, xmax, tolerance,
                               sector_builder::AutomaticSectorBuilder<Matslise<Scalar>>(tolerance)) {}

        MatsliseHalf(std::function<Scalar(Scalar)> V, const Scalar &xmax, const Scalar &tolerance,
                     const sector_builder::SectorBuilder<Matslise<Scalar>> &sectorBuilder);

        Scalar estimatePotentialMinimum() const override {
            return ms->estimatePotentialMinimum();
        }

        using AbstractMatslise<Scalar>::eigenvalues;

        std::vector<std::pair<int, Scalar>>
        eigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &side) const override;

        using AbstractMatslise<Scalar>::eigenvaluesByIndex;

        std::vector<std::pair<int, Scalar>>
        eigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &side) const override;

        using AbstractMatslise<Scalar>::eigenvalueError;

        Scalar eigenvalueError(const Scalar &E, const matslise::Y<Scalar> &side, int index = -1) const override;

        using AbstractMatslise<Scalar>::eigenfunction;

        std::unique_ptr<Eigenfunction>
        eigenfunction(const Scalar &E, const matslise::Y<Scalar> &side, int index = -1) const override;

        virtual ~MatsliseHalf() = default;
    };

    template<typename Scalar=double>
    class PeriodicMatslise {
    public:
        Matslise<Scalar> matslise;
        Eigen::Matrix<Scalar, 2, 2> K;

        using Eigenfunction = typename AbstractMatslise<Scalar>::Eigenfunction;

        PeriodicMatslise(std::function<Scalar(const Scalar &)> V, const Scalar &xmin, const Scalar &xmax,
                         const Scalar &tolerance = 1e-8,
                         const Eigen::Matrix<Scalar, 2, 2> &K = Eigen::Matrix<Scalar, 2, 2>::Identity())
                : matslise{V, xmin, xmax, tolerance}, K(K) {}

        PeriodicMatslise(std::function<Scalar(const Scalar &)> V, const Scalar &xmin, const Scalar &xmax,
                         const Scalar &tolerance, const sector_builder::SectorBuilder<Matslise<Scalar>> &sb,
                         const Eigen::Matrix<Scalar, 2, 2> &K = Eigen::Matrix<Scalar, 2, 2>::Identity())
                : matslise{V, xmin, xmax, tolerance, sb}, K(K) {}


        std::pair<matslise::Y<Scalar, 1, 2>, Eigen::Array<Scalar, 2, 1>>
        propagate(const Scalar &E, const matslise::Y<Scalar, 1, 2> &y0,
                  const Scalar &a, const Scalar &b, bool use_h = true) const;

        std::pair<matslise::Y<Scalar, 1, 2>, Eigen::Array<Scalar, 2, 1>>
        matchingError(const Scalar &E, bool use_h = true) const;

        std::vector<std::unique_ptr<Eigenfunction>> eigenfunction(const Scalar &E) const;

        std::vector<std::tuple<int, Scalar, std::vector<std::unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction>>>>
        eigenpairsByIndex(int iMin, int iMax) const;

        std::vector<std::tuple<int, Scalar, int>> eigenvaluesByIndex(int iMin, int iMax) const;

        std::vector<std::tuple<int, Scalar, int>> eigenvalues(const Scalar &eMin, const Scalar &eMax) const;

        std::vector<std::tuple<int, Scalar, std::vector<std::unique_ptr<typename PeriodicMatslise<Scalar>::Eigenfunction>>>>
        eigenpairs(const Scalar &eMin, const Scalar &eMax) const;
    };
}

#endif //SCHRODINGER_MATSLISE_H
