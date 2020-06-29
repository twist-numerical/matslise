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
#include "util/rectangle.h"

#define MATSLISE_HMAX_delta 17
#define MATSLISE_ETA_delta 8
#define MATSLISE_ETA_h 10
#define MATSLISE_N 16

namespace matslise {
    template<typename Scalar>
    class AbstractMatslise {
    public:
        static void checkSymmetry(const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right) {
            if (left != right)
                throw std::runtime_error("Halfrange::checkSymmetry(), left and right sides have to be identical");
        }

    public:
        const std::function<Scalar(Scalar)> potential;
        const Rectangle<1, Scalar> domain;

        AbstractMatslise(const std::function<Scalar(Scalar)> &potential, const Rectangle<1, Scalar> &domain)
                : potential(potential), domain(domain) {
        }

        struct Eigenfunction {
        public:
            std::function<Y<Scalar>(const Scalar &)> f_scalar;
            std::function<Eigen::Array<Y<Scalar>, Eigen::Dynamic, 1>(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &)>
                    f_array;

            Y<Scalar> operator()(const Scalar &x) const {
                return f_scalar(x);
            }

            Eigen::Array<Y<Scalar>, Eigen::Dynamic, 1>
            operator()(const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) const {
                return f_array(x);
            }
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

        virtual Scalar
        eigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left,
                        const matslise::Y<Scalar> &right, int index = -1) const {
            checkSymmetry(left, right);
            return eigenvalueError(E, right, index);
        }

        virtual Scalar eigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left, int index = -1) const {
            return eigenvalueError(E, left, left, index);
        };

        virtual Eigenfunction eigenfunction(
                const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                int index = -1) const {
            checkSymmetry(left, right);
            return eigenfunction(E, right, index);
        };

        virtual Eigenfunction eigenfunction(
                const Scalar &E, const matslise::Y<Scalar> &left, int index = -1) const {
            return eigenfunction(E, left, left, index);
        };

        virtual ~AbstractMatslise() = default;
    };

    template<typename _Scalar=double>
    class Matslise : public AbstractMatslise<_Scalar> {
    public:
        typedef _Scalar Scalar;
        static const int order = 16;

        class Sector;

        using AbstractMatslise<Scalar>::domain;
        using AbstractMatslise<Scalar>::potential;
        typedef typename AbstractMatslise<Scalar>::Eigenfunction Eigenfunction;

        int sectorCount;
        int matchIndex;
        std::vector<Matslise::Sector *> sectors;

        Scalar tolerance;
    public:
        Matslise(std::function<Scalar(const Scalar &)> V, const Scalar &xmin, const Scalar &xmax,
                 const Scalar &tolerance = 1e-8)
                : Matslise(V, xmin, xmax, tolerance, sector_builder::automatic<Matslise<Scalar>>(tolerance)) {
        }

        Matslise(std::function<Scalar(const Scalar &)> V, const Scalar &xmin, const Scalar &xmax,
                 const Scalar &tolerance, SectorBuilder<Matslise<Scalar>> sectorBuilder)
                : AbstractMatslise<Scalar>(V, {xmin, xmax}), tolerance(tolerance) {
            auto buildSectors = sectorBuilder(this, domain.min, domain.max);
            sectors = std::move(buildSectors.sectors);
            matchIndex = buildSectors.matchIndex;
            sectorCount = sectors.size();
        }

        std::pair<matslise::Y<Scalar>, Scalar>
        propagate(const Scalar &E, const matslise::Y<Scalar> &y, const Scalar &a, const Scalar &b,
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

        Eigenfunction eigenfunction(const Scalar &E, const matslise::Y<Scalar> &left,
                                    const matslise::Y<Scalar> &right, int index = -1) const override;

        virtual ~Matslise();

    public:
        class Sector {
        public:
            static const bool expensive = false;

            Eigen::Array<Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign>, MATSLISE_ETA_delta, MATSLISE_HMAX_delta, Eigen::DontAlign> t_coeff;
            Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign> t_coeff_h[MATSLISE_ETA_h];
            const Matslise<Scalar> *s;
            std::array<Scalar, MATSLISE_N> vs;
            Scalar min, max, h;
            bool backward;

            Sector(const Matslise *problem, const Scalar &min, const Scalar &max, bool backward);

            void calculateTCoeffs();

            bool contains(const Scalar &point) const {
                return point <= max && point >= min;
            }

            T<Scalar> calculateT(const Scalar &E, bool use_h = true) const;

            T<Scalar> calculateT(const Scalar &E, const Scalar &delta, bool use_h = true) const;

            template<bool withPrufer>
            typename std::conditional<withPrufer, std::pair<matslise::Y<Scalar>, Scalar>, matslise::Y<Scalar>>::type
            propagate(const Scalar &E, const Y<Scalar> &y0, const Scalar &a, const Scalar &b, bool use_h = true) const;

            Scalar theta0(const Scalar &E, const Y<Scalar> &y0) const;

            Scalar prufer(const Scalar &E, const Scalar &delta, const Y<Scalar> &y0, const Y<Scalar> &y1) const;

            Scalar error() const;

            ~Sector() = default;

            static bool compare(const Sector &a, const Sector &b) {
                return a.vs[0] < b.vs[0];
            }
        };
    };

    template<typename _Scalar = double>
    class MatsliseHalf : public AbstractMatslise<_Scalar> {
        typedef _Scalar Scalar;
        typedef typename AbstractMatslise<Scalar>::Eigenfunction Eigenfunction;

    public:
        const Matslise<Scalar> *ms;
    public:
        MatsliseHalf(std::function<Scalar(Scalar)> V, const Scalar &xmax, const Scalar &tolerance)
                : MatsliseHalf(V, xmax, tolerance, sector_builder::automatic<Matslise<Scalar>>(tolerance)) {}

        MatsliseHalf(std::function<Scalar(Scalar)> V, const Scalar &xmax, const Scalar &tolerance,
                     matslise::SectorBuilder<Matslise<Scalar>> sectorBuilder);

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

        Eigenfunction eigenfunction(const Scalar &E, const matslise::Y<Scalar> &side, int index = -1) const override;

        virtual ~MatsliseHalf();
    };
}

#include "matscs.h"
#include "matslise2d.h"
#include "matslise3d.h"

#endif //SCHRODINGER_MATSLISE_H
