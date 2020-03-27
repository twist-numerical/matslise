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
#include "util/SectorBuilder_header.h"

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
        virtual Scalar estimatePotentialMinimum() const = 0;

        virtual std::vector<std::pair<int, Scalar>>
        computeEigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left,
                           const matslise::Y<Scalar> &right) const {
            checkSymmetry(left, right);
            return computeEigenvalues(Emin, Emax, right);
        }


        virtual std::vector<std::pair<int, Scalar>>
        computeEigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left) const {
            return computeEigenvalues(Emin, Emax, left, left);
        }

        virtual std::vector<std::pair<int, Scalar>>
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left,
                                  const matslise::Y<Scalar> &right) const {
            checkSymmetry(left, right);
            return computeEigenvaluesByIndex(Imin, Imax, right);
        }

        virtual std::vector<std::pair<int, Scalar>>
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &side) const {
            return computeEigenvaluesByIndex(Imin, Imax, side, side);
        };

        virtual Scalar
        computeEigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left,
                               const matslise::Y<Scalar> &right, int index = -1) const {
            checkSymmetry(left, right);
            return computeEigenvalueError(E, right, index);
        }

        virtual Scalar computeEigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left, int index = -1) const {
            return computeEigenvalueError(E, left, left, index);
        };

        virtual Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1>
        computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                             const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x, int index = -1) const {
            checkSymmetry(left, right);
            return computeEigenfunction(E, right, x, index);
        }

        virtual Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1>
        computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &side,
                             const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x, int index = -1) const {
            return computeEigenfunction(E, side, side, x, index);
        }

        virtual std::function<Y<Scalar>(Scalar)> eigenfunctionCalculator(
                const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                int index = -1) const {
            checkSymmetry(left, right);
            return eigenfunctionCalculator(E, right, index);
        };

        virtual std::function<Y<Scalar>(Scalar)> eigenfunctionCalculator(
                const Scalar &E, const matslise::Y<Scalar> &left, int index = -1) const {
            return eigenfunctionCalculator(E, left, left, index);
        };

        virtual ~AbstractMatslise() = default;
    };

    template<typename _Scalar=double>
    class Matslise : public AbstractMatslise<_Scalar> {
    public:
        typedef _Scalar Scalar;
        static const int order = 16;

        class Sector;

        std::function<Scalar(Scalar)> V;
        Scalar xmin, xmax;
        int sectorCount;
        Scalar match;
        std::vector<Matslise::Sector *> sectors;
    public:
        static std::shared_ptr<matslise::SectorBuilder<Matslise<Scalar>>> UNIFORM(int sectorCount) {
            return std::shared_ptr<matslise::SectorBuilder<Matslise<Scalar>>>(
                    new matslise::sectorbuilder::Uniform<Matslise<Scalar>>(sectorCount));
        }

        static std::shared_ptr<matslise::SectorBuilder<Matslise<Scalar>>> AUTO(const Scalar &tolerance) {
            return std::shared_ptr<matslise::SectorBuilder<Matslise<Scalar>>>(
                    new matslise::sectorbuilder::Auto<Matslise<Scalar>>(tolerance));
        }

        Matslise(std::function<Scalar(const Scalar &)> V, const Scalar &xmin, const Scalar &xmax,
                 std::shared_ptr<matslise::SectorBuilder<Matslise<Scalar>>> sectorBuilder) : V(V), xmin(xmin),
                                                                                             xmax(xmax) {
            sectorBuilder->build(this, xmin, xmax);
            sectorCount = sectors.size();
        }

        Matslise(std::function<Scalar(const Scalar &)> V, const Scalar &xmin, const Scalar &xmax, int sectorCount)
                : Matslise(V, xmin, xmax, UNIFORM(sectorCount)) {};

        bool contains(const Scalar &point) const {
            return point <= xmax && point >= xmin;
        }

        std::pair<matslise::Y<Scalar>, Scalar>
        propagate(const Scalar &E, const matslise::Y<Scalar> &y, const Scalar &a, const Scalar &b,
                  bool use_h = true) const;

        std::tuple<Scalar, Scalar, Scalar>
        calculateError(const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                       bool use_h = true) const;

    public: // Override
        Scalar estimatePotentialMinimum() const override;

        using AbstractMatslise<Scalar>::computeEigenvalues;

        std::vector<std::pair<int, Scalar>>
        computeEigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left,
                           const matslise::Y<Scalar> &right) const override;

        using AbstractMatslise<Scalar>::computeEigenvaluesByIndex;

        std::vector<std::pair<int, Scalar>>
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left,
                                  const matslise::Y<Scalar> &right) const override;

        using AbstractMatslise<Scalar>::computeEigenvalueError;

        Scalar computeEigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left,
                                      const matslise::Y<Scalar> &right, int index = -1) const override;

        using AbstractMatslise<Scalar>::computeEigenfunction;

        Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1>
        computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                             const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x, int index = -1) const override;

        using AbstractMatslise<Scalar>::eigenfunctionCalculator;

        std::function<Y<Scalar>(Scalar)>
        eigenfunctionCalculator(const Scalar &E, const matslise::Y<Scalar> &left,
                                const matslise::Y<Scalar> &right, int index = -1) const override;

        virtual ~Matslise();

    public:
        class Sector {
        public:
            Eigen::Array<Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign>, MATSLISE_ETA_delta, MATSLISE_HMAX_delta, Eigen::DontAlign> t_coeff;
            Eigen::Matrix<Scalar, 2, 2, Eigen::DontAlign> t_coeff_h[MATSLISE_ETA_h];
            const Matslise<Scalar> *s;
            Scalar *vs;
            Scalar min, max, h;
            bool backward;

            Sector(const Matslise *problem, const Scalar &min, const Scalar &max, bool backward);

            void calculateTCoeffs();

            bool contains(const Scalar &point) const {
                return point <= max && point >= min;
            }

            T<Scalar> calculateT(const Scalar &E, bool use_h = true) const;

            T<Scalar> calculateT(const Scalar &E, const Scalar &delta, bool use_h = true) const;

            std::pair<matslise::Y<Scalar>, Scalar>
            propagate(const Scalar &E, const Y<Scalar> &y0, const Scalar &a, const Scalar &b, bool use_h = true) const;

            std::pair<matslise::Y<Scalar>, Scalar>
            propagateForward(const Scalar &E, const matslise::Y<Scalar> &y, bool use_h = true) const {
                return propagate(E, y, min, max, use_h);
            }

            std::pair<matslise::Y<Scalar>, Scalar>
            propagateBackward(const Scalar &E, const matslise::Y<Scalar> &y, bool use_h = true) const {
                return propagate(E, y, max, min, use_h);
            }

            Scalar prufer(const Scalar &E, const Scalar &delta, const Y<Scalar> &y0, const Y<Scalar> &y1) const;

            Scalar error() const;

            virtual ~Sector();

        private :
            std::pair<Y<Scalar>, Scalar>
            propagateDelta(
                    const Scalar &E, const Y<Scalar> &y0, Scalar delta, bool use_h) const;

        };
    };

    template<typename _Scalar = double>
    class HalfRange : public AbstractMatslise<_Scalar> {
        typedef _Scalar Scalar;

    public:
        const Matslise<Scalar> *ms;
    public:
        HalfRange(std::function<Scalar(Scalar)> V, const Scalar &xmax,
                  std::shared_ptr<matslise::SectorBuilder<Matslise<Scalar>>> sectorBuilder);

        Scalar estimatePotentialMinimum() const override {
            return ms->estimatePotentialMinimum();
        }

        using AbstractMatslise<Scalar>::computeEigenvalues;

        std::vector<std::pair<int, Scalar>>
        computeEigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &side) const override;

        using AbstractMatslise<Scalar>::computeEigenvaluesByIndex;

        std::vector<std::pair<int, Scalar>>
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &side) const override;

        using AbstractMatslise<Scalar>::computeEigenvalueError;

        Scalar computeEigenvalueError(const Scalar &E, const matslise::Y<Scalar> &side, int index = -1) const override;

        using AbstractMatslise<Scalar>::computeEigenfunction;

        Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1>
        computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &side,
                             const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x, int index = -1) const override;

        using AbstractMatslise<Scalar>::eigenfunctionCalculator;

        std::function<Y<Scalar>(Scalar)>
        eigenfunctionCalculator(const Scalar &E, const matslise::Y<Scalar> &side, int index = -1) const override;

        virtual ~HalfRange();
    };
}

#include "matscs.h"
#include "se2d.h"
#include "util/SectorBuilder.h"

#endif //SCHRODINGER_MATSLISE_H