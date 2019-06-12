//
// Created by toon on 5/16/18.
//

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

#define MATSLISE_HMAX_delta 15
#define MATSLISE_ETA_delta 7
#define MATSLISE_ETA_h 9
#define MATSLISE_N 15

namespace matslise {
    template<typename Scalar>
    class AbstractMatslise {
    public:
        virtual std::vector<std::pair<int, Scalar>> *
        computeEigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left,
                           const matslise::Y<Scalar> &right) const = 0;

        virtual std::vector<std::pair<int, Scalar>> *
        computeEigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left) const {
            return computeEigenvalues(Emin, Emax, left, left);
        }

        virtual Scalar
        computeEigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left,
                               const matslise::Y<Scalar> &right) const = 0;

        virtual Scalar computeEigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left) const {
            return computeEigenvalueError(E, left, left);
        };

        virtual std::vector<std::pair<int, Scalar>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left,
                                  const matslise::Y<Scalar> &right) const = 0;

        virtual std::vector<std::pair<int, Scalar>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left) const {
            return computeEigenvaluesByIndex(Imin, Imax, left, left);
        };


        virtual Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1>
        computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                             const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) const = 0;

        virtual Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1>
        computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &left,
                             const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) const {
            return computeEigenfunction(E, left, left, x);
        }

        virtual std::function<Y<Scalar>(Scalar)> eigenfunctionCalculator(
                const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right) const = 0;

        virtual std::function<Y<Scalar>(Scalar)> eigenfunctionCalculator(
                const Scalar &E, const matslise::Y<Scalar> &left) const {
            return eigenfunctionCalculator(E, left, left);
        };

        virtual ~AbstractMatslise() {}
    };

    template<typename _Scalar=double>
    class Matslise : public AbstractMatslise<_Scalar> {
    public:
        typedef _Scalar Scalar;

        class Sector;

        std::function<Scalar(Scalar)> V;
        Scalar xmin, xmax;
        int sectorCount;
        Scalar match;
        Matslise::Sector **sectors;
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
        std::vector<std::pair<int, Scalar>> *
        computeEigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left,
                           const matslise::Y<Scalar> &right) const override;

        std::vector<std::pair<int, Scalar>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left,
                                  const matslise::Y<Scalar> &right) const override;

        Scalar computeEigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left,
                                      const matslise::Y<Scalar> &right) const override;

        Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1>
        computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                             const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) const override;

        std::function<Y<Scalar>(Scalar)>
        eigenfunctionCalculator(const Scalar &E, const matslise::Y<Scalar> &left,
                                const matslise::Y<Scalar> &right) const override;

        virtual ~Matslise();

    public:
        class Sector {
        public:
            Eigen::Array<Eigen::Matrix<Scalar, 2, 2>, MATSLISE_ETA_delta, MATSLISE_HMAX_delta>
                    t_coeff;
            Eigen::Matrix<Scalar, 2, 2> t_coeff_h[MATSLISE_ETA_h];
            Matslise<Scalar> *s;
            Scalar *vs;
            Scalar min, max, h;
            bool backward;

            Sector(Matslise *problem, const Scalar &min, const Scalar &max, bool backward);

            void calculateTCoeffs();

            bool contains(const Scalar &point) const {
                return point <= max && point >= min;
            }

            T<Scalar> calculateT(const Scalar &E, bool use_h = true) const;

            T<Scalar> calculateT(const Scalar &E, const Scalar &delta, bool use_h = true) const;

            matslise::Y<Scalar>
            propagate(const Scalar &E, const matslise::Y<Scalar> &y, bool forward, bool use_h = true) const;

            matslise::Y<Scalar>
            propagate(const Scalar &E, const Y<Scalar> &y0, const Scalar &a, const Scalar &b, Scalar &theta,
                      bool use_h = true) const;

            Scalar prufer(const Scalar &E, const Scalar &delta, const Y<Scalar> &y0, const Y<Scalar> &y1) const;

            Scalar calculateError() const;

            virtual ~Sector();
        };
    };

    template<typename _Scalar>
    class HalfRange : public AbstractMatslise<_Scalar> {
        typedef _Scalar Scalar;
    public:
        static void checkSymmetry(const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right) {
            if (left != right)
                throw std::runtime_error("Halfrange::checkSymmetry(), left and right sides have to be identical");
        }

    public:
        const Matslise<Scalar> *ms;
        static const int AUTO = -1;
        static const int ODD = 0;
        static const int EVEN = 1;
    public:
        HalfRange(std::function<Scalar(Scalar)> V, const Scalar &xmax,
                  std::shared_ptr<matslise::SectorBuilder<Matslise<Scalar>>> sectorBuilder);

        using AbstractMatslise<Scalar>::computeEigenvalues;

        std::vector<std::pair<int, Scalar>> *
        computeEigenvalues(const Scalar &Emin, const Scalar &Emax, const matslise::Y<Scalar> &left,
                           const matslise::Y<Scalar> &right) const override;

        using AbstractMatslise<Scalar>::computeEigenvaluesByIndex;

        std::vector<std::pair<int, Scalar>> *
        computeEigenvaluesByIndex(int Imin, int Imax, const matslise::Y<Scalar> &left,
                                  const matslise::Y<Scalar> &right) const override;


        using AbstractMatslise<Scalar>::computeEigenvalueError;

        Scalar computeEigenvalueError(const Scalar &E, const matslise::Y<Scalar> &side, int even) const;

        Scalar computeEigenvalueError(const Scalar &E, const matslise::Y<Scalar> &left,
                                      const matslise::Y<Scalar> &right) const override {
            checkSymmetry(left, right);
            return computeEigenvalueError(E, left, AUTO);
        }

        using AbstractMatslise<Scalar>::computeEigenfunction;

        Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1>
        computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &left, const matslise::Y<Scalar> &right,
                             const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) const override {
            checkSymmetry(left, right);
            return computeEigenfunction(E, left, x, AUTO);
        }

        Eigen::Array<matslise::Y<Scalar>, Eigen::Dynamic, 1>
        computeEigenfunction(const Scalar &E, const matslise::Y<Scalar> &side,
                             const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x,
                             int even) const;

        using AbstractMatslise<Scalar>::eigenfunctionCalculator;

        std::function<Y<Scalar>(Scalar)>
        eigenfunctionCalculator(const Scalar &E, const matslise::Y<Scalar> &left,
                                const matslise::Y<Scalar> &right) const override {
            checkSymmetry(left, right);
            return eigenfunctionCalculator(E, left, AUTO);
        }

        std::function<Y<Scalar>(Scalar)>
        eigenfunctionCalculator(const Scalar &E, const matslise::Y<Scalar> &side, int even) const;

        virtual ~HalfRange();
    };
}

#include "matscs.h"
#include "se2d.h"
#include "util/SectorBuilder.h"

#endif //SCHRODINGER_MATSLISE_H