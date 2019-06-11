//
// Created by toon on 6/13/18.
//

#ifndef MATSLISE_SE2D_H
#define MATSLISE_SE2D_H

#include "schrodinger.h"

namespace matslise {
    namespace SEnD_util {
        template<typename Scalar>
        bool NEWTON_RAPHSON_SORTER(const std::pair<Scalar, Scalar> &a, const std::pair<Scalar, Scalar> &b) {
            if (abs(a.first) > 100 || abs(b.first) > 100)
                return abs(a.first) < abs(b.first);
            return abs(a.first / a.second) < abs(b.first / b.second);
        }

        template<typename Scalar>
        bool ABS_SORTER(const std::pair<Scalar, Scalar> &a, const std::pair<Scalar, Scalar> &b) {
            return abs(a.first) < abs(b.first);
        }
    }

    template<typename _Scalar>
    class SE2D {
    public:
        typedef _Scalar Scalar;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> ArrayXXs;

        class Sector;

        std::function<Scalar(Scalar, Scalar)> V;
        MatrixXs *M;
        Rectangle<2, Scalar> domain;
        int sectorCount;
        SE2D<Scalar>::Sector **sectors;
        ArrayXs grid;
        int N;
        Scalar match;
        Options2<Scalar> options;
    public:
        SE2D(const std::function<Scalar(Scalar, Scalar)> &V, const Rectangle<2, Scalar> &domain,
             const Options2<Scalar> &options);

        std::pair<MatrixXs, MatrixXs> calculateErrorMatrix(Scalar E) const;

        std::pair<Scalar, Scalar> calculateError(Scalar E, const std::function<bool(
                std::pair<Scalar, Scalar>,
                std::pair<Scalar, Scalar>)> &sorter = SEnD_util::NEWTON_RAPHSON_SORTER<Scalar>) const;

        std::vector<std::pair<Scalar, Scalar>> *calculateErrors(Scalar E) const;

        std::vector<std::pair<Scalar, Scalar>> *sortedErrors(Scalar E, const std::function<bool(
                std::pair<Scalar, Scalar>,
                std::pair<Scalar, Scalar>)> &sorter = SEnD_util::NEWTON_RAPHSON_SORTER<Scalar>) const;

        Y<Scalar, Eigen::Dynamic>
        propagate(Scalar E, const Y<Scalar, Eigen::Dynamic> &y0, Scalar a, Scalar b, bool use_h = true) const;

        Scalar findEigenvalue(Scalar Eguess, Scalar tolerance = 1e-9, int maxIterations = 30,
                              Scalar minTolerance = 1e-5) const;

        std::vector<ArrayXXs> *
        computeEigenfunction(Scalar E, const ArrayXs &x, const ArrayXs &y) const;

        std::vector<Scalar> *computeEigenvaluesByIndex(int Imin, int Imax) const;
        // std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

        bool contains(Scalar point) const {
            return point <= domain.max && point >= domain.min;
        }

        virtual ~SE2D();

    protected:
        Y<Scalar, Eigen::Dynamic> *computeEigenfunctionSteps(Scalar E) const;

        MatrixXs calculateM(int k) const;

    public:
        class Sector {
        public:
            SE2D<Scalar> *se2d;
            AbstractMatslise<Scalar> *matslise;
            matslise::Matscs<Scalar> *matscs;
            ArrayXs vbar;
            Scalar min, max;

            Scalar *eigenvalues;
            ArrayXs *eigenfunctions;

            Sector(SE2D<Scalar> *se2d, Scalar min, Scalar max, bool backward);

            Y<Scalar, Eigen::Dynamic> propagate(Scalar E, const Y<Scalar, Eigen::Dynamic> &y0, Scalar a, Scalar b,
                                                bool use_h = true) const;

            bool contains(Scalar point) const {
                return point <= max && point >= min;
            }

            ArrayXs computeEigenfunction(int index, const ArrayXs &x) const;

            Scalar calculateError() const;

            virtual ~Sector();

        private:
            MatrixXs calculateDeltaV(Scalar y) const;
        };
    };
}

#endif //MATSLISE_SE2D_H
