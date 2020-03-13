#ifndef MATSLISE_SE2D_H
#define MATSLISE_SE2D_H

#include "schrodinger.h"

namespace matslise {
    namespace SEnD_util {
        template<typename Scalar=double>
        bool NEWTON_RAPHSON_SORTER(const std::pair<Scalar, Scalar> &a, const std::pair<Scalar, Scalar> &b) {
            if (abs(a.first) > 100 || abs(b.first) > 100)
                return abs(a.first) < abs(b.first);
            return abs(a.first / a.second) < abs(b.first / b.second);
        }

        template<typename Scalar=double>
        bool ABS_SORTER(const std::pair<Scalar, Scalar> &a, const std::pair<Scalar, Scalar> &b) {
            return abs(a.first) < abs(b.first);
        }
    }

    template<typename _Scalar>
    class SE2D {
    public:
        typedef _Scalar Scalar;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> ArrayXXs;
        static const int order = matslise::Matscs<Scalar>::order;

        class Sector;

        std::function<Scalar(const Scalar &, const Scalar &)> V;
        MatrixXs *M;
        Rectangle<2, Scalar> domain;
        int sectorCount;
        std::vector<SE2D<Scalar>::Sector *> sectors;
        ArrayXs grid;
        int N;
        Scalar match;
        Options2<Scalar> options;
        Y<Scalar, Eigen::Dynamic> y0Left;
        Y<Scalar, Eigen::Dynamic> y0Right;
    public:
        SE2D(const std::function<Scalar(const Scalar &, const Scalar &)> &V, const Rectangle<2, Scalar> &domain,
             const Options2<Scalar> &options);

        std::pair<MatrixXs, MatrixXs> calculateErrorMatrix(const Scalar &E) const;

        std::pair<Scalar, Scalar> calculateError(const Scalar &E, const std::function<bool(
                std::pair<Scalar, Scalar>,
                std::pair<Scalar, Scalar>)> &sorter = SEnD_util::NEWTON_RAPHSON_SORTER<Scalar>) const;

        std::vector<std::pair<Scalar, Scalar>> calculateErrors(const Scalar &E) const;

        std::vector<std::pair<Scalar, Scalar>> sortedErrors(
                const Scalar &E,
                const std::function<bool(std::pair<Scalar, Scalar>, std::pair<Scalar, Scalar>)> &sorter
                = SEnD_util::NEWTON_RAPHSON_SORTER<Scalar>) const;

        Y<Scalar, Eigen::Dynamic>
        propagate(const Scalar &E, const Y<Scalar, Eigen::Dynamic> &y0, const Scalar &a, const Scalar &b,
                  bool use_h = true) const;

        Scalar findEigenvalue(const Scalar &Eguess, const Scalar &tolerance = 1e-9, int maxIterations = 30,
                              const Scalar &minTolerance = 1e-5) const;

        std::vector<Scalar> findEigenvalues(
                const Scalar &Emin, const Scalar &Emax,
                const int &initialSteps = 16) const;

        Scalar findFirstEigenvalue() const;

        std::vector<ArrayXXs>
        computeEigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const;

        std::vector<Scalar> computeEigenvaluesByIndex(int Imin, int Imax) const;
        // std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunctionCalculator(const Scalar &E) const;

        bool contains(const Scalar &point) const {
            return point <= domain.max && point >= domain.min;
        }

        virtual ~SE2D();

    protected:
        std::vector<Y<Scalar, Eigen::Dynamic>> computeEigenfunctionSteps(const Scalar &E) const;

        MatrixXs calculateM(int k) const;

        MatrixXs conditionY(Y<Scalar, Eigen::Dynamic> &y) const;

    public:
        class Sector {
        public:
            const SE2D<Scalar> *se2d;
            AbstractMatslise<Scalar> *matslise;
            matslise::Matscs<Scalar> *matscs;
            ArrayXs vbar;
            Scalar min, max;

            Scalar *eigenvalues;
            ArrayXs *eigenfunctions;

            Sector(const SE2D<Scalar> *se2d, const Scalar &min, const Scalar &max, bool backward);

            template<int r>
            Y<Scalar, Eigen::Dynamic, r> propagate(
                    const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b,
                    bool use_h = true) const;

            bool contains(const Scalar &point) const {
                return point <= max && point >= min;
            }

            ArrayXs computeEigenfunction(int index, const ArrayXs &x) const;

            Scalar calculateError() const;

            std::function<ArrayXs(Scalar)> basisCalculator() const;

            virtual ~Sector();

        public:
            MatrixXs calculateDeltaV(const Scalar &y) const;
        };
    };

    template<typename _Scalar=double>
    class SE2DHalf {
    public:
        typedef _Scalar Scalar;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> ArrayXXs;
    private:
        Y<Scalar, Eigen::Dynamic, Eigen::Dynamic> oddBoundary;
        Y<Scalar, Eigen::Dynamic, Eigen::Dynamic> evenBoundary;

        void setParity(bool even);

    public:
        SE2D<Scalar> *se2d;
        Rectangle<2, Scalar> domain;

    public:
        SE2DHalf(const std::function<Scalar(const Scalar &, const Scalar &)> &V, const Rectangle<2, Scalar> &domain,
                 const Options2<Scalar> &options);

        Scalar findEigenvalue(const Scalar &Eguess, const Scalar &tolerance = 1e-9, int maxIterations = 30,
                              const Scalar &minTolerance = 1e-5);

        std::vector<Scalar> findEigenvalues(
                const Scalar &Emin, const Scalar &Emax,
                const int &initialSteps = 16);

        Scalar findFirstEigenvalue();

        std::vector<ArrayXXs>
        computeEigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y);

        std::pair<Scalar, Scalar> calculateError(const Scalar &E, const std::function<bool(
                std::pair<Scalar, Scalar>,
                std::pair<Scalar, Scalar>)> &sorter = SEnD_util::NEWTON_RAPHSON_SORTER<Scalar>);

        std::vector<Scalar> computeEigenvaluesByIndex(int Imin, int Imax);
        // std::vector<double> *computeEigenvalues(double Emin, double Emax) const;

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunctionCalculator(const Scalar &E);

        virtual ~SE2DHalf();
    };
}

#endif //MATSLISE_SE2D_H
