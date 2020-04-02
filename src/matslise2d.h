#ifndef MATSLISE_MATSLISE2D_H
#define MATSLISE_MATSLISE2D_H

#include "schrodinger.h"
#include "util/rectangle.h"

namespace matslise {
    template<typename _Scalar>
    class Matslise2D {
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
        std::vector<Matslise2D<Scalar>::Sector *> sectors;
        ArrayXs grid;
        int N;
        Scalar match;
        int matchIndex;
        Options2<Scalar> options;
        Y<Scalar, Eigen::Dynamic> y0Left;
        Y<Scalar, Eigen::Dynamic> y0Right;
    public:
        Matslise2D(const std::function<Scalar(const Scalar &, const Scalar &)> &V, const matslise::Rectangle<2, Scalar> &domain,
                   const Options2<Scalar> &options);

        std::pair<MatrixXs, MatrixXs> matchingErrorMatrix(const Scalar &E) const;

        std::pair<Scalar, Scalar> matchingError(const Scalar &E) const;

        std::vector<std::pair<Scalar, Scalar>> matchingErrors(const Scalar &E) const;

        std::vector<std::pair<Scalar, Scalar>> sortedErrors(const Scalar &E) const;

        Y<Scalar, Eigen::Dynamic>
        propagate(const Scalar &E, const Y<Scalar, Eigen::Dynamic> &y0, const Scalar &a, const Scalar &b,
                  bool use_h = true) const;

        Scalar eigenvalue(const Scalar &Eguess, const Scalar &tolerance = 1e-9, int maxIterations = 30,
                          const Scalar &minTolerance = 1e-5) const;

        std::vector<Scalar> eigenvalues(const Scalar &Emin, const Scalar &Emax, const int &initialSteps = 16) const;

        std::vector<Scalar> eigenvaluesByIndex(int Imin, int Imax) const;

        Scalar firstEigenvalue() const;

        std::vector<Scalar> firstEigenvalues(int n) const;

        std::vector<ArrayXXs>
        eigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const;

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunctionCalculator(const Scalar &E) const;

        bool contains(const Scalar &point) const {
            return point <= domain.max && point >= domain.min;
        }

        virtual ~Matslise2D();

    protected:
        std::vector<Y<Scalar, Eigen::Dynamic>> eigenfunctionSteps(const Scalar &E) const;

        MatrixXs calculateM(int k) const;

        MatrixXs conditionY(Y<Scalar, Eigen::Dynamic> &y) const;

    public:
        class Sector {
        public:
            const Matslise2D<Scalar> *se2d;
            AbstractMatslise<Scalar> *matslise;
            matslise::Matscs<Scalar> *matscs;
            ArrayXs vbar;
            Scalar min, max;

            Scalar *eigenvalues;
            ArrayXs *eigenfunctions;

            Sector(const Matslise2D<Scalar> *se2d, const Scalar &min, const Scalar &max, bool backward);

            template<int r>
            Y<Scalar, Eigen::Dynamic, r> propagate(
                    const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b,
                    bool use_h = true) const;

            bool contains(const Scalar &point) const {
                return point <= max && point >= min;
            }

            ArrayXs eigenfunction(int index, const ArrayXs &x) const;

            Scalar error() const;

            std::function<ArrayXs(Scalar)> basisCalculator() const;

            virtual ~Sector();

        public:
            MatrixXs calculateDeltaV(const Scalar &y) const;
        };
    };

    template<typename _Scalar=double>
    class Matslise2DHalf {
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
        Matslise2D<Scalar> *se2d;
        Rectangle<2, Scalar> domain;

    public:
        Matslise2DHalf(const std::function<Scalar(const Scalar &, const Scalar &)> &V,
                       const Rectangle<2, Scalar> &domain,
                       const Options2<Scalar> &options);

        Scalar eigenvalue(const Scalar &Eguess, const Scalar &tolerance = 1e-9, int maxIterations = 30,
                          const Scalar &minTolerance = 1e-5);

        std::vector<Scalar> eigenvalues(const Scalar &Emin, const Scalar &Emax, const int &initialSteps = 16);

        std::vector<Scalar> eigenvaluesByIndex(int Imin, int Imax);

        Scalar firstEigenvalue();

        std::vector<ArrayXXs> eigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y);

        std::pair<Scalar, Scalar> matchingError(const Scalar &E);

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunctionCalculator(const Scalar &E);

        virtual ~Matslise2DHalf();
    };
}

#endif //MATSLISE_MATSLISE2D_H
