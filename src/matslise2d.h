#ifndef MATSLISE_MATSLISE2D_H
#define MATSLISE_MATSLISE2D_H

#include "schrodinger.h"
#include "util/rectangle.h"

namespace matslise {
    template<typename Scalar>
    class AbstractMatslise2D {
    public:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> ArrayXXs;

        virtual Scalar firstEigenvalue() const = 0;

        virtual Scalar eigenvalue(const Scalar &Eguess) const = 0;

        virtual std::vector<Scalar> eigenvalues(const Scalar &Emin, const Scalar &Emax) const = 0;

        virtual std::vector<Scalar> eigenvaluesByIndex(int Imin, int Imax) const = 0;

        virtual std::vector<ArrayXXs> eigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const = 0;

        virtual std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunctionCalculator(const Scalar &E) const = 0;

        virtual ~AbstractMatslise2D() = default;
    };

    template<typename _Scalar>
    class Matslise2D : public AbstractMatslise2D<_Scalar> {
    public:
        typedef _Scalar Scalar;
        using typename AbstractMatslise2D<Scalar>::VectorXs;
        using typename AbstractMatslise2D<Scalar>::MatrixXs;
        using typename AbstractMatslise2D<Scalar>::ArrayXs;
        using typename AbstractMatslise2D<Scalar>::ArrayXXs;
        static const int order = matslise::Matscs<Scalar>::order;

        class Sector;

        std::function<Scalar(const Scalar &, const Scalar &)> V;
        MatrixXs *M;
        Rectangle<2, Scalar> domain;
        int sectorCount;
        std::vector<typename Matslise2D<Scalar>::Sector *> sectors;
        ArrayXs grid;
        int N;
        int matchIndex;
        Options2<Scalar> options;
        Y<Scalar, Eigen::Dynamic> dirichletBoundary;
    public:
        Matslise2D(const std::function<Scalar(const Scalar &, const Scalar &)> &V,
                   const matslise::Rectangle<2, Scalar> &domain,
                   const Options2<Scalar> &options);

        virtual ~Matslise2D();

        Y<Scalar, Eigen::Dynamic>
        propagate(const Scalar &E, const Y<Scalar, Eigen::Dynamic> &y0, const Scalar &a, const Scalar &b,
                  bool use_h = true) const;

        std::pair<MatrixXs, MatrixXs> matchingErrorMatrix(const Scalar &E) const {
            return matchingErrorMatrix(dirichletBoundary, E);
        }

        std::vector<std::pair<Scalar, Scalar>> matchingErrors(const Scalar &E) const {
            return matchingErrors(dirichletBoundary, E);
        }

        std::pair<Scalar, Scalar> matchingError(const Scalar &E) const {
            return matchingError(dirichletBoundary, E);
        }

        std::vector<Scalar> firstEigenvalues(int n) const {
            return firstEigenvalues(dirichletBoundary, n);
        }

    public: // overrides
        Scalar firstEigenvalue() const override {
            return firstEigenvalue(dirichletBoundary);
        }

        Scalar eigenvalue(const Scalar &Eguess) const override {
            return eigenvalue(dirichletBoundary, Eguess);
        }

        std::vector<Scalar> eigenvalues(const Scalar &Emin, const Scalar &Emax) const override {
            return eigenvalues(dirichletBoundary, Emin, Emax);
        }

        std::vector<Scalar> eigenvaluesByIndex(int Imin, int Imax) const override {
            return eigenvaluesByIndex(dirichletBoundary, Imin, Imax);
        }

        std::vector<ArrayXXs> eigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const override {
            return eigenfunction(dirichletBoundary, E, x, y);
        }

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunctionCalculator(const Scalar &E) const override {
            return eigenfunctionCalculator(dirichletBoundary, E);
        }

    public: // left boundary conditions
        std::pair<MatrixXs, MatrixXs> matchingErrorMatrix(const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

        std::vector<std::pair<Scalar, Scalar>> matchingErrors(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

        std::pair<Scalar, Scalar> matchingError(const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

        std::vector<Scalar> firstEigenvalues(const Y<Scalar, Eigen::Dynamic> &left, int n) const;

        Scalar firstEigenvalue(const Y<Scalar, Eigen::Dynamic> &left) const;

        Scalar eigenvalue(const Y<Scalar, Eigen::Dynamic> &left, const Scalar &Eguess) const;

        std::vector<Scalar> eigenvalues(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &Emin, const Scalar &Emax) const;

        std::vector<Scalar> eigenvaluesByIndex(const Y<Scalar, Eigen::Dynamic> &left, int Imin, int Imax) const;

        std::vector<ArrayXXs> eigenfunction(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E, const ArrayXs &x, const ArrayXs &y) const;

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunctionCalculator(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

    protected:
        std::vector<Y<Scalar, Eigen::Dynamic>> eigenfunctionSteps(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

        MatrixXs calculateM(int k) const;

        MatrixXs conditionY(Y<Scalar, Eigen::Dynamic> &y) const;

    public:
        class Sector {
        public:
            const Matslise2D<Scalar> *se2d;
            AbstractMatslise <Scalar> *matslise;
            matslise::Matscs<Scalar> *matscs;
            ArrayXs vbar;
            Scalar min, max;

            Scalar *eigenvalues;
            ArrayXs *eigenfunctions;

            Sector(const Matslise2D<Scalar> *se2d, const Scalar &min, const Scalar &max, bool backward);

            virtual ~Sector();

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


            static bool compare(const Sector &a, const Sector &b) {
                return a.vbar.minCoeff() < b.vbar.minCoeff();
            }

        public:
            MatrixXs calculateDeltaV(const Scalar &y) const;
        };
    };

    template<typename _Scalar=double>
    class Matslise2DHalf : public AbstractMatslise2D<_Scalar> {
    public:
        typedef _Scalar Scalar;
        using typename AbstractMatslise2D<_Scalar>::VectorXs;
        using typename AbstractMatslise2D<_Scalar>::MatrixXs;
        using typename AbstractMatslise2D<_Scalar>::ArrayXs;
        using typename AbstractMatslise2D<_Scalar>::ArrayXXs;
    private:
        Y<Scalar, Eigen::Dynamic, Eigen::Dynamic> neumannBoundary;
        Y<Scalar, Eigen::Dynamic, Eigen::Dynamic> dirichletBoundary;
    public:
        Matslise2D<Scalar> *se2d;
        Rectangle<2, Scalar> domain;

    public:
        Matslise2DHalf(const std::function<Scalar(const Scalar &, const Scalar &)> &V,
                       const Rectangle<2, Scalar> &domain,
                       const Options2<Scalar> &options);

        virtual ~Matslise2DHalf();

        Scalar firstEigenvalue() const override;

        Scalar eigenvalue(const Scalar &Eguess) const override;

        std::vector<Scalar> eigenvalues(const Scalar &Emin, const Scalar &Emax) const override;

        std::vector<Scalar> eigenvaluesByIndex(int Imin, int Imax) const override;

        std::vector<ArrayXXs> eigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const override;

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunctionCalculator(const Scalar &E) const override;
    };
}

#endif //MATSLISE_MATSLISE2D_H
