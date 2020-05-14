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

        const std::function<Scalar(Scalar, Scalar)> potential;
        const Rectangle<2, Scalar> domain;

        AbstractMatslise2D(const std::function<Scalar(Scalar, Scalar)> &potential, const Rectangle<2, Scalar> &domain)
                : potential(potential), domain(domain) {
        }

        virtual Scalar firstEigenvalue() const = 0;

        virtual Scalar eigenvalue(const Scalar &Eguess) const = 0;

        virtual Scalar eigenvalueError(const Scalar &E) const = 0;

        virtual std::vector<Scalar> eigenvalues(const Scalar &Emin, const Scalar &Emax) const = 0;

        virtual std::vector<Scalar> eigenvaluesByIndex(int Imin, int Imax) const = 0;

        virtual std::vector<ArrayXXs> eigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const = 0;

        virtual std::vector<std::tuple<ArrayXXs, ArrayXXs, ArrayXXs>>
        eigenfunctionDerivatives(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const = 0;

        virtual std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunction(const Scalar &E) const = 0;

        virtual std::vector<std::function<std::tuple<Scalar, Scalar, Scalar>(Scalar, Scalar)>>
        eigenfunctionDerivatives(const Scalar &E) const = 0;

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

        using AbstractMatslise2D<Scalar>::domain;
        using AbstractMatslise2D<Scalar>::potential;
        MatrixXs *M;
        int sectorCount;
        std::vector<typename Matslise2D<Scalar>::Sector *> sectors;
        ArrayXs grid;
        int N;
        int matchIndex;
        Options2<Scalar> options;
        Y<Scalar, Eigen::Dynamic> dirichletBoundary;
    public:
        Matslise2D(const std::function<Scalar(Scalar, Scalar)> &potential,
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

        Scalar eigenvalueError(const Scalar &E) const override {
            return eigenvalueError(dirichletBoundary, E);
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

        std::vector<std::tuple<ArrayXXs, ArrayXXs, ArrayXXs>>
        eigenfunctionDerivatives(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const override {
            return eigenfunctionDerivatives(dirichletBoundary, E, x, y);
        }

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunction(const Scalar &E) const override {
            return eigenfunction(dirichletBoundary, E);
        }

        std::vector<std::function<std::tuple<Scalar, Scalar, Scalar>(Scalar, Scalar)>>
        eigenfunctionDerivatives(const Scalar &E) const override {
            return eigenfunctionDerivatives(dirichletBoundary, E);
        }

    public: // left boundary conditions
        std::pair<MatrixXs, MatrixXs> matchingErrorMatrix(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E, bool use_h = true) const;

        std::vector<std::pair<Scalar, Scalar>> matchingErrors(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E, bool use_h = true) const;

        std::pair<Scalar, Scalar> matchingError(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E, bool use_h = true) const;

        std::vector<Scalar> firstEigenvalues(const Y<Scalar, Eigen::Dynamic> &left, int n) const;

        Scalar firstEigenvalue(const Y<Scalar, Eigen::Dynamic> &left) const;

        Scalar eigenvalue(const Y<Scalar, Eigen::Dynamic> &left, const Scalar &Eguess, bool use_h = true) const;

        Scalar eigenvalueError(const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

        std::vector<Scalar> eigenvalues(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &Emin, const Scalar &Emax) const;

        std::vector<Scalar> eigenvaluesByIndex(const Y<Scalar, Eigen::Dynamic> &left, int Imin, int Imax) const;

        std::vector<ArrayXXs> eigenfunction(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E, const ArrayXs &x, const ArrayXs &y) const {
            return eigenfunctionHelper<false>(left, E, x, y);
        }

        std::vector<std::tuple<ArrayXXs, ArrayXXs, ArrayXXs>> eigenfunctionDerivatives(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E, const ArrayXs &x, const ArrayXs &y) const {
            return eigenfunctionHelper<true>(left, E, x, y);
        }

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunction(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const {
            return eigenfunctionHelper<false>(left, E);
        }

        std::vector<std::function<std::tuple<Scalar, Scalar, Scalar>(Scalar, Scalar)>> eigenfunctionDerivatives(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const {
            return eigenfunctionHelper<true>(left, E);
        }

        template<bool withDerivative, typename returnType=std::vector<std::function<
                typename std::conditional<withDerivative, std::tuple<Scalar, Scalar, Scalar>, Scalar>::type
                        (Scalar, Scalar)>>>
        returnType eigenfunctionHelper(const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

        template<bool withDerivative, typename returnType=std::vector<
                typename std::conditional<withDerivative, std::tuple<ArrayXXs, ArrayXXs, ArrayXXs>, ArrayXXs>::type>>
        returnType eigenfunctionHelper(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E, const ArrayXs &x, const ArrayXs &y) const;

    protected:
        std::vector<Y<Scalar, Eigen::Dynamic>> eigenfunctionSteps(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

        MatrixXs calculateM(int k) const;

        MatrixXs conditionY(Y<Scalar, Eigen::Dynamic> &y) const;

    public:
        class Sector {
        public:
            const Matslise2D<Scalar> *se2d;
            AbstractMatslise<Scalar> *matslise;
            typename matslise::Matscs<Scalar>::Sector *matscs;
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

            template<bool withDerivative, typename diffType=typename std::conditional<
                    withDerivative, std::pair<ArrayXXs, ArrayXXs>, ArrayXXs>::type>
            diffType basis(const ArrayXs &x) const;

            template<bool withDerivative, typename diffType=typename std::conditional<
                    withDerivative, std::pair<ArrayXs, ArrayXs>, ArrayXs>::type>
            std::function<diffType(Scalar)> basis() const;

            Scalar error() const;

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

        using AbstractMatslise2D<Scalar>::domain;
        using AbstractMatslise2D<Scalar>::potential;

    public:
        Matslise2DHalf(const std::function<Scalar(const Scalar &, const Scalar &)> &potential,
                       const Rectangle<2, Scalar> &domain,
                       const Options2<Scalar> &options);

        virtual ~Matslise2DHalf();

        Scalar firstEigenvalue() const override;

        Scalar eigenvalue(const Scalar &Eguess) const override;

        Scalar eigenvalueError(const Scalar &E) const override;

        std::vector<Scalar> eigenvalues(const Scalar &Emin, const Scalar &Emax) const override;

        std::vector<Scalar> eigenvaluesByIndex(int Imin, int Imax) const override;

        std::vector<ArrayXXs> eigenfunction(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const override {
            return eigenfunctionHelper<false>(E, x, y);
        }

        std::vector<std::tuple<ArrayXXs, ArrayXXs, ArrayXXs>>
        eigenfunctionDerivatives(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const override {
            return eigenfunctionHelper<true>(E, x, y);
        }

        std::vector<std::function<Scalar(Scalar, Scalar)>> eigenfunction(const Scalar &E) const override {
            return eigenfunctionHelper<false>(E);
        }

        std::vector<std::function<std::tuple<Scalar, Scalar, Scalar>(Scalar, Scalar)>>
        eigenfunctionDerivatives(const Scalar &E) const override {
            return eigenfunctionHelper<true>(E);
        }

        template<bool withDerivative, typename returnType=std::vector<std::function<
                typename std::conditional<withDerivative, std::tuple<Scalar, Scalar, Scalar>, Scalar>::type
                        (Scalar, Scalar)>>>
        returnType eigenfunctionHelper(const Scalar &E) const;

        template<bool withDerivative, typename returnType=std::vector<
                typename std::conditional<withDerivative, std::tuple<ArrayXXs, ArrayXXs, ArrayXXs>, ArrayXXs>::type>>
        returnType eigenfunctionHelper(const Scalar &E, const ArrayXs &x, const ArrayXs &y) const;
    };
}

#endif //MATSLISE_MATSLISE2D_H
