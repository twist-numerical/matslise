#ifndef MATSLISE_MATSLISE2D_H
#define MATSLISE_MATSLISE2D_H

#include "schrodinger.h"
#include "util/rectangle.h"
#include "matslise2d/basisQuadrature.h"

namespace matslise {
    template<typename Scalar=double, bool withDerivatives = false>
    struct Eigenfunction2D {
    public:
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> ArrayXXs;
        typedef typename std::conditional<withDerivatives, std::tuple<Scalar, Scalar, Scalar>, Scalar>::type ScalarReturn;
        typedef typename std::conditional<withDerivatives, std::tuple<ArrayXXs, ArrayXXs, ArrayXXs>, ArrayXXs>::type ArrayReturn;

        std::function<ScalarReturn(const Scalar &, const Scalar &)> f_scalar;
        std::function<ArrayReturn(const ArrayXs &, const ArrayXs &)> f_array;

        ScalarReturn operator()(const Scalar &x, const Scalar &y) const {
            return f_scalar(x, y);
        }

        ArrayReturn operator()(const ArrayXs &x, const ArrayXs &y) const {
            return f_array(x, y);
        }
    };

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

        virtual std::vector<Eigenfunction2D<Scalar, false>> eigenfunction(const Scalar &E) const = 0;

        virtual std::vector<Eigenfunction2D<Scalar, true>> eigenfunctionWithDerivatives(const Scalar &E) const = 0;

        virtual Scalar estimatePotentialMinimum() const = 0;

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
        std::vector<MatrixXs> M;
        int sectorCount;
        std::vector<typename Matslise2D<Scalar>::Sector *> sectors;
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
        Scalar estimatePotentialMinimum() const override;

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

        std::vector<Eigenfunction2D<Scalar>> eigenfunction(const Scalar &E) const override {
            return eigenfunction < false > (E);
        }

        std::vector<Eigenfunction2D<Scalar, true>> eigenfunctionWithDerivatives(const Scalar &E) const override {
            return eigenfunction < true > (E);
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

        template<bool withDerivatives = false>
        std::vector<Eigenfunction2D<Scalar, withDerivatives>> eigenfunction(const Scalar &E) const {
            return eigenfunction<withDerivatives>(dirichletBoundary, E);
        }

        template<bool withDerivatives = false>
        std::vector<Eigenfunction2D<Scalar, withDerivatives>>
        eigenfunction(const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

    protected:
        std::vector<Y<Scalar, Eigen::Dynamic>> eigenfunctionSteps(
                const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;

        MatrixXs conditionY(Y<Scalar, Eigen::Dynamic> &y) const;

    public:
        class Sector {
        public:
            static const bool expensive = true;

            const Matslise2D<Scalar> *se2d;
            std::shared_ptr<AbstractMatslise<Scalar>> matslise;
            std::shared_ptr<AbstractBasisQuadrature<Scalar>> quadratures;
            typename matslise::Matscs<Scalar>::Sector *matscs;
            Scalar min, max;
            Scalar ybar;

            std::vector<Scalar> eigenvalues;
            std::vector<typename Matslise<Scalar>::Eigenfunction> eigenfunctions;
            Direction direction = none;;

            Sector(const Matslise2D<Scalar> *se2d) : se2d(se2d) {}

            Sector(const Matslise2D<Scalar> *se2d, const Scalar &min, const Scalar &max, Direction);

            virtual ~Sector();

            void setDirection(Direction);

            template<int r>
            Y<Scalar, Eigen::Dynamic, r> propagate(
                    const Scalar &E, const Y<Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b,
                    bool use_h = true) const;

            bool contains(const Scalar &point) const {
                return point <= max && point >= min;
            }


            template<bool withDerivative>
            typename std::conditional<withDerivative, std::pair<ArrayXXs, ArrayXXs>, ArrayXXs>::type
            basis(const ArrayXs &x) const;

            template<bool withDerivative>
            typename std::conditional<withDerivative, std::pair<ArrayXs, ArrayXs>, ArrayXs>::type
            basis(const Scalar &x) const;

            Scalar error() const;

            static bool compare(const Sector &a, const Sector &b) {
                return a.matslise->estimatePotentialMinimum() < b.matslise->estimatePotentialMinimum();
            }

            Sector *
            refine(const Matslise2D<Scalar> *problem, const Scalar &min, const Scalar &max, Direction) const;
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

        Scalar estimatePotentialMinimum() const override {
            return se2d->estimatePotentialMinimum();
        }

        Scalar firstEigenvalue() const override;

        Scalar eigenvalue(const Scalar &Eguess) const override;

        Scalar eigenvalueError(const Scalar &E) const override;

        std::vector<Scalar> eigenvalues(const Scalar &Emin, const Scalar &Emax) const override;

        std::vector<Scalar> eigenvaluesByIndex(int Imin, int Imax) const override;

        std::vector<Eigenfunction2D<Scalar>> eigenfunction(const Scalar &E) const override {
            return eigenfunction < false > (E);
        }

        std::vector<Eigenfunction2D<Scalar, true>>
        eigenfunctionWithDerivatives(const Scalar &E) const override {
            return eigenfunction < true > (E);
        }

        template<bool withDerivatives>
        std::vector<Eigenfunction2D<Scalar, withDerivatives>> eigenfunction(const Scalar &E) const;
    };
}

#endif //MATSLISE_MATSLISE2D_H
