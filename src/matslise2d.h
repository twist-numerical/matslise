#ifndef MATSLISE_MATSLISE2D_H
#define MATSLISE_MATSLISE2D_H

#include "util/rectangle.h"
#include "matslise2d/basisQuadrature.h"
#include "matsliseNd.h"

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
        const Rectangle<Scalar, 2> domain;

        AbstractMatslise2D(const std::function<Scalar(Scalar, Scalar)> &potential, const Rectangle<Scalar, 2> &domain)
                : potential(potential), domain(domain) {
        }

        virtual std::pair<Scalar, Eigen::Index> eigenvalue(const Scalar &Eguess) const = 0;

        virtual Scalar eigenvalueError(const Scalar &E) const = 0;

        virtual std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvalues(const Scalar &Emin, const Scalar &Emax) const = 0;

        virtual std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvaluesByIndex(Eigen::Index Imin, Eigen::Index Imax) const = 0;

        virtual std::vector<Eigenfunction2D<Scalar, false>> eigenfunction(const Scalar &E) const = 0;

        virtual std::vector<Eigenfunction2D<Scalar, true>> eigenfunctionWithDerivatives(const Scalar &E) const = 0;

        virtual Scalar estimatePotentialMinimum() const = 0;

        virtual Eigen::Index estimateIndex(const Scalar &E) const = 0;

        virtual ~AbstractMatslise2D() = default;
    };

    template<typename Scalar>
    class Matslise2DSector : public MatsliseNDSector<Scalar> {
    public:
        static const bool expensive = true;
        typedef typename AbstractMatslise2D<Scalar>::VectorXs VectorXs;
        typedef typename AbstractMatslise2D<Scalar>::MatrixXs MatrixXs;
        typedef typename AbstractMatslise2D<Scalar>::ArrayXs ArrayXs;
        typedef typename AbstractMatslise2D<Scalar>::ArrayXXs ArrayXXs;

        const Matslise2D<Scalar> *se2d;
        std::shared_ptr<AbstractMatslise<Scalar>> matslise;
        std::shared_ptr<AbstractBasisQuadrature<Scalar>> quadratures;
        using MatsliseNDSector<Scalar>::matscs;
        using MatsliseNDSector<Scalar>::min;
        using MatsliseNDSector<Scalar>::max;
        Scalar ybar;

        std::vector<Scalar> eigenvalues;
        std::vector<typename Matslise<Scalar>::Eigenfunction> eigenfunctions;
        Direction direction = none;

        Matslise2DSector(const Matslise2D<Scalar> *_se2d) : se2d(_se2d) {};

        Matslise2DSector(const Matslise2D<Scalar> *se2d, const Scalar &min, const Scalar &max, Direction);

        void setDirection(Direction);

        template<bool withDerivative>
        typename std::conditional<withDerivative, std::pair<ArrayXXs, ArrayXXs>, ArrayXXs>::type
        basis(const ArrayXs &x) const;

        template<bool withDerivative>
        typename std::conditional<withDerivative, std::pair<ArrayXs, ArrayXs>, ArrayXs>::type
        basis(const Scalar &x) const;

        static bool compare(const Matslise2DSector &a, const Matslise2DSector &b) {
            return a.matslise->estimatePotentialMinimum() < b.matslise->estimatePotentialMinimum();
        }

        Matslise2DSector *
        refine(const Matslise2D<Scalar> *problem, const Scalar &min, const Scalar &max, Direction) const;
    };

    template<typename _Scalar=double>
    class Matslise2D : public AbstractMatslise2D<_Scalar>, private MatsliseND<_Scalar, Matslise2DSector<_Scalar>> {
    public:
        typedef _Scalar Scalar;
        using typename AbstractMatslise2D<Scalar>::VectorXs;
        using typename AbstractMatslise2D<Scalar>::MatrixXs;
        using typename AbstractMatslise2D<Scalar>::ArrayXs;
        using typename AbstractMatslise2D<Scalar>::ArrayXXs;
        static constexpr int order = matslise::Matscs<Scalar>::order;
        static constexpr bool refineSectors = true;

        struct Config {
            Scalar tolerance = 1e-6;
            bool xSymmetric = false;
            Eigen::Index basisSize = 12;
            Eigen::Index stepsPerSector = 2;
            std::optional<SectorBuilder<Matslise<Scalar>, Scalar>> xSectorBuilder;
            std::optional<SectorBuilder<Matslise2D<Scalar>, Scalar>> ySectorBuilder;
        };

        using Sector = Matslise2DSector<Scalar>;
        using AbstractMatslise2D<Scalar>::domain;
        using AbstractMatslise2D<Scalar>::potential;
        using MatsliseND<Scalar, Matslise2DSector<Scalar>>::M;
        using MatsliseND<Scalar, Matslise2DSector<Scalar>>::sectors;
        using MatsliseND<Scalar, Matslise2DSector<Scalar>>::matchIndex;
        using MatsliseND<Scalar, Matslise2DSector<Scalar>>::dirichletBoundary;
        Config config;
    public:
        Matslise2D(const std::function<Scalar(Scalar, Scalar)> &potential,
                   const matslise::Rectangle<Scalar, 2> &domain, const Config &config = Config());

        virtual ~Matslise2D();

        using MatsliseND<Scalar, Sector>::propagate;
        using MatsliseND<Scalar, Sector>::matchingErrorMatrix;
        using MatsliseND<Scalar, Sector>::matchingErrors;
        using MatsliseND<Scalar, Sector>::matchingError;
        using MatsliseND<Scalar, Sector>::eigenvalue;
        using MatsliseND<Scalar, Sector>::eigenvalueError;
        using MatsliseND<Scalar, Sector>::eigenvalues;
        using MatsliseND<Scalar, Sector>::eigenvaluesByIndex;
        using MatsliseND<Scalar, Sector>::estimateIndex;

    public: // overrides
        Scalar estimatePotentialMinimum() const override;

        std::pair<Scalar, Eigen::Index> eigenvalue(const Scalar &Eguess) const override {
            return MatsliseND<Scalar, Sector>::eigenvalue(Eguess);
        }

        Scalar eigenvalueError(const Scalar &E) const override {
            return MatsliseND<Scalar, Sector>::eigenvalueError(E);
        }

        std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvalues(const Scalar &Emin, const Scalar &Emax) const override {
            return MatsliseND<Scalar, Sector>::eigenvalues(Emin, Emax);
        }

        std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvaluesByIndex(Eigen::Index Imin, Eigen::Index Imax) const override {
            return MatsliseND<Scalar, Sector>::eigenvaluesByIndex(Imin, Imax);
        }

        Eigen::Index estimateIndex(const Scalar &E) const override {
            return MatsliseND<Scalar, Sector>::estimateIndex(E);
        }

        std::vector<Eigenfunction2D<Scalar>> eigenfunction(const Scalar &E) const override {
            return eigenfunction < false > (E);
        }

        std::vector<Eigenfunction2D<Scalar, true>> eigenfunctionWithDerivatives(const Scalar &E) const override {
            return eigenfunction < true > (E);
        }

    public: // left boundary conditions
        template<bool withDerivatives = false>
        std::vector<Eigenfunction2D<Scalar, withDerivatives>> eigenfunction(const Scalar &E) const {
            return eigenfunction<withDerivatives>(dirichletBoundary, E);
        }

        template<bool withDerivatives = false>
        std::vector<Eigenfunction2D<Scalar, withDerivatives>>
        eigenfunction(const Y<Scalar, Eigen::Dynamic> &left, const Scalar &E) const;
    };

    template<typename _Scalar=double>
    class Matslise2DHalf : public AbstractMatslise2D<_Scalar> {
    public:
        typedef _Scalar Scalar;
        using typename AbstractMatslise2D<_Scalar>::VectorXs;
        using typename AbstractMatslise2D<_Scalar>::MatrixXs;
        using typename AbstractMatslise2D<_Scalar>::ArrayXs;
        using typename AbstractMatslise2D<_Scalar>::ArrayXXs;
        typedef typename Matslise2D<_Scalar>::Config Config;
    public:
        Y<Scalar, Eigen::Dynamic, Eigen::Dynamic> neumannBoundary;
        Y<Scalar, Eigen::Dynamic, Eigen::Dynamic> dirichletBoundary;

        Matslise2D<Scalar> *se2d;

        using AbstractMatslise2D<Scalar>::domain;
        using AbstractMatslise2D<Scalar>::potential;

    public:
        Matslise2DHalf(const std::function<Scalar(const Scalar &, const Scalar &)> &potential,
                       const Rectangle<Scalar, 2> &domain, const Config &config);

        virtual ~Matslise2DHalf();

        Scalar estimatePotentialMinimum() const override {
            return se2d->estimatePotentialMinimum();
        }

        std::pair<Scalar, Eigen::Index> eigenvalue(const Scalar &Eguess) const override;

        Scalar eigenvalueError(const Scalar &E) const override;

        std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvalues(const Scalar &Emin, const Scalar &Emax) const override;

        std::vector<std::tuple<Eigen::Index, Scalar, Eigen::Index>>
        eigenvaluesByIndex(Eigen::Index Imin, Eigen::Index Imax) const override;

        std::vector<Eigenfunction2D<Scalar>> eigenfunction(const Scalar &E) const override {
            return eigenfunction < false > (E);
        }

        std::vector<Eigenfunction2D<Scalar, true>>
        eigenfunctionWithDerivatives(const Scalar &E) const override {
            return eigenfunction < true > (E);
        }

        Eigen::Index estimateIndex(const Scalar &) const override;

        template<bool withDerivatives>
        std::vector<Eigenfunction2D<Scalar, withDerivatives>> eigenfunction(const Scalar &E) const;
    };
}

#endif //MATSLISE_MATSLISE2D_H
