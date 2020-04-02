#ifndef SCHRODINGER_MATSCS_H
#define SCHRODINGER_MATSCS_H

#include <ostream>
#include <array>
#include <vector>
#include <functional>

#define MATSCS_HMAX_delta 7
#define MATSCS_ETA_delta 3
#define MATSCS_ETA_h 6
#define MATSCS_N 9

namespace matslise {
    template<typename _Scalar=double>
    class Matscs {
    public:
        typedef _Scalar Scalar;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
        typedef Eigen::Matrix<std::complex<Scalar>, Eigen::Dynamic, Eigen::Dynamic> MatrixXcs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ArrayXs;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> ArrayXXs;
        static const int order = 12;

        class Sector;

        std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(Scalar)> V;
        int n;
        Scalar xmin, xmax;
        int sectorCount;
        std::vector<Matscs::Sector *> sectors;
        Scalar match;
        int matchIndex;
    public:
        static std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>> UNIFORM(int sectorCount) {
            return std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>>(
                    new matslise::sectorbuilder::Uniform<Matscs<Scalar>>(sectorCount));
        }

        static std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>> AUTO(const Scalar &tolerance) {
            return std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>>(
                    new matslise::sectorbuilder::Auto<Matscs<Scalar>>(tolerance));
        }

        Matscs(std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(Scalar)> V,
               int n, const Scalar &xmin, const Scalar &xmax,
               std::shared_ptr<matslise::SectorBuilder<Matscs<Scalar>>> sectorBuilder) :
                V(V), n(n), xmin(xmin), xmax(xmax) {
            sectorBuilder->build(this, xmin, xmax);
            sectorCount = sectors.size();
        }

        Matscs(std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(Scalar)> V,
               int n, const Scalar &xmin, const Scalar &xmax, int sectorCount)
                : Matscs(V, n, xmin, xmax, UNIFORM(sectorCount)) {};

        bool contains(const Scalar &point) const {
            return point <= xmax && point >= xmin;
        }

        template<int r>
        matslise::Y<Scalar, Eigen::Dynamic, r>
        propagateColumn(const Scalar &E, const matslise::Y<Scalar, Eigen::Dynamic, r> &y, const Scalar &a,
                        const Scalar &b,
                        bool use_h = true) const;

        std::pair<matslise::Y<Scalar, Eigen::Dynamic>, Scalar>
        propagate(const Scalar &E, const matslise::Y<Scalar, Eigen::Dynamic> &y, const Scalar &a, const Scalar &b,
                  bool use_h = true) const;

        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
        propagatePsi(const Scalar &E, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &psi,
                     const Scalar &a, const Scalar &b) const;

        std::vector<matslise::Y<Scalar, Eigen::Dynamic>> *
        eigenfunction(const Scalar &E, std::vector<Scalar> &x) const;

        std::function<Y<Scalar, Eigen::Dynamic, 1>(Scalar)> eigenfunctionCalculator(
                const Scalar &E, const matslise::Y<Scalar, Eigen::Dynamic, 1> &left,
                const matslise::Y<Scalar, Eigen::Dynamic, 1> &right);

        ~Matscs();


    public:

        class Sector {
        private:
            const Matscs *s;
        public:
            Eigen::Array<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, MATSCS_ETA_delta, MATSCS_HMAX_delta>
                    t_coeff;
            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> t_coeff_h[MATSCS_ETA_h];
            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> diagonalize;
            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *vs;
            Scalar min, max, h;
            bool backward;

            Sector(const Matscs *problem, const Scalar &min, const Scalar &max, bool backward = false);

            void calculateTCoeffs();

            bool contains(Scalar point) const {
                return point <= max && point >= min;
            }

            T <Scalar, Eigen::Dynamic> calculateT(const Scalar &E, bool use_h = true) const;

            T <Scalar, Eigen::Dynamic> calculateT(const Scalar &E, const Scalar &delta, bool use_h = true) const;

            template<int r = Eigen::Dynamic>
            Y <Scalar, Eigen::Dynamic, r>
            propagateColumn(const Scalar &E, const Y <Scalar, Eigen::Dynamic, r> &y0, const Scalar &a, const Scalar &b,
                            bool use_h = true) const;

            std::pair<Y < Scalar, Eigen::Dynamic>, Scalar>

            propagate(const Scalar &E, const Y <Scalar, Eigen::Dynamic> &y0, const Scalar &a, const Scalar &b,
                      bool use_h = true) const;

            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
            propagatePsi(const Scalar &E, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &psi,
                         const Scalar &a, const Scalar &b) const;

            Scalar error() const;

            ~Sector();

        private :
            MatrixXcs theta(const Y <Scalar, Eigen::Dynamic> &) const;

            template<int r>
            Y <Scalar, Eigen::Dynamic, r> propagateDeltaColumn(
                    const Scalar &E, const Y <Scalar, Eigen::Dynamic, r> &y0, const Scalar &_delta, bool use_h) const;

            std::pair<Y < Scalar, Eigen::Dynamic>, Scalar>

            propagateDelta(
                    const Scalar &E, const Y <Scalar, Eigen::Dynamic> &y0, const Scalar &_delta, bool use_h) const;

        };
    };
}

#endif