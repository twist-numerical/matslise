#ifndef MATSLISE_BASISQUADRATURE_H
#define MATSLISE_BASISQUADRATURE_H

namespace matslise {
    template<typename Scalar>
    class AbstractBasisQuadrature {
    public:
        virtual Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> dV(
                const typename matslise::Matslise2D<Scalar>::Sector &sector2d, const Scalar &y) = 0;

        virtual ~AbstractBasisQuadrature() {}
    };

    template<typename Scalar, int hmax, bool halfrange = false>
    class BasisQuadrature : public AbstractBasisQuadrature<Scalar> {
    public:
        const matslise::Matslise<Scalar> *matslise;
        std::vector<Eigen::Array<Scalar, hmax, 1>> quadData;

        BasisQuadrature(const matslise::Matslise<Scalar> *matslise) : matslise(matslise) {
        }

        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
        dV(const typename matslise::Matslise2D<Scalar>::Sector &sector2d, const Scalar &y) override;

    public:
        void calculateQuadData(const typename matslise::Matslise2D<Scalar>::Sector &sector2d);
    };
}

#endif //MATSLISE_BASISQUADRATURE_H
