#include "./quadrature.h"

using namespace Eigen;
using namespace quadrature;

namespace quadrature::gauss_konrod {
    template<typename Scalar>
    const Array<Scalar, 31, 1> abscissae;
    template<typename Scalar>
    const Array<Scalar, 15, 1> weightsGauss;
    template<typename Scalar>
    const Eigen::Array<Scalar, 31, 1> weightsKonrod;

    template<>
    const Array<double, 31, 1> abscissae<double> = (Array<double, 31, 1>()
            << 0.998002298693397060285172840152271, 0.987992518020485428489565718586613, 0.967739075679139134257347978784337, 0.937273392400705904307758947710209, 0.897264532344081900882509656454496, 0.848206583410427216200648320774217, 0.790418501442465932967649294817947, 0.724417731360170047416186054613938, 0.650996741297416970533735895313275, 0.570972172608538847537226737253911, 0.485081863640239680693655740232351, 0.394151347077563369897207370981045, 0.299180007153168812166780024266389, 0.201194093997434522300628303394596, 0.101142066918717499027074231447392, 0.000000000000000000000000000000000, -0.101142066918717499027074231447392, -0.201194093997434522300628303394596, -0.299180007153168812166780024266389, -0.394151347077563369897207370981045, -0.485081863640239680693655740232351, -0.570972172608538847537226737253911, -0.650996741297416970533735895313275, -0.724417731360170047416186054613938, -0.790418501442465932967649294817947, -0.848206583410427216200648320774217, -0.897264532344081900882509656454496, -0.937273392400705904307758947710209, -0.967739075679139134257347978784337, -0.987992518020485428489565718586613, -0.998002298693397060285172840152271
    ).finished();

    template<>
    const Array<double, 31, 1> weightsKonrod<double> = (Array<double, 31, 1>()
            << 0.005377479872923348987792051430128, 0.015007947329316122538374763075807, 0.025460847326715320186874001019653, 0.035346360791375846222037948478360, 0.044589751324764876608227299373280, 0.053481524690928087265343147239430, 0.062009567800670640285139230960803, 0.069854121318728258709520077099147, 0.076849680757720378894432777482659, 0.083080502823133021038289247286104, 0.088564443056211770647275443693774, 0.093126598170825321225486872747346, 0.096642726983623678505179907627589, 0.099173598721791959332393173484603, 0.100769845523875595044946662617570, 0.101330007014791549017374792767493, 0.100769845523875595044946662617570, 0.099173598721791959332393173484603, 0.096642726983623678505179907627589, 0.093126598170825321225486872747346, 0.088564443056211770647275443693774, 0.083080502823133021038289247286104, 0.076849680757720378894432777482659, 0.069854121318728258709520077099147, 0.062009567800670640285139230960803, 0.053481524690928087265343147239430, 0.044589751324764876608227299373280, 0.035346360791375846222037948478360, 0.025460847326715320186874001019653, 0.015007947329316122538374763075807, 0.005377479872923348987792051430128
    ).finished();

    template<>
    const Array<double, 15, 1> weightsGauss<double> = (Array<double, 15, 1>()
            << 0.030753241996117268354628393577204, 0.070366047488108124709267416450667, 0.107159220467171935011869546685869, 0.139570677926154314447804794511028, 0.166269205816993933553200860481209, 0.186161000015562211026800561866423, 0.198431485327111576456118326443839, 0.202578241925561272880620199967519, 0.198431485327111576456118326443839, 0.186161000015562211026800561866423, 0.166269205816993933553200860481209, 0.139570677926154314447804794511028, 0.107159220467171935011869546685869, 0.070366047488108124709267416450667, 0.030753241996117268354628393577204
    ).finished();

    template<typename Scalar, typename Value=Scalar>
    inline std::pair<Value, Scalar> applyGaussKonrod(
            const std::function<Eigen::Array<Value, Eigen::Dynamic, 1>(
                    const Eigen::Array<Scalar, Eigen::Dynamic, 1> &)> &f, Scalar a, Scalar b,
            const std::function<Scalar(const Value &)> &error) {
        Scalar h = (b - a) / 2;
        Eigen::Array<Scalar, Eigen::Dynamic, 1> x = a + h + h * abscissae<Scalar>;
        Eigen::Array<Value, Eigen::Dynamic, 1> y = f(x);
        Value valueKonrad;
        Value valueGauss;
        if constexpr(std::is_same<Scalar, Value>::value) {
            valueKonrad = (weightsKonrod<Scalar> * y).sum();
            valueGauss = (weightsGauss<Scalar> *
                          Eigen::Map<Eigen::Array<Value, 15, 1>, 0, Eigen::InnerStride<2>>(y.data() + 1)).sum();
        } else {
            valueKonrad = weightsKonrod<Scalar>(0) * y(0);
            for (Index i = 1; i < 31; ++i)
                valueKonrad += weightsKonrod<Scalar>(i) * y(i);

            valueGauss = weightsGauss<Scalar>(1) * y(1);
            for (Index i = 3; i < 31; i += 2)
                valueKonrad += weightsGauss<Scalar>(i) * y(i);
        }
        valueGauss *= h;
        valueKonrad *= h;
        return {valueKonrad, error((valueKonrad - valueGauss))};
    }

    template<typename Scalar, typename Value=Scalar>
    inline std::pair<Value, Scalar> applyGaussKonrod(
            const std::function<Value(const Scalar &)> &f, Scalar a, Scalar b,
            const std::function<Scalar(const Value &)> &error) {
        return applyGaussKonrod((std::function<Eigen::Array<Value, Eigen::Dynamic, 1>(
                const Eigen::Array<Scalar, Eigen::Dynamic, 1> &)>) [&](
                const Eigen::Array<Scalar, Eigen::Dynamic, 1> &x) -> Eigen::Array<Value, Eigen::Dynamic, 1> {
            return x.unaryExpr(f);
        }, a, b, error);
    }

    template<typename Scalar, typename Value, bool bulk>
    Value adaptive(
            const std::function<typename std::conditional<bulk, Eigen::Array<Value, Eigen::Dynamic, 1>, Value>::type(
                    const typename std::conditional<bulk, Eigen::Array<Scalar, Eigen::Dynamic, 1>, Scalar>::type &)> &f,
            Scalar a, Scalar b, const Scalar &tolerance,
            const std::function<Scalar(const Value &)> &error) {
        typedef std::tuple<Scalar, Scalar, std::pair<Value, Scalar >> Item;
        static std::function<bool(const Item &, const Item &)> comp = [](const Item &a, const Item &b) {
            return std::get<2>(a).second < std::get<2>(b).second;
        };
        std::vector<Item> heap;
        std::pair<Value, Scalar> onInterval = applyGaussKonrod<Scalar, Value>(f, a, b, error);
        heap.push_back({a, b, onInterval});
        Scalar totalError = onInterval.second;
        while (totalError > tolerance) {
            std::tie(a, b, onInterval) = heap.front();
            std::pop_heap(heap.begin(), heap.end(), comp);
            heap.pop_back();
            totalError -= onInterval.second;
            Scalar c = (a + b) / 2;

            onInterval = applyGaussKonrod(f, a, c, error);
            totalError += onInterval.second;
            heap.push_back({a, c, onInterval});
            std::push_heap(heap.begin(), heap.end(), comp);

            onInterval = applyGaussKonrod(f, c, b, error);
            totalError += onInterval.second;
            heap.push_back({c, b, onInterval});
            std::push_heap(heap.begin(), heap.end(), comp);
        }

        Value total = std::get<2>(heap.front()).first;
        for (auto i = heap.begin() + 1; i != heap.end(); ++i)
            total += std::get<2>(*i).first;
        return total;
    }
}


#define INSTANTIATE_GK_BULK(Scalar, Value, bulk) \
template Value quadrature::gauss_konrod::adaptive<Scalar, Value, bulk>( \
    const std::function<std::conditional<bulk, Eigen::Array<Value, Eigen::Dynamic, 1>, Value>::type( \
        const std::conditional<bulk, Eigen::Array<Scalar, Eigen::Dynamic, 1>, Scalar>::type &)> &, \
    Scalar, Scalar, const Scalar &, const std::function<Scalar(const Value &)> &);\

#define INSTANTIATE_GK(Scalar, Value) \
INSTANTIATE_GK_BULK(Scalar, Value, true) \
INSTANTIATE_GK_BULK(Scalar, Value, false)

#define INSTANTIATE_QUADRATURES(Scalar, Value) \
INSTANTIATE_GK(Scalar, Scalar) \
INSTANTIATE_GK(Scalar, Value)

INSTANTIATE_QUADRATURES(double, MatrixXd)

#ifdef MATSLISE_long_double

typedef Eigen::Matrix<long double, Dynamic, Dynamic> MatrixXld;
INSTANTIATE_QUADRATURES(long double, MatrixXld)

#endif

#ifdef MATSLISE_float128

#include <boost/multiprecision/float128.hpp>


typedef Eigen::Matrix<boost::multiprecision::float128, Dynamic, Dynamic> MatrixXq;
INSTANTIATE_QUADRATURES(boost::multiprecision::float128)

#endif