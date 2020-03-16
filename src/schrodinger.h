#ifndef MATSLISE_SCHRODINGER_H
#define MATSLISE_SCHRODINGER_H

namespace matslise {
    template<typename Scalar=double>
    class SE2D;

    template<typename Scalar=double>
    struct Options1 {
        std::shared_ptr<matslise::SectorBuilder<matslise::Matslise<Scalar>>> _builder = Matslise<Scalar>::UNIFORM(26);
        bool _symmetric = false;

        Options1<Scalar> &symmetric(bool symmetric) {
            _symmetric = symmetric;
            return *this;
        }

        Options1<Scalar> &sectorCount(int count) {
            _builder = Matslise<Scalar>::UNIFORM(count);
            return *this;
        }

        Options1<Scalar> &tolerance(Scalar tol) {
            _builder = Matslise<Scalar>::AUTO(tol);
            return *this;
        }
    };

    template<typename Scalar=double>
    struct Options2 {
        Options1<Scalar> nestedOptions;
        std::shared_ptr<matslise::SectorBuilder<matslise::SE2D<Scalar>>> _builder
                = std::shared_ptr<matslise::SectorBuilder<matslise::SE2D<Scalar>>>(
                        new matslise::sectorbuilder::Uniform<SE2D<Scalar>>(17));
        int _stepsPerSector = 1;
        int _N = 12;
        int _gridPoints = 52;

        Options2<Scalar> &nested(Options1<Scalar> &nested) {
            nestedOptions = nested;
            return *this;
        }

        Options2<Scalar> &sectorCount(int count) {
            _builder.reset(new matslise::sectorbuilder::Uniform<SE2D<Scalar>>(count));
            return *this;
        }

        Options2<Scalar> &tolerance(Scalar tol) {
            _builder.reset(new matslise::sectorbuilder::Auto<SE2D<Scalar>>(tol));
            return *this;
        }

        Options2<Scalar> &stepsPerSector(int count) {
            _stepsPerSector = count;
            return *this;
        }

        Options2<Scalar> &N(int value) {
            _N = value;
            return *this;
        }

        Options2<Scalar> &gridPoints(int value) {
            _gridPoints = value;
            return *this;
        }
    };


    template<int n, typename Scalar>
    class Rectangle {
    public:
        Rectangle<n - 1, Scalar> sub;
        Scalar min, max;

        Scalar getMin(int axis) const {
            if (axis + 1 == n)
                return min;
            return sub.getMin(axis);
        }

        Scalar getMax(int axis) const {
            if (axis + 1 == n)
                return max;
            return sub.getMax(axis);
        }

        Scalar diameter() const {
            return hypot(max - min, sub.diameter());
        }
    };

    template<typename Scalar>
    class Rectangle<1, Scalar> {
    public:
        Scalar min, max;

        Scalar getMin(int axis) const {
            assert(axis == 0);
            (void) (axis); // UNUSED
            return min;
        }

        Scalar getMax(int axis) const {
            assert(axis == 0);
            (void) (axis); // UNUSED
            return max;
        }

        Scalar diameter() const {
            return abs(max - min);
        }
    };
}

#endif //MATSLISE_SCHRODINGER_H
