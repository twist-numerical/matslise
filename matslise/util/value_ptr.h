#ifndef MATSLISE_VALUE_PTR_H
#define MATSLISE_VALUE_PTR_H

#include <memory>

namespace matslise {
    template<typename T>
    class value_ptr : public std::unique_ptr<T> {
    public:
        using std::unique_ptr<T>::unique_ptr;

        // copy constructor
        value_ptr(value_ptr<T> const &other) : std::unique_ptr<T>() {
            *this = other;
        }

        // Copy assignment
        value_ptr<T> &operator=(value_ptr<T> const &v) {
            this->reset(new T(*v));
            return *this;
        }
    };
}

#endif //MATSLISE_VALUE_PTR_H
