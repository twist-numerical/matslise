#ifndef MATSLISE_HEAP_H
#define MATSLISE_HEAP_H

#include <vector>
#include <algorithm>
#include <functional>

namespace matslise {
    template<typename T>
    class Heap {
    private:
        std::vector<T> _data;
        std::function<bool(const T &, const T &)> comp;
    public:
        typedef typename std::vector<T>::size_type size_type;

        explicit Heap(const std::function<bool(const T &, const T &)> &lessThen = std::less<T>()) : comp(lessThen) {
        }

        bool empty() const {
            return _data.empty();
        }

        size_type size() const {
            return _data.size();
        }

        void push(const T &item) {
            _data.push_back(item);
            std::push_heap(_data.begin(), _data.end(), comp);
        }

        void push(T &&item) {
            _data.push_back(std::forward<T>(item));
            std::push_heap(_data.begin(), _data.end(), comp);
        }

        template<class... Args>
        void emplace(Args &&... args) {
            _data.emplace_back(std::forward<Args>(args)...);
            std::push_heap(_data.begin(), _data.end(), comp);
        }

        void pop() {
            std::pop_heap(_data.begin(), _data.end(), comp);
            _data.pop_back();
        }

        const T &front() const {
            return _data.front();
        }

        const std::vector<T> data() const {
            return _data;
        }

        std::vector<T> data() {
            return _data;
        }
    };
}

#endif //MATSLISE_HEAP_H
