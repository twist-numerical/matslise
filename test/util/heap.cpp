#include "../catch.hpp"
#include "../../matslise/util/heap.h"

using namespace matslise;

struct IntWrapper {
    int data;

    explicit IntWrapper(int value) : data(value) {
    }

    bool operator<(const IntWrapper &rhs) const {
        return data < rhs.data;
    }

    bool operator>(const IntWrapper &rhs) const {
        return data > rhs.data;
    }
};

TEST_CASE("Max heap", "[util][rectangle]") {
    Heap<IntWrapper> heap;

    REQUIRE(heap.empty());
    REQUIRE(heap.size() == 0);

    IntWrapper m(4);
    heap.push(std::move(m));

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 1);
    REQUIRE(heap.front().data == 4);

    heap.push(IntWrapper(5));

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 2);
    REQUIRE(heap.front().data == 5);

    heap.emplace(3);

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 3);
    REQUIRE(heap.front().data == 5);

    heap.pop();

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 2);
    REQUIRE(heap.front().data == 4);

    heap.pop();

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 1);
    REQUIRE(heap.front().data == 3);

    heap.pop();

    REQUIRE(heap.empty());
    REQUIRE(heap.size() == 0);
}

TEST_CASE("Min heap", "[util][rectangle]") {
    Heap<IntWrapper> heap((std::greater<IntWrapper>()));

    REQUIRE(heap.empty());
    REQUIRE(heap.size() == 0);

    IntWrapper m(4);
    heap.push(std::move(m));

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 1);
    REQUIRE(heap.front().data == 4);

    heap.push(IntWrapper(5));

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 2);
    REQUIRE(heap.front().data == 4);

    heap.emplace(3);

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 3);
    REQUIRE(heap.front().data == 3);

    heap.pop();

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 2);
    REQUIRE(heap.front().data == 4);

    heap.pop();

    REQUIRE(!heap.empty());
    REQUIRE(heap.size() == 1);
    REQUIRE(heap.front().data == 5);

    heap.pop();

    REQUIRE(heap.empty());
    REQUIRE(heap.size() == 0);
}
