//
// Created by toon on 5/16/18.
//

#ifndef SCHRODINGER_ARRAY2D_H
#define SCHRODINGER_ARRAY2D_H


template<class T, int rows, int cols>
struct Array2D {
    T data[rows * cols];
public:
    T *operator[](int i) {
        return data + i * cols;
    }
};


#endif //SCHRODINGER_ARRAY2D_H
