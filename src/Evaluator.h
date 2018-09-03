

#ifndef SCHRODINGER_EVALUATOR_H
#define SCHRODINGER_EVALUATOR_H

template<class R, class... A>
class Evaluator {
public:

    virtual R eval(A...) const = 0;

    R operator()(A... args) const {
        return eval(args...);
    }

    virtual ~Evaluator() {};
};

#endif