
template<typename D, typename F, typename Scalar>
inline D horner(const F &f, Scalar x, int n) {
    D r = f[n - 1];
    for (int i = n - 2; i >= 0; --i)
        r = r * x + f[i];
    return r;
}