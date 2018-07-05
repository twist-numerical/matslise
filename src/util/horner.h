
template<typename D>
inline D horner(const D *f, double x, int n) {
    D r = f[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        r *= x;
        r += f[i];
    }
    return r;
}