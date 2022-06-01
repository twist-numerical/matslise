module.exports = function approx(t, exact, measured, relEps = 1e-5) {
    const absError = Math.abs(exact - measured);
    const error = Math.abs(exact) < 1 ? absError : absError / Math.abs(exact);
    t.true(error < relEps, `${exact} ≉ ${measured}`);
};
