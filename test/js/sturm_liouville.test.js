const test = require('ava');
const Module = require('../matslise/matslise');
const approx = require('./approx');

(new Module()).then(({SturmLiouville}) => {

    test('klotter', t => {
        const slp = new SturmLiouville(
            (x) => 1,
            (x) => 3 / (4 * x * x),
            (x) => {
                const r = 8 * Math.PI / (3 * x * x * x);
                return r * r;
            },
            8 / 7, 8, 1e-6);

        const n = 20;
        const eigenvalues = slp.eigenvaluesByIndex(0, n, [0, 1], [0, 1]);
        t.is(n, eigenvalues.length);
        for (let i = 0; i < n; ++i) {
            t.is(i, eigenvalues[i].index);
            approx(t, (i + 1) * (i + 1), eigenvalues[i].eigenvalue);
        }

        slp.delete();
    });

});
