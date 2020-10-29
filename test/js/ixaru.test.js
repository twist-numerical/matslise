const test = require('ava');
const Module = require('../src/matslise');

const correct = [
    [0, 3.1959181, 1],
    [1, 5.5267439, 2],
    [3, 7.5578033, 1],
    [4, 8.0312723, 1],
    [5, 8.4445814, 1],
    [6, 9.9280611, 2],
    [8, 11.3118171, 2],
    [10, 12.1032536, 1],
    [11, 12.2011790, 1],
    [12, 13.3323313, 1],
];

(new Module()).then(({Matslise2D}) => {
    test('ixaru', t => {
        const matslise = new Matslise2D(
            (x, y) => (1 + x * x) * (1 + y * y),
            -5.5, 5.5, -5.5, 5.5,
            {
                tolerance: 1e-8,
                xSymmetric: true
            });
        const eigenvalues = matslise.eigenvaluesByIndex(0, 13);
        for (let i = 0; i < eigenvalues.length && i < correct.length; ++i) {
            t.is(correct[i][0], eigenvalues[i].index);
            t.true(Math.abs(correct[i][1] - eigenvalues[i].value) < 1e-3);
            t.is(correct[i][2], eigenvalues[i].multiplicity);
        }
    });
});