const test = require('ava');
const Module = require('../src/matslise');

const correct = [
    3.1959181, 5.5267439, 7.5578033, 8.0312723, 8.4445814, 9.9280611, 11.3118171,
    12.1032536, 12.2011790, 13.3323313];

(new Module()).then(({Matslise2D}) => {
    test('ixaru', t => {
        const matslise = new Matslise2D(
            (x, y) => (1 + x * x) * (1 + y * y),
            -5, 5, -5, 5,
            {
                tolerance: 1e-5,
                nested: {
                    tolerance: 1e-5,
                }
            });
        const eigenvalue = matslise.firstEigenvalue();
        t.true(Math.abs(correct[0] - eigenvalue) < 1e-5);
/*
        for (let i = 0; i < correct.length; ++i) {
            t.true(Math.abs(correct[i] - eigenvalues[i]) < 1e-5,
                `${correct[i]} != ${eigenvalues[i]}`);
        }*/
    });
});