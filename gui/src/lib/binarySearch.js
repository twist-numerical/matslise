export default (f, a, b) => {
  let fa = f(a),
    fb = f(b);
  if (fa === 0) return a;
  if (fb === 0) return b;
  // eslint-disable-next-line
  if (fa < 0 === fb < 0)
    throw new Error("f(a) and f(b) must have a different sign");
  let c, fc;
  while (Math.abs(a - b) > 1e-12) {
    c = (a + b) / 2;
    fc = f(c);
    if (fc === 0) return c;
    // eslint-disable-next-line
    if (fa < 0 === fc < 0) {
      fa = fc;
      a = c;
    } else {
      fb = fc;
      b = c;
    }
  }
  return c;
};
