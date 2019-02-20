const indexOf = (arr, x) => {
  let a = 0,
    b = arr.length;
  while (a + 1 < b) {
    const c = (a + b) >> 1;
    if (x < arr[c]) {
      b = c;
    } else {
      a = c;
    }
  }
  return a;
};

const interpolate = (metric, d) => (xs, ys, outside = "interpolate") => {
  const minX = xs[0];
  const maxX = xs[xs.length - 1];
  return x => {
    if ((x >= minX && x <= maxX) || outside === "interpolate") {
      const i = Math.min(
        xs.length - d - 2,
        Math.max(0, indexOf(xs, x) - d + 1)
      );
      const xnorm = [],
        y = [];
      const norm = xs[i + 2 * d - 1] - xs[i];
      for (let j = 0; j < 2 * d; ++j) {
        xnorm.push((xs[i + j] - xs[i]) / norm);
        y.push(ys[i + j]);
      }
      const m = metric(xnorm, y, (x - xs[i]) / norm);
      return m;
    } else if (outside === "constant") {
      return x <= minX ? ys[0] : ys[ys.length - 1];
    } else {
      return outside;
    }
  };
};

export default {
  nearest: interpolate((_, [y1, y2], x) => (x < 0.5 ? y1 : y2), 1),
  linear: interpolate((_, [y1, y2], x) => (1 - x) * y1 + x * y2, 1)
};
