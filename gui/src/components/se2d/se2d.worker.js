import loadMatslise from "../../lib/loadMatslise";
import DF from "../../lib/Differentiable";

const onInit = new class {
  se2d = null;
  settings = {};
  callbacks = [];
  listen(f) {
    if (this.se2d === null) {
      this.callbacks.push(f);
    } else {
      f(this.se2d, this.settings);
    }
  }

  done(se2d, settings) {
    this.se2d = se2d;
    this.settings = settings;
    while (this.callbacks.length > 0) {
      this.callbacks.pop()(se2d, settings);
    }
  }
}();

const actions = {
  init: ({ V: rawV, x, y, E }, updateState) => {
    loadMatslise.then(({ SE2D }) => {
      const V = DF.parse(rawV, ["x", "y"]).toFunction("x", "y");
      const settings = { V, x, y };
      const se2d = new SE2D(V, ...x, ...y, {
        tolerance: 1e-5,
        stepsPerSector: 3,
        nested: {
          tolerance: 1e-8
        }
      });
      updateState("init", {
        calculating: true,
        init: true
      });
      actions.errors({ E }, updateState);
      onInit.done(se2d, settings);
    });
  },
  errors: (() => {
    let callCount = 0;
    return ({ E: [Emin, Emax] }, updateState) => {
      const callIndex = ++callCount;
      onInit.listen((se2d, settings) => {
        const Es = [];
        const n = 50;
        for (let i = 0; i < n; ++i) {
          Es.push(Emin + ((Emax - Emin) * (i + Math.random())) / n);
        }
        let errors,
          active,
          max = 0;
        const calculate = i => {
          if (callIndex !== callCount) return;
          if (i === n) {
            updateState("errors", {
              calculating: false,
              errors: {
                functions: errors,
                max: max
              }
            });
          } else {
            const values = se2d
              .calculateErrors(Es[i])
              .map(v => {
                return {
                  ...v,
                  order: Math.hypot(v.first, v.first / v.second)
                };
              })
              .sort((a, b) => a.order - b.order);
            //values.splice(4);
            max = Math.max(
              max,
              values.map(v => Math.abs(v.first)).sort((a, b) => a - b)[1]
            );
            const newPoint = v => {
              return {
                E: [Es[i]],
                errors: [v.first]
              };
            };
            const addPoint = ({ E, errors }, v) => {
              E.push(Es[i]);
              errors.push(v);
            };
            if (i === 0) {
              errors = values.map(newPoint);
              active = errors.map(v => v);
            } else {
              const h = Es[i] - Es[i - 1];
              values.forEach(v => {
                v.closest = active
                  .map(({ errors, E }, i) => {
                    const current = errors[errors.length - 1];
                    const prev =
                      errors.length >= 2 ? errors[errors.length - 2] : current;
                    const hp =
                      errors.length >= 2
                        ? E[E.length - 1] - E[E.length - 2]
                        : h;
                    const diff = (current - prev) / hp;
                    const d = Math.hypot(
                      current + h * diff - v.first,
                      0.1 * Math.atan(v.second - diff)
                    );
                    return { i, d };
                  })
                  .sort((a, b) => a.d - b.d);
                v.distance = v.closest[0].d;
                v.closest = v.closest.map(({ i }) => i);
              });
              values.sort((a, b) => a.distance - b.distance);
              values.forEach(v => {
                for (let j = 0; j < v.closest.length; ++j) {
                  const Ej = active[v.closest[j]].E;
                  if (Ej[Ej.length - 1] < Es[i]) {
                    if (j < 0.3 * v.closest.length) {
                      addPoint(active[v.closest[j]], v.first);
                    } else {
                      errors.push((active[v.closest[j]] = newPoint(v)));
                    }
                    break;
                  }
                }
              });
            }
            setTimeout(() => calculate(i + 1));
          }
        };
        calculate(0);
      });
    };
  })(),
  locate: ({ E: Eguess, eigenvalues }, updateState) => {
    onInit.listen((se2d, settings) => {
      const E = se2d.findEigenvalue(Eguess);
      if (
        !isNaN(E) &&
        !eigenvalues.filter(old => Math.abs(old - E) < 1e-5).length
      )
        eigenvalues.push(E);
      eigenvalues.sort((a, b) => a - b);
      updateState("locate", {
        calculating: false,
        eigenvalues
      });
    });
  },
  eigenfunction: ({ E, xn, yn }, updateState) => {
    onInit.listen((se2d, settings) => {
      const {
        x: [xmin, xmax],
        y: [ymin, ymax]
      } = settings;

      const x = [];
      for (let i = 0; i <= xn; ++i) x.push(xmin + ((xmax - xmin) * i) / xn);
      const y = [];
      for (let i = 0; i <= yn; ++i) y.push(ymin + ((ymax - ymin) * i) / yn);

      const zs = se2d.computeEigenfunction(E, x, y);
      updateState("eigenfunction", {
        eigenfunction: { x, y, zs, sectorPoints: se2d.sectorPoints() }
      });
    });
  }
};

self.addEventListener("message", message => {
  const time = +new Date();
  if (message.data.type in actions)
    actions[message.data.type](message.data.data, (type, state) => {
      self.postMessage({ type, data: state });
      self.postMessage({
        type: "time",
        data: { time: +new Date() - time, action: type }
      });
    });
  else console.error(`Unknown action '${message.data.type}'`);
});
