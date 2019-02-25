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
      const se2d = new SE2D(V, ...x, ...y, 27, 27);
      updateState("init", {
        calculating: true,
        init: true
      });
      actions.errors({ E }, updateState);
      onInit.done(se2d, settings);
    });
  },
  errors: ({ E: [Emin, Emax] }, updateState) => {
    onInit.listen((se2d, settings) => {
      const Es = [];
      const n = 50;
      for (let i = 0; i < n; ++i) {
        Es.push(Emin + ((Emax - Emin) * (i + Math.random())) / n);
      }
      updateState("errors", {
        calculating: false,
        errors: {
          E: Es,
          errors: Es.map(E => se2d.calculateError(E).first)
        }
      });
    });
  },
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
        eigenfunction: { x, y, zs }
      });
    });
  }
};

self.addEventListener("message", message => {
  const time = +new Date();
  if (message.data.type in actions)
    actions[message.data.type](message.data.data, (type, state) => {
      self.postMessage({ type, data: state });
      self.postMessage({ type: "time", data: +new Date() - time });
    });
  else console.error(`Unknown action '${message.data.type}'`);
});
