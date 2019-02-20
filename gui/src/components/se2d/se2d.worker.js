import loadMatslise from "../../lib/loadMatslise";
import DF from "../../lib/Differentiable";

let se2d = null;
let settings = {};

const actions = {
  init: ({ f: raw_f, x, y, E }, updateState) => {
    loadMatslise.then(({ SE2D }) => {
      const f = DF.parse(raw_f, ["x", "y"]).toFunction("x", "y");
      settings = { f, x, y };
      se2d = new SE2D(f, ...x, ...y, 27, 27);
      updateState("init", {
        calculating: true,
        init: true
      });
      actions.errors({ E }, updateState);
    });
  },
  errors: ({ E: [Emin, Emax] }, updateState) => {
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
  },
  locate: ({ E: Eguess, eigenvalues }, updateState) => {
    let err;
    let E = Eguess;
    for (let i = 0; i < 20; ++i) {
      err = se2d.calculateError(E);
      E -= err.first / err.second;
    }
    if (
      err.first < 1e-10 &&
      !eigenvalues.filter(old => Math.abs(old - E) < 1e-8).length
    )
      eigenvalues.push(E);
    eigenvalues.sort((a, b) => a - b);
    updateState("locate", {
      calculating: false,
      eigenvalues
    });
  },
  eigenfunction: ({ E, xn, yn }, updateState) => {
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
