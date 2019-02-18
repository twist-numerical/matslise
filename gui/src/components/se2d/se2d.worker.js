import loadMatslise from "../../lib/loadMatslise";
import DF from "../../lib/Differentiable";

let se2d = null;
const actions = {
  init: ({ settings: { f, x, y } }, updateState) => {
    loadMatslise.then(({ SE2D }) => {
      se2d = new SE2D(
        DF.parse(f, ["x", "y"]).toFunction("x", "y"),
        ...x,
        ...y,
        27,
        27
      );
      updateState({ init: true });
    });
  },
  errors: ({ E: [Emin, Emax] }, updateState) => {
    const Es = [];
    const n = 200;
    for (let i = 0; i < n; ++i) {
      Es.push(Emin + ((Emax - Emin) * (i + Math.random())) / n);
    }
    updateState({
      errors: {
        E: Es,
        errors: Es.map(E => se2d.calculateError(E).first)
      }
    });
  }
};

self.addEventListener("message", message => {
  const time = +new Date();
  actions[message.data.type](message.data, state => {
    state.time = +new Date() - time;
    self.postMessage(state);
  });
});
