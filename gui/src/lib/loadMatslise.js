/* global Module */
if (self.importScripts) {
  self.importScripts("matslise.js");
}

class MatsliseLoader {
  Matslise = null;
  callbacks = [
    Matslise => {
      if (typeof window !== "undefined") window.Matslise = Matslise;
    }
  ];

  constructor() {
    if (typeof window !== "undefined") {
      if (!window.Module) window.Module = {};
    }

    if (Module.Matslise !== undefined) this.found();
    else if (!Module.onRuntimeInitialized)
      Module.onRuntimeInitialized = () => this.found();
  }

  found() {
    this.Matslise = {
      Matslise: Module.Matslise,
      SE2D: Module.SE2D,
      HalfRange: Module.HalfRange
    };

    while (this.callbacks.length > 0) this.callbacks.pop()(this.Matslise);
  }

  then(callback) {
    if (this.Matslise != null) callback(this.Matslise);
    else this.callbacks.push(callback);
  }

  catch(callback) {}
}

export default new MatsliseLoader();
