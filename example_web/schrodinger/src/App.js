/* global Module */


import React, { Component } from 'react';
import './App.css';
import ParsedInput from './components/ParsedInput'
import Graph from './components/Graph'

class App extends Component {
  state = {
    f: (x) => 0,
    x: [0, Math.PI],
    steps: 32,
    matslise: null
  };
  toDelete = [];

  constructor(...args) {
    super(...args);

    if(!Module)
      window.Module = {};
    if(!Module.Matslise && !Module.onRuntimeInitialized) {
      Module.onRuntimeInitialized = () => {
        this.updateFunction();
      }
    }
  }

  render() {
    if(!this.state.matslise)
      return "Loading";
    
    const {x, f} = this.state;
    window.f= f;
    let funcs = [f];
    this.eigenvalues(0, 2).forEach(({value}) => 
      funcs.push(this.eigenfunction(value)));

    return (
      <div className="App">
      <ParsedInput
      onParsed={parsed => {
        window.p = parsed;
        this.updateFunction({f: parsed.toFunction("x")});
      }}
      value="2*cos(2*x)" />
      <Graph
      x={x}
      func={funcs} />
      </div>
      );
  }

  updateFunction(newState) {
    let {f, x, steps} = {
      ...this.state,
      ...newState
    };

    while(this.toDelete.length > 0)
      this.toDelete.pop().delete();

    const matslise = new Module.Matslise(f, ...x, steps);

    this.toDelete.push(matslise);
    this.setState({f, x, steps, matslise});
  }

  eigenvalues(imin, imax) {
    return this.state.matslise.eigenvaluesByIndex(imin, imax, [0,1], [0,1]);
  }

  eigenfunction(e) {
    let eigen = this.state.matslise.eigenfunction(e, [0,1], [0,1]);
    this.toDelete.push(eigen);
    return (x) => eigen(x)[0];
  }
}

export default App;
