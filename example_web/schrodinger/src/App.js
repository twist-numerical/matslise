/* global Module */


import React, { Component } from 'react';
import './App.css';
import ParsedInput from './components/ParsedInput'
import Graph from './components/Graph'

class App extends Component {
  state = {
    parsed: 0,
    x: [-5.5, 5.5]
  }
  Matslise = false;

  constructor(...args) {
    super(...args);

    if(!Module)
      window.Module = {};
    if(!Module.Matslise && !Module.onRuntimeInitialized) {
      Module.onRuntimeInitialized = () => {
        this.Matslise = Module.Matslise;
        this.forceUpdate();
      }
    }
  }

  render() {
    if(!this.Matslise) {
      return "Loading";
    }

    let func = this.state.parsed.toFunction("x");

    let allFuncs = [func, this.eigenfunction(func, 3)];

    return (
      <div className="App">
      <ParsedInput onParsed={parsed => this.setState({parsed})} value="x^2" />
      <Graph
      x={this.state.x}
      func={allFuncs} />
      </div>
      );
  }

  eigenfunction(f, e) {
    const ms = new Module.Matslise(f, ...this.state.x, 32);
    let eigen = ms.eigenfunction(e, [0,1], [0,1]);
    return (x) => eigen(x)[0];
  }
}

export default App;
