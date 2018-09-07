/* global Module */


import React, { Component } from 'react';
import './App.css';
import SettingsForm from './components/SettingsForm'
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

    let funcs = [f];
    this.eigenvalues(0, 50).forEach(({value}) => 
      funcs.push(this.eigenfunction(value)));

    return (
      <div className="container-fluid">
        <h1>Schrödinger</h1>
        <div className="row">
          <div className="col-4">
            <SettingsForm onSubmit={data => this.updateFunction(data)} />
          </div>
          <div className="col-8"  style={{minHeight: '300px'}}>
            <Graph
            x={x}
            func={funcs} />
          </div>
        </div>
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
