import React, { Component } from 'react';
import SettingsForm from './SettingsForm'
import Graph from './Graph'
import loadMatslise from '../lib/loadMatslise'

let Matslise = null;

class MatsliseGUI extends Component {
  state = {
    f: (x) => 0,
    x: [0, 1],
    steps: 32,
    matslise: null,
    eigenvalues: [],
    showPotential: true,
  };
  toDelete = [];

  componentDidMount() {
    if(Matslise == null) {
      loadMatslise.then((_Matslise) => {
        Matslise = _Matslise;
        this.updateFunction();
      });
    }
  }

  render() {
    if(!this.state.matslise)
      return "Loading";
    
    const {x, f} = this.state;

    let funcs = [];
    if(this.state.showPotential)
      funcs.push(f);
    this.state.eigenvalues.forEach((eigenvalue) => {
      const {value, show} = eigenvalue;
      if(show) {
        if(eigenvalue.eigenfunction === undefined)
          eigenvalue.eigenfunction = this.eigenfunction(value);
        funcs.push(eigenvalue.eigenfunction);
      }
    });

    return (
      <div className="container-fluid h-100">
        <h1 style={{position: 'absolute', top:0, width: '100%', height: "50px"}}>Schrödinger</h1>
        <div className="row h-100">
          <div className="col-5 col-md-4 col-xl-3" style={{height: '100vh', overflowX: 'auto', borderTop: "solid transparent 60px"}}>
            <SettingsForm onSubmit={data => this.updateFunction(data)} onInit={data => this.updateFunction(data)} />
            <div className="table-responsive">
              <table className="table table-sm table-striped">
                <tbody>
                  <tr onClick={_ => this.setState({showPotential: !this.state.showPotential})}>
                    <th><input type="checkbox" checked={this.state.showPotential} onChange={_=>_} /></th>
                    <td></td>
                    <td>V(x) = {f.value}</td>
                  </tr>
                  { this.state.eigenvalues.map(({index, show, value}, i) => 
                    <tr key={i} onClick={_ => this.toggleShowEigenvalue(i)}>
                      <th><input type="checkbox" checked={!!show} onChange={_=>_} /></th>
                      <td>{index}</td>
                      <td>{value}</td>
                    </tr>)}
                  <tr>
                    <td colSpan="3">
                      <button onClick={e => this.moreEigenvalues()} className="btn btn-link">
                        + more eigenvalues
                      </button>
                    </td>
                  </tr>
                  </tbody>
              </table>
            </div>
          </div>
          <div className="col-7 col-md-8 col-9-xl"  style={{minHeight: '300px', maxHeight: '100vh'}}>
            <Graph
              symmetricY={true}
              x={x}
              func={funcs} />
          </div>
        </div>
      </div>
      );
  }

  moreEigenvalues() {
    let eigenvalues = this.state.eigenvalues;
    let n = eigenvalues.length;
    eigenvalues = eigenvalues.concat(
      this.eigenvalues(n, n + 10));
    this.setState({eigenvalues});
  }

  toggleShowEigenvalue(i) {
    let e = this.state.eigenvalues[i];
    e.show = !e.show;
    this.setState({eigenvalues: this.state.eigenvalues});
  }

  updateFunction(newState) {
    let {f, x, steps} = {
      ...this.state,
      ...newState
    };
    
    while(this.toDelete.length > 0)
      this.toDelete.pop().delete();

    const matslise = new Matslise.Matslise(f, ...x, steps);

    this.toDelete.push(matslise);
    this.setState({f, x, steps, matslise,
      eigenvalues: this.eigenvalues(0, 10, matslise),
      showPotential: true,
    });
  }

  eigenvalues(imin, imax, matslise=this.state.matslise) {
    return matslise.eigenvaluesByIndex(imin, imax, [0,1], [0,1]);
  }

  eigenfunction(e) {
    let eigen = this.state.matslise.eigenfunction(e, [0,1], [0,1]);
    this.toDelete.push(eigen);
    return (x) => eigen(x)[0];
  }
}

export default MatsliseGUI;
