import React, { Component } from "react";
import SettingsForm from "./SettingsForm";
import Graph from "./Graph";
import { ReferenceLine } from "recharts";
import loadMatslise from "../lib/loadMatslise";
import Color from "color";

let Matslise = null;

class MatsliseGUI extends Component {
  state = {
    f: x => 0,
    x: [0, 1],
    ymin: [1, 0],
    ymax: [1, 0],
    matslise: null,
    eigenvalues: [],
    showPotential: true,
    symmetric: false,
    sectorPoints: [],
    sectors: {
      type: "auto",
      //count: 31
      tolerance: 1e-3
    }
  };
  toDelete = [];

  componentDidMount() {
    if (Matslise == null) {
      loadMatslise.then(_Matslise => {
        Matslise = _Matslise;
        this.updateFunction();
      });
    }
  }

  render() {
    if (!this.state.matslise) return "Loading";

    const { x, f } = this.state;

    let funcs = [];
    if (this.state.showPotential) {
      const potential = this.state.symmetric ? x => f(x < 0 ? -x : x) : f;
      potential.color = "#000000";
      funcs.push(potential);
    }
    this.state.eigenvalues.forEach(eigenvalue => {
      const { value, show, index } = eigenvalue;
      if (show) {
        if (eigenvalue.eigenfunction === undefined) {
          eigenvalue.eigenfunction = this.eigenfunction(value, index);
          eigenvalue.eigenfunction.color = Color.hsl(
            ((index * 1.618) % 1) * 360,
            80,
            50
          ).hex();
        }
        funcs.push(eigenvalue.eigenfunction);
      }
    });

    return (
      <div className="container-fluid h-100">
        <h1
          style={{
            position: "absolute",
            top: 0,
            width: "100%",
            height: "50px"
          }}
        >
          Schrödinger
        </h1>
        <div className="row h-100">
          <div
            className="col-5 col-md-4 col-xl-3"
            style={{
              height: "100vh",
              overflowX: "auto",
              borderTop: "solid transparent 60px"
            }}
          >
            <SettingsForm
              onSubmit={data => this.updateFunction(data)}
              onInit={data => this.updateFunction(data)}
            />
            <div className="table-responsive">
              <table className="table table-sm table-striped">
                <tbody>
                  <tr
                    onClick={_ =>
                      this.setState({
                        showPotential: !this.state.showPotential
                      })
                    }
                  >
                    <th>
                      <input
                        type="checkbox"
                        checked={this.state.showPotential}
                        onChange={_ => _}
                      />
                    </th>
                    <td />
                    <td colSpan="2">V(x) = {f.value}</td>
                  </tr>
                  {this.state.eigenvalues.map(
                    ({ index, show, value, error }, i) => (
                      <tr key={i} onClick={_ => this.toggleShowEigenvalue(i)}>
                        <th>
                          <input
                            type="checkbox"
                            checked={!!show}
                            onChange={_ => _}
                          />
                        </th>
                        <td>{index}</td>
                        <td>{value.toFixed(10)}</td>
                        <td>{error.toExponential(1)}</td>
                      </tr>
                    )
                  )}
                  <tr>
                    <td colSpan="4">
                      <button
                        onClick={e => this.moreEigenvalues()}
                        className="btn btn-link"
                      >
                        + more eigenvalues
                      </button>
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>
          <div
            className="col-7 col-md-8 col-9-xl"
            style={{ minHeight: "300px", maxHeight: "100vh" }}
          >
            <Graph
              symmetricY={true}
              x={x}
              func={funcs}
              grid={false}
              strokeWidth={1.25}
            >
              {this.state.sectorPoints.map(p => (
                <ReferenceLine
                  key={"ref_" + p}
                  x={p}
                  stroke="#f99"
                  strokeDasharray="3 3"
                />
              ))}
            </Graph>
          </div>
        </div>
      </div>
    );
  }

  moreEigenvalues() {
    let eigenvalues = this.state.eigenvalues;
    let n = eigenvalues.length;
    eigenvalues = eigenvalues.concat(this.eigenvalues(n, n + 10));
    this.setState({ eigenvalues });
  }

  toggleShowEigenvalue(i) {
    let e = this.state.eigenvalues[i];
    e.show = !e.show;
    this.setState({ eigenvalues: this.state.eigenvalues });
  }

  updateFunction(newState) {
    let { f, x, ymin, ymax, symmetric, sectors } = {
      ...this.state,
      ...newState
    };

    while (this.toDelete.length > 0) this.toDelete.pop().delete();

    const matslise = symmetric
      ? new Matslise.HalfRange(f, x[1], sectors)
      : new Matslise.Matslise(f, ...x, sectors);

    this.toDelete.push(matslise);
    this.setState({
      f,
      x,
      ymin,
      ymax,
      matslise,
      symmetric,
      sectors,
      eigenvalues: this.eigenvalues(0, 10, { matslise, ymin, ymax, symmetric }),
      showPotential: true,
      sectorPoints: matslise.sectorPoints()
    });
  }

  eigenvalues(imin, imax, options = {}) {
    const { matslise, ymin, ymax, symmetric } = {
      ...this.state,
      ...options
    };
    const [ya, dya] = ymin,
      [yb, dyb] = ymax;
    return matslise
      .eigenvaluesByIndex(imin, imax, ...(symmetric ? [] : [[dya, -ya]]), [
        dyb,
        -yb
      ])
      .map(({ first: index, second: E }) => {
        return {
          index: index,
          value: E,
          error: symmetric
            ? matslise.eigenvalueError(E, [dyb, -yb], 1 - (index % 2))
            : matslise.eigenvalueError(E, [dya, -ya], [dyb, -yb])
        };
      });
  }

  eigenfunction(e, index) {
    const {
      ymin: [ya, dya],
      ymax: [yb, dyb],
      symmetric
    } = this.state;
    let eigen = symmetric
      ? this.state.matslise.eigenfunction(e, [dyb, -yb], 1 - (index % 2))
      : this.state.matslise.eigenfunction(e, [dya, -ya], [dyb, -yb]);
    this.toDelete.push(eigen);
    return x => eigen(x)[0];
  }
}

export default MatsliseGUI;
