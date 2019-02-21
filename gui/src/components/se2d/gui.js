import React, { Component } from "react";
import Settings from "./Settings";
import Eigenfunction from "./Eigenfunction";
import interpolate from "../../lib/interpolate";
import Graph from "../Graph";
import { ReferenceDot } from "recharts";
import WorkerDelegator from "../../lib/WorkerDelegator";

class MatsliseGUI extends Component {
  state = {
    calculating: true,
    settings: null,
    E: [0, 10],
    drawErrors: false,
    errors: null,
    time: 0,
    eigenvalues: []
  };
  worker = null;
  errorHistory = [];

  constructor() {
    super();
    this.worker = new WorkerDelegator(
      new Worker("./se2d.worker.js", { type: "module" })
    );
    const setState = state => this.setState(state);
    this.worker.addListener("init", setState);
    this.worker.addListener("errors", setState);
    this.worker.addListener("locate", setState);
    this.worker.addListener("time", time => this.setState({ time }));
  }

  render() {
    const [Emin, Emax] = this.state.E;
    const Ed = (Emax - Emin) / 2;
    if (this.state.errors) {
      if (this.errorHistory.length > 5) {
        this.errorHistory.pop();
      }
      this.errorHistory.unshift(
        interpolate.linear(this.state.errors.E, this.state.errors.errors, NaN)
      );
    }

    return (
      <div className="container-fluid h-100">
        <div
          style={{
            position: "absolute",
            bottom: "5px",
            right: "5px",
            zIndex: 10,
            textAlign: "right"
          }}
        >
          {this.state.calculating ? ["Loading...", <br key="br" />] : ""}
          Last calculation: {this.state.time}ms
        </div>
        <h1
          style={{
            position: "absolute",
            top: 0,
            width: "100%",
            height: "50px"
          }}
        >
          2D Schrödinger
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
            <Settings
              onSubmit={data => this.updateFunction(data)}
              onInit={data => this.updateFunction(data)}
            />
            <div style={{ height: "200px" }}>
              <Graph
                selectable={true}
                onClick={e => {
                  if (e && isFinite(e.activeLabel))
                    this.findEigenvalue(+e.activeLabel);
                }}
                height={400}
                x={this.state.E}
                func={
                  this.state.errors
                    ? [
                        (() => {
                          const f = x => 0;
                          f.color = "#666";
                          return f;
                        })(),
                        (() => {
                          const f = x => {
                            for (let i = 0; i < this.errorHistory.length; ++i) {
                              const v = this.errorHistory[i](x);
                              if (isFinite(v)) return v;
                            }
                          };
                          f.color = "#ff0000";
                          f.selectable = true;
                          return f;
                        })()
                      ]
                    : []
                }
              >
                {this.state.eigenvalues
                  .filter(E => Emin <= E && E <= Emax)
                  .map(E => (
                    <ReferenceDot
                      key={E}
                      x={E}
                      y={0}
                      r={3}
                      fill={"#30a"}
                      ifOverflow={"extendDomain"}
                    />
                  ))}
              </Graph>
            </div>
            {this.state.errors ? (
              <div className="text-center">
                <div
                  className="btn-group"
                  role="group"
                  aria-label="Basic example"
                >
                  <button
                    className="btn"
                    onClick={() => this.setE(Emin - Ed, Emax - Ed)}
                  >
                    &lt;&lt;
                  </button>
                  <button
                    className="btn"
                    onClick={() => this.setE(Emin - Ed, Emax + Ed)}
                  >
                    -
                  </button>
                  <button className="btn" onClick={() => this.setE(0, 10)}>
                    0
                  </button>
                  <button
                    className="btn"
                    onClick={() => this.setE(Emin + Ed / 2, Emax - Ed / 2)}
                  >
                    +
                  </button>
                  <button
                    className="btn"
                    onClick={() => this.setE(Emin + Ed, Emax + Ed)}
                  >
                    &gt;&gt;
                  </button>
                </div>
              </div>
            ) : null}
            <div>
              <ul>
                {this.state.eigenvalues.map(E => (
                  <li key={E}>
                    <button
                      className="btn btn-link"
                      onClick={e => {
                        e.preventDefault();
                        this.viewEigenfunctions(E);
                      }}
                    >
                      {E}
                    </button>
                  </li>
                ))}
              </ul>
            </div>
          </div>
          <div className="col-7 col-md-8 col-9-xl">
            <Eigenfunction
              E={this.state.eigenfunction}
              worker={this.worker}
              xn={60}
              yn={60}
            />
          </div>
        </div>
      </div>
    );
  }

  viewEigenfunctions(E) {
    this.setState({ eigenfunction: E });
  }

  findEigenvalue(Eguess) {
    this.setState({
      calculating: true
    });
    this.worker.send("locate", {
      E: Eguess,
      eigenvalues: this.state.eigenvalues
    });
  }

  setE(Emin, Emax) {
    this.setState({
      calculating: true,
      E: [Emin, Emax]
    });
    this.worker.send("errors", {
      E: [Emin, Emax]
    });
  }

  updateFunction(data) {
    const E = [0, 10];
    this.setState({
      init: false,
      settings: data,
      errors: null,
      drawErrors: false,
      calculating: true,
      eigenvalues: [],
      E
    });
    this.errorHistory = [];
    this.worker.send("init", {
      ...data,
      f: data.f.value,
      E
    });
  }
}

export default MatsliseGUI;
