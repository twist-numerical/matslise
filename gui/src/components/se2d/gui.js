import React, { Component } from "react";
import Settings from "./Settings";
import Eigenfunction from "./Eigenfunction";
import interpolate from "../../lib/interpolate";
import Graph from "../Graph";
import { ReferenceDot } from "recharts";
import WorkerDelegator from "../../lib/WorkerDelegator";
import { addUrlProps } from "react-url-query";

const urlPropsQueryConfig = {
  eigenfunction: {
    type: {
      encode: v => (isNaN(v) ? undefined : v),
      decode: v => +v
    },
    addChangeHandlers: false
  },
  settings: {
    type: {
      encode: ({ x, y, V }) => [x.join(","), y.join(","), V].join(","),
      decode: v => {
        if (!v) return v;
        const [xmin, xmax, ymin, ymax, V] = v.split(",");
        return {
          x: [xmin, xmax],
          y: [ymin, ymax],
          V
        };
      }
    },
    addChangeHandlers: false
  }
};

class MatsliseGUI extends Component {
  static defaultProps = {
    eigenfunction: NaN,
    settings: {
      x: [-5.5, 5.5],
      y: [-5.5, 5.5],
      V: "(1+x^2)*(1+y^2)"
    }
  };

  worker = null;
  errorHistory = [];
  state = {
    calculating: true,
    settings: this.props.settings,
    eigenfunction: this.props.eigenfunction,
    E: [0, 10],
    drawErrors: false,
    errors: null,
    last: {
      time: 0,
      action: ""
    },
    eigenvalues: isNaN(this.props.eigenfunction)
      ? []
      : [this.props.eigenfunction]
  };

  constructor(...args) {
    super(...args);
    this.worker = new WorkerDelegator(
      new Worker("./se2d.worker.js", { type: "module" })
    );
    const setState = state => this.setState(state);
    this.worker.addListener("init", setState);
    this.worker.addListener("errors", setState);
    this.worker.addListener("locate", setState);
    this.worker.addListener("time", last => this.setState({ last }));
  }

  render() {
    const [Emin, Emax] = this.state.E;
    const Ed = (Emax - Emin) / 2;
    let maxError = 1;

    const functions = [x => 0];
    functions[0].color = "#666";
    if (this.state.errors) {
      maxError = Graph.nextPower(2., this.state.errors.max);
      this.state.errors.functions.forEach(({ E, errors }) => {
        const f = interpolate.linear(E, errors, NaN);
        f.color = "#ff0000";
        f.selectable = true;
        functions.push(f);
      });
    }
    window.functions = functions;

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
          Last calculation ({this.state.last.action}): {this.state.last.time}ms
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
              onInit={data => this.updateFunction(data, false)}
              {...this.state.settings}
            />
            <div style={{ height: "200px" }}>
              <Graph
                selectable={true}
                onClick={e => {
                  if (e && isFinite(e.activeLabel))
                    this.findEigenvalue(+e.activeLabel);
                }}
                height={400}
                y={[-maxError, maxError]}
                x={this.state.E}
                func={functions}
              >
                {this.state.eigenvalues
                  .filter(E => Emin <= E && E <= Emax)
                  .map(E => (
                    <ReferenceDot
                      key={E}
                      x={E}
                      y={0}
                      r={E === this.state.eigenfunction ? 4 : 3}
                      fill={E === this.state.eigenfunction ? "#03a" : "#30a"}
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
                      className={
                        "btn btn-link" +
                        (E === this.state.eigenfunction
                          ? " font-weight-bold"
                          : "")
                      }
                      onClick={e => this.viewEigenfunctions(E)}
                    >
                      {E}
                    </button>
                  </li>
                ))}
              </ul>
            </div>
          </div>
          <div className="col-7 col-md-8 col-xl-9" style={{ paddingRight: 0 }}>
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

  setStateAndUrl(state) {
    this.props.onChangeUrlQueryParams(state);
    this.setState(state);
  }

  viewEigenfunctions(E) {
    this.setStateAndUrl({ eigenfunction: E });
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

  updateFunction(data, reset = true) {
    const E = [0, 10];
    const settings = {
      ...data,
      V: data.V.value
    };
    this.setStateAndUrl({
      settings,
      ...(reset ? { eigenfunction: NaN } : {})
    });
    this.setState({
      init: false,
      errors: null,
      drawErrors: false,
      calculating: true,
      E,
      ...(reset
        ? {
            eigenvalues: [],
            eigenfunction: NaN
          }
        : {})
    });
    this.errorHistory = [];
    this.worker.send("init", {
      ...settings,
      E
    });
  }
}

export default addUrlProps({ urlPropsQueryConfig })(MatsliseGUI);
