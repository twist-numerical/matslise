import React, { Component } from "react";
import Settings from "./Settings";
import ParsedInput from "../ParsedInput";
import DF from "../../lib/Differentiable";
import interpolate from "../../lib/interpolate";
import Graph from "../Graph";

class MatsliseGUI extends Component {
  state = {
    init: false,
    settings: null,
    E: [DF.parse("0"), DF.parse("10")],
    drawErrors: false,
    errors: null,
    time: 0
  };
  worker = null;

  constructor() {
    super();
    this.worker = new Worker("./se2d.worker.js", { type: "module" });
    this.worker.onmessage = ({ data: stateUpdate }) =>
      this.setState(stateUpdate);
  }

  render() {
    return (
      <div className="container-fluid h-100">
        <div style={{ position: "absolute", bottom: "5px", left: "5px" }}>
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
            <div>
              {this.state.init
                ? "Loaded"
                : this.state.settings
                ? "Loading..."
                : "Waiting for input"}
            </div>
          </div>
          <div className="col-7 col-md-8 col-9-xl">
            {this.state.drawErrors ? (
              this.state.errors ? (
                <Graph
                  x={this.state.E.map(v => v.toValue())}
                  func={interpolate.linear(
                    this.state.errors.E,
                    this.state.errors.errors
                  )}
                />
              ) : (
                <div>Loading</div>
              )
            ) : (
              <div style={{ maxWidth: "400px", margin: "10vh auto" }}>
                <div className="row">
                  <label className="input-group col-6">
                    <div className="input-group-prepend">
                      <div className="input-group-text">
                        E<sub>min</sub>&nbsp;=
                      </div>
                    </div>
                    <ParsedInput
                      className="form-control"
                      onParsed={val =>
                        this.setState({ E: [val, this.state.E[1]] })
                      }
                      parsed={this.state.E[0]}
                    />
                  </label>
                  <label className="input-group col-6">
                    <div className="input-group-prepend">
                      <div className="input-group-text">
                        E<sub>max</sub>&nbsp;=
                      </div>
                    </div>
                    <ParsedInput
                      className="form-control"
                      onParsed={val =>
                        this.setState({ E: [this.state.E[0], val] })
                      }
                      parsed={this.state.E[1]}
                    />
                  </label>
                </div>
                <div className="row">
                  <button
                    className="btn btn-primary"
                    style={{ margin: "2vh auto" }}
                    onClick={() => {
                      this.setState({ drawErrors: true });
                      this.worker.postMessage({
                        type: "errors",
                        E: this.state.E.map(v => v.toValue())
                      });
                    }}
                  >
                    Calculate errors
                  </button>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    );
  }

  updateFunction(data) {
    this.setState({
      init: false,
      settings: data,
      errors: null,
      drawErrors: false
    });
    this.worker.postMessage({
      type: "init",
      settings: {
        ...data,
        f: data.f.value
      }
    });
  }
}

export default MatsliseGUI;
