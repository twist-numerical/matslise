import React, { Component } from "react";
import ParsedInput from "./ParsedInput";
import DF from "../lib/Differentiable";
import problems from "./problems";

class SettingsForm extends Component {
  static defaultProps = {
    onSubmit: () => {},
    onInit: () => {}
  };
  state = {
    x: [DF.parse("0"), DF.parse("1")],
    ymin: [DF.parse("1"), DF.parse("0")],
    ymax: [DF.parse("1"), DF.parse("0")],
    parsed: DF.parse("0"), //"(1-cos(2*pi*x))/2*1000", ["x"])
    symmetric: false,
    loadDialog: false,
    sectors: {
      type: "auto",
      tolerance: DF.parse("1e-3"),
      count: DF.parse("30")
    }
  };

  componentDidMount() {
    this.props.onInit(this.calculateData());
  }

  render() {
    window.DF = DF;
    const { symmetric, ymin, ymax, x } = this.state;
    const sectorTypeValue = {
      auto: "tolerance",
      uniform: "count"
    }[this.state.sectors.type];
    const sectorTypeValueName = {
      auto: "Tolerance",
      uniform: "Count"
    }[this.state.sectors.type];

    return (
      <form
        onSubmit={e => {
          e.preventDefault();
          this.onSubmit();
        }}
        className={this.state.loadDialog ? "modal-open" : ""}
      >
        <div
          className="modal show fade"
          style={this.state.loadDialog ? { display: "block" } : {}}
          tabIndex="-1"
          role="dialog"
        >
          <div className="modal-dialog" role="document">
            <div className="modal-content">
              <div className="modal-header">
                <h5 className="modal-title" id="exampleModalLabel">
                  Load a problem
                </h5>
                <button
                  className="close"
                  onClick={() =>
                    this.setState({
                      loadDialog: false
                    })
                  }
                >
                  &times;
                </button>
              </div>
              <div className="modal-body">
                {Object.keys(problems)
                  .sort()
                  .map(name => {
                    const problem = problems[name];
                    return (
                      <button
                        key={name}
                        onClick={() => {
                          const newState = {
                            loadDialog: false,
                            symmetric: !!problem.symmetric,
                            parsed: DF.parse(problem.f, ["x"])
                          };
                          ["x", "ymin", "ymax"].forEach(k => {
                            newState[k] = problem[k].map(v => DF.parse(v));
                          });
                          this.setState(newState);
                        }}
                        className="btn btn-link"
                      >
                        {name}
                      </button>
                    );
                  })}
              </div>
              <div className="modal-footer">
                <button
                  className="btn btn-secondary"
                  onClick={() =>
                    this.setState({
                      loadDialog: false
                    })
                  }
                >
                  Cancel
                </button>
              </div>
            </div>
          </div>
        </div>
        <div
          className="modal-backdrop show"
          style={{ display: this.state.loadDialog ? "block" : "none" }}
        />
        <div className="form-row">
          <label className="input-group col-12 mb-2">
            <div className="input-group-prepend">
              <div className="input-group-text">V(x) = </div>
            </div>
            <ParsedInput
              className="form-control"
              onParsed={parsed => {
                this.setState({ parsed });
              }}
              parsed={this.state.parsed}
              variables={["x"]}
            />
          </label>
        </div>

        <div className="form-row">
          <label
            className={
              "input-group col-6 mb-2 " + (symmetric ? "disabled" : "enabled")
            }
          >
            <div className="input-group-prepend">
              <div className="input-group-text">
                x<sub>min</sub>&nbsp;=
              </div>
            </div>
            <ParsedInput
              disabled={symmetric}
              className="form-control"
              onParsed={val => this.setState({ x: [val, x[1]] })}
              parsed={symmetric ? DF.parse("-(" + x[1].value + ")") : x[0]}
            />
          </label>
          <label className="input-group col-6 mb-2">
            <div className="input-group-prepend">
              <div className="input-group-text">
                x<sub>max</sub>&nbsp;=
              </div>
            </div>
            <ParsedInput
              className="form-control"
              onParsed={val => this.setState({ x: [x[0], val] })}
              parsed={x[1]}
            />
          </label>
        </div>

        <div className={"form-row " + (symmetric ? "disabled" : "enabled")}>
          <div className="input-group col-12 mb-2">
            <ParsedInput
              disabled={symmetric}
              className="form-control"
              onParsed={val => this.setState({ ymin: [val, ymin[1]] })}
              parsed={symmetric ? ymax[0] : ymin[0]}
            />
            <div className="input-group-append input-group-prepend">
              <span className="input-group-text">
                y(x<sub>min</sub>) +{" "}
              </span>
            </div>
            <ParsedInput
              disabled={symmetric}
              className="form-control"
              onParsed={val => this.setState({ ymin: [ymin[0], val] })}
              parsed={symmetric ? ymax[1] : ymin[1]}
            />
            <div className="input-group-append">
              <span className="input-group-text">
                y'(x<sub>min</sub>)&nbsp;=&nbsp;0
              </span>
            </div>
          </div>
        </div>

        <div className="form-row">
          <div className="input-group col-12 mb-2">
            <ParsedInput
              className="form-control"
              onParsed={val => this.setState({ ymax: [val, ymax[1]] })}
              parsed={ymax[0]}
            />
            <div className="input-group-append input-group-prepend">
              <span className="input-group-text">
                y(x<sub>max</sub>) +{" "}
              </span>
            </div>
            <ParsedInput
              className="form-control"
              onParsed={val => this.setState({ ymax: [ymax[0], val] })}
              parsed={ymax[1]}
            />
            <div className="input-group-append">
              <span className="input-group-text">
                y'(x<sub>max</sub>)&nbsp;=&nbsp;0
              </span>
            </div>
          </div>
        </div>
        <div className="form-row">
          <label className="input-group col-6 mb-2">
            <div className="input-group-prepend">
              <span className="input-group-text">Mesh</span>
            </div>
            <select
              className="custom-select"
              value={this.state.sectors.type}
              onChange={e => {
                this.setState(
                  {
                    sectors: {
                      ...this.state.sectors,
                      type: e.target.value
                    }
                  },
                  () => this.onSubmit()
                );
              }}
            >
              {["auto", "uniform"].map(t => (
                <option value={t} key={t}>
                  {t}
                </option>
              ))}
            </select>
          </label>
          <label className="input-group col-6 mb-2">
            <div className="input-group-prepend">
              <span className="input-group-text">{sectorTypeValueName}</span>
            </div>
            <ParsedInput
              className="form-control"
              onParsed={val =>
                this.setState({
                  sectors: {
                    ...this.state.sectors,
                    [sectorTypeValue]: val
                  }
                })
              }
              parsed={this.state.sectors[sectorTypeValue]}
            />
          </label>
        </div>
        <div className="form-row">
          <label className="input-group col-6">
            <div className="input-group">
              <div className="input-group-prepend">
                <span className="input-group-text">
                  <input
                    type="checkbox"
                    checked={symmetric}
                    onChange={e => {
                      this.setState({ symmetric: !!e.target.checked }, () =>
                        this.onSubmit()
                      );
                    }}
                  />
                </span>
              </div>
              <div className="input-group-append">
                <span className="input-group-text">Symmetric</span>
              </div>
            </div>
          </label>
          <div className="input-group col-6 justify-content-end">
            <button
              className="btn btn-link"
              onClick={() =>
                this.setState({
                  loadDialog: true
                })
              }
            >
              Load a problem
            </button>
          </div>
        </div>
      </form>
    );
  }

  onSubmit() {
    this.props.onSubmit(this.calculateData());
  }

  calculateData() {
    const {
      parsed,
      x,
      ymin,
      ymax,
      symmetric,
      sectors: { type, tolerance, count }
    } = this.state;
    const f = parsed.toFunction("x");
    f.value = parsed.value;
    const x_val = x.map(v => v.toValue());
    const ymax_val = ymax.map(v => v.toValue());
    return {
      ymin: symmetric ? ymax_val : ymin.map(v => v.toValue()),
      ymax: ymax_val,
      x: symmetric ? [-x_val[1], x_val[1]] : x_val,
      f: f,
      sectors: {
        type,
        tolerance: tolerance.toValue(),
        count: count.toValue()
      },
      symmetric
    };
  }
}

export default SettingsForm;
