import React, { Component } from "react";
import ParsedInput from "./ParsedInput";
import DF from "../lib/Differentiable";

class SettingsForm extends Component {
  static defaultProps = {
    onSubmit: () => {},
    onInit: () => {}
  };
  state = {
    x: [DF.parse("0"), DF.parse("1")],
    ymin: [DF.parse("1"), DF.parse("0")],
    ymax: [DF.parse("0"), DF.parse("1")],
    parsed: DF.parse("(1-cos(2*pi*x))/2*1000", ["x"]), //"0"),
    symmetric: false
  };

  componentDidMount() {
    this.props.onInit(this.calculateData());
  }

  render() {
    window.DF = DF;
    const { symmetric, ymin, ymax, x } = this.state;

    return (
      <form
        onSubmit={e => {
          e.preventDefault();
          this.onSubmit();
        }}
      >
        <div className="form-row">
          <label className="input-group col-12">
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
              "input-group col-6 " + (symmetric ? "disabled" : "enabled")
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
          <label className="input-group col-6">
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
          <label className="input-group col-12">
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
          </label>
        </div>

        <div className="form-row">
          <label className="input-group col-12">
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
          </label>
          <label>
            <input
              type="checkbox"
              checked={symmetric}
              onChange={e => {
                this.setState({ symmetric: !!e.target.checked }, () =>
                  this.onSubmit()
                );
              }}
            />{" "}
            Symmetric
          </label>
        </div>
      </form>
    );
  }

  onSubmit() {
    this.props.onSubmit(this.calculateData());
  }

  calculateData() {
    const { parsed, x, ymin, ymax, symmetric } = this.state;
    const f = parsed.toFunction("x");
    f.value = parsed.value;
    const x_val = x.map(v => v.toValue());
    const ymax_val = ymax.map(v => v.toValue());
    return {
      ymin: symmetric ? ymax_val : ymin.map(v => v.toValue()),
      ymax: ymax_val,
      x: symmetric ? [-x_val[1], x_val[1]] : x_val,
      f: f,
      symmetric
    };
  }
}

export default SettingsForm;
