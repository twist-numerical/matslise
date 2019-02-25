import React, { Component } from "react";
import ParsedInput from "../ParsedInput";
import DF from "../../lib/Differentiable";

const zero = DF.parse("0");
const parse = vars => t => {
  try {
    return DF.parse("" + t, vars);
  } catch (e) {
    console.error(e);
    return zero;
  }
};

class Settings extends Component {
  static defaultProps = {
    onSubmit: () => {},
    onInit: () => {},
    x: ["-1", "1"],
    y: ["-1", "1"],
    V: "0"
  };
  state = {
    x: this.props.x.map(parse([])),
    y: this.props.y.map(parse([])),
    V: parse(["x", "y"])(this.props.V)
  };

  componentDidMount() {
    this.props.onInit(this.calculateData());
  }

  render() {
    window.DF = DF;
    const { x, y } = this.state;

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
              onParsed={V => {
                this.setState({ V });
              }}
              parsed={this.state.V}
              variables={["x", "y"]}
            />
          </label>
        </div>

        <div className="form-row">
          <label className="input-group col-6">
            <div className="input-group-prepend">
              <div className="input-group-text">
                x<sub>min</sub>&nbsp;=
              </div>
            </div>
            <ParsedInput
              className="form-control"
              onParsed={val => this.setState({ x: [val, x[1]] })}
              parsed={x[0]}
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

        <div className="form-row">
          <label className="input-group col-6">
            <div className="input-group-prepend">
              <div className="input-group-text">
                y<sub>min</sub>&nbsp;=
              </div>
            </div>
            <ParsedInput
              className="form-control"
              onParsed={val => this.setState({ y: [val, y[1]] })}
              parsed={y[0]}
            />
          </label>
          <label className="input-group col-6">
            <div className="input-group-prepend">
              <div className="input-group-text">
                y<sub>max</sub>&nbsp;=
              </div>
            </div>
            <ParsedInput
              className="form-control"
              onParsed={val => this.setState({ y: [y[0], val] })}
              parsed={y[1]}
            />
          </label>
        </div>
      </form>
    );
  }

  onSubmit() {
    this.props.onSubmit(this.calculateData());
  }

  calculateData() {
    const { V, x, y } = this.state;
    const f = V.toFunction("x", "y");
    f.value = V.value;
    return {
      x: x.map(v => v.toValue()),
      y: y.map(v => v.toValue()),
      V: f
    };
  }
}

export default Settings;
