import React, { Component } from "react";
import ParsedInput from "../ParsedInput";
import DF from "../../lib/Differentiable";

class Settings extends Component {
  static defaultProps = {
    onSubmit: () => {},
    onInit: () => {}
  };
  state = {
    x: [DF.parse("-5.5"), DF.parse("5.5")],
    y: [DF.parse("-5.5"), DF.parse("5.5")],
    parsed: DF.parse("(1+x^2)*(1+y^2)", ["x", "y"])
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
              onParsed={parsed => {
                this.setState({ parsed });
              }}
              parsed={this.state.parsed}
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
    const { parsed, x, y } = this.state;
    const f = parsed.toFunction("x");
    f.value = parsed.value;
    return {
      x: x.map(v => v.toValue()),
      y: y.map(v => v.toValue()),
      f: f
    };
  }
}

export default Settings;
