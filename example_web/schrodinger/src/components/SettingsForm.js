import React, { Component } from 'react';
import ParsedInput from './ParsedInput';
import DF from '../lib/Differentiable';

class SettingsForm extends Component {
  static defaultProps = {
    onSubmit: () => {},
    onInit: () => {},
  };
  state = {
    x: [DF.parse('0'), DF.parse('pi')],
    ymin: [DF.parse('1'),DF.parse('0')],
    ymax: [DF.parse('1'),DF.parse('0')],
    parsed: DF.parse('0'),
  };

  componentDidMount() {
    this.props.onInit(this.calculateData());
  }

  render() {
    window.DF = DF;
    const {ymin, ymax, x} = this.state;

    return <form onSubmit={e => {
      e.preventDefault();
      this.onSubmit();
    }}>

    <div className="form-row">
      <label className="input-group col-12">
        <div className="input-group-prepend">
          <div className="input-group-text">V(x) = </div>
        </div>
        <ParsedInput
          className="form-control"
          onParsed={parsed => {this.setState({parsed})}}
          parsed={this.state.parsed}
          variables={['x']} />
      </label>
    </div>

    <div className="form-row">
      <label className="input-group col-6">
        <div className="input-group-prepend">
          <div className="input-group-text">x<sub>min</sub>&nbsp;=</div>
        </div>
        <ParsedInput className="form-control"
        onParsed={(val) => this.setState({x: [val, x[1]]})}
        parsed={x[0]} />
      </label>
      <label className="input-group col-6">
        <div className="input-group-prepend">
          <div className="input-group-text">x<sub>max</sub>&nbsp;=</div>
        </div>
        <ParsedInput className="form-control"
        onParsed={(val) => this.setState({x: [x[0], val]})}
        parsed={x[1]} />
      </label>
    </div>

    <div className="form-row">
      <label className="input-group col-12">
        <ParsedInput className="form-control"
        onParsed={(val) => this.setState({ymin: [val, ymin[1]]})}
        parsed={ymin[0]} />
        <div className="input-group-append input-group-prepend">
          <span className="input-group-text">y(x<sub>min</sub>) + </span>
        </div>
        <ParsedInput className="form-control"
        onParsed={(val) => this.setState({ymin: [ymin[0], val]})}
        parsed={ymin[1]} />
        <div className="input-group-append">
          <span className="input-group-text">y'(x<sub>min</sub>)&nbsp;=&nbsp;0</span>
        </div>
      </label>
    </div>

    <div className="form-row">
      <label className="input-group col-12">
        <ParsedInput className="form-control"
        onParsed={(val) => this.setState({ymax: [val, ymax[1]]})}
        parsed={ymax[0]} />
        <div className="input-group-append input-group-prepend">
          <span className="input-group-text">y(x<sub>max</sub>) + </span>
        </div>
        <ParsedInput className="form-control"
        onParsed={(val) => this.setState({ymax: [ymax[0], val]})}
        parsed={ymax[1]} />
        <div className="input-group-append">
          <span className="input-group-text">y'(x<sub>max</sub>)&nbsp;=&nbsp;0</span>
        </div>
      </label>
    </div>

    </form>
  }

  onSubmit() {
    this.props.onSubmit(this.calculateData());
  }

  calculateData() {
    const {parsed, x, ymin, ymax} = this.state;
    const f = parsed.toFunction('x');
    f.value = parsed.value;
    return {
      ymin: ymin.map(v => v.toValue()),
      ymax: ymax.map(v => v.toValue()),
      x: x.map(v => v.toValue()),
      f: f,
    };
  }
}

export default SettingsForm;