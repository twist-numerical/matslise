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
    parsed: DF.parse('0'),
  };

  componentDidMount() {
    this.props.onInit(this.calculateData());
  }

  render() {
    window.DF = DF;
    const {x} = this.state;

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

    </form>
  }

  onSubmit() {
    this.props.onSubmit(this.calculateData());
  }

  calculateData() {
    const {parsed, x} = this.state;
    window.parsed = parsed;
    const f = parsed.toFunction('x');
    f.value = parsed.value;
    return {
      x: x.map(v => v.toFunction()()),
      f: f,
    };
  }
}

export default SettingsForm;