import React, { Component } from 'react';
import { ComposedChart, CartesianGrid, ResponsiveContainer, XAxis, YAxis, Line } from 'recharts';


class Graph extends Component {
  static defaultProps = {
    numPoints: 200,
    func: x => x,
    x: [-1, 1],
  };

  constructor(props) {
    super(props);

    this.state = this.calculateState(props);
  }

  calculateState(props) {
    return {
      x: props.x,
      f: props.func ? Array.isArray(props.func) ? props.func : [props.func] : [],
    };
  }

  getFunc() {
    if(this.props.func === undefined || this.props.func === null)
      return [];
    if(Array.isArray(this.props.func))
      return this.props.func;
    return [this.props.func];
  }

  componentWillReceiveProps(newProps) {
    this.setState(this.calculateState(newProps));
  }

  render() {
    const {x: [xmin, xmax]} = this.state;
    let ymin = 'auto', ymax = 'auto';
    let func = this.getFunc();

    const data = [];
    if(func.length > 0) {
      for(let i = 0; i < this.props.numPoints; ++i) {
        const x = xmin + (xmax-xmin)*(i+.5+Math.random()/4)/this.props.numPoints;
        let obj = {
          x: x,
        };
        func.forEach((f, i) => {
          obj['y'+i] = f(x);
        });
        data.push(obj);
      }

      ymin = Number.POSITIVE_INFINITY;
      ymax = Number.NEGATIVE_INFINITY;

      data.forEach((point) => {
        this.getFunc().forEach((_, i) => {
          let y = point['y'+i];
          if(y < ymin)
            ymin = y;
          if(y > ymax)
            ymax = y;
        });
      })

      let d = Math.max(ymax-ymin, 1e-13);
      let step = 1;
      while(step < d)
        step *= 1e1;
      while(step > d)
        step /= 1e1;

      if(step < 1) {
        const istep = Math.round(1/step);
        ymin = Math.floor(ymin*istep)/istep;
        ymax = Math.ceil(ymax*istep)/istep;
      } else {
        ymin = Math.floor(ymin/step)*step;
        ymax = Math.ceil(ymax/step)*step;
      }
    }

    return (
      <ResponsiveContainer>
      <ComposedChart
      data={data}>
      <CartesianGrid strokeDasharray="3 3"/>
      <XAxis 
      allowDataOverflow={true}
      dataKey="x"
      type="number"
      domain={[xmin, xmax]}
      />
      <YAxis 
      allowDataOverflow={true}
      type="number"
      domain={[ymin, ymax]}
      />
      { func.map((_, i) => 
        <Line key={i} isAnimationActive={false}
        dataKey={'y'+i}
        stroke={'#' + ('000000'+Math.floor(Math.random()*0x10000007).toString(16)).substr(-6)}
        dot={false} />
        )}
      </ComposedChart>
      </ResponsiveContainer>
      );
  }
}

export default Graph;
