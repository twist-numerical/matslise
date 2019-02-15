import React, { Component } from 'react';
import { LineChart, CartesianGrid, ResponsiveContainer, XAxis, YAxis, Line } from 'recharts';


class Graph extends Component {
  static defaultProps = {
    numPoints: 200,
    func: x => x,
    x: [-1, 1],
    symmetricY: false,
    y: ['auto', 'auto'],
    xAxis: true,
    yAxis: true,
    grid: true,
    strokeWidth: 1,
  };

  constructor(props) {
    super(props);

    this.state = this.calculateState(props);
  }

  calculateState(props) {
    return {
      x: props.x,
      y: props.y,
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
    let [ymin, ymax] = this.state.y;
    let func = this.getFunc();

    const data = [];
    const color = func.map(f => f.color || '#' + ('000000'+Math.floor(Math.random()*0x1000000).toString(16)).substr(-6));
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

      let newYmin = Number.POSITIVE_INFINITY;
      let newYmax = Number.NEGATIVE_INFINITY;

      data.forEach((point) => {
        this.getFunc().forEach((_, i) => {
          let y = point['y'+i];
          if(y < newYmin)
            newYmin = y;
          if(y > newYmax)
            newYmax = y;
        });
      })

      let d = Math.max(newYmax-newYmin, 1e-13);
      let step = 1;
      while(step < d)
        step *= 1e1;
      while(step > d)
        step /= 1e1;

      if(step < 1) {
        const istep = Math.round(1/step);
        newYmin = Math.floor(newYmin*istep)/istep;
        newYmax = Math.ceil(newYmax*istep)/istep;
      } else {
        newYmin = Math.floor(newYmin/step)*step;
        newYmax = Math.ceil(newYmax/step)*step;
      }

      if(this.props.symmetricY) {
        newYmax = Math.max(newYmax, -newYmin);
        newYmin = -newYmax;
      }

      if(ymin === 'auto')
        ymin = newYmin;
      if(ymax === 'auto')
        ymax = newYmax;
    }

    return (
      <ResponsiveContainer>
      <LineChart
      data={data}>
      {this.props.grid?<CartesianGrid strokeDasharray="3 3"/>:""}
      <XAxis 
      allowDataOverflow={true}
      dataKey="x"
      type="number"
      domain={[xmin, xmax]}
      hide={!this.props.xAxis}
      />
      <YAxis 
      allowDataOverflow={true}
      type="number"
      domain={[ymin, ymax]}
      hide={!this.props.yAxis}
      />
      { func.map((_, i) => 
        <Line key={i} isAnimationActive={false}
        dataKey={'y'+i}
        strokeWidth={this.props.strokeWidth}
        stroke={color[i]}
        dot={false} />
        )}
      {this.props.children}
      </LineChart>
      </ResponsiveContainer>
      );
  }
}

export default Graph;
