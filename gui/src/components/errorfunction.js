import React, { Component } from 'react';
import Graph from './Graph'
import loadMatslise from '../lib/loadMatslise'
import { ReferenceDot } from 'recharts';

const Es = [3.20, 5.53, 7.55, 8.04, 8.46, 9.94, 11.34];

const binarySearch = (f, a, b) => {
  let fa = f(a), fb = f(b);
  if(fa === 0)
    return a;
  if(fb === 0)
    return b;
  if((fa < 0) === (fb < 0))
    throw new Error("f(a) and f(b) must have a different sign");
  let c, fc; 
  while(Math.abs(a - b) > 1e-12) {
    c = (a+b)/2;
    fc = f(c);
    if(fc === 0)
      return c;
    if((fa < 0) === (fc < 0)) {
      fa = fc;
      a = c;
    } else {
      fb = fc;
      b = c;
    }
  }
  return c
};

class Errorfunction extends Component {
  state = {
    loaded: false,
    x: [-5.5, 5.5],
    y: [-5.5, 5.5],
    V: (x, y) => (1+x*x)*(1+y*y),
    se2d: null,
  };
  
  SE2D = null;

  componentDidMount() {
    loadMatslise.then(({Matslise, SE2D}) => {
      this.SE2D = SE2D;

      this.setState({
        loaded: true,
        se2d: new SE2D(this.state.V, ...this.state.x, ...this.state.y, 32, 32),
      });
    })
  }

  render() {
    if(!this.state.loaded)
      return "Loading";
    const {se2d} = this.state;
    window.se2d = se2d;

    const xmin = 2, xmax = 12;
    const numPoints = 500;

    const eigenvalues = Es
    .filter(E => E > xmin && E < xmax)
    .map(E => binarySearch(
      se2d.calculateError.bind(se2d), E-.01, E+.01));
    console.log(eigenvalues);

    const f = window.f = (e) => {
      const err = se2d.calculateError(e);
      const dErr = se2d.calculateError(e+1.1*(xmax-xmin)/numPoints) - err;
      if(dErr > .5) {
        return Number.NaN;
      }
      return err;
    }

    return <Graph
    yAxis={false}
    xAxis={false}
    grid={false}
    strokeWidth={4}
    x={[xmin,xmax]}
    numPoints={numPoints}
    func={f}>

    {eigenvalues.map((e, i) => 
      <ReferenceDot key={"eigen"+i} x={e} y={0} r={7} stroke="none" fill="#999" />
    )}</Graph>
  }
}

export default Errorfunction;