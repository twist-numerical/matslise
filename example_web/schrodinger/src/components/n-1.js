import React, { Component } from 'react';
import Graph from './Graph'
import loadMatslise from '../lib/loadMatslise'
import { ReferenceDot, Label } from 'recharts';

let Matslise = null;

class Nm1 extends Component {
  state = {
    f: (x) => 1+x*x,
    x: [-5.5, 5.5],
    steps: 32,
    matsliseLoaded: false,
  };
  toDelete = [];

  componentDidMount() {
    if(!this.state.matsliseLoaded) {
      loadMatslise.then((_Matslise) => {
        Matslise = _Matslise;
        this.setState({matsliseLoaded: true});
      });
    }
  }

  normalize(f) {
    const [xmin, xmax] = this.state.x;
    const n = 200;
    const scale = Math.max(...Array.from({length: n}, (x,i) => {
      return f(xmin + (i+.2+.6*Math.random())/n*(xmax - xmin))
    }).map(Math.abs));

    const r = (x) => {
      return f(x)/scale;
    };
    r.delete = () => f.delete();
    return r;
  }

  render() {
    if(!this.state.matsliseLoaded)
      return "Loading";
    while(this.toDelete.length)
      this.toDelete.pop().delete();

    const {x, f, steps} = this.state;

    const matslise = new Matslise.Matslise(f, ...x, steps);
    this.toDelete.push(matslise);

    const eigenvalues = matslise.eigenvaluesByIndex(0, 6, [0,1], [0,1])
    const eigenfunctions = eigenvalues.map(e => {
      const f = matslise.eigenfunction(e.value, [0,1], [0,1]);
      this.toDelete.push(f)
      const normalized = this.normalize(x => f(x)[0]);
      return (x) => normalized(x) + e.value;
    });

    return (
      <div className="container-fluid h-100">
      <Graph
      yAxis={false}
      xAxis={false}
      grid={false}
      strokeWidth={6}
      x={[-5,5]}
      y={[
        eigenvalues[0].value-3,
        eigenvalues[eigenvalues.length-1].value
        ]}
        func={eigenfunctions.slice(0,-1)}>

        {eigenvalues.map(e => 
          <ReferenceDot key={e.index} x={-4.5} y={e.value} stroke="none" fill="none">
          <Label x={-5} y={e.value} content={({viewBox: {x, y}}) => {
            return (e.index === eigenvalues.length-1?
              <text x={x+30} y={y+100} fill="#666666" fontSize="50px">&#x2807;</text> 
              : 
              <text x={x} y={y+65} fill="#666666" fontSize="35px">
              b
              <tspan baselineShift="sub">{e.index}</tspan>
              <tspan dx="-15px" baselineShift="super">(k)</tspan>
              (<tspan fontWeight="bold">y</tspan>)
              </text>);
          }} />
          </ReferenceDot>)};
        </Graph>
        </div>
        );
  }

  eigenvalues(imin, imax, matslise=this.state.matslise) {
    return 
  }

  eigenfunction(e) {
  }
}

export default Nm1;