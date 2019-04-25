import React, { Component } from "react";
import Graph from "../util/Graph";
import loadMatslise from "../../lib/loadMatslise";
import { ReferenceDot } from "recharts";

const Es = []; // [3.20, 5.53, 7.55, 8.04, 8.46, 9.94, 11.34];

class Errorfunction extends Component {
  state = {
    loaded: false,
    x: [-5.5, 5.5],
    y: [-5.5, 5.5],
    V: (x, y) => (1 + x * x) * (1 + y * y),
    se2d: null
  };

  SE2D = null;

  componentDidMount() {
    loadMatslise.then(({ Matslise, SE2D }) => {
      this.SE2D = SE2D;

      this.setState({
        loaded: true,
        se2d: new SE2D(this.state.V, ...this.state.x, ...this.state.y, {
          sectorCount: 33,
          nested: {
            sectorCount: 33
          }
        })
      });
    });
  }

  render() {
    if (!this.state.loaded) return "Loading";
    const { se2d } = this.state;
    window.se2d = se2d;

    const xmin = 2,
      xmax = 12;
    const numPoints = 500;

    const eigenvalues = Es.filter(E => E > xmin && E < xmax).map(E =>
      se2d.findEigenvalue(E)
    );
    console.log(eigenvalues);

    const f = e => {
      return se2d.calculateError(e).first;
    };

    return (
      <Graph
        yAxis={true}
        xAxis={true}
        grid={false}
        strokeWidth={4}
        x={[xmin, xmax]}
        numPoints={numPoints}
        func={f}
        verbose={true}
      >
        {eigenvalues.map((e, i) => (
          <ReferenceDot
            key={"eigen" + i}
            x={e}
            y={0}
            r={7}
            stroke="none"
            fill="#999"
          />
        ))}
      </Graph>
    );
  }
}

export default Errorfunction;
