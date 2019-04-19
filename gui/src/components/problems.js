const dirichlet = {
  ymin: ["1", "0"],
  ymax: ["1", "0"]
};

export default {
  Mathieu: {
    f: "2*cos(2*x)",
    x: ["0", "pi"],
    ...dirichlet
  },
  "Coffey Evans": {
    f: "-2*20*cos(2*x)+20^2*sin(2*x)^2",
    x: ["-pi/2", "pi/2"],
    ...dirichlet,
    symmetric: true
  },
  "Quartic Anharm. Osc.": {
    f: "x^4+x^2",
    x: ["-7", "7"],
    ...dirichlet,
    symmetric: true
  },
  Airy: {
    f: "x",
    x: ["0", "20"],
    ...dirichlet
  },
  Bessel: {
    f: ".25/x^2",
    x: ["0", "1"],
    ...dirichlet
  },
  Hydrogen: {
    f: "-1/x+2/x^2",
    x: ["0", "150"],
    ...dirichlet
  },
  Marletta: {
    f: "3*(x-31)/(4*(x+1)*(x+4)^2)",
    x: ["0", "12"],
    ymin: ["5", "8"],
    ymax: ["1", "0"]
  }
};
