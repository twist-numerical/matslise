import { Component } from "react";

class MatsliseGUI extends Component {
  state = {};

  componentDidMount() {
    const worker = new Worker("./se2d.worker.js", { type: "module" });
    worker.onmessage = event => {
      console.log("worker: " + event.data);
    };
    worker.postMessage(42);
  }

  render() {
    return "";
  }
}

export default MatsliseGUI;
