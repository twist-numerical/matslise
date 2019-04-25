import React, { Component, lazy } from "react";
import "./App.css";
import { Switch, Route } from "react-router-dom";

const MatsliseGUI = lazy(() => import("./components/matslise/MatsliseGUI"));
const Poster = lazy(() => import("./components/misc/Poster"));
const SE2DGUI = lazy(() => import("./components/se2d/gui.js"));
const Nm1 = lazy(() => import("./components/misc/n-1"));
const Errorfunction = lazy(() => import("./components/misc/errorfunction"));
const Eigenfunctions = lazy(() => import("./components/misc/eigenfunctions"));

class Loading extends Component {
  render() {
    return "Loading...";
  }
}

class App extends Component {
  render() {
    return (
      <React.Suspense fallback={<Loading />}>
        <Switch>
          <Route path="/poster" exact component={Poster} />
          <Route
            path="/eigen/:index/:multiplicity"
            exact
            component={Eigenfunctions}
          />
          <Route path="/n-1" exact component={Nm1} />
          <Route path="/error" exact component={Errorfunction} />
          <Route path="/se2d" component={SE2DGUI} />
          <Route path="/" exact component={MatsliseGUI} />
        </Switch>
      </React.Suspense>
    );
  }
}

export default App;
