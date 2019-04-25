import React from "react";
import ReactDOM from "react-dom";
import "./index.css";
import App from "./App";
import { HashRouter as Router } from "react-router-dom";
import { RouterToUrlQuery } from "react-url-query";

ReactDOM.render(
  <Router>
    <RouterToUrlQuery>
      <App />
    </RouterToUrlQuery>
  </Router>,
  document.getElementById("root")
);
