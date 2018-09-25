import React, { Component } from 'react';
import './App.css';
import MatsliseGUI from './components/MatsliseGUI'
import { Switch, Route } from 'react-router-dom'
import Poster from './components/Poster'

class App extends Component {
  render() {
    return <Switch>
      <Route path='/' exact component={Poster} />
      <Route path='/matslise' component={MatsliseGUI} />
    </Switch>
  }
}

export default App;
