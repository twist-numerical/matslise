import React, { Component } from 'react';
import './App.css';
import MatsliseGUI from './components/MatsliseGUI'
import SE2DGUI from './components/se2d/gui.js'
import { Switch, Route } from 'react-router-dom'
import Poster from './components/Poster'
import Nm1 from './components/n-1'
import Errorfunction from './components/errorfunction'
import Eigenfunctions from './components/eigenfunctions'

class App extends Component {
  render() {
    return <Switch>
      <Route path='/poster' exact component={Poster} />
      <Route path='/eigen/:index/:multiplicity' exact component={Eigenfunctions} />
      <Route path='/n-1' exact component={Nm1} />
      <Route path='/error' exact component={Errorfunction} />
      <Route path='/se2d' component={SE2DGUI} />
      <Route path='/' component={MatsliseGUI} />
    </Switch>
  }
}

export default App;
