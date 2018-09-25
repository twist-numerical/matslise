import React, { Component } from 'react';
import loadMatslise from '../lib/loadMatslise'
import * as THREE from 'three';
import OrbitControls from 'orbit-controls-es6';

const height = 3;

const binarySearch = (f, a, b) => {
  let fa = f(a), fb = f(b);
  if(fa === 0)
    return a;
  if(fb === 0)
    return b;
  if((fa < 0) == (fb < 0))
    throw new Error("f(a) and f(b) must have a different sign");
  let c, fc; 
  while(Math.abs(a - b) > 1e-12) {
    c = (a+b)/2;
    fc = f(c);
    if(fc === 0)
      return c;
    if((fa < 0) == (fc < 0)) {
      fa = fc;
      a = c;
    } else {
      fb = fc;
      b = c;
    }
  }
  return c
};

class Poster extends Component {
  state = {
    loaded: false,
    x: [-5.5, 5.5],
    y: [-5.5, 5.5],
    V: (x, y) => (1+x*x)*(1+y*y),
    se2d: null,
  };

  xn = 100;
  yn = 160;
  x = [];
  y = [];
  
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
    this.x = [];
    this.y = [];
    const {x: [xmin, xmax], y: [ymin, ymax]} = this.state;
    for(let i = 0; i <= this.xn; ++i)
      this.x.push(xmin + i*(xmax-xmin)/this.xn);
    for(let i = 0; i <= this.yn; ++i)
      this.y.push(ymin + i*(ymax-ymin)/this.yn);

    if(!this.state.loaded)
      return <div />;
    return <div>
    <div ref={e => this.injectThree(e)} />
    </div>
  }

  computeEigenfunctions(E) {
    let zs = this.state.se2d.computeEigenfunction(E, this.x, this.y);
    console.log(zs);
    return zs.map(z => {
      let m = Math.max(...z.map(r => Math.max(...r.map(Math.abs))));
      m /= height;
      return z.map(r => r.map(v => v/m));
    });
  }

  buildGeometries(E) {
    return this.computeEigenfunctions(E).map(z => {
      return new THREE.ParametricGeometry((u, v, p) => {
        const vx = u*this.xn, vy = v*this.yn;
        const ix = Math.floor(vx), iy = Math.floor(vy);
        const sx = vx - ix, sy = vy - iy;

        const get = (r, s, i) => s === 0 ? r[i] : (1-s) * r[i] + s * r[i+1]

        p.x = get(this.x, sx, ix);

        if(sx === 0)
          p.y = get(z[ix], sy, iy);
        else
          p.y = (1-sx) * get(z[ix], sy, iy) + sx * get(z[ix+1], sy, iy)

        p.z = get(this.y, sy, iy);
      }, this.xn, this.yn);
    });
  }

  buildScene() {
    this.camera = new THREE.PerspectiveCamera(70, window.innerWidth / window.innerHeight, 0.01, 1000);
    this.camera.position.x = 10;
    this.camera.position.y = 10;
    this.camera.position.z = 10;

    this.scene = new THREE.Scene();

    const colors = [0xff0000, 0xff00, 0xff, 0x00ffff];
    [3.2, 5.5, 7.5, 8].forEach((approxE, i) => {
      const E = binarySearch((E) => this.state.se2d.calculateError(E), approxE-.1, approxE+.1);
      this.buildGeometries(E).forEach(geometry => {
        const material = new THREE.MeshLambertMaterial({
          color: colors[i%colors.length],
          opacity: .5,
          transparent: true,
          side: THREE.DoubleSide,
        });

        const mesh = new THREE.Mesh(geometry, material);
        mesh.position.y = i*height
        this.scene.add(mesh);
      });
    });

    this.scene.add(new THREE.AmbientLight(0x404040));
    [[0,20,20], [0,-50,50]].forEach(([x, y, z]) => {
      const light = new THREE.PointLight(0x808080, 1, 0);
      light.position.x = x;
      light.position.y = y;
      light.position.z = z;
      this.scene.add(light);
    });

    this.renderer = new THREE.WebGLRenderer({antialias: true});
    this.renderer.setSize(window.innerWidth, window.innerHeight);

    const controls = new OrbitControls(this.camera, this.renderer.domElement);
    controls.enabled = true;
    controls.maxDistance = 1500;
    controls.minDistance = 0;
  }

  animate() {   
    requestAnimationFrame(() => this.animate());

    this.renderer.render(this.scene, this.camera);
  }

  injectThree(element) {
    this.buildScene();
    element.appendChild(this.renderer.domElement);
    this.animate();
  }
}

export default Poster;