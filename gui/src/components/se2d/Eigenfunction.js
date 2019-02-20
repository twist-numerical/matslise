import React, { Component } from "react";
import * as THREE from "three";
import OrbitControls from "orbit-controls-es6";

class Eigenfunction extends Component {
  static defaultProps = {
    E: false,
    worker: null,
    xn: 50,
    yn: 50
  };

  state = {
    value: false
  };
  currentGroup = null;

  updateWorker(oldWorker, currentWorker) {
    if (oldWorker !== currentWorker) {
      if (oldWorker) oldWorker.removeListener("eigenfunction");
      if (currentWorker)
        currentWorker.addListener("eigenfunction", state => {
          this.setState(state);
        });
    }
  }

  componentDidMount() {
    this.updateWorker(null, this.props.worker);
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    this.updateWorker(prevProps.worker, this.props.worker);

    if (prevProps.E !== this.props.E) {
      const E = this.props.E;
      const isValue = E === 0 || (!!E && isFinite(E));
      this.setState({
        value: isValue,
        loading: isValue,
        eigenfunction: null
      });
      if (isValue) {
        this.props.worker.send("eigenfunction", {
          E,
          xn: this.props.xn,
          yn: this.props.yn
        });
      }
    }
  }

  render() {
    if (this.currentGroup) {
      this.scene.remove(this.currentGroup);
      this.currentGroup = null;
    }
    if (this.state.eigenfunction) {
      const { x, y, zs } = this.state.eigenfunction;
      this.currentGroup = this.buildObjects(x, y, zs);
      this.scene.add(this.currentGroup);
    }
    return <div ref={div => this.injectThree(div)} />;
  }

  calculateGeometry(x, y, z) {
    const get = (r, s, i) => (s === 0 ? r[i] : (1 - s) * r[i] + s * r[i + 1]);

    return new THREE.ParametricGeometry(
      (u, v, p) => {
        const vx = u * (x.length - 1),
          vy = v * (y.length - 1);
        const ix = Math.floor(vx),
          iy = Math.floor(vy);
        const sx = vx - ix,
          sy = vy - iy;
        p.x = get(x, sx, ix);

        if (sx === 0) p.y = get(z[ix], sy, iy);
        else p.y = (1 - sx) * get(z[ix], sy, iy) + sx * get(z[ix + 1], sy, iy);

        p.z = get(y, sy, iy);
      },
      x.length,
      y.length
    );
  }

  buildGraph(x, y, z) {
    const scale = 3 / Math.max(...z.map(r => Math.max(...r.map(Math.abs))));
    z = z.map(r => r.map(v => v * scale));

    const group = new THREE.Group();
    const materialSettings = {
      color: 0x26b6ff,
      side: THREE.DoubleSide
    };

    group.add(
      new THREE.Mesh(
        this.calculateGeometry(x, y, z),
        new THREE.MeshStandardMaterial(materialSettings)
      )
    );

    return group;
  }

  buildObjects(x, y, zs) {
    const group = new THREE.Group();

    zs.forEach(z => {
      group.add(this.buildGraph(x, y, z));
    });
    return group;
  }

  buildScene() {
    this.camera = new THREE.PerspectiveCamera(
      70,
      window.innerWidth / window.innerHeight,
      0.01,
      1000
    );
    this.camera.position.x = 10;
    this.camera.position.y = 10;
    this.camera.position.z = 10;

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0xffffff);

    this.scene.add(new THREE.AmbientLight(0x404040));
    [[0, 20, 20], [0, -50, 50]].forEach(([x, y, z]) => {
      const light = new THREE.PointLight(0x808080, 1, 0);
      light.position.x = x;
      light.position.y = y;
      light.position.z = z;
      this.scene.add(light);
    });

    this.renderer = new THREE.WebGLRenderer({ antialias: true });
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
    if (!element) return;
    if (!this.renderer) {
      this.buildScene();
      this.animate();
    }
    element.appendChild(this.renderer.domElement);
  }
}

export default Eigenfunction;
