class WorkerDelegator {
  listeners = {};

  constructor(worker) {
    this.worker = worker;
    this.worker.onmessage = e => {
      const { type, data } = e.data;
      if (!(type in this.listeners)) {
        console.error(`No one is listening to '${type}'`);
      } else {
        this.listeners[type](data, e);
      }
    };
  }

  send(type, data) {
    this.worker.postMessage({
      type: type,
      data: data
    });
  }

  addListener(type, callback) {
    if (type in this.listeners) {
      console.error(`Someone is already listening to '${type}'`);
    } else {
      this.listeners[type] = callback;
    }
  }
}

export default WorkerDelegator;
