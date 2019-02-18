import loadMatslise from "../../lib/loadMatslise";

self.addEventListener("message", message => {
  loadMatslise.then(Matslise => {
    console.log(message);
    console.log(Matslise);
  });
});
