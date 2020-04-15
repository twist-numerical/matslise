<template lang="pug">
  div
    .modal-backdrop.show(v-show="value")
    .modal(
        :style="`display:${value?'block':'none'}`")
      .modal-dialog
        .modal-content
          .modal-header
            h4.modal-title Select an example problem
          .modal-body
            a.btn.btn-link(
              v-for="problem of problems",
              @click.prevent="selectProblem(problem)",
              href="#",) {{problem.name}}
          .modal-footer
            button.btn.btn-link(@click.prevent="hideModal()") Cancel
</template>

<script lang="ts">
import Vue from "vue";

type Problem = {
  name: string;
  symmetric: boolean;
  potential: string;
  x: [string, string];
  ymin: [string, string];
  ymax: [string, string];
  tolerance: string;
};

const dirichlet = {
  ymin: ["1", "0"],
  ymax: ["1", "0"]
};

const problems: Problem[] = [
  <Problem>{
    name: "Airy",
    symmetric: false,
    potential: "x",
    x: ["0", "20"],
    tolerance: "1e-6",
    ...dirichlet
  },
  <Problem>{
    name: "Coffey Evans",
    symmetric: true,
    potential: "-2*20*cos(2*x)+20^2*sin(2*x)^2",
    x: ["-pi/2", "pi/2"],
    tolerance: "1e-5",
    ...dirichlet
  },
  <Problem>{
    name: "Bessel",
    symmetric: false,
    potential: ".25/x^2",
    x: ["0", "1"],
    tolerance: "1e-6",
    ...dirichlet
  },
  <Problem>{
    name: "Hydrogen",
    symmetric: false,
    potential: "-1/x+2/x^2",
    x: ["0", "200"],
    tolerance: "1e-6",
    ...dirichlet
  },
  <Problem>{
    name: "Marletta",
    symmetric: false,
    potential: "3*(x-31)/(4*(x+1)*(x+4)^2)",
    x: ["0", "12"],
    tolerance: "1e-6",
    ymin: ["5", "8"],
    ymax: ["1", "0"]
  },
  <Problem>{
    name: "Mathieu",
    symmetric: false,
    potential: "2*cos(2*x)",
    x: ["0", "pi"],
    tolerance: "1e-5",
    ...dirichlet
  },
  <Problem>{
    name: "Quartic Anharm. Osc.",
    symmetric: true,
    potential: "x^4+x^2",
    x: ["-7", "7"],
    tolerance: "1e-6",
    ...dirichlet
  }
];

export default Vue.extend({
  props: ["value", "problem"],
  data() {
    return { problems };
  },
  methods: {
    hideModal() {
      this.$emit("input", false);
    },
    selectProblem(problem: Problem) {
      // Non reactive copy
      problem = JSON.parse(JSON.stringify(problem));
      this.problem.symmetric = problem.symmetric;
      this.problem.potential = problem.potential;
      this.problem.x = problem.x.map(s => "" + s);
      this.problem.ymin = problem.ymin;
      this.problem.ymax = problem.ymax;
      this.problem.tolerance = problem.tolerance;
      this.hideModal();
      this.$emit("problem-updated");
    }
  }
});
</script>