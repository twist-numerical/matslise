<template lang="pug">
div.container
  problem-selector(
      v-model="selectProblem",
      :problem="controller.problem",
      @problem-updated="controller.calculate()")
  .row
    .col-12
      h1 Matslise
  .row
    .col-6
      h2 Problem
      form(action="" @submit.prevent="controller.calculate")
        .form-row
          label.input-group.col-12
            .input-group-prepend 
              .input-group-text V(x) =
            input.form-control(
                @input="controller.reset()",
                v-model="controller.problem.potential")
            .input-group-append
              .input-group-text
                label.m-0(for="toggle-show-potential") show&nbsp;
                toggle-switch(v-model="showPotential" name="toggle-show-potential")
        .form-row
          label.input-group.col-6
            .input-group-prepend
              .input-group-text x<sub>min</sub>&nbsp;=
            input.form-control(
                v-bind:disabled="controller.problem.xSymmetric",
                @input="controller.reset()",
                v-model="controller.problem.xSymmetric?`-(${controller.problem.x[1]})`:controller.problem.x[0]")
          label.input-group.col-6
            .input-group-prepend
              .input-group-text x<sub>max</sub>&nbsp;=
            input.form-control(
                @input="controller.reset()",
                v-model="controller.problem.x[1]")
        .form-row
          label.input-group.col-6
            .input-group-prepend
              .input-group-text y<sub>min</sub>&nbsp;=
            input.form-control(
                v-bind:disabled="controller.problem.ySymmetric",
                @input="controller.reset()",
                v-model="controller.problem.ySymmetric?`-(${controller.problem.y[1]})`:controller.problem.y[0]")
          label.input-group.col-6
            .input-group-prepend
              .input-group-text y<sub>max</sub>&nbsp;=
            input.form-control(
                @input="controller.reset()",
                v-model="controller.problem.y[1]")
        .form-row
          label.input-group.col-12.col-lg-6
            .input-group-prepend
              .input-group-text Tolerance
            input.form-control(
                @input="controller.reset()",
                v-model="controller.problem.tolerance")
        .form-row
          label.input-group.col-12.col-lg-6.justify-content-end
            .input-group-prepend
              .input-group-text x-symmetric
            .input-group-append
              .input-group-text
                toggle-switch(
                    v-model="controller.problem.xSymmetric" 
                    @change="controller.reset()")
          label.input-group.col-12.col-lg-6.justify-content-end
            .input-group-prepend
              .input-group-text y-symmetric
            .input-group-append
              .input-group-text
                toggle-switch(
                    v-model="controller.problem.ySymmetric" 
                    @change="controller.reset()")
        .form-row
          .col-6
            a.btn.btn-link(
              href="#",
              @click.prevent="selectProblem=true") Select a problem
          .col-6.text-right
            input.btn.btn-primary(
              type="submit",
              v-bind:disabled="disabledSubmit",
              value="Calculate")

      div(v-if="showPotential")
        h2 Potential
        function-graph(
            :f="controller.problem.potential"
            :x="controller.problem.x",
            n="210"
            :symmetric="controller.problem.symmetric")

      div(v-if="controller.eigenvalues !== null")
        h2 Eigenfunctions
        eigenfunctions-graph(      
            :eigenvalues="visibleEigenvalues"
            :x="controller.xValues")
    .col-6
      h2 Eigenvalues
      eigenvalues-table(
          v-if="controller.eigenvalues !== null"
          :eigenvalues="controller.eigenvalues"
          @more-eigenvalues="() => controller.moreEigenvalues()")
      div(v-else)
        p Press 'Calculate' to compute the eigenvalues of this problem.
</template>

<script lang="ts">
import Vue from "vue";
import Controller from "./Controller";
import EigenvaluesTable from "./EigenvaluesTable.vue";
import EigenfunctionsGraph from "./EigenfunctionsGraph.vue";
import FunctionGraph from "./FunctionGraph.vue";
import ProblemSelector from "./ProblemSelector.vue";
import ToggleSwitch from "../util/ToggleSwitch.vue";

export default Vue.extend({
  data() {
    return {
      controller: new Controller(),
      selectProblem: false,
      showPotential: false
    };
  },
  components: {
    "eigenvalues-table": EigenvaluesTable,
    "eigenfunctions-graph": EigenfunctionsGraph,
    "function-graph": FunctionGraph,
    "problem-selector": ProblemSelector,
    "toggle-switch": ToggleSwitch
  },
  computed: {
    disabledSubmit() {
      return !this.controller.problem.ready;
    },
    visibleEigenvalues() {
      const eigenvalues = (this.controller as Controller).eigenvalues;
      if (eigenvalues === null) return [];
      return eigenvalues.filter(e => e.visible);
    }
  }
});
</script>