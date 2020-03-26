<template lang="pug">
div.container
  problem-selector(
      v-model="selectProblem",
      :problem="matslise.problem",
      @problem-updated="matslise.calculate()")
  .row
    .col-12
      h1 Matslise
  .row
    .col-6
      h2 Problem
      form(action="" @submit.prevent="matslise.calculate")
        .form-row
          label.input-group.col-12
            .input-group-prepend 
              .input-group-text V(x) =
            input.form-control(
                @input="matslise.reset()",
                v-model="matslise.problem.potential")
            .input-group-append
              .input-group-text
                label.m-0(for="toggle-show-potential") show&nbsp;
                toggle-switch(v-model="showPotential" name="toggle-show-potential")
        .form-row
          label.input-group.col-6
            .input-group-prepend
              .input-group-text x<sub>min</sub>&nbsp;=
            input.form-control(
                v-bind:disabled="matslise.problem.symmetric",
                @input="matslise.reset()",
                v-model="matslise.problem.symmetric?`-(${matslise.problem.x[1]})`:matslise.problem.x[0]")
          label.input-group.col-6
            .input-group-prepend
              .input-group-text x<sub>max</sub>&nbsp;=
            input.form-control(
                @input="matslise.reset()",
                v-model="matslise.problem.x[1]")
        .form-row
          .input-group.col-12
            input.form-control#problem-ymin0(
                v-bind:disabled="matslise.problem.symmetric",
                @input="matslise.reset()",
                v-model="matslise.problem.symmetric?matslise.problem.ymax[0]:matslise.problem.ymin[0]")
            label.input-group-append.input-group-prepend(for="problem-ymin0")
              span.input-group-text y(x<sub>min</sub>)
            input.form-control#problem-ymin1(
                v-bind:disabled="matslise.problem.symmetric",
                @input="matslise.reset()",
                v-model="matslise.problem.symmetric?matslise.problem.ymax[1]:matslise.problem.ymin[1]")
            label.input-group-append(for="problem-ymin1")
              span.input-group-text y'(x<sub>min</sub>) = 0
        .form-row
          .input-group.col-12
            input.form-control#problem-ymax0(
                @input="matslise.reset()",
                v-model="matslise.problem.ymax[0]")
            label.input-group-append.input-group-prepend(for="problem-ymax0")
              span.input-group-text y(x<sub>max</sub>)
            input.form-control#problem-ymax1(
                @input="matslise.reset()",
                v-model="matslise.problem.ymax[1]")
            label.input-group-append(for="problem-ymax1")
              span.input-group-text y'(x<sub>max</sub>) = 0
        .form-row
          label.input-group.col-12.col-lg-6
            .input-group-prepend
              .input-group-text Tolerance
            input.form-control(
                @input="matslise.reset()",
                v-model="matslise.problem.tolerance")
          label.input-group.col-12.col-lg-6.justify-content-end
            .input-group-prepend
              .input-group-text Symmetric
            .input-group-append
              .input-group-text
                toggle-switch(
                    v-model="matslise.problem.symmetric" 
                    @change="matslise.reset()")
        .form-row
          .col-6
            a.btn.btn-link(
              href="#",
              @click.prevent="selectProblem=true") Select a problem
          .col-6.text-right
            input.btn.btn-primary(
              type="submit",
              v-bind:disabled="matslise.problem.Matslise === null",
              value="Calculate")

      div(v-if="showPotential")
        h2 Potential
        function-graph(
            :f="matslise.problem.potential"
            :x="matslise.problem.x",
            n="210"
            :symmetric="matslise.problem.symmetric")

      div(v-if="matslise.eigenvalues !== null")
        h2 Eigenfunctions
        eigenfunctions-graph(      
            :eigenvalues="visibleEigenvalues"
            :x="matslise.xValues")
    .col-6
      h2 Eigenvalues
      eigenvalues-table(
          v-if="matslise.eigenvalues !== null"
          :eigenvalues="matslise.eigenvalues"
          @more-eigenvalues="() => matslise.moreEigenvalues()")
      div(v-else)
        p Press 'Calculate' to compute the eigenvalues of this problem.
</template>

<script lang="ts">
import Vue from "vue";
import MatsliseController from "./MatsliseController";
import EigenvaluesTable from "./EigenvaluesTable.vue";
import EigenfunctionsGraph from "./EigenfunctionsGraph.vue";
import FunctionGraph from "./FunctionGraph.vue";
import ProblemSelector from "./ProblemSelector.vue";
import ToggleSwitch from "./ToggleSwitch.vue";

const matsliseController = new MatsliseController();
export default Vue.extend({
  data() {
    return {
      matslise: matsliseController,
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
    visibleEigenvalues() {
      const eigenvalues = (this.matslise as MatsliseController).eigenvalues;
      if (eigenvalues === null) return [];
      return eigenvalues.filter(e => e.visible);
    }
  }
});
(window as any).mc = matsliseController;
</script>