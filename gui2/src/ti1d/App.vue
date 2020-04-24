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
                label.m-0.pr-1(for="toggle-show-potential") show
                toggle-switch(v-model="showPotential" name="toggle-show-potential")
        .form-row
          label.input-group.col-6
            .input-group-prepend
              .input-group-text x<sub class="pr-1">min</sub>=
            input.form-control(
                v-bind:disabled="controller.problem.symmetric",
                @input="controller.reset()",
                v-model="controller.problem.symmetric?`-(${controller.problem.x[1]})`:controller.problem.x[0]")
          label.input-group.col-6
            .input-group-prepend
              .input-group-text x<sub class="pr-1">max</sub>=
            input.form-control(
                @input="controller.reset()",
                v-model="controller.problem.x[1]")
        .form-row
          .input-group.col-12
            input.form-control#problem-ymin0(
                v-bind:disabled="controller.problem.symmetric",
                @input="controller.reset()",
                v-model="controller.problem.symmetric?controller.problem.ymax[0]:controller.problem.ymin[0]")
            label.input-group-append.input-group-prepend(for="problem-ymin0")
              span.input-group-text y(x<sub>min</sub>)
            input.form-control#problem-ymin1(
                v-bind:disabled="controller.problem.symmetric",
                @input="controller.reset()",
                v-model="controller.problem.symmetric?controller.problem.ymax[1]:controller.problem.ymin[1]")
            label.input-group-append(for="problem-ymin1")
              span.input-group-text y'(x<sub>min</sub>) = 0
        .form-row
          .input-group.col-12
            input.form-control#problem-ymax0(
                @input="controller.reset()",
                v-model="controller.problem.ymax[0]")
            label.input-group-append.input-group-prepend(for="problem-ymax0")
              span.input-group-text y(x<sub>max</sub>)
            input.form-control#problem-ymax1(
                @input="controller.reset()",
                v-model="controller.problem.ymax[1]")
            label.input-group-append(for="problem-ymax1")
              span.input-group-text y'(x<sub>max</sub>) = 0
        .form-row
          label.input-group.col-12.col-lg-6
            .input-group-prepend
              .input-group-text Tolerance
            input.form-control(
                @input="controller.reset()",
                v-model="controller.problem.tolerance")
          label.input-group.col-12.col-lg-6.justify-content-end
            .input-group-prepend
              .input-group-text Symmetric
            .input-group-append
              .input-group-text
                toggle-switch(
                    v-model="controller.problem.symmetric" 
                    @change="controller.reset()")
        .form-row
          .col-6
            a.btn.btn-link(
              href="#",
              @click.prevent="selectProblem=true") Select a problem
          .col-6.text-right
            input.btn.btn-primary(
              type="submit",
              v-bind:disabled="controller.problem.Matslise === null",
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
          :errors="controller.errors"
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
    visibleEigenvalues() {
      const eigenvalues = (this.controller as Controller).eigenvalues;
      if (eigenvalues === null) return [];
      return eigenvalues.filter(e => e.visible);
    }
  }
});
</script>