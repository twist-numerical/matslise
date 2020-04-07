<template lang="pug">
div.container
  problem-selector(
      v-model="selectProblem",
      :problem="problem",
      @problem-updated="calculate()")
  .row
    .col-12
      h1 Matslise
  .row
    .col-6
      h2 Problem
      form(action="" @submit.prevent="calculate")
        .form-row
          label.input-group.col-12
            .input-group-prepend 
              .input-group-text V(x) =
            input.form-control(
                @input="reset()",
                v-model="problem.potential")
            .input-group-append
              .input-group-text
                label.m-0(for="toggle-show-potential") show&nbsp;
                toggle-switch(v-model="showPotentialSwitch" name="toggle-show-potential")
        .form-row
          label.input-group.col-6
            .input-group-prepend
              .input-group-text x<sub>min</sub>&nbsp;=
            input.form-control(
                v-bind:disabled="problem.xSymmetric",
                @input="reset()",
                v-model="problem.xSymmetric?`-(${problem.x[1]})`:problem.x[0]")
          label.input-group.col-6
            .input-group-prepend
              .input-group-text x<sub>max</sub>&nbsp;=
            input.form-control(
                @input="reset()",
                v-model="problem.x[1]")
        .form-row
          label.input-group.col-6
            .input-group-prepend
              .input-group-text y<sub>min</sub>&nbsp;=
            input.form-control(
                v-bind:disabled="problem.ySymmetric",
                @input="reset()",
                v-model="problem.ySymmetric?`-(${problem.y[1]})`:problem.y[0]")
          label.input-group.col-6
            .input-group-prepend
              .input-group-text y<sub>max</sub>&nbsp;=
            input.form-control(
                @input="reset()",
                v-model="problem.y[1]")
        .form-row
          label.input-group.col-12.col-lg-6
            .input-group-prepend
              .input-group-text Tolerance
            input.form-control(
                @input="reset()",
                v-model="problem.tolerance")
        .form-row
          label.input-group.col-12.col-lg-6.justify-content-end
            .input-group-prepend
              .input-group-text x-symmetric
            .input-group-append
              .input-group-text
                toggle-switch(
                    v-model="problem.xSymmetric" 
                    @change="reset()")
        //-
          label.input-group.col-12.col-lg-6.justify-content-end
            .input-group-prepend
              .input-group-text y-symmetric
            .input-group-append
              .input-group-text
                toggle-switch(
                    v-model="problem.ySymmetric" 
                    @change="reset()")
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
            v-if="showPotential"
            :f="problem.potentialData"
            :x="problem.parsed.xPoints"
            :y="problem.parsed.yPoints"
            :xLimit="problem.parsed.x"
            :yLimit="problem.parsed.y")

      div(v-if="eigenvalues !== null")
        h2 Eigenfunctions
        div(v-for="eigenvalue of visibleEigenvalues")
          function-graph(
            v-if="eigenvalue.eigenfunctions !== null"
            v-for="eigenfunction of eigenvalue.eigenfunctions"
            :f="eigenfunction"
            :x="problem.parsed.xPoints"
            :y="problem.parsed.yPoints"
            :xLimit="problem.parsed.x"
            :yLimit="problem.parsed.y"
            :symmetric="true")
    .col-6
      h2 Eigenvalues
      eigenvalues-table(
          v-if="eigenvalues !== null"
          :eigenvalues="eigenvalues"
          @more-eigenvalues="() => moreEigenvalues()")
      div(v-else)
        p Press 'Calculate' to compute the eigenvalues of this problem.
      h3 History
      history-list(:history="problem.history")
</template>

<script lang="ts">
import Vue from "vue";
import Color from "color";
import EigenvaluesTable from "./EigenvaluesTable.vue";
import FunctionGraph from "./FunctionGraph.vue";
import ProblemSelector from "./ProblemSelector.vue";
import ToggleSwitch from "../util/ToggleSwitch.vue";
import HistoryList from "./HistoryList.vue";
import Problem from "./Problem";
import Eigenvalue from "./Eigenvalue";

export default Vue.extend({
  data() {
    return {
      eigenvalues: null,
      problem: new Problem(),
      selectProblem: false,
      showPotentialSwitch: false
    };
  },
  components: {
    "eigenvalues-table": EigenvaluesTable,
    "function-graph": FunctionGraph,
    "problem-selector": ProblemSelector,
    "toggle-switch": ToggleSwitch,
    "history-list": HistoryList
  },
  computed: {
    disabledSubmit() {
      return !this.problem.ready;
    },
    visibleEigenvalues() {
      const eigenvalues = this.eigenvalues;
      if (eigenvalues === null) return [];
      return eigenvalues
        .filter((e: any) => e.visible)
        .map((e: Eigenvalue) => {
          e.ensureEigenfunctions();
          return e;
        });
    },
    showPotential() {
      return this.showPotentialSwitch && this.problem.parsed !== null;
    }
  },
  methods: {
    reset() {
      this.problem.reset();
      this.eigenvalues = null;
    },
    async calculate() {
      this.reset();
      await this.problem.parse();
      await this.moreEigenvalues();
    },
    containsEigenvalue(E: number) {
      if (this.eigenvalues === null) return false;
      for (const { value } of this.eigenvalues) {
        if (Math.abs(E - value) < 10 * this.problem.parsed!.tolerance)
          return true;
      }
      return false;
    },
    async moreEigenvalues() {
      let eigenvaluesFound;
      if (this.eigenvalues === null) {
        eigenvaluesFound = await this.problem.eigenvaluesByIndex(0, 2);
        this.eigenvalues = [];
      } else {
        const max = this.eigenvalues.length - 1;
        const min = Math.max(0, max - 5);
        eigenvaluesFound = await this.problem.eigenvalues(
          this.eigenvalues[max].value,
          2 * this.eigenvalues[max].value - this.eigenvalues[min].value
        );
      }
      const self = this;
      for (const [value, error] of eigenvaluesFound) {
        if (!this.containsEigenvalue(value)) {
          const index = this.eigenvalues!.length;
          this.eigenvalues.push(
            new Eigenvalue(index, value, error, this.problem)
          );
        }
      }
    }
  }
});
</script>