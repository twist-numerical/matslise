<template lang="pug">
div.container
  h1 Matslise
  form(action="" @submit.prevent="matslise.calculate")
    .form-row
      label.input-group.col-12
        .input-group-prepend 
          .input-group-text V(x) =
        input.form-control(
            @input="matslise.reset()",
            v-model="matslise.problem.potential")
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
      label.input-group.col-6
        .input-group-prepend
          .input-group-text Tolerance
        input.form-control(
            @input="matslise.reset()",
            v-model="matslise.problem.tolerance")
      label.input-group.col-6.justify-content-end
        .input-group-prepend
          .input-group-text Symmetric
        .input-group-append
          .input-group-text
            input(
                type="checkbox"
                @input="matslise.reset()",
                v-model="matslise.problem.symmetric")
    .form-row.justify-content-end
      input.btn.btn-primary(type="submit" v-bind:disabled="matslise.problem.Matslise === null", value="Calculate")

  eigenfunctions-graph(
    v-if="matslise.eigenvalues !== null"
    :eigenfunctions="visibleEigenfunctions"
    :x="matslise.xValues")

  eigenvalues-table(
      :eigenvalues="matslise.eigenvalues"
      @more-eigenvalues="() => matslise.moreEigenvalues()")
</template>

<script lang="ts">
import Vue from "vue";
import MatsliseController from "./MatsliseController";
import EigenvaluesTable from "./EigenvaluesTable.vue";
import EigenfunctionsGraph from "./EigenfunctionsGraph.vue";

const matsliseController = new MatsliseController();
export default Vue.extend({
  data() {
    return { matslise: matsliseController };
  },
  components: {
    "eigenvalues-table": EigenvaluesTable,
    "eigenfunctions-graph": EigenfunctionsGraph
  },
  computed: {
    visibleEigenfunctions() {
      console.log("test");
      const eigenvalues = (this.matslise as MatsliseController).eigenvalues;
      if (eigenvalues === null) return [];
      const r = eigenvalues.filter(e => e.visible).map(e => e.eigenfunction);
      console.log(r);
      return r;
    }
  }
});
window.mc = matsliseController;
</script>