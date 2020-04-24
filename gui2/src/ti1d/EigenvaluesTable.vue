<template lang="pug">
div
  div.alert.alert-danger(v-if="errors.length > 0")
    p Some errors were encountered.
    ul
      li(v-for="error of errors") {{ error }}
    p A possible solution could be to change the tolerance.
  table.table.table-striped
    thead
      tr
        th
        th #
        th Eigenvalue
        th Error
    tbody
      tr(v-for="(eigenvalue, index) of eigenvalues")
        th.legend-label
          input(:id="`eigenvalue-${index}`",type="checkbox",v-model:checked="eigenvalue.visible")
          label(
              :for="`eigenvalue-${index}`",
              :style="`background:${eigenvalue.color}`")
        td {{ eigenvalue.index }}
        td {{ eigenvalue.value.toPrecision(13) }}
        td {{ eigenvalue.error == 0 ? 0 : eigenvalue.error.toExponential(2) }}
  div.mb-4.text-center(v-if="eigenvalues !== null")
    a(
        @click="(e)=>{$emit('more-eigenvalues');e.preventDefault();}"
        href="#") + more eigenvalues
</template>

<script lang="ts">
import Vue from "vue";

export default Vue.extend({
  props: ["eigenvalues", "errors"]
});
</script>

<style scoped lang="scss">
.legend-label {
  label {
    padding: 3px;
    margin: 0;
    cursor: pointer;

    &:before {
      display: block;
      width: 1rem;
      height: 1rem;
      content: " ";
      background: white;
    }
  }
  input {
    display: none;
  }
  input:checked + label:before {
    background: transparent;
  }
}
</style>