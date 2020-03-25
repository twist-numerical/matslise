<template lang="pug">
div
  h2 Eigenvalues
  table.table.table-striped
    thead
      tr
        th Index
        th Eigenvalue
        th Error
    tbody
      tr(v-for="eigenvalue of eigenvalues")
        th 
          input(type="checkbox",v-model:checked="eigenvalue.visible")
          span {{ eigenvalue.index }}
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
  props: ["eigenvalues"]
});
</script>