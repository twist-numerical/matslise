<template lang="pug">
div
  h2 Eigenvalues
  table.table
    thead
      tr
        th Index
        th Eigenvalue
        th Error
    tbody(v-if="eigenvalues !== null")
      tr(v-for="[index, value, error] of eigenvalues")
        td {{ index }}
        td {{ value.toPrecision(13) }}
        td {{ error == 0 ? 0 :error.toExponential(2) }}
      tr
        td(colspan=3)
          a(
            @click="(e)=>{$emit('more-eigenvalues');e.preventDefault();}"
            href="#") More eigenvalues
</template>

<script lang="ts">
import Vue from "vue";

export default Vue.extend({
  props: ["eigenvalues"]
});
</script>