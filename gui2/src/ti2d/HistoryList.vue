<template lang="pug">
ul.history
  li(v-for="{index, time, name} of reversedHistory" :key="index")
    span.time {{time === null ? '' : `${(time/1000).toFixed(1)}s`}}
    span.name {{name}}
</template>

<script lang="ts">
import Vue from "vue";
export default Vue.extend({
  props: ["history"],
  computed: {
    reversedHistory: function() {
      return this.history
        .map((v: any, i: number) => ({
          index: i,
          ...v
        }))
        .reverse();
    }
  }
});
</script>

<style lang="scss" scoped>
ul.history {
  margin: 0;
  padding: 0;
  li {
    margin: 0;
    padding: 0 0 0 3em;
    display: block;
    border-top: 1px black solid;

    .time {
      display: inline-block;
      margin-left: -3em;
      width: 3em;
    }
  }
}
</style>