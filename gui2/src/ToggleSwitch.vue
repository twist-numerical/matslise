<template lang="pug">
  .toggle-switch
    input(type="checkbox" :id="realName" :checked="checked" @change="(e) => {$emit('change', e.target.checked)}" )
    label(:for="realName")
</template>

<script lang="ts">
import Vue from "vue";
export default Vue.extend({
  model: {
    prop: "checked",
    event: "change"
  },
  data() {
    return {
      realName: (this as any).getName()
    };
  },
  methods: {
    getName() {
      if (this.name === undefined)
        return "toggle-switch-" + Math.floor(1000000 * Math.random());
      return this.name;
    }
  },
  watch: {
    name() {
      this.realName = this.getName();
    }
  },
  props: ["checked", "name"]
});
</script>

<style lang="scss" scoped>
.toggle-switch {
  $blue: #007bff;
  $white: #ffffff;
  $gray-400: #ced4da;
  $gray: #6c757d;

  input {
    display: none;
  }

  input + label {
    background-color: $gray;
    width: 2.5rem;
    box-sizing: content-box;
    margin: 0;
    color: $white;
    border-radius: 0.25rem;
    position: relative;
    border: 1px solid $gray-400;
    overflow: hidden;
    height: 1rem;
    cursor: pointer;

    &:after {
      content: " ";
      position: absolute;
      left: 0;
      top: 0;
      bottom: 0;
      background-color: $white;
      width: 1rem;
      height: 100%;
      transition: left 0.2s;
    }
  }

  input:checked + label {
    background-color: $blue;

    &:after {
      left: 1.5rem;
    }
  }
}
</style>