<template lang="pug">
  div(ref="graph-container")
</template>

<script lang="ts">
import Vue from "vue";
import Chart from "chart.js";

function ensureLength<T>(arr: T[], length: number, newElement: () => T) {
  while (arr.length < length) arr.push(newElement());
  while (arr.length > length) arr.pop();
}

export default Vue.extend({
  props: ["eigenfunctions", "x"],
  data() {
    const canvas = document.createElement("canvas");
    const context = canvas.getContext("2d");

    return {
      canvas,
      chart: Chart.Line(context, {
        data: {
          datasets: []
        },
        options: {
          scales: {
            xAxes: [
              {
                type: "linear",
                position: "bottom"
              }
            ]
          }
        }
      })
    };
  },
  mounted() {
    const container = this.$refs["graph-container"];
    container.appendChild(this.canvas);
    this.updateCanvas();
  },
  watch: {
    eigenfunctions: function() {
      this.updateCanvas();
    }
  },
  methods: {
    updateCanvas() {
      this.chart.options.scales.xAxes[0].ticks.min = this.x[0];
      this.chart.options.scales.xAxes[0].ticks.max = this.x[this.x.length - 1];

      ensureLength(
        this.chart.data.datasets,
        this.eigenfunctions.length,
        () => ({
          label: "",
          data: []
        })
      );

      this.eigenfunctions.forEach((eigenfunction: number[], i: number) => {
        const dataset = this.chart.data.datasets[i];
        ensureLength(dataset.data, eigenfunction.length, () => ({
          x: 0,
          y: 0
        }));
        eigenfunction.forEach((y: number, j: number) => {
          const point = dataset.data[j];
          point.x = this.x[j];
          point.y = y;
        });
      });
      console.log(this.chart.data);
      this.chart.update();
    }
  }
});
</script>