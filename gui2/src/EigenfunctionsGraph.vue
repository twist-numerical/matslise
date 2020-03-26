<template lang="pug">
  div(ref="graph-container")
</template>

<script lang="ts">
import Vue from "vue";
import Chart from "chart.js";
import { Eigenvalue } from "./MatsliseController";

function ensureLength<T>(arr: T[], length: number, newElement: () => T) {
  while (arr.length < length) arr.push(newElement());
  while (arr.length > length) arr.pop();
}

export default Vue.extend({
  props: ["eigenvalues", "x"],
  data() {
    const canvas = document.createElement("canvas");

    return {
      canvas,
      chart: new Chart(canvas, {
        type: "line",
        data: {
          datasets: []
        },
        options: {
          tooltips: {
            enabled: false
          },
          legend: {
            display: false
          },
          elements: {
            point: {
              radius: 0,
              hoverRadius: 0,
              hitRadius: 0
            }
          },
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
    eigenvalues() {
      this.updateCanvas();
    },
    x() {
      this.updateCanvas();
    }
  },
  methods: {
    updateCanvas() {
      this.chart.options.scales.xAxes[0].ticks.min = this.x[0];
      this.chart.options.scales.xAxes[0].ticks.max = this.x[this.x.length - 1];

      const eigenvalueMap: Map<number, Eigenvalue> = new Map(
        this.eigenvalues.map((eigenvalue: Eigenvalue) => [
          eigenvalue.index,
          eigenvalue
        ])
      );

      const eigenvalueSet: Set<number> = new Set(eigenvalueMap.keys());
      const eigenvaluesSorted: number[] = [];
      const remove: number[] = [];
      this.chart.data.datasets.forEach(({ index }: any, i: number) => {
        if (eigenvalueSet.has(index)) {
          eigenvaluesSorted.push(index);
          eigenvalueSet.delete(index);
        } else {
          remove.push(i);
        }
      });
      for (const i of remove) {
        this.chart.data.datasets.splice(i, 1);
      }

      ensureLength(this.chart.data.datasets, this.eigenvalues.length, () => ({
        index: -1,
        label: "",
        data: [],
        fill: false,
        cubicInterpolationMode: "monotone"
      }));
      eigenvaluesSorted.push(...eigenvalueSet);

      eigenvaluesSorted
        .map(i => eigenvalueMap.get(i))
        .forEach(({ index, eigenfunction, color }: Eigenvalue, i: number) => {
          const dataset = this.chart.data.datasets[i];
          dataset.borderColor = color;
          dataset.index = index;
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
      this.chart.update();
    }
  }
});
</script>