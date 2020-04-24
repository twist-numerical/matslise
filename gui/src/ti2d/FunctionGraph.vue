<template lang="pug">
  div(ref="graph-container")
</template>

<script lang="ts">
import Vue from "vue";
import Chart from "chart.js";
import math from "mathjs-expression-parser";
import Color from "color";

interface DataPoint {
  x: number;
  y: number;
  f: number;
}

const red = Color("rgb(255,99,132)");
const blue = Color("rgb(54,162,235)");
const white = Color("rgb(255,255,255)");

export default Vue.extend({
  props: {
    x: Array,
    y: Array,
    f: Array,
    xLimit: Array,
    yLimit: Array,
    symmetric: {
      type: Boolean,
      default: false
    },
    colors: {
      type: Array,
      default() {
        return [blue, white, red];
      }
    }
  },
  canvas: null,
  chart: null,

  data() {
    const chart = this.$options.chart;
    return {
      xPixelsPerUnit: 1,
      yPixelsPerUnit: 1
    };
  },

  beforeCreate() {
    const onAxisUpdate = () => {
      if (!this.$options.chart || this.$options.updatingRadius)
        return;

      const scales = this.$options.chart.scales;

      this.xPixelsPerUnit =
        scales.xAxis.getPixelForValue(1) - scales.xAxis.getPixelForValue(0);

      this.yPixelsPerUnit =
        scales.yAxis.getPixelForValue(1) - scales.yAxis.getPixelForValue(0);
    };

    this.$options.updatingRadius = false;
    this.$options.canvas = document.createElement("canvas");
    this.$options.chart = new Chart(this.$options.canvas, <any>{
      type: "scatter",
      data: {
        datasets: [
          {
            data: []
          }
        ]
      },
      options: {
        aspectRatio: 1,
        tooltips: {
          enabled: false
        },
        hover: false,
        legend: {
          display: false
        },
        elements: {
          point: {
            pointStyle: "rect",
            radius: 6,
            hoverRadius: 0,
            hitRadius: 0,
            backgroundColor: "white"
          }
        },
        scales: {
          xAxes: [
            {
              id: "xAxis",
              type: "linear",
              position: "bottom",
              ticks: {},
              afterUpdate: onAxisUpdate
            }
          ],
          yAxes: [
            {
              id: "yAxis",
              type: "linear",
              position: "left",
              ticks: {},
              afterUpdate: onAxisUpdate
            }
          ]
        }
      }
    });
    (window as any).chart = this.$options.chart;
  },
  mounted() {
    const container = this.$refs["graph-container"];
    container.appendChild(this.$options.canvas);
    this.updateAspectRatio();
    this.updateCanvas();
  },
  computed: {
    toWatchAspect() {
      return [this.xLimit, this.yLimit];
    },
    toWatchGraph() {
      return [this.dataPoints, this.f];
    },
    fMin() {
      return Math.min(...this.f.map((r: number[]) => Math.min(...r)));
    },
    fMax() {
      return Math.max(...this.f.map((r: number[]) => Math.max(...r)));
    },
    radius() {
      const r = Math.ceil(
        0.5 *
        1.5 * // 1.5 ~= sqrt(2)
          Math.max(
            this.xStep * this.xPixelsPerUnit,
            this.yStep * this.yPixelsPerUnit
          )
      );
      return r;
    },
    xStep() {
      return Math.max(
        ...this.x.map((v: number, i: number) =>
          i == 0 ? 0 : Math.abs(v - this.x[i - 1])
        )
      );
    },
    yStep() {
      return Math.max(
        ...this.y.map((v: number, i: number) =>
          i == 0 ? 0 : Math.abs(v - this.y[i - 1])
        )
      );
    },
    dataPoints() {
      const data: DataPoint[] = (this.$options.chart.data.datasets[0].data = []);

      for (const y of this.y)
        for (const x of this.x) {
          data.push({ x, y, f: 0 });
        }

      return data;
    }
  },
  watch: {
    toWatchAspect: {
      handler: function() {
        this.updateAspectRatio();
      },
      deep: true
    },
    toWatchGraph: {
      handler: function() {
        this.updateCanvas();
      },
      deep: true
    },
    radius: function() {
      const chart = this.$options.chart;
      chart.options.elements.point.radius = this.radius;
      this.$options.updatingRadius = true;
      chart.update();
      this.$options.updatingRadius = false;
    }
  },
  methods: {
    updateAspectRatio() {
      const chart = this.$options.chart;

      const x = this.xLimit;
      const y = this.yLimit;

      const xAxis = chart.options.scales.xAxes[0];
      const yAxis = chart.options.scales.yAxes[0];

      xAxis.ticks.min = x[0];
      xAxis.ticks.max = x[1];
      yAxis.ticks.min = y[0];
      yAxis.ticks.max = y[1];

      chart.aspectRatio =
        (this.xLimit[1] - this.xLimit[0]) / (this.yLimit[1] - this.yLimit[0]);
      chart.update();
    },
    getColor(value: number) {
      const n = this.colors.length;
      value *= n - 1;
      for (let i = 1; i < n; ++i, value -= 1)
        if (value < 1 || i + 1 == n)
          return this.colors[i - 1].mix(this.colors[i], value);
    },
    updateCanvas() {
      const chart = this.$options.chart;

      const max = this.symmetric ? Math.max(this.fMax, -this.fMin) : this.fMax;
      const min = this.symmetric ? -max : this.fMin;
      chart.options.elements.point.backgroundColor = ({
        dataset,
        dataIndex
      }: any) => {
        const value = dataset.data[dataIndex].f;
        return this.getColor((value - min) / (max - min));
      };

      const data: DataPoint[] = this.dataPoints;
      let i = 0;
      for (const row of this.f)
        for (const value of row) {
          data[i++].f = value;
        }
      chart.update();
    }
  }
} as any);
</script>