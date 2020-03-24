import Vue from "vue";
// @ts-ignore
import App from "./App.vue";
import 'bootstrap-css-only/css/bootstrap.css'

new Vue({ render: createElement => createElement(App) }).$mount("#matslise");
