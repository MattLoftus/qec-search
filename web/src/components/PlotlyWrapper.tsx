// @ts-nocheck
import plotComponentFactory from 'react-plotly.js/factory';
import Plotly from 'plotly.js/dist/plotly';

const factory = plotComponentFactory.default || plotComponentFactory;
const Plot = factory(Plotly.default || Plotly);
export default Plot;
