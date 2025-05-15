---
title: AIM Overview
---

# My Software Modules

<style>
.tab {
  display: none;
}
.tab-button {
  background-color: #eee;
  border: none;
  padding: 10px;
  cursor: pointer;
}
.tab-button.active {
  background-color: #ddd;
  font-weight: bold;
}
.tab-content {
  display: none;
  padding: 10px;
  border: 1px solid #ccc;
}
</style>

<div>
  <button class="tab-button" onclick="showTab('mod1')">IsoFit</button>
  <button class="tab-button" onclick="showTab('mod2')">HeatFit</button>
  <button class="tab-button" onclick="showTab('mod3')">MixPred</button>
  <button class="tab-button" onclick="showTab('mod4')">BreakLab</button>
</div>

<div id="mod1" class="tab-content">
  <h3>Isotherm Fitting</h3>
  <p>This module fits adsorption isotherms to single temperture isotherm data</p>
</div>
<div id="mod2" class="tab-content">
  <h3>Multi-temperature Isotherm Fitting and Isosteric Heat of Prediction</h3>
  <p>This module fits adsorption isotherm to multiple temperature isotherm data and also predicts isosteric heat of adsorption using Clausius-Clapeyron and Virial equation</p>
</div>
<div id="mod3" class="tab-content">
  <h3>Mixture Isotherm Prediction</h3>
  <p>This module predicts mixture adsorption isotherm for the given compositions using Ideal Adsorption Solution Theory (IAST) and Extended Dual-site Langmuir (EDSL) model</p>
</div>
<div id="mod4" class="tab-content">
  <h3>Breakthrough Simulation</h3>
  <p>This module simulates Isothermal/Non-isothermal, Non-isobaric fixed bed adsorption breakthrough for up to 5 components</p>
</div>

<script>
function showTab(id) {
  var tabs = document.querySelectorAll('.tab-content');
  var buttons = document.querySelectorAll('.tab-button');
  tabs.forEach(t => t.style.display = 'none');
  buttons.forEach(b => b.classList.remove('active'));
  document.getElementById(id).style.display = 'block';
  event.target.classList.add('active');
}
document.querySelector('.tab-button').click(); // Show first tab on load
</script>
