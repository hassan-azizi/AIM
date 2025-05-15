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
  <button class="tab-button" onclick="showTab('mod1')">Module 1</button>
  <button class="tab-button" onclick="showTab('mod2')">Module 2</button>
  <button class="tab-button" onclick="showTab('mod3')">Module 3</button>
</div>

<div id="mod1" class="tab-content">
  <h3>Isotherm Fitting</h3>
  <p>This module fits adsorption isotherms to experimental data...</p>
</div>
<div id="mod2" class="tab-content">
  <h3>Mixture Isotherm</h3>
  <p>This module calculates multicomponent adsorption using IAST...</p>
</div>
<div id="mod3" class="tab-content">
  <h3>Breakthrough Simulation</h3>
  <p>This module simulates packed bed performance for mixtures...</p>
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
