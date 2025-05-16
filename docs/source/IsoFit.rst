.. AIM Documentation documentation master file, created by
   sphinx-quickstart on Fri May 16 14:38:34 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

IsoFit
===============================

The isotherm models suported in **IsoFit** are:

.. list-table:: IsoFit Isotherm Models
   :header-rows: 1
   :widths: auto

   * - Isotherm Model
     - Isotherm Equation
     - Fitting Parameters
   * - Langmuir
     - :math:`q = \frac{q_\mathrm{max} b c}{1 + b c}`
     - :math:`q_\mathrm{max},\ b`
   * - Freundlich
     - :math:`q = K c^{1/n}`
     - :math:`K,\ n`
   * - Toth
     - :math:`q = \frac{q_\mathrm{max} b c}{\left(1 + (b c)^t\right)^{1/t}}`
     - :math:`q_\mathrm{max},\ b,\ t`
   * - Sips
     - :math:`q = \frac{q_\mathrm{max} (b c)^n}{1 + (b c)^n}`
     - :math:`q_\mathrm{max},\ b,\ n`

---------------------------------------
.. note:: The isotherm models are fitted using the Levenberg-Marquardt algorithm.
   The fitting parameters are estimated using the initial guess values provided by the user.
   The initial guess values are set to 1.0 for all parameters if not specified.

The isotherm fitting results from **IsoFit** and **HeatFit** can be directly loaded into **MixPred** and **BreakLab** modules,
providing seamless integration.



