.. geppy documentation master file, created by
   sphinx-quickstart on Mon Jul 23 23:21:14 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. image:: _static/geppy-icon.png
   :width: 350 px
   :align: left

=====================================
Welcome to *geppy*'s documentation!
=====================================

*geppy* is an evolutionary algorithm framework specially designed for gene expression programming (GEP) in Python. *geppy* is built on top of the more general evolutionary computation framework `DEAP <https://github.com/deap/deap>`_ , which lacks support for GEP by itself. *geppy* conforms to DEAP's design philosophy that it seeks to make algorithms explicit and data structures transparent in GEP. In the following, please first check the **Get Started** part to learn the core concepts of GEP and the related important data structures, with which you can build your own GEP application easily. Then, you can go through the **Tutorials & Examples** to get yourself familiar with *geppy*.  All the APIs of *geppy* are well documented in details and can be found in the **Library Reference** part. Since *geppy* depends on DEAP, do remember to check the documentation of DEAP if you get stuck, especially the `Overview <http://deap.readthedocs.io/en/master/overview.html>`_ of how a DEAP program is composed. 

Get started
===============
* :doc:`Installation <installation>`
* :doc:`Introduction to GEP theory <intro_GEP>`
* :doc:`Overview of geppy for GEP implementation <overview>`
* :doc:`Conventions of genetic operator design and registration in geppy <convention>`

.. _tutorial_example:

Tutorials & Examples
================================

Simple symbolic regression
--------------------------------
1. `Boolean model identification <https://github.com/ShuhuaGao/geppy/blob/master/examples/sr/Boolean_function_identification.ipynb>`_ (Getting started)
2. `Simple numerical expression inference 1 <https://github.com/ShuhuaGao/geppy/blob/master/examples/sr/numerical_expression_inference-ENC.ipynb>`_  (Using ephemeral numerical constants (ENC))
3. `Simple numerical expression inference 2 <https://github.com/ShuhuaGao/geppy/blob/master/examples/sr/numerical_expression_inference-RNC.ipynb>`_ (The GEP-RNC algorithm for random numerical constant evolution) 

Advanced symbolic regression
------------------------------------------
1. `Improving symbolic regression with linear scaling <https://github.com/ShuhuaGao/geppy/blob/master/examples/sr/numerical_expression_inference-Linear_scaling.ipynb>`_ (Paper: *Improving    Symbolic Regression with Interval Arithmetic and Linear Scaling*)
2. `Apply symbolic regression to teh UCI Power Plant dataset <https://github.com/ShuhuaGao/geppy/blob/master/examples/sr/GEP_RNC_for_ML_with_UCI_Power_Plant_dataset.ipynb>`_

.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. _lib_ref:

Library Reference
================================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
